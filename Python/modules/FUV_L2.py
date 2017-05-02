# coding: utf-8
'''
This module contains the algorithms that are used for FUV level 2.5 nighttime data product. 
Includes:
    - the reading and writing of NetCDF files
    - the regularization algorithms to estimate the volume emission rate for the 135.6nm O+ emission during nighttime
    - the calculation of the electron density content assuming radiative recombination and mutual neutralization
'''

# Basic numerical python and math modules
import numpy as np
import math

# Pyglow is needed in order to call MSIS to get the [O] density that is needed for the electron density calculation
import pyglow
# From the datetime module we need the datetime and timedelta functions
import datetime

# From the scipy module we need the optimize, linalg.sqrtm, io.netcdf, interpolate.interp1d
import scipy
from scipy import interpolate

# From the time modeule we need the gmtime, strftime functions
import time
import calendar

# Module that contains all the geometry calculations
import ICON as ic

from scipy import io

#scipy.io.netcdf for the reader
import netCDF4

# For Monte Carlo processing:
from numpy.random import multivariate_normal

from scipy.linalg import sqrtm # matrix square root

'''
From the multiprocssing tool we need Pool. 
This might not be used since we dont know how many cores we are going to have in our disposal.
It is needed when performing Tikhonov regularization and when we calculate the prior for the MAP estimation
'''
import multiprocessing

# From the sys module we need the exit function
import sys


# Function that contains all the instrument parameters for ICON_FUV
def get_FUV_instrument_constants():
    '''
    Notes:
        12-Dec-2014: Last known parameters for ICON FUV   
    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu)
        01-Apr-2015: Last known parameters for ICON FUV - Rescell Dimensions
        28-Aug-2015: Last known parameters for ICON FUV - Sensitivity = 0.083 -> 0.0873
    '''

    instrument = {    'npixelx': 6,            # number of pixels per slice of the interferogram in a rescell(horizontal) [512px total]
                      'npixely': 256,           # number of rescells in vertical (altitude) direction on CCD.
                'aperture_area': 0.006*0.032,   # [cm^2]
                   'coneangle1': 18,            # fov of instrument in horizontal direction, [deg]
                   'coneangle2': 24,            # fov of instrument in vertical direction, [deg]
                 'coneangle1_r': 0.314159265359,# fov of instrument in horizontal direction, [rad]
                 'coneangle2_r': 0.418879020479,# fov of instrument in vertical direction, [rad]
                    'exposure' : 12,            # exposure [sec]
                  'sensitivity': 0.0873,         # combined transmittance of all optics. [counts/res_cell/s/R]
                        'fov_l': 98,            # Field of View (Vertical) Lower Bound (degrees)
                        'fov_u': 123,           # Field of View (Vertical) Upper Bound (degrees)
                       'fovr_l': 1.71042266695, # Field of View (Vertical) Lower Bound (radians)
                       'fovr_u': 2.14675497995, # Field of View (Vertical) Upper Bound (radians)
                 'stripes_used': 1,             # Number of stripes used in the ccd the CCD
                   'Dark Noise': 4000,          # Dark noise [Not to be of our concern, only for reference]
                   'Read Noise': 60,            # Read noise [e-][Not to be of our concern, only for reference]
                  }

    return instrument
 
'''

DISTANCE MATRIX CALCULATIONS

''' 
 
# Calculation of the Distance matrix assuming Spherical Earth    
def create_cells_Matrix_spherical_symmetry(theta,Horbit,RE=6371.):
    '''
    Returns the boundaries for the created cells of the atmosphere assuming sperical symmetry of earth.
    INPUTS:
        theta      - zenith angle vector, [degrees]
        Horbit     - satellite altitude, [km]
        RE         - radius of Earth, [km]
    OUTPUTS:
        S          - Weight Matrix
    NOTES:
        - The number of shells is defined to be the same as the number of observation zenith angles. 
          This will make the distance matrix, square, lower-triangular, invertible.
    HISTORY:
        16-Jun-2016: Written by Dimitrios Iliou (iliou2@illinois.edu) and Brian J. Harding
                     (bhardin2@illinois.edu) 
    '''
    
    # Calculate the number of observation zenith angles
    Ntheta = len(theta)
    
    # Convert zenith angles to radians
    theta = np.deg2rad(theta)
    
    #for angle in theta:
    rbot= ic.ze_to_tang_alt(np.rad2deg(theta), Horbit, RE)
    # Define top of each layer.
    
    rtop = rbot.copy()

    # Create the top boundary of the shell
    rtop[1:] = rbot[:-1]
    
    # Top shell is at satellite altitude. This can be changed if need it. 
    # Top shell will be significanlty bigger than the second to top due to this choice. 
    rtop[0] = Horbit-1
    
    # Define middle bound of each layer
    rmid = (rbot + rtop)/2

    # Build observation matrix
    S = np.zeros((Ntheta, len(rmid)))
    
    # For each observation zenith angle
    for i in range(Ntheta):
        # For each defined altitude shell
        for j in range(len(rmid)):
            
            th = theta[i]
            rb = rbot[j]
            rt = rtop[j]
            sb2 = -math.sin(th)**2  + ((RE+rb)/(RE+Horbit))**2
            st2 = -math.sin(th)**2  + ((RE+rt)/(RE+Horbit))**2  
            
            if sb2 < 0: # there is no intersection of LOS with altitude rb. Set term to 0.
                sb2 = 0.
            if st2 < 0: # there is no intersection of LOS with altitude rt. Set term to 0.
                st2 = 0.
                
            path_len_km = 2*(RE+Horbit) * ( math.sqrt(st2) - math.sqrt(sb2) )
            
            # Go to cm
            path_len_cm = path_len_km*1e5
            # Result should be a conversion from VER to Rayleighs. path_length_cm * 10^-6 to match paper
            S[i,j] = path_len_cm * 1e-6
            
    return S

# Calculation of the Distance matrix assuming Non-Spherical Earth  
def Calculate_D_Matrix_WGS84(satlatlonalt,az,ze):
    '''
    Returns the Distance matrix calculated given the satellite coordinates and viewing geometry for the WGS84 model.
    INPUTS:
        satlatlonalt -  vector containing the satellite coordinates [lat-lon-alt] (km)
        az           -  azimuth vector (deg)
        ze           -  zenith vector (deg)
    OUTPUT:
        S  - Distance matrix assuming WGS84 (km 10^{-1})
    NOTES:
    HISTORY:
        03-Jun-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    
    # Check the size of the zenith vector            
    if np.size(np.shape(ze))==2:
        ze = ze[:,0]
    elif np.size(np.shape(ze))>2:
        raise Exception('Invalid observation zenith angle vector')
    
    # Check the size of the azimuth vector            
    if np.size(np.shape(az))==2:
        az = az[:,0]
    elif np.size(np.shape(az))>2:
        raise Exception('Invalid observation zenith angle vector')
        
    # Create as many shells as the distinct observation zenith angles. 
    rbot = np.zeros(len(ze))
    
    # For each observation zenith and azimuth angles calculate the lower bound tangent altitudes
    rbot = [ic.tangent_point(satlatlonalt,az_i,ze_i)[2] for az_i,ze_i in zip(az,ze)]

    # Calculate the top boundary for each shell
    rtop = rbot.copy()
    rtop[1:] = rbot[:-1]
    
    # Top shell top boundary is defined at satellite altitude.
    rtop[0] = satlatlonalt[2] - 1
    
    # Calculate the mid boundary for every shell
    rmid = (rbot + rtop)/2
    
    # Initiallize the distance matrix 
    S = np.zeros((np.size(ze),np.size(rbot)))

    # For every observation zenith angle
    for i,ze_i in enumerate(ze):
        # For every distinct altitude shell
        for j,rbot_j in enumerate(rbot):
            
            # Calculate the distance from satellite to the top of the shell
            ub = ic.distance_to_shell(satlatlonalt, az[i], ze_i, rtop[j])
            # Calculate the distance from satellite to the bottom of the shell
            lb = ic.distance_to_shell(satlatlonalt, az[i], ze_i, rbot_j)
            
            # If results is NaN for bith  the upper and lower bound then that line-of-sight doesnt go throught the specific shell
            if np.isnan(ub) and np.isnan(lb):
                S[i,j] = 0
            # In case we dont have a lower bound is because we are the top shell and the entry and exit point of the raypath will be the same shell
            elif np.isnan(lb):
                lb = ic.distance_to_shell(satlatlonalt, az[i], ze_i, rtop[j],'second')
                S[i,j] = abs(ub - lb)
            
            #If both upper and lower bounds are defined then the total distance travelled is the distance between the entry and the exit point for that shell. That values is multiplied by 2 since we assume spherically symmetric atmosphere. This means that the distances for the same shells are symmetric in regards to the tangent altitude point
            else:
                S[i,j] = 2*abs(ub - lb)

    # Scaling factor due to the conversion from VER to Rayleigh
    S = S*1e-1
    
    return S
    
'''

TIKHONOV REGULARIZATION

'''

# Tikhonov Regularization
def Tikhonov(A, Bright, reg_deg, reg_param=0., Sig_Bright=None, weight_resid=False):
    '''
    Tikhonov Regularization function. Solves the minimization problem.
    The regularization parameter is calculated using the L-Curve and finding the maximum curvature point. 
    INPUTS:
        A            - Distance matrix
        Bright       - Brightness Profile (Rayleigh)
        reg_deg      - The degree of Tikhonov regularization [valid is 0,1,2]
        reg_param    - regularization parameter. [Input 0: calls create_alpha_values to create the reguraliation
                       parameter vector. Input, any scalar or matrix: perform regularization given that input]
        Sig_Bright   - Covariance matrix of Bright. If None, ver is the only output. If specified, it is propagated
                       through, and both ver and Sig_ver are outputs.
        weight_resid - Whether to weight the data by its covariance. If this is True, Sig_Bright must be specified.
    OUTPUT:
        ver          - Estimated Volume emission rate [ph/cm^{-3}/s]
        Sig_ver      - Covariance matrix of ver (only outputted if Sig_Bright is not None)
    NOTES:
    HISTORY:
        15-May-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        20-Apr-2016: Uncertainty propagation added by Brian Harding (bhardin2@illinois.edu)
    CALLS:
        -create_alpha_values
        -calc_solution
        -Maximum_Curvature_gradiens
    '''
    
    # Make copy of inputs so we don't accidentally overwrite them
    A = A.copy()
    Bright = Bright.copy()

            
    # Check if weight inputs make sense
    if weight_resid and Sig_Bright is None:
        raise Exception('weight_resid==True does not make sense if Sig_Bright is None. Either specify Sig_Bright or ' +\
                        'set weight_resid==False')
        
    # If we want to use a weighted least-squares formulation, we can "whiten" the errors and use
    # the exact same code.
    if weight_resid:
        W = sqrtm(np.linalg.inv(Sig_Bright)) # whitening matrix
        A = W.dot(A)
        Bright = W.dot(Bright)
        Sig_Bright = np.eye(len(Bright)) # whitening diagonalizes the covariance matrix
        

    # Check if reg_param is a vector or not
    if np.size(reg_param) == 1:
        # if its not a vector only if zero it will calculate the vector. Otherwise it will use the single value given. 
        if reg_param == 0:
            reg_param = create_alpha_values(A)
    
    residual = np.zeros(len(reg_param))
    seminorm = np.zeros(len(reg_param))
    
    # create the matrix that defines the order of the regularization.
    L = get_rough_matrix(len(Bright),reg_deg)

    # For every regularization parameter, estimate solution
    for i in range(0,len(reg_param)):
        sol = calc_solution(A,Bright,reg_param[i],L) # for speed, we can omit the uncertainty prop here 
        r = A.dot(sol) - Bright
        residual[i] = np.linalg.norm(r)
        seminorm[i] = np.linalg.norm(L.dot(sol))
        
    # Find the optimal regularization parameter using the maximum second derivative method
    reg_corner = Maximum_Curvature_gradiens(residual,seminorm,reg_param)
    
    # Calculate the solution with the optimal parameter (and, if desired, also the uncertainty)
    if Sig_Bright is None:
        ver = calc_solution(A,Bright,reg_corner,L)
        return ver
    else:
        # Note that if weight_resid was True, we already whitened the errors so don't need to weight
        # again in the sub-function
        ver,Sig_ver = calc_solution(A,Bright,reg_corner,L,Sig_Bright=Sig_Bright,weight_resid=False)
        return ver,Sig_ver
    

# Calculate regularized inverse rolution
def calc_solution(A,Bright,alpha,L,Sig_Bright=None,weight_resid=False):
    '''
    Calculates the solution using non-negative least square for the regularization problem.
    Optionally propagates uncertainty from input brightness to output VER.
    INPUTS:
        A          - Distance matrix
        Bright     - Brightness Profile (Rayleigh)
        lamda      - The selected regularization parameter
        L          - Roughening matrix
    OPTIONAL INPUTS:
        Sig_Bright   - Covariance matrix of Bright. If not None, the uncertainty propagation calculation will be performed, and
                       output uncertainty will be returned. Default is None, for speed.
        weight_resid - Whether to weight the data by its covariance. If this is True, Sig_Bright must be specified.
    OUTPUT:
        ver        - Estimated electron density (cm^-3)
    OPTIONAL OUTPUT:
        Sig_ver    - Covariance matrix of ver.
                     
    NOTES:
    HISTORY:
        15-May-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        20-Apr-2016: Uncertainty propagation added by Brian Harding (bhardin2@illinois.edu)
    CALLS:
    '''
    
    # Check if weight inputs make sense
    if weight_resid and Sig_Bright is None:
        raise Exception('weight_resid==True does not make sense if Sig_Bright is None. Either specify Sig_Bright or ' +\
                        'set weight_resid==False')
        
    # If we want to use a weighted least-squares formulation, we can "whiten" the errors and use
    # the exact same code.
    if weight_resid:
        W = sqrtm(np.linalg.inv(Sig_Bright)) # whitening matrix
        A = W.dot(A)
        Bright = W.dot(Bright)
        Sig_Bright = np.eye(len(Bright)) # whitening diagonalizes the covariance matrix
    
    # Create the augmented matrix C
    C = np.concatenate((A, alpha*L))
    
    # Check of the size of the brightness is correct
    if np.size(np.shape(Bright))==2:
        Bright = Bright[:,0]
    elif np.size(np.shape(Bright))>2:
        raise Exception('Invalid data vector')
        
    # Constuct vector d = [Bright;0]
    d_temp = np.zeros(len(L))
    d = np.concatenate((Bright, d_temp))

    # Solve the non-negative constrained least-squares
    ver,_ = scipy.optimize.nnls(C,d)
    
    # Solve using normal least squares
    #ver = np.linalg.lstsq(C,d)[0]
    
    if Sig_Bright is None:
        return ver
    else:
        # Use linear error propagation, ignoring the non-negativity constrant.
        R = np.linalg.inv(A.T.dot(A) + alpha**2 * L.T.dot(L))
        Sig_ver = R.dot(A.T).dot(Sig_Bright).dot(A).dot(R.T)
        return ver, Sig_ver
    
# Calculate the optimal regularization parameter
def Maximum_Curvature_gradiens(residual,x_lamda,reg_param,method='derivative'):
    '''
    Given the residual and the seminorm of the L2 norm problem the corner of the L-Curve is calculated by
    finding the maximum curvature by finding the maximum second derivative
    INPUTS:
        residual   - The residual norm of our minimization problem
        x_lamda    - The seminorm of ||Lx||
        reg_param  - Vector containing all the regularization parameters
        method     - either 'curvature' or 'derivative':
                        'curvature': find the point in the L-curve with the minimum radius of curvature
                        'derivative': find the point in the L-curve with the maximum 2nd derivative
    OUTPUT:
        reg_corner - Optimal regularization parameter 
    NOTES:
    HISTORY:
        17-Sep-2015: Written by Dimitrios Iliou (iliou2@illinois.edu) and Jianqi Qin (jianqi@illinois.edu)
        26-Apr-2017: Added option for maximum curvature (instead of maximum 2nd derivative) Brian Harding (bhardin2@illinois.edu) 
    CALLS:
    '''
    
    #transform rho and eta into log-log space
    Xvec=np.log(residual);
    Yvec=np.log(x_lamda);
    
    if method == 'derivative': # Compute second derivative at each point
        grad1 =  np.gradient(Yvec)/np.gradient(Xvec)
        grad2 = np.gradient(grad1)/np.gradient(Xvec)
    
    elif method == 'curvature': # Compute radius of curvature at each point
        grad2 = np.zeros(len(reg_param))
        for i in range(1,len(reg_param)-1):
            da = reg_param[i+1]-reg_param[i-1]
            dxda = (Xvec[i+1]-Xvec[i-1])/da
            dyda = (Yvec[i+1]-Yvec[i-1])/da
            da2 = ((reg_param[i+1]-reg_param[i-1])/2)**2
            d2yda2 = (Yvec[i-1] - 2*Yvec[i] + Yvec[i+1])/da2
            d2xda2 = (Xvec[i-1] - 2*Xvec[i] + Xvec[i+1])/da2
            grad2[i] = (dxda*d2yda2 - dyda*d2xda2)/(dxda**2 + dyda**2)**(1.5)
    else:
        raise Exception('Input "method" should be "curvature" or "derivative"')
    
    ireg_corner=np.argmax(grad2)
    reg_corner=reg_param[ireg_corner]
            
    return reg_corner

# Create regularization parameter vector
def create_alpha_values(A,npoints = 100):
    '''
    Given a distance matrix A a number of points the regularization parameter vector is created. 
    INPUTS:
        A       -  Distance matrix
        npoints -  Number of distinct values for the regularization parameter
    OUTPUT:
        reg_param - Vector containing the regularization parameters
    NOTES:
    HISTORY:
        17-Sep-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:
    '''
    # SVD Decomposition of matrix A (Distance matrix)
    U, s, V = np.linalg.svd(A, full_matrices=True)
    
    # multiplication ratio
    smin_ratio = 16*np.spacing(1)
    #smin_ratio = 16*np.finfo(np.float32).eps
    reg_param = np.zeros(npoints)
    reg_param[npoints-1] = max([s[np.size(s,0)-1],s[1]*smin_ratio])
    ratio = (s[0]*100/reg_param[npoints-1])**(1./(npoints-1));
    
    # Put regularization parameters in descending order
    for i in np.arange(npoints-2,-1,-1):
        reg_param[i] = ratio*reg_param[i+1]
        
    return reg_param

# Create matrix to define regularization order 
def get_rough_matrix(n,deg):
    '''
    Creates the matrix needed for higher order Tikhonov regularization 
    INPUTS:
        n   - Size of the matrix [Depends on the size of the Distance matrix]
        deg - The degree of Tikhonov regularization [valid is 0,1,2]
    OUTPUT:
        L   - Roughening matrix
    NOTES:
    HISTORY:
        15-May-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:
    '''
    
    if deg==0:
        L = np.eye(n)
        return L
    elif deg==1:
            # First Order
        L = np.zeros((n-1,n))
        for i in range(0,n-1):
            L[i][i] = -1
            if i<(n-1):
                L[i][i+1] = 1
        return L
    elif deg==2:
        L = np.zeros((n-2,n))
        for i in range(0,n-2):
            L[i][i] = 1
            if i<(n-1):
                L[i][i+1] = -2
            if i<(n-2):
                L[i][i+2] = 1
        return L
    else:
        raise Exception('Invalid degree. Degree can be 0,1 or 2')


'''

MAP ESTIMATION

'''
# MAP ESTIMATION FUNCTION
def MAP_Estimation(S,Bright,VER_mean=0.,VER_profiles=0.,weight_resid=False):
    '''
    Bayesian MAP Estimation Method for regularization. This method takes as input the prior distribution mean and ver profiles
    INPUTS:
        S            - Distance matrix
        Bright       - Brightness Profile (Rayleigh)  
        VER_mean     - Prior mean
        VER_profiles - Prior variance matrix  
        weight_resid - Whether to weight the data
    OUTPUT:
        VER          - Estimated VER (ph/cm^-3/s)
    NOTES:
        ==> VER Prior is loaded for default path. This can be changed at will inside the function
    HISTORY:
        18-Apr-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:

    '''
    # Check whether or not the prior is there. If not raise exception
    if not hasattr(VER_profiles, "__len__"):
        if (VER_profiles==0):
            print 'PRIOR IS MISSING'
            raise Exception
            
    if weight_resid:
        raise Exception('Weighted residual not yet implemented for MAP_Estimation')
    
    # Interpolate brightness to remove zeros and calculate the covariance matrix
    # This is needed so the covariance matrix for the data is not singular
    x = range(0,len(Bright))
    xp = np.where(Bright!=0)[0]
    Bright_interp= np.interp(x, xp, Bright[xp])
    
    # Get instrument parameters
    instr_params = get_FUV_instrument_constants()
    
    # Get sensitivity of the instrument and exposure time
    sensitivity =  instr_params['sensitivity']
    exposure = instr_params['exposure']
    
    # To calculate the variance of the noise we just take a single Brightness profile and we find the variance. More accurate estimation can be performed if we sum the 6 stripes before feeding them here. We also need to have the sensitivity of the instrument to calculate the number of counts. 
    W_var = Bright_interp * sensitivity * exposure
    # Create Diagonal matrix of the VER variances.
    W_var = W_var/((exposure*sensitivity)**2)
    
    # Create Diagonal matrix of the Brightness variances.
    Rw = np.diag(W_var)

    # Prior VER profiles are created for the full CCD
    # Truncate them to the limb altitudes
    VER_profiles = VER_profiles[:len(Bright),:]
    VER_mean =  VER_mean[:len(Bright)]
    
    # Calculate the covariance matrix 
    Rx_ver = np.cov(VER_profiles);

    # Calculate the square root of the data and prior covariance matrices
    covd12=sqrtm(np.linalg.inv(Rw))
    covm12=sqrtm(np.linalg.inv(Rx_ver))

    # Create the augmented matrices to solve the least-squares problem 
    A = np.concatenate((covd12.dot(S), covm12), axis=0)
    rhs= np.concatenate((covd12.dot(Bright), covm12.dot(VER_mean)), axis=0)

    # Solve the least-squares problem
    VER,_,_,_ = np.linalg.lstsq(A,rhs)
    
    # Check if the values of the VER are negative. 
    for i in range(0,np.size(VER)):
            if VER[i]<0:
                VER[i] = 0
    
    # This is the non-negative least-squares solver
    #VER,_ = scipy.optimize.nnls(A,rhs)
    
    return VER





'''

CALCULATE ELECTRON DENSITY

'''

# Calculate the electron density given VER
def calculate_electron_density(VER,satlatlonalt,tang_altitude,dt,Sig_VER=None,contribution='RR'):
    '''
    Given a VER profile the electron density is calculated. The physical process for the VER must be defined.
    INPUTS:
        VER             - Volume Emission Rate profile (ph/cm^-3)
        satlatlonalt    - satellite coordinates [latitude(degrees),longitude(degrees),altitude(km)]
        tang_altitude   - tangent altitude [km]
        dt              - datetime (python datetime)
        Sig_VER         - Covariance matrix of VER. If None, Ne will be the only output. [(ph/cm^-3)**2]
        contribution    - contributions to the 135.6nm emission, RR: Radiative Recombination, RRMN: Radiative Recombination and Mutual Neutralization 
    OUTPUTS:
        Ne              - The estimated Electron Density profile (1/cm^3/s)
        Sig_Ne          - (OPTIONAL) covariance matrix of Ne. If Sig_VER is None, this is not returned.
    NOTES:
        For the calculation, Ne is approximately equal O+. This is assumption is valid for the altitudes that the FUV limb measurement are looking.
    TODO:
        (BJH) For MN, we should use the location of the tangent point, not the satellite location (though it probably won't matter much).
    HISTORY:
        17-Sep-2015: Written by Dimitrios Iliou (iliou2@illinois.edu) 
        24-Apr-2017: Uncertainty propagation added by Brian Harding (bhardin2@illinois.edu)
    CALLS:
        Pyglow
    
    '''
    
    # Check the size of the VER vector            
    if np.size(np.shape(VER))==2:
        VER = VER[:,0]
    elif np.size(np.shape(VER))>2:
        raise Exception('Invalid data vector')

    # Check if uncertainties should be returned
    ret_cov = True
    if Sig_VER is None:
        Sig_VER = np.eye(len(VER)) # dummy variable
        ret_cov = False

    b1356 = 0.54    # yield parameter (unitless)
    a1356 = 7.3e-13 # radiative recombination rate (cm^3/s)
    
    # consider that to be constant through altitude, normally it should change with temperature
    k1 = 1.3e-15    # radiative attachment rate (cm^3/s)
    k2 = 1e-7       # ion-ion neutralization rate (cm^3/s)
    k3 = 1.4e-10    # ion-atom neutralization rate (cm^3/2)
    
    if contribution=='RRMN':
        O = np.zeros(len(VER))
        for i,height in enumerate(tang_altitude):
            
            # Create pyglow point
            pt = pyglow.pyglow.Point(dt, satlatlonalt[0], satlatlonalt[1], height)
            # Run MSIS-00 to get O density
            pt.run_msis()
    
            # Store the Oxygen density for every line of sight - tangent altitude
            O[i] = pt.nn['O']  
        
        # Solve the third order equation to calculate the electron density
        a0 = a1356/O
        b0 = a1356*k3/k2+b1356*k1
        c0 = -VER/O
        d0 = -VER*k3/k2

        a1 = b0*c0/(6.0*a0**2)-b0**3./a0**3/27.0-d0/a0/2.0
        
        b1 = c0/a0/3.0-b0**2./a0**2/9.0

        c1 = -b0/a0/3.0;
        
        d1 = a1/np.sqrt((-b1)**3)
        
        # Numerical correction
        d1[np.where(d1<-1.)]=-1.
        d1[np.where(d1>1.)]=1.
            
        # used arccos instead of MATLAB acos
        Ne = c1+2.0*np.sqrt(-b1)*np.cos(np.arccos(d1)/3.0)
        
        # Propagation of uncertainty through this process
        # TODO: use the correct formulation for RRMN
        # CURRENT IMPLEMENTATION: copy of RR implementation
        dVER = np.sqrt(np.diag(Sig_VER))
        dNe_dVER = 1/np.sqrt(a1356) * (np.sqrt(VER+dVER) - np.sqrt(VER))/dVER
        J = np.diag(dNe_dVER) # Ne at altitude i only depends on VER at altitude i
        Sig_Ne = J.dot(Sig_VER).dot(J.T)
        
    elif contribution == 'RR': 
       
        # VER to ne conversion
        Ne = np.sqrt(VER/a1356)
        
        # Propagate uncertainty through this conversion, using the Jacobian formula:
        # if X is a random vector, and Y=f(X), then SigY = J.dot(SigX).dot(J.T) where
        # J is the Jacobian of f.
        
        # The only problem is that when VER is really small, the Jacobian grows to infinity.
        # For these samples, simple linearized error propagation is not valid. Instead of
        # using the analytical dy/dx, we will use dy/dx estimated from the secant method.
        # This gives a better indication of the true confidence intervals.
        dVER = np.sqrt(np.diag(Sig_VER))
        dNe_dVER = 1/np.sqrt(a1356) * (np.sqrt(VER+dVER) - np.sqrt(VER))/dVER
        J = np.diag(dNe_dVER) # Ne at altitude i only depends on VER at altitude i
        Sig_Ne = J.dot(Sig_VER).dot(J.T)
        
    else:
        raise Exception('Invalid input for the process selection')

    #Ne[np.isnan(Ne)] = 0
    
    if ret_cov:
        return Ne,Sig_Ne
    else:
        return Ne
    
    
'''

EXTRACT F2 PEAK DENSITY AND HEIGHT

'''
def find_hm_Nm_F2(NE,alt_vector,interpolate_F2=True,interval_increase = 20.,kind ='cubic', Sig_NE=None, Nmc=100):
    '''
    Calculates the NmF2 and hmF2 of a given electron density profile
    INPUTS:
        NE             - Electron density profile [cm^{-3}]
        alt_vector     - altitude vector [km]
        interpolate_F2 - Flag indicating if we interpolation in the F2-peak region will be performed 
        interval_increase  - Number indicating the amount of points to be created. => len(hmf2Region)* int_increase
        kind           - Goes as input to the interp1d function and declares the type of interpolation
        Sig_NE         - covariance matrix of NE. If this is not None, then the sigma_hmF2 and sigma_NmF2 will be outputted.
        Nmc            - number of Monte Carlo trials to calculate uncertainty (Only used if Sig_NE is not None).
    OUTPUT:
        hmF2       - altitude of the maximum intesity value of the Ne profile
        NmF2       - peak intensity value of the altitude profile
        sigma_hmF2 - (OPTIONAL) propagated uncertainty of hmF2 (only outputted if Sig_NE is not None)
        sigma_NmF2 - (OPTIONAL) propagated uncertainty of NmF2 (only outputted if Sig_NE is not None)
    NOTES:

    HISTORY:
        06-Jun-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        16-Jun-2016: Added F2-peak region interpolation option
        25-Apr-2017: Uncertainty propagation (Monte Carlo) added by Brian Harding (bhardin2@illinois.edu)
    '''
    
    # Save copy of the inputs in case interpolate_F2 is True.
    NE_orig = NE.copy()
    alt_vector_orig = alt_vector.copy()
    
    # If flagged do interpolation in the F2-peak region
    if interpolate_F2:
        NE, alt_vector = hmF2_region_interpolate(NE,alt_vector,interval_increase,kind)
        
    # Locate the index that has the maximum electron density values
    indexf2= NE.argmax()

    # Return the value in the altitude vector that corresponds to the above index
    hmF2 = alt_vector[indexf2]
        
    # Return the value in the electron density vector that corresponds to the above index
    NmF2 = NE[indexf2]

    
    # If desired, propagate uncertainty using bootstrapped Monte Carlo
    if Sig_NE is not None:
        # For this we only need to consider the region near hmF2. By restricting the domain:
        #     - it runs faster
        #     - it avoids problems with huge uncertainties on bottomside
        idx = np.where(abs(alt_vector_orig-hmF2) < 100)[0] # focus on +/- 100 km from peak
        NE2 = NE_orig[idx]
        alt2 = alt_vector_orig[idx]
        Sig_NE2 = Sig_NE[np.ix_(idx,idx)] # awkward indexing to get covariance-matrix of subset
        hmF2_vec = np.zeros(Nmc)
        NmF2_vec = np.zeros(Nmc)
        for n in range(Nmc):
            # Generate a new random profile by sampling from a multivariate Gaussian
            NE_n = multivariate_normal(NE2, Sig_NE2)
            # Calculate and store the hmF2 and NmF2 of this random profile.
            h,N = find_hm_Nm_F2(NE_n, alt2, interpolate_F2, interval_increase, kind, Sig_NE=None)
            hmF2_vec[n] = h
            NmF2_vec[n] = N
        
        sigma_hmF2 = np.std(hmF2_vec)
        sigma_NmF2 = np.std(NmF2_vec)
        
        return hmF2, NmF2, sigma_hmF2, sigma_NmF2
    
    else:
        return hmF2, NmF2

# Perform interpolation in the F2-region
def hmF2_region_interpolate(Ne,alt_vector,interval_increase = 20.,kind ='cubic'):
    '''
    Interpolates the values around the hmF2 region in order to increase altitude vector interval which helps in minimizing the hmF2 error
    INPUTS:
        Ne                 - Original Electron density profile [cm^{-3}]
        alt_vector         - altitude vector for original electron density profile [km]
        interval_increase  - Number indicating the amount of points to be created. => len(hmf2Region)* int_increase
        kind               - Goes as input to the interp1d function and declares the type of interpolation
    OUTPUT:
        NE_inter           -  Interpolated altitude profile for the estimated electron density
        alt_vector_interp  -  Interpolated altitude vector for the true electron density
    NOTES:
    HISTORY:
        29-Apr-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    
    # Locate the F2-peak index
    nm_true = np.argmax(Ne)
    
    # Check that the profile is in altitudes that make sense.
    # Noise can create maximum values that are located at the top or the bottom of the profile
    # In that case interpolation won't be applied
    if nm_true>5 and nm_true<len(Ne)-5:
        
        # Create the region of electron density values to be interpolated
        # It will be +- 5 indeces from the F2-peak index
        y_true = Ne[nm_true-5:nm_true+5]
        # Create the corresponding altitude vector
        x_true = alt_vector[nm_true-5:nm_true+5]

        # Create the interplator
        f2 = interpolate.interp1d(x_true, y_true, kind)
        
        # Create the new altitude vector by determining how many new points we want in the F2-peak region interval
        alt_vector_interp = np.linspace(x_true[0], x_true[-1], num=len(x_true)*interval_increase, endpoint=True)
        
        # Calculate the interpolated electron density values
        NE_inter = f2(alt_vector_interp)
        
        return NE_inter,alt_vector_interp
    else:
        print 'Ne maximum ill defined - No interpolation performed'
        # In this case the values remain the same. No exception is required since the F2-peak code can run
        return Ne,alt_vector
  
        
'''

L2 Calculation top-level function

'''

def FUV_Level_2_Density_Calculation(Bright,alt_vector,satlatlonalt,az,ze, Sig_Bright = None, weight_resid = False, limb = 150.,Spherical = True,reg_method='Tikhonov',regu_order = 2,reg_param=0,contribution='RRMN', VER_mean=0.,VER_profiles=0,dn = datetime.datetime(2012, 9, 26, 20, 0, 0)):                                  
    '''
    Within this function, given data from LVL1 Input file, the VER and Ne profiles for the tangent altitude point are
    calculated
    INPUTS:
        Bright      - Brightness Profile [Rayleigh]
        alt_vector  - Tangent altitudes vector [km]
        satlatlonalt- Satellite Coordinates [degrees,degrees,km]
        az          - Azimuth [degrees]
        ze          - Zenith [degrees]
        Sig_Bright  - Covariance matrix of Bright [Rayleigh**2]. If not None, the covariance of VER and Ne are also returned.
        weight_resid- Whether to weight the data by the uncertainties. Passed to Tikhonov or MAP_Estimation function.
        limb        - Defines the lower bound that defines the limb [km]
        Spherical   - Flag indicating whether Sperical or Non-Spherical Earth is assumed [True, False]
        reg_method  - Regularization method selection [Tikhonov, MAP]
        regu_order  - Regularization Order [int] (0,1,2 Possible Values)
        reg_param   - Regularization parameter [0: regularization parameter vector is calculated for the corresponding distance matrix; any scalar or vector: Use this value instead for the regularization pamaremeter]
        contribution- contributions to the 135.6nm emission, RR: Radiative Recombination, RRMN: Radiative Recombination and Mutual Neutralization 
        VER_mean    - Prior mean [len(altitudes)](ph/cm^{-3}/s)
        VER_profiles- Prior ver profiles [len(altitudes)xnum_of_profiles](ph/cm^{-3}/s) 
        dn          - Datetime object of the measurement; is needed to calculate the [O] density from MSIS to calculate the MN effect
    OUTPUT:
        VER         - Volume Emission Rate tangent altitude profile (ph/cm^{-3}/s)
        Ne          - Electron Density tangent altitude profile (cm^{-3})
        h           - Truncated tangent altitude vector [km]
        Sig_VER     - (OPTIONAL) Returned if Sig_Bright is not None. Covariance matrix of VER.
        Sig_Bright  - (OPTIONAL) Returned if Sig_Bright is not None. Covariance matrix of Ne.
    NOTES:

    HISTORY:
        16-Jun-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)
        24-Apr-2017: Uncertainty propagation added by Brian Harding (bhardin2@illinois.edu)
    '''
    
    
    # Determine if covariances should be returned
    ret_cov = True
    if Sig_Bright is None:
        ret_cov = False
        Sig_Bright = np.eye(len(Bright)) # dummy variable
    
    # Truncate the tangent altitude vector to include only limb measurements
    h = alt_vector[np.where(alt_vector>limb)]
    Bright= Bright[0:len(h)]
    Sig_Bright = Sig_Bright[:len(h),:len(h)]
    ze = ze[:len(h)]
    az = az[:len(h)]
    
            
    # Check of we are estimating the distance matrix for a spherical or ellipsoid earth
    # Spherical earth calculations need ~0.26 secs
    # Non-Spherical earth calculations need ~240 secs. Using Non-Spherical earth will create a bottleneck for the data pipeline
    if Spherical:
        S = create_cells_Matrix_spherical_symmetry(ze,satlatlonalt[2])
    else:
        S = Calculate_D_Matrix_WGS84(satlatlonalt,az,ze)    
  

    if reg_method=='Tikhonov':
        VER, Sig_VER = Tikhonov(S,Bright,regu_order,reg_param=reg_param,Sig_Bright=Sig_Bright,weight_resid=weight_resid)
    elif reg_method =='MAP':
        # BJH: I don't think MAP will be used in practice, but if so, uncertainty propagation
        # will need to be added before we can use it. 
        raise Exception('Uncertainty propagation not yet implemented for MAP Estimation.')
        VER = MAP_Estimation(S,Bright,VER_mean=VER_mean,VER_profiles=VER_profiles,weight_resid=weight_resid)
    else:
        raise Exception('Incorrect regularization method chosen. Choices are: Tikhonov or MAP')
        
    Ne, Sig_Ne = calculate_electron_density(VER=VER, satlatlonalt=satlatlonalt, tang_altitude=h, dt=dn, Sig_VER=Sig_VER, contribution=contribution)         
   
    if ret_cov:
        return VER, Ne, h, Sig_VER, Sig_Ne
    else:
        return VER,Ne,h
    
'''

CREATE NETCDF FILE

'''
def FUV_Level_2_OutputProduct_NetCDF(dn,satlatlonalt,az,ze,tanlatlonalt,Bright,VER,Ne,NmF2,hmF2,path='/home/dimitris/public_html/Datafiles/LVL2TEST/',description =  '.NOMINAL' ,year = str(datetime.datetime.now().year) ,doy = datetime.datetime.now().timetuple().tm_yday ,version = 'v0001'):
    '''
    This function takes as input all the necessary outputs from LVL2.5 processing and writes them on a NetCDF file
    INPUTS:
        dn          - Universal date and time
        satlatlonalt- Satellite Coordinates [degrees,degrees,km]
        az          - Azimuth [degrees]
        ze          - Zenith [degrees]
        tatlatlonalt- Tangent point coordinates Coordinates [degrees,degrees,km]
        Bright      - Brightness Profile [Rayleigh]
        VER         - Volume Emission Rate profile at the tangent altitude [ph/cm^3/sec]
        NE          - Electron Density profile at the tangent altitude [cm^(-3)]
        NmF2        - NmF2 peak [cm^(-3)]
        hmF2        - hmF2 peak [km]
        path        - Pathname where the NetCDF file is saved
    OUTPUT:
        Creates a NetCDF file on the desired path
    NOTES:
        
    HISTORY:
        24-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    # Open the NetCDF file 
    TOD = '.' + datetime.datetime.strftime(datetime.datetime.utcnow(), '%H.%M.%S') 
    fn = 'ICON.L2_5.FUV%s.%s.%s%s.%s.nc' %(description,year,doy,TOD,version)
    datagrp = netCDF4.Dataset(path+fn, 'w', format='NETCDF4')
    
    
    ## The details for the creation of the NetCDF file can be found in the corresponding file by Tory Fae
    # Global atributes
    # 3.1.1 
    datagrp.setncattr('Acknowledgement','[TBD] Awaiting input from PI')

    # 3.1.2 
    datagrp.setncattr('ADID_Ref','NASA Contract > NNG12FA45C')

    # 3.1.3 - This should be empty if no calibration files were used
    cal_file_name = 'ICON.L1.EUV.Calibration.v0001.NC'
    datagrp.setncattr('Calibration_File',cal_file_name)

    # 3.1.4
    datagrp.setncattr('Conventions','SPDF ISTP/IACG Modified for NetCDF')

    # 3.1.5
    datagrp.setncattr('Data_Level','L2.5')

    # 3.1.6
    datagrp.setncattr('Data_Type','DP25 > Data Product 2.5:FUV Nighttime O+ profile')

    # 3.1.7
    datagrp.setncattr('Data_Version','1')

    # 3.1.8 

    #- Date_end: Awaiting input file
    datagrp.setncattr('Date_End','Awaiting Input File')

    # Date_start - Date_end: Awaiting input file
    datagrp.setncattr('Date_Start','Awaiting Input File')

    # File_Date
    File_date = datetime.date.strftime(datetime.datetime.utcnow(),'%a, %d %b %Y, ') + datetime.datetime.strftime(datetime.datetime.utcnow(), '%Y-%m-%dT%H:%M:%S.%f')[:-3] +' UTC'
    datagrp.setncattr('File_date',File_date)

    # Generation Date
    Generation_date = datetime.date.strftime(datetime.datetime.now(),'%Y%m%d')
    datagrp.setncattr('Generation_Date',Generation_date)

    # 3.1.9
    datagrp.setncattr('Description','ICON Level 2 Nominal FUV Daily Downlink')

    # 3.1.10
    datagrp.setncattr('Descriptor','FUV > Intensified Far Ultraviolet Imager')

    # 3.1.11
    datagrp.setncattr('Discipline','Space Physics > Ionospheric Science')

    # 3.1.12
    #doy = datetime.datetime.utcnow().timetuple().tm_yday
    #year = str(datetime.datetime.utcnow().year)
    #version = 'v0001'

    # Optional descriptions - Need fullstop at the beginning
    #description = '.NOMINAL' 
    #TOD = '.' + datetime.datetime.strftime(datetime.datetime.utcnow(), '%H.%M.%S') 

    fn = 'ICON.L2_5.FUV%s.%s.%s%s.%s.nc' %(description,year,doy,TOD,version)
    datagrp.setncattr('File',fn)

    # 3.1.13 - Need to specify processor here
    datagrp.setncattr('Generated_By','ICON SDC > ICON UOI FUV Processor v1.0, D. Iliou')

    # 3.1.14
    File_date_2 = datetime.datetime.strftime(datetime.datetime.utcnow(), '%Y-%m-%dT%H:%M:%S')
    datagrp.setncattr('History','Version 1, D. Iliou, '+ File_date_2 +', FUV L2 Processor v1')
    datagrp.setncattr('MODS','Version 1, D. Iliou, '+ File_date_2 +', FUV L2 Processor v1')

    # 3.1.15
    datagrp.setncattr('Instrument','FUV')

    # 3.1.16
    datagrp.setncattr('Instrument_Type','Imagers (space)')

    # 3.1.17
    datagrp.setncattr('HTTP_LINK','http://icon.ssl.berkeley.edu/FUV')
    datagrp.setncattr('Link_Text','12-second FUV SW altitude profiles')
    datagrp.setncattr('Link_Title','ICON FUV')

    # 3.1.18
    fn2 = 'ICON.L2_5.FUV.SWP.%s.%s.%s' %(year,doy,version)
    datagrp.setncattr('Logical_File_ID',fn2)

    datagrp.setncattr('Logical_Source','ICON.L2_5.FUV.SWP.')

    datagrp.setncattr('Logical_Source_Description','FUV Short Wavelength Channel - 135.6 Altitude Profiles')

    # 3.1.19 
    datagrp.setncattr('Mission_Group','Ionospheric Investigations')

    # 3.1.20
    parent_file = 'ICON.L1.FUV.SWP.2015.017.v0001' # This is an example. It should change dynamically
    datagrp.setncattr('Parents','NC > '+parent_file)

    # 3.1.21
    datagrp.setncattr('PI_Affiliation','UC Berkeley > SSL')

    # 3.1.22
    datagrp.setncattr('PI_Name','T. J. Immel')

    # 3.1.23
    datagrp.setncattr('Project','NASA > ICON')

    # 3.1.24
    datagrp.setncattr('Rules_of_Use','Public Data for Scientific Use')

    # 3.1.25
    datagrp.setncattr('Software_Version','ICON SDC > FUV Profile Generator v1.0.0')

    # 3.1.26
    datagrp.setncattr('Source_Name','ICON > Ionospheric Connection Explorer')

    # 3.1.27
    datagrp.setncattr('Spacecraft_ID','NASA > ICON - 493')

    # 3.1.28
    datagrp.setncattr('Text','ICON explores the boundary between Earth and space – the ionosphere – to understand the physical connection between our world and the immediate space environment around us. Visit ‘http://icon.ssl.berkeley.edu’ for more details.')

    # 3.1.29 - This Should change when we have an average runtime on the server
    datagrp.setncattr('Time_Resolution','7 seconds')

    # 3.1.30
    datagrp.setncattr('Title','ICON FUV O+ Altitude Profiles (DP 2.5)')

    
    
    # ------------- # 
    # 3.2 -- Need to add epoch

    Epoch_dim = datagrp.createDimension('Epoch')
    Epoch = datagrp.createVariable('UTC_EPOCH','i8',('Epoch',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=-999999999999999L )

    Epoch.CatDesc = 'Milliseconds since 1970-01-01 00:00:00 UTC at middle of image integration.'
    Epoch.Long_Name = 'Milliseconds since 1970-01-01 00:00:00 UTC at middle of image integration.'

    Epoch.Depend_0 = 'Epoch' # This needs fixing in regards to the epoch

    Epoch.Display_Type = 'Time_Series'

    Epoch.FieldNam = 'UTC Time'

    Epoch.Format = 'I8'

    Epoch.LablAxis = 'UTC Time'

    Epoch.MonoTon = 'Increase' # ISTP optional

    Epoch.ScaleTyp = 'Linear'

    Epoch.Time_Base = '1970-01-01 00:00:00.000 UTC'
    Epoch.Time_Scale = 'UTC'
    # 3.3.3.10
    Epoch.Units = 'milliseconds'

    # 3.3.3.11 
    Epoch.ValidMax = 6000000000000L
    Epoch.ValidMin = 0L
    Epoch.Valid_Max = 6000000000000L
    Epoch.Valid_Min = 0L
    Epoch.Valid_Range = [0L,6000000000000L]

    # 3.3.3.12 
    Epoch.Var_Notes = ["Milliseconds since 1970-01-01 00:00:00 UTC at middle of image integration.","Derived from original GPS values reported from spacecraft (Time_GPS_Seconds and Time_GPS_Subseconds).", "Time calculation is offset by 615ms (flush time) for the first image in the series and for all other images are adjusted by subtracting (integration time + 308 milliseconds) from the reported GPS time then adding the difference between the readout FRT and the header FRT.","Time may be delayed by up to 10 ms due to FSW polling delay.", "Maximum time is ~2150 UTC and minimum time is ~1970 UTC."]
    
    Epoch.Var_Type = "Support_Data" 

    Epoch[:] = calendar.timegm(time.gmtime())
    
    # Create Dimensions for FUV 2.5 output values

    # Single sized variable for the satellite location
    sat_location = datagrp.createDimension("sat_location", 1)

    # Size of the altitude vector 
    alt_profile = datagrp.createDimension("alt_profile", len(Ne))

    # Size of the ionosphere parameters
    F_peak_params = datagrp.createDimension("F_peak_params", 1)

    # Size of the Datetime - This is a string
    date_time = datagrp.createDimension("date_time", len(dn))

    # Datetime 
    # This might be more complicated than just one variable. Look at MIGHTI example. 
    # Datetime can be string, UTC, GPS, end, start time. Since no input provided for now only one is used and the rest can be 
    # written in future versions
    ut_dt = datagrp.createVariable('ICON_L2_FUV_Time_UTC','S1',('date_time',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = ut_dt, catdesc='ISO 9601 formatted UTC timestamp',
                            depend = 'Epoch',displayType = 'Time_Series', Fieldname = 'UTC Time',
                            form = 'A27', labelAxis = 'UTC Time', scaletype = 'Linear', units = '',
                            minmax = ['1970-01-01 00:00:00 UTC','2150-01-01 00:00:00 UTC'],
                            Notes = 'ISO 9601 formatted UTC timestamp [[ VALUES CAN BE PASSED FROM INPUT]]', varType='Support_Data',
                            monotone = 'Increase')
    ut_dt[:] = dn

    # Satellite coordinates
    satalt = datagrp.createVariable('ICON_WGS84_ALTITUDE','f8',('sat_location',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = satalt, catdesc='ICON WGS84 altitude coordinate at center time',
                            depend = 'Epoch',displayType = 'Coordinates', Fieldname = 'ICON Altitude (km)',
                            form = 'F8', labelAxis = 'Altitude', scaletype = 'Linear', units = 'km',
                            minmax = [130,580],Notes = 'A float that determines the satellite altitude coordinate', varType='Data')

    satalt[:] = satlatlonalt[2]

    satlat = datagrp.createVariable('ICON_WGS84_LATITUDE','f8',('sat_location',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = satlat, catdesc='ICON WGS84 latitude coordinate at center time',
                            depend = 'Epoch',displayType = 'Coordinates', Fieldname = 'ICON Latitude (degrees)',
                            form = 'F8', labelAxis = 'Latitude', scaletype = 'Linear', units = 'degrees',
                            minmax = [-90,90],Notes = 'A float that determines the satellite latitude coordinate', varType='Data')

    satlat[:] = satlatlonalt[0]

    satlon = datagrp.createVariable('ICON_WGS84_LONGITUDE','f8',('sat_location',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = satlon, catdesc='ICON WGS84 longitude coordinate at center time',
                            depend = 'Epoch',displayType = 'Coordinates', Fieldname = 'ICON Longitude (degrees)',
                            form = 'F8', labelAxis = 'Longitude', scaletype = 'Linear', units = 'degrees',
                            minmax = [-90,90],Notes = 'A float that determines the satellite longitude coordinate', varType='Data')

    satlon[:] = satlatlonalt[1]

    # Observation Geometry - This might need to be changed if 1 or 6 stripes are used on the output product
    zenith = datagrp.createVariable('ICON_L2_FUV_ZENITH','f8',('alt_profile',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = zenith, catdesc='FUV Zenith angle',
                            depend = 'Epoch',displayType = 'Time_Series', Fieldname = 'Zenith angles (degrees)',
                            form = 'F8', labelAxis = 'Zenith angle', scaletype = 'Linear', units = 'degrees',
                            minmax = [0,180],Notes = 'A float vector containing the zenith angles that correspond to the viewing geometry of the FUV instrument',
                            varType='data')

    zenith[:] = ze

    azimuth = datagrp.createVariable('ICON_L2_FUV_AZIMUTH','f8',('alt_profile',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = azimuth, catdesc='FUV Azimuth angle',
                            depend = 'Epoch',displayType = 'Time_Series', Fieldname = 'Azimuth angles (degrees)',
                            form = 'F8', labelAxis = 'Azimuth angle', scaletype = 'Linear', units = 'degrees',
                            minmax = [0,180],Notes = 'A float vector containing the azimuth angles that correspond to the viewing geometry of the FUV instrument',
                            varType='data')

    azimuth[:] = az


    # Tangent altitude coordinates
    tanalt = datagrp.createVariable('FUV_TANGENT_ALTITUDES','f8',('alt_profile',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = tanalt, catdesc='FUV tangent altitude coordinate at center time',
                            depend = 'Epoch',displayType = 'Coordinates', Fieldname = 'FUV tangent altitudes (km)',
                            form = 'F8', labelAxis = 'Tangent Altitude', scaletype = 'Linear', units = 'km',
                            minmax = [130,580],Notes = 'A float that determines the tangent altitude coordinates [[ THIS CAN COME FROM INPUT ]]', varType='Data')

    tanalt[:] = tanlatlonalt[:,2]

    tanlat = datagrp.createVariable('FUV_TANGENT_LATITUDES','f8',('alt_profile',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = tanlat, catdesc='FUV tangent latitude coordinates at center time',
                            depend = 'Epoch',displayType = 'Coordinates', Fieldname = 'FUV tangent latitude (degrees)',
                            form = 'F8', labelAxis = 'Tangent Latitude', scaletype = 'Linear', units = 'degrees',
                            minmax = [-90,90],Notes = 'A float that determines the tangent latitude coordinates [[ THIS CAN COME FROM INPUT ]]', varType='Data')

    tanlat[:] = tanlatlonalt[:,0]

    tanlon = datagrp.createVariable('FUV_TANGENT_LONGITUDES','f8',('alt_profile',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = tanlon, catdesc='FUV tangent longitude coordinates at center time',
                            depend = 'Epoch',displayType = 'Coordinates', Fieldname = 'FUV tangent longitude (degrees)',
                            form = 'F8', labelAxis = 'Tangent Longitude', scaletype = 'Linear', units = 'degrees',
                            minmax = [-90,90],Notes = 'A float that determines the tangent longitude coordinates [[ THIS CAN COME FROM INPUT ]]', varType='Data')

    tanlon[:] = tanlatlonalt[:,1]

    # Intensity emissions
    brightness = datagrp.createVariable('FUV_1356_INTENSITY_ALTITUDE_PROFILE','f8',('alt_profile',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = brightness, catdesc='Intensity of 135.6-nm emission at each observing angle',
                            depend = 'Epoch',displayType = 'Altitude_profile', Fieldname = 'Brightness (Rayleigh)',
                            form = 'F8', labelAxis = 'Brightness', scaletype = 'Linear', units = 'Rayleigh',
                            minmax = [0,1e7],Notes = 'A float vector containing the FUV nighttime intensity altitude profile for the 135.6-nm emissions ',
                            varType='data')

    brightness[:] = Bright

    # VER
    VolumeEmissionRate = datagrp.createVariable('FUV_1356_VER_ALTITUDE_PROFILE','f8',('alt_profile',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = VolumeEmissionRate, catdesc='VER of 135.6-nm emission at each observing angle',
                            depend = 'Epoch',displayType = 'Altitude_profile', Fieldname = 'VER (ph/cm^3/s)',
                            form = 'F8', labelAxis = 'Volume Emission Rate', scaletype = 'Linear', units = 'ph/cm^3/s',
                            minmax = [0,100],Notes = 'A float vector containing the FUV nighttime VER altitude profile for the 135.6-nm emissions ',
                            varType='data')

    VolumeEmissionRate[:] = VER

    # Electron density
    ElectronDensity = datagrp.createVariable('FUV_ELECTRON_DENSITY_ALTITUDE_PROFILE','f8',('alt_profile',),complevel=9,shuffle=True,chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = ElectronDensity, catdesc='Electron density at each observing angle',
                            depend = 'Epoch',displayType = 'Altitude_profile', Fieldname = 'Ne (1/cm^3)',
                            form = 'F8', labelAxis = 'Electron density', scaletype = 'Linear', units = '1/cm^3',
                            minmax = [0,1e10],Notes = 'A float vector containing the FUV nighttime Electron density altitude profile',
                            varType='data')

    ElectronDensity[:] = Ne


    # Satellite coordinates
    F2_Height = datagrp.createVariable('FUV_HMF2','f8',('F_peak_params',) ,complevel=9, shuffle=True, chunksizes=None, endian='little',fill_value=None)
    set_variable_attributes(var = F2_Height, catdesc='FUV F2-peak height at UTC time',
                            depend = 'Epoch',displayType = 'Altitude', Fieldname = 'hmF2(km)',
                            form = 'F8', labelAxis = 'hmF2', scaletype = 'Linear', units = 'km',
                            minmax = [130,580],Notes = 'A float that determines the electron density peak tangent altitude', varType='Data')

    F2_Height[:] = hmF2

    F2_Density = datagrp.createVariable('FUV_NMF2','f8',('F_peak_params',), complevel=9, shuffle=True, chunksizes=None, endian='little', fill_value=None)
    set_variable_attributes(var = F2_Density, catdesc='FUV F2-peak density at UTC time',
                            depend = 'Epoch',displayType = 'Density', Fieldname = 'NmF2 (1/cm^3)',
                            form = 'F8', labelAxis = 'NmF2', scaletype = 'Linear', units = '1/cm^3',
                            minmax = [0,1e10],Notes = 'A float that determines the electron density peak density', varType='Data')

    F2_Density[:] = NmF2

    
    datagrp.close()
  

'''
NETCDF ATTRIBUTE COMPLETE FUNCTION 
  
'''

def set_variable_attributes(var, catdesc, depend, displayType, Fieldname, form, labelAxis, scaletype, units, minmax,Notes, varType,monotone = ''):

      # 3.3.3.1
    var.CatDesc = catdesc
    var.Long_Name = catdesc

    # 3.3.3.2
    var.Depend_0 = depend # This needs fixing in regards to the epoch

    # 3.3.3.3
    var.Display_Type = displayType

    # 3.3.3.4
    var.FieldNam = Fieldname

    # 3.3.3.5 
    #satlat._FillValue = np.nan

    # 3.3.3.6
    var.Format = form

    # 3.3.3.7
    var.LablAxis = labelAxis

    # 3.3.3.8 - Monotony
    var.MonoTon = monotone # ISTP optional

    # 3.3.3.9
    var.ScaleTyp = scaletype

    # 3.3.3.10
    var.Units = units

    # 3.3.3.11 
    var.ValidMax = minmax[1]
    var.ValidMin = minmax[0]
    var.Valid_Max = minmax[1]
    var.Valid_Min = minmax[0]
    var.Valid_Range = minmax

    # 3.3.3.12 
    var.Var_Notes = Notes

    # 3.3.3.13 
    var.Var_Type = varType



'''

FUV LEVEL 2 TOP LEVEL FUNCTION

'''
  
def Get_lvl2_5_product(path_input='/home/dimitris/Data_Files/ICON_FUV_ray_UT_15sec_night.nc',path_output='/home/dimitris/public_html/Datafiles/LVL2TEST/',limb = 150.,Spherical = True,reg_method = 'Tikhonov',regu_order = 2, contribution ='RRMN'):
    '''
    Operational Code that reads Lvl1 file and creates the corresponding Lvl2.5
    INPUTS:
        path_input  - Input file path 
        path_output - Output file path
        limb        - Defines the lower bound that defines the limb [km]
        Spherical   - Flag indicating whether Sperical or Non-Spherical Earth is assumed [True, False]
        reg_method  - Regularization method selection [Tikhonov, MAP]
        regu_order  - Regularization Order [int] (0,1,2 Possible Values)
        contribution- contributions to the 135.6nm emission, RR: Radiative Recombination, RRMN: Radiative Recombination and Mutual Neutralization 
    OUTPUT:
        Creates a NetCDF file on the desired output_path
    NOTES:
        This versions uses paths as input and output. That can change in case needed. Also, the lvl1 file that is used
        is elementary. No actual LVL1 file has been given yet.
    HISTORY:
        24-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    TODO:
        Add uncertainty propagation, once the input uncertainties are specified
    '''
    # Open input NetCDF file
    data = io.netcdf.netcdf_file(path_input,mode='r')
    
    # Get Data from file.
    # For now we use only one stripe from the image. This process must be done for the whole CCD
    FUV_1356_IMAGE = data.variables['FUV_1356_IMAGE'][0][:][:]

    FUV_TANGENT_ALTITUDES = data.variables['FUV_TANGENT_ALTITUDES'][0][:][:]
    FUV_TANGENT_ALTITUDES_END = data.variables['FUV_TANGENT_ALTITUDES_END'][0][:][:]
    FUV_TANGENT_ALTITUDES_START = data.variables['FUV_TANGENT_ALTITUDES_START'][0][:][:]

    FUV_TANGENT_LATITUDES = data.variables['FUV_TANGENT_LATITUDES'][0][:][:]
    FUV_TANGENT_LATITUDES_END = data.variables['FUV_TANGENT_LATITUDES_END'][0][:][:]
    FUV_TANGENT_LATITUDES_START = data.variables['FUV_TANGENT_LATITUDES_START'][0][:][:]

    FUV_TANGENT_LONGITUDES = data.variables['FUV_TANGENT_LONGITUDES'][0][:][:]
    FUV_TANGENT_LONGITUDES_END = data.variables['FUV_TANGENT_LONGITUDES_END'][0][:][:]
    FUV_TANGENT_LONGITUDES_START = data.variables['FUV_TANGENT_LONGITUDES_START'][0][:][:]

    FUV_TANGENT_N2 = data.variables['FUV_TANGENT_N2'][0][:][:]
    FUV_TANGENT_O1 = data.variables['FUV_TANGENT_O1'][0][:][:]
    FUV_TANGENT_O2 = data.variables['FUV_TANGENT_O2'][0][:][:]
    FUV_TANGENT_OP = data.variables['FUV_TANGENT_OP'][0][:][:]

    FUV_TANGENT_POINT_INDEX = data.variables['FUV_TANGENT_POINT_INDEX'][0][:][:]

    FUV_ECEF_VECTORS_START = data.variables['FUV_ECEF_VECTORS_START'][0][:][:]

    ICON_WGS84_LATITUDE_START =  data.variables['ICON_WGS84_LATITUDE_START'][0]
    ICON_WGS84_LONGITUDE_START =  data.variables['ICON_WGS84_LONGITUDE_START'][0]
    ICON_WGS84_ALTITUDE_START =  data.variables['ICON_WGS84_ALTITUDE_START'][0]

    satlatlonalt = [ICON_WGS84_LATITUDE_START,ICON_WGS84_LONGITUDE_START,ICON_WGS84_ALTITUDE_START]

    ICON_UT_START = data.variables['ICON_UT_START'].data
    
    # Stripe 4 used for now. We know that this stripe contains all the information without any loss of data
    STRIPE = 4
    
    # Calculate viewing geometry vectors
    FUV_ZE_VECTORS_START = np.zeros((np.size(FUV_ECEF_VECTORS_START,0),np.size(FUV_ECEF_VECTORS_START,1)));
    FUV_AZ_VECTORS_START = np.zeros((np.size(FUV_ECEF_VECTORS_START,0),np.size(FUV_ECEF_VECTORS_START,1)));
    for i in range(0,np.size(FUV_ECEF_VECTORS_START,0)):
        for j in range(0,np.size(FUV_ECEF_VECTORS_START,1)):
            [FUV_AZ_VECTORS_START[i][j],FUV_ZE_VECTORS_START[i][j]]= ic.ecef_to_azze(satlatlonalt,FUV_ECEF_VECTORS_START[i,j,:])
    
    limb = np.where(FUV_TANGENT_ALTITUDES[:,STRIPE]>=150)
    # Call the function that calculates the solution
    dn = ICON_UT_START[limb][::-1]
    
    az = FUV_AZ_VECTORS_START[limb,STRIPE][0][::-1]
    ze = FUV_ZE_VECTORS_START[limb,STRIPE][0][::-1]
    
    #tanlatlonalt = [FUV_TANGENT_LATITUDES[limb,STRIPE][0][::-1],FUV_TANGENT_LONGITUDES[limb,STRIPE][0][::-1],FUV_TANGENT_ALTITUDES[limb,STRIPE][0][::-1]]
    tanlatlonalt = np.column_stack((FUV_TANGENT_LATITUDES[limb,STRIPE][0][::-1],FUV_TANGENT_LONGITUDES[limb,STRIPE][0][::-1],FUV_TANGENT_ALTITUDES[limb,STRIPE][0][::-1]))

    bright = FUV_1356_IMAGE[limb,:][0][::-1]
    bright = np.mean(bright,axis=1)
    # h needs to be changed when non-symmetric earth is assumed. For now we have symmetry
    h = FUV_TANGENT_ALTITUDES[limb,STRIPE][0][::-1]
    # That might not work well now if non-symmetric earth is assumed from Scott
    O = FUV_TANGENT_O1[limb,STRIPE][0][::-1]

    data.close()
    
    VER_mean = 0
    VER_profiles = 0
    
    ver,Ne,h = FUV_Level_2_Density_Calculation(Bright = bright,alt_vector=h,satlatlonalt=satlatlonalt,az=az,
                                               ze = ze, limb = limb,Spherical = Spherical, reg_method = reg_method,
                                               regu_order = regu_order, contribution =contribution,VER_mean= VER_mean,
                                               VER_profiles = VER_profiles,dn = dn)
  
    hmF2,NmF2 = find_hm_Nm_F2(Ne,h)
    
    # Filename needs to be added here, or the parameters to set the fn
    FUV_Level_2_OutputProduct_NetCDF(dn,satlatlonalt,az,ze,tanlatlonalt,bright,ver,Ne,NmF2,hmF2,path_output)
    
    print 'LVL2.5 Processing Terminated. File created!'

    
    
  
    
    
