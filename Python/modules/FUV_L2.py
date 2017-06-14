# coding: utf-8
'''
This module contains the algorithms that are used for FUV level 2.5 nighttime data product.
Includes:
    - the reading and writing of NetCDF files
    - the regularization algorithms to estimate the volume emission rate for the 135.6nm O+ emission during nighttime
    - the calculation of the electron density content assuming radiative recombination and mutual neutralization
Todo:
    - error flags (e.g., chi^2, ill-defined hmF2-finding)
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

# Modules needed for parsing dates
from dateutil import parser

import getpass

import os

# For Monte Carlo processing:
from numpy.random import multivariate_normal

from scipy.linalg import sqrtm # matrix square root

# For plotting
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

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
        31-May-2017: Recast calclation of S using np matrix math, rather than nested
                     for loops (Jonathan J. Makela; jmakela@illinois.edu)
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

    # Build observation matrix and other required matrix variables
    S = np.zeros((Ntheta, len(rmid)))
    th = np.tile(theta,(len(theta),1)).T
    rb = np.tile(rbot,(len(rbot),1))
    rt = np.tile(rtop,(len(rtop),1))
    sb2 = -np.sin(th)**2 + ((RE+rb)/(RE+Horbit))**2
    st2 = -np.sin(th)**2 + ((RE+rt)/(RE+Horbit))**2

    sb2[sb2<0] = 0. # there is no intersection of LOS with altitude rb. Set term to 0.
    st2[st2<0] = 0. # there is no intersection of LOS with altitude rt. Set term to 0.

    # Calculate path length
    path_len_km = 2.*(RE+Horbit) * (np.sqrt(st2) - np.sqrt(sb2))

    # Go to cm
    path_len_cm = path_len_km*1e5

    # Result should be a conversion from VER to Rayleighs. path_length_cm * 10^-6 to match paper
    S = path_len_cm * 1e-6

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
        03-May-2017: Simplified by Brian Harding (bhardin2@illinois.edu)
    CALLS:
    '''



    # Old code from Dimitris. It scales nicely with A, but doesn't scale
    # with changes in brightness.
    ## SVD Decomposition of matrix A (Distance matrix)
    #U, s, V = np.linalg.svd(A, full_matrices=True)

    # multiplication ratio
    #smin_ratio = 16*np.spacing(1)
    ##smin_ratio = 16*np.finfo(np.float32).eps
    #reg_param = np.zeros(npoints)
    #reg_param[npoints-1] = max([s[np.size(s,0)-1],s[1]*smin_ratio])
    #ratio = (s[0]*100/reg_param[npoints-1])**(1./(npoints-1));

    ## Put regularization parameters in descending order
    #for i in np.arange(npoints-2,-1,-1):
    #    reg_param[i] = ratio*reg_param[i+1]



    # Updated by Brian. Will we ever need a larger range?
    reg_param = np.logspace(5.5,1,npoints)

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
def calculate_electron_density(VER,satlatlonalt,tang_altitude,dt,Sig_VER=None,contribution='RR',f107=None, f107a=None,apmsis=None):
    '''
    Given a VER profile the electron density is calculated. The physical process for the VER must be defined.
    INPUTS:
        VER             - Volume Emission Rate profile (ph/cm^-3)
        satlatlonalt    - satellite coordinates [latitude(degrees),longitude(degrees),altitude(km)]
        tang_altitude   - tangent altitude [km]
        dt              - datetime (python datetime)
        Sig_VER         - Covariance matrix of VER. If None, Ne will be the only output. [(ph/cm^-3)**2]
        contribution    - contributions to the 135.6nm emission, RR: Radiative Recombination, RRMN: Radiative Recombination and Mutual Neutralization
        f107        - f107 value for the date of interest (None is use values in PyGlow)
        f107a       - f107a value for the date of interest (None is use values in PyGlow)
        apmsis      - apmsis vector for the date of interest (None is use values in PyGlow)
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
        02-Jun-2017: Added GPI values by Jonathan Makela (jmakela@illinois.edu)
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

            # Create pyglow point, either use the default GPI or the passed in GPI
            if apmsis is None:
                pt = pyglow.pyglow.Point(dt, satlatlonalt[0], satlatlonalt[1], height)
                # Run MSIS-00 to get O density
            else:
                pt = pyglow.pyglow.Point(dt, satlatlonalt[0], satlatlonalt[1], height,user_ind=True)
                pt.f107 = f107
                pt.f107a = f107a
                pt.apmsis = apmsis

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
def hmF2_region_interpolate(Ne,alt_vector,interval_increase = 20.,kind ='cubic',debug=False):
    '''
    Interpolates the values around the hmF2 region in order to increase altitude vector interval which helps in minimizing the hmF2 error
    INPUTS:
        Ne                 - Original Electron density profile [cm^{-3}]
        alt_vector         - altitude vector for original electron density profile [km]
        interval_increase  - Number indicating the amount of points to be created. => len(hmf2Region)* int_increase
        kind               - Goes as input to the interp1d function and declares the type of interpolation
        debug              - Flag to print debug warnings
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
        if debug:
            print 'Ne maximum ill defined - No interpolation performed'
        # In this case the values remain the same. No exception is required since the F2-peak code can run
        return Ne,alt_vector


'''

L2 Calculation top-level function

'''

def FUV_Level_2_Density_Calculation(Bright,alt_vector,satlatlonalt,az,ze, Sig_Bright = None, weight_resid = False, limb = 150.,Spherical = True,reg_method='Tikhonov',regu_order = 2,reg_param=0,contribution='RRMN', VER_mean=0.,VER_profiles=0,dn = datetime.datetime(2012, 9, 26, 20, 0, 0), f107=None, f107a=None, apmsis=None):
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
        f107        - f107 value for the date of interest (None is use values in PyGlow)
        f107a       - f107a value for the date of interest (None is use values in PyGlow)
        apmsis      - apmsis vector for the date of interest (None is use values in PyGlow)
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
        02-Jun-2017: Added GPI values by Jonathan Makela (jmakela@illinois.edu)
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

    Ne, Sig_Ne = calculate_electron_density(VER=VER, satlatlonalt=satlatlonalt, tang_altitude=h, dt=dn, Sig_VER=Sig_VER, contribution=contribution,f107=f107, f107a=f107a, apmsis=apmsis)

    if ret_cov:
        return VER, Ne, h, Sig_VER, Sig_Ne
    else:
        return VER,Ne,h

'''

HELPER FUNCTION FOR CREATING NETCDF4 variables

'''

def _create_variable(ncfile, name, value, format_nc='f8', format_fortran='F', dimensions=(), zlib=True, complevel=6,
                    shuffle=True,  depend_0=None, depend_1=None, depend_2=None, depend_3=None, chunk_sizes=None, desc='',
                    display_type='scalar', field_name='', fill_value=None,label_axis='', bin_location=0.5,
                    time_base='FIXED: 1970 (POSIX)', time_scale='UTC', units='', valid_min=None, valid_max=None,
                    notes='', var_type='data'):
    '''
    A helper function to write a variable to a netCDF file.

    INPUTS:

      *  Self evident from the parameters above. Notes:

            * fill_value = None --> default fill values will be used, if they exist. See netCDF4.default_fillvals
            * display_type: for now, 'scalar', 'altitude_profile', or 'image' will be used
            * var_type: one of 'data', 'support_data', 'metadata', 'ignore_data'
            * format_fortran: Used by ISTP. See http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
            * except as specified above, if a variable attribute is left as the default None, it will not be written to the file

    OUTPUT:

      *  The netCDF4._netCDF4.Variable object that was created and to which was written

    '''

    # Rudimentary error-checking:
    valid_var_types = ['data','support_data','metadata','ignore_data']
    if var_type not in valid_var_types:
        raise Exception('var_type="%s" is not valid. Try one of: %s' % (var_type, valid_var_types) )
    if len(desc) > 80:
        raise Exception('"desc" is too long (%i chars). Shorten to 80 characters:\n"%s"' % (len(desc),desc))
    if len(field_name) > 30:
        raise Exception('field_name="%s" is too long (%i chars). Shorten to 30 characters:\n' % (field_name,len(field_name)))
    if len(label_axis) > 10:
        raise Exception('label_axis="%s" is too long (%i chars). Shorten to 10 characters.' % (label_axis,len(label_axis)))

    # If fill value was not specified, use the default value, if it exists.
    # It will not exist for strings, for example, for which fill values
    # cannot be set. (TODO: is this right?)
    if fill_value is None and format_nc in netCDF4.default_fillvals.keys():
        fill_value = netCDF4.default_fillvals[format_nc]

    var = ncfile.createVariable(name, format_nc, dimensions=dimensions, zlib=zlib, complevel=complevel,
                               shuffle=shuffle, chunksizes=chunk_sizes, fill_value=fill_value)
    var.CatDesc            = desc
    var.Long_Name          = desc
    if chunk_sizes is not None:
        var._ChunkingSizes = chunk_sizes
    var._DeflateLevel      = complevel
    var._Shuffle           = str(shuffle).lower()
    if depend_0 is not None:
        var.Depend_0       = depend_0
    if depend_1 is not None:
        var.Depend_1       = depend_1
    if depend_2 is not None:
        var.Depend_2       = depend_2
    if depend_3 is not None:
        var.Depend_3       = depend_3
    var.Display_Type       = display_type
    var.FieldNam           = field_name
    # Note: t_var._FillValue not expliclity needed since that is set by the createVariable function argument "fill_value"
    #var._FillValue         = fill_value
    if fill_value is not None:
        var.FillVal        = var._FillValue
    elif fill_value is None and format_nc == str:
        # Special case for strings. Make sure to set FillVal even thought _FillValue can't be set
        var.FillVal        = ''

    var.Format             = format_fortran
    var.LablAxis           = label_axis
    var.Bin_Location       = bin_location
    var.Time_Base          = time_base
    var.Time_Scale         = time_scale
    var.Units              = units
    if valid_min is not None:
        var.ValidMin       = valid_min
        var.Valid_Min      = valid_min
    if valid_max is not None:
        var.ValidMax       = valid_max
        var.Valid_Max      = valid_max
    var.Var_Notes          = notes
    var.Var_Type           = var_type

    # If a fill_value was specified, and if there are any np.nan values in
    # the variable, replace them with the fill value.
    if fill_value is not None:
        # For sequences that are not strings:
        if hasattr(value,'__len__') and not isinstance(value,(str,unicode)):
            value[np.isnan(value)] = var._FillValue
        # For non-sequences and strings:
        elif np.isnan(value):
            value = var._FillValue

    # Assign value
    var[...] = value

    return var

'''

CREATE NETCDF FILE

'''
def FUV_Level_2_OutputProduct_NetCDF(L25_full_fn, L25_dict):
    '''
    This function takes as input all the necessary outputs from LVL2.5 processing and writes them on a NetCDF file
    INPUTS:
        L25_full_fn        - Full file name where the NetCDF file is saved
        L25_dict           - Dictionary containing data to be written (see Get_lvl2_5_product for dictionary definition)
    OUTPUT:
        Creates a NetCDF file on the desired path
    NOTES:

    HISTORY:
        24-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        4-May-2017: Rewritten by Jonathan Makela to conform to L25 format (jmakela@illinois.edu)
        13-Jun-2017: Send in full file path, rather than generating on own (jmakela@illinois.edu)
        14-Jun-2017: Pull version and revision information from the provided filename
    '''

    # Parse the filename for version and revision. Assumes the provided filename conforms to
    # the ICON conventions: ICON_<LEVEL>_<INSTRUMENT>[_<DESCRIPTION>]_<DATE>[_<TOD>]_v<VERSION>r<REVISION>.NC
    data_versionmajor = np.int(L25_full_fn[-9:-7])
    data_revsision = np.int(L25_full_fn[-6:-3])
    data_version =  np.float(data_versionmajor)+data_revsision/1000.# TODO: how will this be calculated? It goes into global attr Data_Version
    software_versionmajor = data_versionmajor
    software_revision = 0   # Increment this when making a revision!
    software_version =  np.float(software_versionmajor)+software_revision/1000.

    # TODO: How will sensor be determined? Will it be in L1 file?
    source_files = ' '.join(L25_dict['source_files'])
    sensor = 'FUV'

    #################### Compile variables to write in file ######################
    t_file  = datetime.datetime.now()   # time this file was created
    ### Who's running this process
    user_name = getpass.getuser()
    ### Parent files
    parents = '' # This will go in global attr Parents
    for source_fn in L25_dict['source_files']:
        s = source_fn.split('/')[-1].split('.')
        pre = '.'.join(s[:-1])
        post = s[-1].upper()
        parents += '%s > %s, ' % (post, pre)
    if parents: parents = parents[:-2] # trim trailing comma
    ### Parameters file, if there is one (note, that this will go in Calibration_File global attribute)
    param_fn = L25_dict['param_fn']
    if param_fn is None:
        param_fn = ''

    ########################## Open file for writing #############################
    ncfile = netCDF4.Dataset(L25_full_fn,mode='w',format='NETCDF4')

    # Parse out the short filename
    L25_fn = os.path.basename(L25_full_fn)

    ########################## Global Attributes #################################
    ncfile.Acknowledgement =       ''.join(("This is a data product from the NASA Ionospheric Connection Explorer mission, ",
                                            "an Explorer launched in June 2017.\n",
                                            "\n",
                                            "Responsibility of the mission science falls to the Principal Investigator, ",
                                            "Dr. Thomas Immel at UC Berkeley.\n",
                                            "\n",
                                            "Validation of the L1 data products falls to the instrument lead ",
                                            "investigators/scientists.\n",
                                            "  * EUV  Dr. Eric Korpela\n",
                                            "  * FUV  Dr. Harald Frey\n",
                                            "  * MIGHTI  Dr. Chris Englert\n",
                                            "  * IVM  Dr. Roderick Heelis\n",
                                            "\n",
                                            "Validation of the L2 data products falls to those responsible for those products.\n",
                                            "  * O/N2  Dr. Andrew Stephan\n",
                                            "  * Daytime (EUV) O+ profiles  Dr. Andrew Stephan\n",
                                            "  * Nighttime (FUV) O+ profiles  Dr. Farzad Kamalabadi\n",
                                            "  * Neutral Wind profiles  Dr. Jonathan Makela\n",
                                            "  * Neutral Temperature profiles  Dr. Chris Englert\n",
                                            "\n",
                                            "Responsibility for Level 4 products are detailed on the ICON website ",
                                            "(http://icon.ssl.berkeley.edu).\n",
                                            "\n",
                                            "Overall validation of the products is overseen by the ICON Project Scientist ",
                                            "Dr. Scott England.\n",
                                            "\n",
                                            "NASA oversight for all products is provided by the Mission Scientist ",
                                            "Dr. Douglas Rowland.\n",
                                            "\n",
                                            "Users of these data should contact and acknowledge the Principal Investigator ",
                                            "Dr. Immel and the party directly responsible for the data product and the NASA ",
                                            "Contract Number NNG12FA45C from the Explorers Project Office." ))

    ncfile.ADID_Ref =                       'NASA Contract > NNG12FA45C'
    ncfile.Calibration_File =               param_fn
    ncfile.Conventions =                    'SPDF ISTP/IACG Modified for NetCDF'
    ncfile.Data_Level =                     'L2.5'
    ncfile.Data_Type =                      'DP25 > Data Product 2.5: FUV Nighttime O+ profile'
    ncfile.Data_Version =                   np.float32(data_version)
    ncfile.Data_VersionMajor =              np.ubyte(data_versionmajor)
    ncfile.Data_Revision =                  np.ushort(data_revsision)
    ncfile.Date_Stop =                      L25_dict['FUV_dn'][0].strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC' # single measurement: use midpoint
    ncfile.Date_Start =                     L25_dict['FUV_dn'][-1].strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC' # single measurement: use midpoint
    ncfile.Description =                    'ICON FUV Nighttime O+ profiles (DP 2.5)'
    ncfile.Descriptor =                     'FUV > Intensified Far Ultraviolet Imager'
    ncfile.Discipline =                     'Space Physics > Ionospheric Science'
    ncfile.File =                           L25_fn
    ncfile.File_Date =                      t_file.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC'
    ncfile.Generated_By =                   'ICON SDC > ICON UIUC FUV L2.5 Processor v%.2f, J. J. Makela, D. Iliou' % software_version
    ncfile.Generation_Date =                t_file.strftime('%Y%m%d')
    ncfile.History =                        'Version %i, %s, %s, ' % (data_version, user_name, t_file.strftime('%Y-%m-%dT%H:%M:%S')) +\
                                            'FUV L2.5 Processor v%.2f ' % software_version # TODO: Tori suggested we make this a history
                                                                                              # of the software versions instead of data versions
    ncfile.HTTP_LINK =                      'http://icon.ssl.berkeley.edu/Instruments/FUV'
    ncfile.Instrument =                     'FUV'
    ncfile.Instrument_Type =                'Imagers (space)'
    ncfile.Link_Text =                      '12-second FUV SW altitude profiles'
    ncfile.Link_Title =                     'ICON FUV'
    ncfile.Logical_File_ID =                L25_fn[:-3]
    ncfile.Logical_Source =                 'ICON_L2_5_FUV_SWP'
    ncfile.Logical_Source_Description =     'FUV Short Wavelength Channel - 135.6 Altitude Profiles'
    ncfile.Mission_Group =                  'Ionospheric Investigations'
    ncfile.MODS =                           ncfile.History
    ncfile.Parents =                        parents
    ncfile.PI_Affiliation =                 'UC Berkeley > SSL'
    ncfile.PI_Name =                        'T. J. Immel'
    ncfile.Project =                        'NASA > ICON'
    ncfile.Rules_of_Use =                   'Public Data for Scientific Use'
    ncfile.Software_Version =               'ICON SDC > ICON UIUC FUV L2.5 Processor v%.3f' % software_version
    ncfile.Source_Name =                    'ICON > Ionospheric Connection Explorer'
    ncfile.Spacecraft_ID =                  'NASA > ICON - 493'
    ncfile.Text =                           'ICON explores the boundary between Earth and space - the ionosphere - ' +\
                                            'to understand the physical connection between our world and the immediate '+\
                                            'space environment around us. Visit \'http://icon.ssl.berkeley.edu\' for more details.'
    ncfile.Text_Supplement =                '((TODO: Insert reference to Kamalabadi et al [2017] paper))'
    ncfile.Time_Resolution =                '12 seconds'
    ncfile.Title =                          'ICON FUV O+ Altitude Profiles (DP 2.5)'

    ################################## Dimensions ########################################
    n = np.shape(L25_dict['FUV_TANGENT_ALT'])[1]
    t = np.shape(L25_dict['FUV_TANGENT_ALT'])[0]
    ncfile.createDimension('Epoch',t)
    ncfile.createDimension('Vertical', n)
    ncfile.createDimension('Horizontal',6)
    ncfile.createDimension('Vector',3)

    # Get instrument name (FUVA or FUVB)
    inst = 'FUVA' # TODO: This needs to be automatically calculated

    ################################## Variables #########################################

    ######### Timing Variables #########

    # Time midpoint (the official required "Epoch" variable)
    var_epoch = _create_variable(ncfile, 'EPOCH', L25_dict['FUV_EPOCH'],
                          dimensions=('Epoch'),
                          format_nc='i8', format_fortran='I', desc='Milliseconds since 1970-01-01 00:00:00 UTC at middle of measurement integration.',
                          display_type='Time_Series', field_name='ms epoch', fill_value=-1, label_axis='Time', bin_location=0.5,
                          units='ms', valid_min=np.int64(0), valid_max=np.int64(1000*365*86400e3), var_type='support_data', chunk_sizes=[1],
                          depend_0 = 'Epoch', notes='')

    var_time = _create_variable(ncfile, 'ICON_L2_%s_SWP_CENTER_TIMES' % inst, L25_dict['FUV_CENTER_TIMES'],
                          dimensions=('Epoch'),
                          format_nc=str, format_fortran='A23', desc='Center time of 12-second profile integration',
                          display_type='Time_Series', field_name='Center time', fill_value=None, label_axis='Time', bin_location=0.5,
                          units='s', valid_min='1970-01-01/00:00:00', valid_max='2100-12-31/23:59:59', var_type='support_data', chunk_sizes=[1],
                          depend_0 = 'Epoch', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_START_TIMES' % inst, L25_dict['FUV_START_TIMES'],
                          dimensions=('Epoch'),
                          format_nc=str, format_fortran='A23', desc='Start time of 12-second profile integration',
                          display_type='Time_Series', field_name='Start time', fill_value=None, label_axis='Time', bin_location=0.5,
                          units='s', valid_min='1970-01-01/00:00:00', valid_max='2100-12-31/23:59:59', var_type='support_data', chunk_sizes=[1],
                          depend_0 = 'Epoch', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_STOP_TIMES' % inst, L25_dict['FUV_STOP_TIMES'],
                          dimensions=('Epoch'),
                          format_nc=str, format_fortran='A23', desc='Stop time of 12-second profile integration',
                          display_type='Time_Series', field_name='Stop time', fill_value=None, label_axis='Time', bin_location=0.5,
                          units='s', valid_min='1970-01-01/00:00:00', valid_max='2100-12-31/23:59:59', var_type='support_data', chunk_sizes=[1],
                          depend_0 = 'Epoch', notes='')

    ######### Geometry Variables #########

    # ICON Spacecraft location in WGS
    var = _create_variable(ncfile, 'ICON_L2_FUV_SC_WGS', L25_dict['ICON_WGS'],
                          dimensions=('Epoch','Vector'),
                          format_nc='f8', format_fortran='F', desc='Spacecraft location in WGS84 (Latitude, Longitude, Altitude).',
                          display_type='Time_Series', field_name='Spacecraft location in WGS84', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='', valid_min=-180., valid_max=1000., var_type='support_data', chunk_sizes=[1,1],
                          depend_0 = 'Epoch', depend_1 = 'Vector', notes='')

    # FUV tangent point locations in WGS
    var = _create_variable(ncfile, 'ICON_L2_%s_TANGENT_LAT' % inst, L25_dict['FUV_TANGENT_LAT'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='Tangent latitude in WGS84.',
                          display_type='Time_Series', field_name='Tangent latitude in WGS84', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='degrees North', valid_min=-90., valid_max=90., var_type='support_data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_TANGENT_LON' % inst, L25_dict['FUV_TANGENT_LON'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='Tangent longitude in WGS84.',
                          display_type='Time_Series', field_name='Tangent longitude in WGS84', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='degrees East', valid_min=0., valid_max=360., var_type='support_data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_TANGENT_ALT' % inst, L25_dict['FUV_TANGENT_ALT'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='Tangent altitude in WGS84.',
                          display_type='Time_Series', field_name='Tangent altitude in WGS84', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='km', valid_min=-90., valid_max=90., var_type='support_data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    # FUV look direction AZ and ZE
    var = _create_variable(ncfile, 'ICON_L2_%s_FOV_AZIMUTH_ANGLE' % inst, L25_dict['FUV_AZ'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='FOV Celestial Azimuth',
                          display_type='Time_Series', field_name='FOV Celestial Azimuth', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='degrees', valid_min=-180., valid_max=180., var_type='support_data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_FOV_ZENITH_ANGLE' % inst, L25_dict['FUV_ZE'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='FOV Celestial Zenith',
                          display_type='Time_Series', field_name='FOV Celestial Azimuth', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='degrees', valid_min=-90., valid_max=90., var_type='support_data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    ######### Results Variables #########

    # FUV inverted volume emission rate
    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_VER_ALTITUDE_PROFILE' % inst, L25_dict['FUV_ver'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='VER of 135.6-nm emission at each observing angle',
                          display_type='Time_Series', field_name='VER', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='ph/cm^3/s', valid_min=0., valid_max=100., var_type='data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_VER_ERROR_ALTITUDE_PROFILE' % inst, L25_dict['FUV_ver_error'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='Error in VER of 135.6-nm emission at each observing angle',
                          display_type='Time_Series', field_name='VER', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='ph/cm^3/s', valid_min=0., valid_max=100., var_type='data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    # FUV inverted electron density
    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_ELECTRON_DENSITY_ALTITUDE_PROFILE' % inst, L25_dict['FUV_Ne'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='Electron density at each observing angle',
                          display_type='Time_Series', field_name='Ne', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='1/cm^3', valid_min=0., valid_max=1e10, var_type='data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_ELECTRON_DENSITY_ERROR_ALTITUDE_PROFILE' % inst, L25_dict['FUV_Ne_error'],
                          dimensions=('Epoch','Vertical','Horizontal'),
                          format_nc='f8', format_fortran='F', desc='Error in electron density at each observing angle',
                          display_type='Time_Series', field_name='Ne', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='1/cm^3', valid_min=0., valid_max=1e10, var_type='data',
                          chunk_sizes=[1,ncfile.dimensions['Vertical'].size,ncfile.dimensions['Horizontal'].size],
                          depend_0 = 'Epoch', depend_1='Vertical',depend_2='Horizontal', notes='')

    # FUV inverted hmF2
    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_HMF2' % inst, L25_dict['FUV_hmF2'],
                          dimensions=('Epoch','Horizontal'),
                          format_nc='f4', format_fortran='F', desc='FUV F2-peak height estimate',
                          display_type='Time_Series', field_name='hmF2', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='km', valid_min=130., valid_max=580., var_type='data', chunk_sizes=[1,1],
                          depend_0 = 'Epoch',depend_1='Horizontal', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_HMF2_ERROR' % inst, L25_dict['FUV_hmF2_error'],
                          dimensions=('Epoch','Horizontal'),
                          format_nc='f4', format_fortran='F', desc='FUV F2-peak height estimate',
                          display_type='Time_Series', field_name='hmF2', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='km', valid_min=130., valid_max=580., var_type='data', chunk_sizes=[1,1],
                          depend_0 = 'Epoch',depend_1='Horizontal', notes='')

    # FUV inverted NmF2
    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_NMF2' % inst, L25_dict['FUV_NmF2'],
                          dimensions=('Epoch','Horizontal'),
                          format_nc='f4', format_fortran='F', desc='FUV F2-peak density estimate',
                          display_type='Time_Series', field_name='NmF2', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='1/cm^3', valid_min=0., valid_max=1e10, var_type='data', chunk_sizes=[1,1],
                          depend_0 = 'Epoch',depend_1='Horizontal', notes='')

    var = _create_variable(ncfile, 'ICON_L2_%s_SWP_NMF2_ERROR' % inst, L25_dict['FUV_NmF2_error'],
                          dimensions=('Epoch','Horizontal'),
                          format_nc='f4', format_fortran='F', desc='FUV F2-peak density estimate',
                          display_type='Time_Series', field_name='NmF2', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='1/cm^3', valid_min=0., valid_max=1e10, var_type='data', chunk_sizes=[1,1],
                          depend_0 = 'Epoch',depend_1='Horizontal', notes='')

    # FUV inversion quality flag
    var = _create_variable(ncfile, 'ICON_L2_%s_QUALITY' % inst, L25_dict['FUV_Quality'],
                          dimensions=('Epoch','Horizontal'),
                          format_nc='i1', format_fortran='I', desc='FUV inversion quality flag',
                          display_type='Time_Series', field_name='Quality', fill_value=-999, label_axis='Time', bin_location=0.5,
                          units='', valid_min=0., valid_max=8, var_type='data', chunk_sizes=[1,1],
                          depend_0 = 'Epoch',depend_1='Horizontal', notes='Flag is 0 if the inversion is valid.')

    # Regularization type
    var = _create_variable(ncfile, 'ICON_L2_REGULARIZATION_METHOD', L25_dict['reg_method'],
                          dimensions=(),
                          format_nc=str, format_fortran='A', desc='Regularization method: Tikhonov or MAP',
                          display_type='scalar', field_name='Regularzation Method', fill_value=None, label_axis='Method', bin_location=0.5,
                          units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=1,
                          notes='')

    # Regularization order
    var = _create_variable(ncfile, 'ICON_L2_REGULARIZATION_ORDER', L25_dict['regu_order'],
                          dimensions=(),
                          format_nc='i1', format_fortran='I', desc='Regularization order',
                          display_type='scalar', field_name='Regularzation Order', fill_value=None, label_axis='Order', bin_location=0.5,
                          units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=1,
                          notes='')

    # Contribution type
    var = _create_variable(ncfile, 'ICON_L2_1356_CONTRIBUTION', L25_dict['contribution'],
                          dimensions=(),
                          format_nc=str, format_fortran='A', desc='Chemestry considered in 135.6-nm emission: RR or RRMN',
                          display_type='scalar', field_name='1356 Contribution', fill_value=None, label_axis='Contrib', bin_location=0.5,
                          units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=1,
                          notes='')

    # Limb value used in inversion
    var = _create_variable(ncfile, 'ICON_L2_INVERSION_LIMB', L25_dict['limb'],
                          dimensions=(),
                          format_nc='f4', format_fortran='F', desc='Limb altitude used in inversion',
                          display_type='scalar', field_name='Inversion Limb Alt', fill_value=None, label_axis='Alt', bin_location=0.5,
                          units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=1,
                          notes='')

    ncfile.close()

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
def Get_lvl2_5_product(file_input='/home/jmakela/ICON_L1_FUV_SWP_20090320_v01r006.NC',
                       file_ancillary='/home/jmakela/ICON_L0P_FUV_ANCILLARY_2009-03-20_v001_r003.nc',
                       file_GPI=None,
                       file_output='/home/jmakela/ICON_L2-5_FUV-2009-03-20_00.18.54_v01_r00.nc',
                       limb = 150.,Spherical = True,
                       reg_method = 'Tikhonov',regu_order = 2, contribution ='RR'):
    '''
    Operational Code that reads Lvl1 file and creates the corresponding Lvl2.5
    INPUTS:
        file_input  - Input file path
        file_ancillary - Ancillary file path
        file_GPI - GPI file path; default to None which uses GPIs saved in PyGlow
        file_output - Full output file path for the resultant netCDF4 file
        limb        - Defines the lower bound that defines the limb [km]
        Spherical   - Flag indicating whether Sperical or Non-Spherical Earth is assumed [True, False]
        reg_method  - Regularization method selection [Tikhonov, MAP]
        regu_order  - Regularization Order [int] (0,1,2 Possible Values)
        contribution- contributions to the 135.6nm emission, RR: Radiative Recombination, RRMN: Radiative Recombination and Mutual Neutralization
    OUTPUT:
        Creates a NetCDF file on the desired file_outpu
        Returns 0 on success
    NOTES:
        This versions uses paths as input and output. That can change in case needed. Also, the lvl1 file that is used
        is elementary. No actual LVL1 file has been given yet.
    HISTORY:
        24-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        04-May-2017: Revised by Jonathan Makela to use L1 and ancillary data files (jmakela@illinois.edu)
        02-Jun-2017: Revised by Jonathan Makela to read GPI files (jmakela@illinois.edu)
        13-Jun-2017: Revised by Jonathan Makela to pass in full file output path (jmakela@illinois.edu)
    TODO:
        1) Variable names in L1 and ancillary data files may change
    '''

    try:
        # Open input Level 1 and ancillary NetCDF files
        data = netCDF4.Dataset(file_input,mode='r')
        ancillary = netCDF4.Dataset(file_ancillary,mode='r')

        if file_GPI is not None:
            gpi = netCDF4.Dataset(file_GPI,mode='r')

            # Read the geophysical indeces
            ap3 = gpi['ap3'][:]
            ap = gpi['ap'][:]
            year_day = gpi['year_day'][:]
            f107 = gpi['f107d'][:]
            f107a = gpi['f107a'][:]
        else:
            ap3 = None
            ap = None
            year_day = None
            f107 = None
            f107a = None

        # The tangent point WGS-84 coordinates at the center of the integration time
        FUV_TANGENT_WGS = ancillary.variables['ICON_ANCILLARY_FUV_TANGENTPOINTS_LATLONALT'][:,:,:,:]
        FUV_TANGENT_LATITUDES = ancillary.variables['ICON_ANCILLARY_FUV_TANGENTPOINTS_LATLONALT'][:,:,:,0]
        FUV_TANGENT_LONGITUDES = ancillary.variables['ICON_ANCILLARY_FUV_TANGENTPOINTS_LATLONALT'][:,:,:,1]
        FUV_TANGENT_ALTITUDES = ancillary.variables['ICON_ANCILLARY_FUV_TANGENTPOINTS_LATLONALT'][:,:,:,2]

        # The az/el of the look vector
        FUV_AZ = ancillary.variables['ICON_ANCILLARY_FUV_FOV_AZIMUTH_ANGLE'][:,:,:]
        FUV_ZE = ancillary.variables['ICON_ANCILLARY_FUV_FOV_ZENITH_ANGLE'][:,:,:]

        # The ICON WGS-84 location at the center of the integration time
        ICON_WGS = ancillary.variables['ICON_ANCILLARY_FUV_SC_WGS'][:,:]
        ICON_WGS84_LATITUDE = ancillary.variables['ICON_ANCILLARY_FUV_SC_WGS'][:,0]
        ICON_WGS84_LONGITUDE = ancillary.variables['ICON_ANCILLARY_FUV_SC_WGS'][:,1]
        ICON_WGS84_ALTITUDE = ancillary.variables['ICON_ANCILLARY_FUV_SC_WGS'][:,2]

        # Read the UTC of all measurements and store in a datetime variable
        temp = ancillary.variables['ICON_ANCILLARY_FUV_TIME_UTC']
        ANC_dn = []
        for d in temp:
            ANC_dn.append(parser.parse(d))
        ANC_dn = np.array(ANC_dn)

        # Get Data from file.
        # L1 data stores the individual stripes in different variables. Read them
        # all stripes and combine into a single variable
        mirror_dir = ['M9','M6','M3','P0','P3','P6']
        FUV_1356_IMAGE = np.zeros(np.shape(FUV_AZ))
        FUV_1356_ERROR = np.zeros(np.shape(FUV_AZ))
        for ind, d in enumerate(mirror_dir):
            FUV_1356_IMAGE[:,:,ind] = data.variables['ICON_L1_FUVA_SWP_PROF_%s' % d][:]
            FUV_1356_ERROR[:,:,ind] = data.variables['ICON_L1_FUVA_SWP_PROF_%s_ERROR' % d][:]

        # Get observation times from file and store in a datetime variable
        temp = data.variables['ICON_L1_FUVA_SWP_CENTER_TIMES']
        FUV_dn = []
        for d in temp:
            FUV_dn.append(parser.parse(d))
        FUV_dn = np.array(FUV_dn)

        # The FUV epoch
        FUV_EPOCH = data.variables['EPOCH'][:]

        # Get science mode
        FUV_mode = data.variables['ICON_L1_FUVA_SWP_MODE'][:]

        # Pre-allocate arrays for output variables
        FUV_ver = np.zeros(np.shape(FUV_1356_IMAGE))*np.nan
        FUV_sigma_ver = np.zeros(np.shape(FUV_1356_IMAGE))*np.nan
        FUV_Ne = np.zeros(np.shape(FUV_1356_IMAGE))*np.nan
        FUV_sigma_Ne = np.zeros(np.shape(FUV_1356_IMAGE))*np.nan
        FUV_h = np.zeros(np.shape(FUV_1356_IMAGE))*np.nan
        FUV_b = np.zeros(np.shape(FUV_1356_IMAGE))*np.nan
        FUV_hmF2 = np.zeros((len(FUV_dn),len(mirror_dir)))*np.nan
        FUV_sigma_hmF2 = np.zeros((len(FUV_dn),len(mirror_dir)))*np.nan
        FUV_NmF2 = np.zeros((len(FUV_dn),len(mirror_dir)))*np.nan
        FUV_sigma_NmF2 = np.zeros((len(FUV_dn),len(mirror_dir)))*np.nan
        FUV_quality = np.zeros((len(FUV_dn),len(mirror_dir)))*np.nan
        FUV_tangent_lat = np.zeros(np.shape(FUV_TANGENT_LATITUDES))*np.nan
        FUV_tangent_lon = np.zeros(np.shape(FUV_TANGENT_LONGITUDES))*np.nan
        FUV_tangent_alt = np.zeros(np.shape(FUV_TANGENT_ALTITUDES))*np.nan

        # Variables to keep track of the indices of the valid range of limb pixels
        min_li = 9999
        max_li = -9999
    except Exception as e:
        if hasattr(e, 'message'):
            print 'Error setting up reading files and setting up L2.5 inversion: (%s)' % e.message
        else:
            print 'Error setting up reading files and setting up L2.5 inversion'
        # Close the input netCDF files
        ancillary.close()
        data.close()
        gpi.close()
        return 1

    # Work on each individual stripe
    for stripe, d in enumerate(mirror_dir):
#        print stripe
        night_ind = []
        for ind, mode in enumerate(FUV_mode):
            # Check if we are in night mode
            if mode == 2:
                # We are in nighttime science mode, process the data
                # Save the index of this measurement
                night_ind.append(ind)

                # Vector with the space craft lat, lon, alt at the measurement time
                satlatlonalt = [ICON_WGS84_LATITUDE[ind],ICON_WGS84_LONGITUDE[ind],ICON_WGS84_ALTITUDE[ind]]

                # Only consider values above the limb
                limb_i = np.where(np.squeeze(FUV_TANGENT_ALTITUDES[ind,:,stripe])>=limb)
                bright = np.squeeze(FUV_1356_IMAGE[ind,limb_i,stripe]) # TODO: THIS IS A MASKED ARRAY. HOW TO HANDLE?
                bright = bright[::-1]
                FUV_b[ind,limb_i,stripe] = bright
                err = np.squeeze(FUV_1356_ERROR[ind,limb_i,stripe])
                err = err[::-1]
                h = np.squeeze(FUV_TANGENT_ALTITUDES[ind,limb_i,stripe])
                h = h[::-1]
                az = np.squeeze(FUV_AZ[ind,limb_i,stripe])
                az = az[::-1]
                ze = np.squeeze(FUV_ZE[ind,limb_i,stripe])
                ze = ze[::-1]

                # Check if we have a larger range of valid limb indices
                if max(limb_i[0]) > max_li:
                    max_li = max(limb_i[0])
                if min(limb_i[0]) < min_li:
                    min_li = min(limb_i[0])

                tan_la = np.squeeze(FUV_TANGENT_LATITUDES[ind,limb_i,stripe])
                tan_lo = np.squeeze(FUV_TANGENT_LONGITUDES[ind,limb_i,stripe])
                tan_al = np.squeeze(FUV_TANGENT_ALTITUDES[ind,limb_i,stripe])
                tanlatlonalt = [tan_la[::-1],tan_lo[::-1],tan_al[::-1]]

                # Time of observation
                # TODO: CHECK THAT ANC AND FUV DATES ARE THE SAME
                dn = FUV_dn[ind]

                try:
                    # Get the GPI needed to run MSIS
                    my_f107, my_f107a, my_apmsis = get_msisGPI(dn, year_day, f107, f107a, ap, ap3)

                    # Run the inversion
                    ver,Ne,h,Sig_ver,Sig_Ne = FUV_Level_2_Density_Calculation(bright,h,satlatlonalt,az,ze,
                                                       Sig_Bright = np.diag(err**2), weight_resid=False,
                                                       limb = limb,Spherical = Spherical, reg_method = reg_method,
                                                       regu_order = regu_order, contribution =contribution,dn = dn,
                                                       f107=my_f107, f107a=my_f107a,apmsis=my_apmsis)

                    # Save the values to output arrays
                    FUV_ver[ind,limb_i,stripe] = ver
                    FUV_sigma_ver[ind,limb_i,stripe] = np.sqrt(np.diagonal(Sig_ver))
                    FUV_Ne[ind,limb_i,stripe] = Ne
                    FUV_sigma_Ne[ind,limb_i,stripe] = np.sqrt(np.diagonal(Sig_Ne))
                    FUV_h[ind,limb_i,stripe] = h
                    FUV_tangent_lat[ind,limb_i,stripe] = tan_la[::-1]
                    FUV_tangent_lon[ind,limb_i,stripe] = tan_lo[::-1]
                    FUV_tangent_alt[ind,limb_i,stripe] = tan_al[::-1]

                    # Calculate hmF2 and NmF2
                    hm,Nm,sig_hm,sig_Nm = find_hm_Nm_F2(Ne,h,Sig_NE=Sig_Ne)
                    FUV_hmF2[ind,stripe] = hm
                    FUV_sigma_hmF2[ind,stripe] = sig_hm
                    FUV_NmF2[ind,stripe] = Nm
                    FUV_sigma_NmF2[ind,stripe]= sig_Nm
                    FUV_quality[ind,stripe] = 0
                except ImportError as error:
                    print "You don't have module {0} installed".format(error.message[16:])
                    return 2
                except ValueError as error:
                    # Just an error in the inversion
                    print 'ind: %d, stripe: %d, error: %s' %(ind,stripe,sys.exc_info()[0])
                    FUV_quality[ind,stripe] = 1
                except Exception as e:
                    if hasattr(e, 'message'):
                        print 'Unknown inversion error: (%s)' % e.message
                    else:
                        print "Unknown inversion error"

    # The master index for limb scan indexes
    master_li = np.arange(min_li,max_li+1)

    # Inversion is complete. Form the output dictionary
    L25_dict = {
        'FUV_EPOCH': FUV_EPOCH[night_ind],
        'FUV_CENTER_TIMES': data.variables['ICON_L1_FUVA_SWP_CENTER_TIMES'][:][night_ind],
        'FUV_START_TIMES': data.variables['ICON_L1_FUVA_SWP_START_TIMES'][:][night_ind],
        'FUV_STOP_TIMES': data.variables['ICON_L1_FUVA_SWP_STOP_TIMES'][:][night_ind],
        'FUV_dn': FUV_dn[night_ind],
        'FUV_integration_time': data.variables['ICON_L1_FUVA_SWP_INTEGRATION_TIME'][:][night_ind],
        'FUV_TANGENT_LAT': FUV_tangent_lat[night_ind,min_li:max_li+1,:],
        'FUV_TANGENT_LON': FUV_tangent_lon[night_ind,min_li:max_li+1,:],
        'FUV_TANGENT_ALT': FUV_tangent_alt[night_ind,min_li:max_li+1,:],
        'FUV_TANGENT_WGS': FUV_TANGENT_WGS[night_ind,:,:,:],
        'FUV_AZ': FUV_AZ[night_ind,min_li:max_li+1,:],
        'FUV_ZE': FUV_ZE[night_ind,min_li:max_li+1,:],
        'ICON_WGS': ICON_WGS[night_ind,:],
        'FUV_ver': FUV_ver[night_ind,min_li:max_li+1,:],
        'FUV_ver_error': FUV_sigma_ver[night_ind,min_li:max_li+1,:],
        'FUV_Ne': FUV_Ne[night_ind,min_li:max_li+1,:],
        'FUV_Ne_error': FUV_sigma_Ne[night_ind,min_li:max_li+1,:],
        'FUV_h': FUV_h[night_ind,:,:],
        'FUV_NmF2': FUV_NmF2[night_ind,:],
        'FUV_NmF2_error': FUV_sigma_NmF2[night_ind,:],
        'FUV_hmF2': FUV_hmF2[night_ind,:],
        'FUV_hmF2_error': FUV_sigma_hmF2[night_ind,:],
        'FUV_Quality': FUV_quality[night_ind,:],
        'reg_method': reg_method,
        'contribution': contribution,
        'limb': limb,
        'Spherical': Spherical,
        'regu_order': regu_order,
        'source_files': [file_input,file_ancillary],
        'param_fn': None,
    }

    # Close the input netCDF files
    ancillary.close()
    data.close()
    gpi.close()

    # Write the output
    print 'Writing netcdf %s' % file_output
    try:
        FUV_Level_2_OutputProduct_NetCDF(file_output, L25_dict)
    except Exception as e:
        if hasattr(e, 'message'):
            print 'Error writing netCDF output file: (%s)' % e.message
        else:
            print 'Error writing netCDF output file: %s (%s)' % (file_output, sys.exc_info()[0])

    print 'Successfully wrote %s' % file_output

    return 0

def get_GPI(dn, year_day, f107, f107a, ap, ap3):
    '''
    Operational Code to pull the geophysical indices from the standard GPI
    ancillary file.
    INPUTS:
        dn - datetime to be used
        year_day - vector containing the year/day of the GPIs (yyyddd)
        f107 - vector containing the f107 index
        f107a - vector containing the proxy for f107a average
        ap - vector containing the daily ap index
        ap3 - array containing the 3-hour ap index
    OUTPUT:
        returns the f107, f107a, ap, and ap3 for the requested time (dn)
    NOTES:
        Can be used in conjunction with get_msisGPI to generate the msis ap vector.
    HISTORY:
        02-Jun-2017: Written by Jonathan Makela (jmakela@illinois.edu)
    '''

    yd = dn.year*1000+dn.timetuple().tm_yday

    # Check bounds
    if (yd < year_day.min()) or (yd > year_day.max()):
        return np.nan, np.nan, np.nan, np.nan

    # Get the index into the year_day dimension
    i = np.argwhere(year_day == yd).flatten()[0]

    # Get the index into the 3 hour dimension
    j = dn.hour/3

    return f107[i], f107a[i], ap[i], ap3[j,i]

def get_msisGPI(dn, year_day, f107, f107a, ap, ap3):
    '''
    Operational Code to pull the geophysical indices from the standard GPI
    ancillary file and generate the ap vector required by MSIS.
    INPUTS:
        dn - datetime to be used
        year_day - vector containing the year/day of the GPIs (yyyddd)
        f107 - vector containing the f107 index
        f107a - vector containing the proxy for f107a average
        ap - vector containing the daily ap index
        ap3 - array containing the 3-hour ap index
    OUTPUT:
        returns the f107, f107a, msis ap for the requested time (dn)
    NOTES:
        As defined in the MSIS code the msis ap is a vector containing:
             (1) DAILY AP
             (2) 3 HR AP INDEX FOR CURRENT TIME
             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS
                 PRIOR   TO CURRENT TIME
             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS
                 PRIOR  TO CURRENT TIME
    HISTORY:
        02-Jun-2017: Written by Jonathan Makela (jmakela@illinois.edu)
    '''

    # Where to store the
    my_apmsis = np.zeros(7)

    # Get f107, f107a, daily ap, and 3-hour ap

    my_f107, my_f107a, my_apmsis[0], my_apmsis[1] = get_GPI(dn+datetime.timedelta(hours=0),year_day,f107,f107a,ap,ap3)

    # Get 3-hour ap for -3, -6, and -9 hours
    _, _, _, my_apmsis[2] = get_GPI(dn+datetime.timedelta(hours = -3),year_day,f107,f107a,ap,ap3)
    _, _, _, my_apmsis[3] = get_GPI(dn+datetime.timedelta(hours = -6),year_day,f107,f107a,ap,ap3)
    _, _, _, my_apmsis[4] = get_GPI(dn+datetime.timedelta(hours = -9),year_day,f107,f107a,ap,ap3)

    # Get average 3-hour ap for -12 to -33 hours
    temp_ap = 0
    for delta in range(-12,-33-1,-3):
        _, _, _, temp = get_GPI(dn+datetime.timedelta(hours = delta),year_day,f107,f107a,ap,ap3)
        temp_ap += temp
    my_apmsis[5] = temp_ap/8.

    # Get average 3-hour ap for -36 to -57 hours
    temp_ap = 0
    for delta in range(-36,-57-1,-3):
        _, _, _, temp = get_GPI(dn+datetime.timedelta(hours = delta),year_day,f107,f107a,ap,ap3)
        temp_ap += temp
    my_apmsis[6] = temp_ap/8.

    return my_f107, my_f107a, my_apmsis

def CreateSummaryPlot(file_netcdf, file_png, stripe=2, min_alt=None, max_alt=None,
                      min_dn=None, max_dn=None, min_ne=None, max_ne=None, min_dne=None, max_dne=None):
    '''
    Operational Code to generate the standard summary file for the nighttime FUV 2.5 dataproduct.
    INPUTS:
        file_netcdf - full path to the netCDF4 file
        file_png - full path to the png file to be generated
        stripe - which stripe to plot profiles for (default = 2)
        min_alt, max_alt - min/max range of altitude to plot (km). If no value passed, figure this
                           out from the range of altitudes in the file.
        min_dn, max_dn - min/max date range to plot (datenum). If no value passed, figure this
                         out from the range of altitudes in the file.
        min_ne, max_ne - min/max of electron density to plot (cm^-3). If no value passed, figure this
                         out from the range of altitudes in the file.
        min_dne, max_dne - min/max of electron density error to plot (cm^-3). If no value passed,
                           figure this out from the range of altitudes in the file.
    OUTPUT:
        generates a png file and saves it to disk
    NOTES:
    HISTORY:
        05-Jun-2017: Written by Jonathan Makela (jmakela@illinois.edu)
    '''

    # Use viridis if it exists
    try:
        cm = plt.get_cmap('viridis')
    except:
        cm = plt.get_cmap(plt.rcParams['image.cmap'])

    try:
        # open the netCDF file
        f = netCDF4.Dataset(file_netcdf)

        # Get variables from netCDF file
        dn = []
        for d in f.variables['ICON_L2_FUVA_SWP_CENTER_TIMES']:
            dn.append(parser.parse(d))
        dn = np.array(dn)

        X = np.transpose([dn,]*f.dimensions['Vertical'].size)
        Y = f.variables['ICON_L2_FUVA_TANGENT_ALT'][:,:,stripe]
        Z = f.variables['ICON_L2_FUVA_SWP_ELECTRON_DENSITY_ALTITUDE_PROFILE'][:,:,stripe]
        Ze = f.variables['ICON_L2_FUVA_SWP_ELECTRON_DENSITY_ERROR_ALTITUDE_PROFILE'][:,:,stripe]

        # Min/max ranges for plots
        if min_alt is None:
            min_alt = Y.min()
        if max_alt is None:
            max_alt = Y.max()

        if min_dn is None:
            min_dn = dn[0]
        if max_dn is None:
            max_dn = dn[-1]

        if min_ne is None:
            min_ne = Z.min()
        if max_ne is None:
            max_ne = Z.max()

        if min_dne is None:
            min_dne = np.log10(Ze.min())
        if max_dne is None:
            max_dne = np.log10(Ze.max())

        # generate the plot
        fig, axes = plt.subplots(nrows=5, sharex=True, figsize=(8,11))

        # The electron density estimates
        im = axes[0].pcolormesh(X,Y,Z,vmin=min_ne,vmax=max_ne,cmap=cm)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=15))
        axes[0].set_ylim([min_alt,max_alt])
        axes[0].set_title('Estimated Ne; Stripe #%d' % stripe)
        axes[0].set_ylabel('Altitude [km]')
        axes[0]

        # The electron density error estimates
        im2 = axes[1].pcolormesh(X,Y,np.log10(Ze),vmin=min_dne,vmax=max_dne,cmap=cm)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=15))
        axes[1].set_ylim([min_alt,max_alt])
        axes[1].set_title('Estimated Ne Error; Stripe #%d' % stripe)
        axes[1].set_ylabel('Altitude [km]')

        # The hmF2 estimates
        for stripe in range(0,6):
            Y = f.variables['ICON_L2_FUVA_SWP_HMF2'][:,stripe]
            Ye = f.variables['ICON_L2_FUVA_SWP_HMF2_ERROR'][:,stripe]
            axes[2].plot(dn,Y,'-')
            axes[2].fill_between(dn,Y-Ye,Y+Ye,alpha=.25)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=15))
        plt.xlim([min_dn,max_dn])
        axes[2].set_ylim([min_alt,max_alt])
        axes[2].set_title('Estimated hmF2')
        axes[2].set_ylabel('Altitude [km]')

        # The NmF2 estimates
        for stripe in range(0,6):
            Y = f.variables['ICON_L2_FUVA_SWP_NMF2'][:,stripe]
            Ye = f.variables['ICON_L2_FUVA_SWP_NMF2_ERROR'][:,stripe]
            axes[3].plot(dn,Y,'-')
            axes[3].fill_between(dn,Y-Ye,Y+Ye,alpha=.25)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=15))
        plt.xlim([min_dn,max_dn])
        axes[3].set_title('Estimated NmF2')
        axes[3].set_ylabel('Ne [cm^-3]')

        # The quality flag
        for stripe in range(0,6):
            quality = f.variables['ICON_L2_FUVA_QUALITY'][:]
            axes[4].plot(dn,quality[:,stripe],'.',label=stripe)
        axes[4].set_ylim([-.25,2])
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=15))
        axes[4].set_title('Inversion quality flag')
        axes[4].set_xlabel('UTC')
        lgd = plt.legend(ncol=6,numpoints=1,loc=9, bbox_to_anchor=(0.6, -0.23))

        # labels
        axes[4].annotate('(%s, %s)' % (f.variables['ICON_L2_1356_CONTRIBUTION'][:],
                                      f.variables['ICON_L2_REGULARIZATION_METHOD'][:]),
                         xy=(1, 0.05), xytext=(0, 10),
                         xycoords=('axes fraction', 'figure fraction'),
                         textcoords='offset points',
                         ha='right', va='top')

        an = axes[0].annotate('%s' % (dn[0].strftime('%d %b %Y')),
                         xy=(.5, .95), xytext=(0, 10),
                         xycoords=('axes fraction', 'figure fraction'),
                         textcoords='offset points',
                         size=12,ha='center', va='top')

        # Make some room for the colorbar
        fig.subplots_adjust(left=0.07, right=0.87)

        # Add the colorbar outside...
        box = axes[0].get_position()
        pad, width = 0.02, 0.02
        cax = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
        fig.colorbar(im, cax=cax,format='%.0e',label='Ne [cm^-3]',extend='max')

        box = axes[1].get_position()
        pad, width = 0.02, 0.02
        cax = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
        fig.colorbar(im2, cax=cax, label='log10(Ne_error) [cm^-3]', extend='max')
        #fig.tight_layout()
        fig.savefig(file_png, bbox_extra_artists=(lgd,an), bbox_inches='tight')

        # Close the netCDF file
        f.close()
    except:
        print 'Error generating L2.5 summary file'
        return 1

    return 0
