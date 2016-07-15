import ICON as ic
from pyglow import pyglow
from datetime import datetime, timedelta
from scipy import optimize
from scipy.linalg import sqrtm
import numpy as np
import math
from scipy.io import netcdf
from time import gmtime, strftime
import FUV_Distance_Matrix as dmat
reload(dmat)

import matplotlib.pyplot as plt

import multiprocessing
from multiprocessing import Pool
from sys import exit


from scipy.interpolate import interp1d


def decs(x, pos, sc = 1e-5):
    'The two args are the value and tick position'
    return '%1.f' % (x*sc)

def guassian_elinimation(S,Rayl,pseudo = 0):
    '''
    Given a distance matrix S and a Brightness profile, guassian elimination is used to calculate the VER
    INPUTS:
        S -  Distance matrix
        Rayl -  Brightness profile (Rayleigh)
        pseudo - flag to indicate if pseudoinverse is used [0: Gaussian Elimination, 1:Pseudoinverse || DEFAULT: 0]
    OUTPUT:
        VER - volume emission rate (ph/cm^-3)
    NOTES:
        This works mainly for cases where the Brightness profile doesnt contain noise and matrix S is square
    HISTORY:
        17-Sep-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        16-Apr-2016: Added pseudoinvese option

    CALLS:
        linalg.solve
    '''
    
    if np.size(S,1)==np.size(S,0):
        Ip = np.zeros(np.size(Rayl))
        if pseudo == 1:
            S_pinv = np.linalg.pinv(S)
            Ip = np.dot(S_pinv, Rayl)
        else:
            Ip = np.linalg.solve(S,Rayl)

        # Remove any negative values
        for i in range(0,np.size(Ip)):
            if Ip[i]<0:
                Ip[i] = 0
        return Ip
    else:
        print 'Distance matrix is not square'
        
        return 0

#def calc_electron_density(Ip,O,cont=1):  
def calc_electron_density(Ip,satlatlonalt,h,dn,cont=1):
    '''
    Given a VER profile the electron density is calculated. The physical process for the VER must be defined
    INPUTS:
        Ip   - Volume Emission Rate profile (ph/cm^-3)
        O    - Oxygen profile on the tangent altitude (1/cm^3/s)
        cont - contribution factors, cont=1 -> RR+MN, 2-> RR
    OUTPUTS:
        NE_est - The estimated Electron Density profile (1/cm^3/s)       
    NOTES:
        For the inversion Ne = Op is assumed
    HISTORY:
        17-Sep-2015: Written by Dimitrios Iliou (iliou2@illinois.edu) 
    CALLS:
    
    '''
    
                
    if np.size(np.shape(Ip))==2:
        Ip = Ip[:,0]
    elif np.size(np.shape(Ip))>2:
        raise Exception('Invalid data vector')

    b1356 = 0.54    # yield parameter (unitless)
    a1356 = 7.3e-13 # radiative recombination rate (cm^3/s)
    # consider that to be constant through altitude, normally it should change with temperature
    k1 = 1.3e-15    # radiative attachment rate (cm^3/s)
    k2 = 1e-7       # ion-ion neutralization rate (cm^3/s)
    k3 = 1.4e-10    # ion-atom neutralization rate (cm^3/2)
    
    if cont==1:
        O = np.zeros(len(Ip))
        for i in range(0,len(h)):
            pt = pyglow.Point(dn, satlatlonalt[0], satlatlonalt[1], h[i])

            pt.run_msis()

            O[i] = pt.nn['O']  
        
        a0 = a1356/O
        b0 = a1356*k3/k2+b1356*k1
        c0 = -Ip/O
        d0 = -Ip*k3/k2

        a1 = b0*c0/(6.0*a0**2)-b0**3./a0**3/27.0-d0/a0/2.0
        
        b1 = c0/a0/3.0-b0**2./a0**2/9.0

        c1 = -b0/a0/3.0;
        
        d1 = a1/np.sqrt((-b1)**3)
        
        d1[np.where(d1<-1.)]=-1.
        d1[np.where(d1>1.)]=1.
            
        # used arccos instead of MATLAB acos
        NE_est = c1+2.0*np.sqrt(-b1)*np.cos(np.arccos(d1)/3.0)
    else:    
        NE_est = np.sqrt(Ip/a1356)

    NE_est[np.isnan(NE_est)] = 0
    return NE_est


# Create vector containing regularization parameters
def create_alpha_values(A,npoints = 100.):
    '''
    Given a distance matrix S a number of points the regularization parameter vector is created. 
    INPUTS:
        S       -  Distance matrix
        npoints -  Number of points for the regularization parameter
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

# Find optimum reg param using Maximum Curvature of L-Curve
def Maximum_Curvature(residual,x_lamda,reg_param):
    '''
    Given the residual and the seminorm of the L2 norm problem the corner of the L-Curve is calculated by
    finding the maximum curvature by finding the circle with the smaller radius testing for every three
    consecutive points
    INPUTS:
        residual   - The residual norm of our minimization problem
        x_lamda    - The seminorm of ||Lx||
        reg_param  - Vector containing all the regularization parameters
    OUTPUT:
        reg_corner - Maximum curvature point - Best regularization parameter
        kappa      - Vector containing the calculated curvature [Mainly for testing reasons]
    NOTES:
    HISTORY:
        15-May-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        31-Aug-2015 : Curvature return added
    CALLS:
    '''
    
    # Maximum Curvature
    #transform rho and eta into log-log space
    x=np.log(residual);
    y=np.log(x_lamda);

    # the series of points used for the triangle/circle
    x1 = x[:-2]
    x2 = x[1:-1]
    x3 = x[2:]
    y1 = y[:-2]
    y2 = y[1:-1]
    #y3 = y[:-2]
    y3 = y[2:]

    # the side lengths for each triangle
    a1 = np.sqrt((x3-x2)**2+(y3-y2)**2);
    b1 = np.sqrt((x1-x3)**2+(y1-y3)**2);
    c1 = np.sqrt((x2-x1)**2+(y2-y1)**2);

    # semi-perimeter
    s1=(a1+b1+c1)/2;

    # the radius of each circle
    R=(a1*b1*c1)/(4*np.sqrt((s1*(s1-a1)*(s1-b1)*(s1-c1))));

    # The curvature for each estimate for each value which is
    # the reciprocal of its circumscribed radius. Since there aren't circles for 
    # the end points they have no curvature
    kappa = np.zeros(np.size(R)+2)
    kappa[0] = 0
    kappa[-1] = 0
    kappa[1:-1] = 1/R

    ireg_corner=np.argmax(abs(kappa[1:-1]));
    reg_corner=reg_param[ireg_corner];

    return reg_corner, kappa

def Maximum_Curvature_gradiens(residual,x_lamda,reg_param):
    '''
    Given the residual and the seminorm of the L2 norm problem the corner of the L-Curve is calculated by
    finding the maximum curvature by finding the maximum second derivative
    INPUTS:
        residual   - The residual norm of our minimization problem
        x_lamda    - The seminorm of ||Lx||
        reg_param  - Vector containing all the regularization parameters
    OUTPUT:
        reg_corner - Maximum curvature point - Best regularization parameter
        
    NOTES:
    HISTORY:
        17-Sep-2015: Written by Dimitrios Iliou (iliou2@illinois.edu) and Jianqi Qin (jianqi@illinois.edu)
    CALLS:
    '''
    # Maximum Curvature
    #transform rho and eta into log-log space
    d1=np.log(residual);
    m1=np.log(x_lamda);
    
    grad1 =  np.gradient(m1)/np.gradient(d1)
    grad2 = np.gradient(grad1)/np.gradient(d1)
    
    ireg_corner=np.argmax(grad2)
    reg_corner=reg_param[ireg_corner];
    
    return reg_corner

# Tikhonov Regularization
def Tikhonov(A,b,deg,reg_param=0.,method=1,ireg=1):
    '''
    Tikhonov Regularization function. Solves the minimization problem.
    INPUTS:
        A         - Distance matrix
        b         - Brightness Profile (Rayleigh)
        deg       - The degree of Tikhonov regularization [valid is 0,1,2]
        reg_param - regularization parameter. If zero it calculates the regularization parameters using S. Otherwise uses
                    values given (scalar or matrix)
        method    - Flag indicating which method for the Maximum Curvature is going to be used [0:Maximum Curvature 1: Second order derivative]
        ireg      - Flag indicating if we want to print the parameter chosen when executing [1:print 0:noprint]
    OUTPUT:
        sol       - Estimated electron density (cm^-3)
        reg_corner- Value of the regularization parameter chosen
        residual  - Vector containing the residual norm for all the values of the regularization parameter
        seminorm  - Vector containing the seminorm norm for all the values of the regularization parameter
        curvature - In case the smaller circle maximum curvature is chosen the curvature vector is returned
    NOTES:
    HISTORY:
        15-May-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:
        -create_alpha_values
        -calc_solution
        -Maximum_Curvature
        -Maximum_Curvature_gradiens
    '''
    # Make a change to include the shape of reg_param so we dont have to calculate it each time
    # Check if reg_param is a vector or not
    if np.size(reg_param) == 1:
        # if its not a vector only if zero it will calculate the vector. Otherwise it will use the single value given. 
        if reg_param == 0:
            reg_param = create_alpha_values(A)
    #if reg_param == 0:
        #reg_param = create_alpha_values(A)

    residual = np.zeros(len(reg_param))
    seminorm = np.zeros(len(reg_param))
    
    L = get_rough_matrix(len(b),deg)
    # This is the big bottleneck on my code
    for i in range(0,len(reg_param)):
        sol = calc_solution(A,b,reg_param[i],L) 
        residual[i] = np.linalg.norm(A.dot(sol)-b)
        seminorm[i] = np.linalg.norm(L.dot(sol))
        
    if method == 0:
        reg_corner, curvature = Maximum_Curvature(residual,seminorm,reg_param)
    else:
        reg_corner = Maximum_Curvature_gradiens(residual,seminorm,reg_param)
        curvature = 0
        
    sol = calc_solution(A,b,reg_corner,L) 
    
    if ireg==1:
        print 'Alpha parameter chose:.%4f' %(reg_corner)

    return sol,reg_corner,residual,seminorm,curvature

def Reg_param_repet(A,L,b,reg_param):
    '''
    Takes as input the parameters for the problem and we want to get the residual and seminorm quicker using multiple cores
    '''
    sol = calc_solution(A,b,reg_param,L) 
    residual = np.linalg.norm(A.dot(sol)-b)
    seminorm = np.linalg.norm(L.dot(sol))
    
    return residual,seminorm

#Multiprocessing
def Reg_param_repet_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return Reg_param_repet(*a_b)

def Tikhonov_mp(A,b,deg,reg_param=0.,method=1,ireg=1):
    '''
    Tikhonov Regularization function. Solves the minimization problem.
    INPUTS:
        A         - Distance matrix
        b         - Brightness Profile (Rayleigh)
        deg       - The degree of Tikhonov regularization [valid is 0,1,2]
        reg_param - regularization parameter. If zero it calculates the regularization parameters using S. Otherwise uses
                    values given (scalar or matrix)
        method    - Flag indicating which method for the Maximum Curvature is going to be used [0:Maximum Curvature 1: Second order derivative]
        ireg      - Flag indicating if we want to print the parameter chosen when executing [1:print 0:noprint]
    OUTPUT:
        sol       - Estimated electron density (cm^-3)
        reg_corner- Value of the regularization parameter chosen
        residual  - Vector containing the residual norm for all the values of the regularization parameter
        seminorm  - Vector containing the seminorm norm for all the values of the regularization parameter
        curvature - In case the smaller circle maximum curvature is chosen the curvature vector is returned
    NOTES:
    HISTORY:
        15-May-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:
        -create_alpha_values
        -calc_solution
        -Maximum_Curvature
        -Maximum_Curvature_gradiens
    '''
    try:
        # Make a change to include the shape of reg_param so we dont have to calculate it each time
        # Check if reg_param is a vector or not
        if np.size(reg_param) == 1:
            # if its not a vector only if zero it will calculate the vector. Otherwise it will use the single value given. 
            if reg_param == 0:
                reg_param = create_alpha_values(A)
        #if reg_param == 0:
            #reg_param = create_alpha_values(A)

        residual = np.zeros(len(reg_param))
        seminorm = np.zeros(len(reg_param))

        L = get_rough_matrix(len(b),deg)

        job_args = [(A,L,b,reg_param[i]) for i in range(0,len(reg_param))]
        N = multiprocessing.cpu_count()

        # Create the pool.  Be nice.  Don't use all the cores!
        pool = Pool(processes=16)


        results = pool.map(Reg_param_repet_star,job_args)


        for i in range(0,len(results)):
            residual[i] = results[i][0]
            seminorm[i] = results[i][1]


        pool.close()
        pool.join()


        if method == 0:
            reg_corner, curvature = Maximum_Curvature(residual,seminorm,reg_param)
        else:
            reg_corner = Maximum_Curvature_gradiens(residual,seminorm,reg_param)
            curvature = 0

        sol = calc_solution(A,b,reg_corner,L) 

        if ireg==1:
            print 'Alpha parameter chose:.%4f' %(reg_corner)

    except (KeyboardInterrupt,SystemExit,ZeroDivisionError,BaseException) as inst :
        
        if 'pool' in vars():
            pool.terminate()
            
        #print "You cancelled the program!"
        print type(inst)
        print inst
        exit(1) 

    except Exception:
    
        print "Something Happened :("
        print type(inst)
        print inst
        exit(1)
    
    return sol,reg_corner,residual,seminorm,curvature

# Find optimum regularization parameter using GCV
def gcv(A,b,deg,reg_param=0,ireg=1):
    '''
    Generalized cross-validation method for regularization
    INPUTS:
        A         - Distance matrix
        b         - Brightness Profile (Rayleigh)
        deg       - The degree of Tikhonov regularization [valid is 0,1,2]
        reg_param - regularization parameter. If zero it calculates the regularization parameters using S. Otherwise uses
                    values given (scalar or matrix)
        ireg      - Flag indicating if we want to print the parameter chosen when executing [1:print 0:noprint]
    OUTPUT:
        sol       - Estimated electron density (cm^-3)
    NOTES:
    HISTORY:
        15-May-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:
        -calc_solution
    '''

    # Make a change to include the shape of reg_param so we dont have to calculate it each time
    # Check if reg_param is a vector or not
    if np.size(reg_param) == 1:
        # if its not a vector only if zero it will calculate the vector. Otherwise it will use the single value given. 
        if reg_param == 0:
            reg_param = create_alpha_values(A)
    #if reg_param == 0:
        #reg_param = create_alpha_values(A)

    m = len(b)
    I = np.eye(np.size(A,0))
    g = np.zeros(len(reg_param))
    
    L = get_rough_matrix(len(b),deg)

    for i in range(0,len(reg_param)):
        A_sharp =  np.linalg.inv(((A.T).dot(A)+reg_param[i]**2*(L.T).dot(L))).dot(A.T)
        '''
        C = np.concatenate((A, reg_param[i]*L))
        
        # Constuct vector d = [b;0]
        d_temp = np.zeros(len(b))
        d = np.concatenate((b, d_temp))
       
        # Might need to add augmentation here
        sol,rnorm = optimize.nnls(C,d)       
        '''
        sol = calc_solution(A,b,reg_param[i],L) 
        residual = np.linalg.norm(b - A.dot(sol))**2 
        #residual[i] = rnorm
        denom = np.trace(I - A.dot(A_sharp))**2

        g[i] = m*residual/denom

    ireg_corner = np.argmin(g)
    reg_corner=reg_param[ireg_corner];

    if ireg==1:
        print reg_corner

    sol = calc_solution(A,b,reg_corner,L) 

    return sol
    
# Do regularization using MAP Estimation
def MAP_Estimation(S,Bright,VER_mean=0.,VER_reshape=0.):
    '''
    Bayesian MAP Estimation Method for regularization
    INPUTS:
        S         - Distance matrix
        Bright    - Brightness Profile (Rayleigh)  
        VER_mean  - Prior mean
        VER_reshape  - Prior variance matrix  
        path      - Path where the prior file is located (default prior IRI)  '/home/dimitris/public_html/Datafiles/All_files_simulations/Priors/Prior_mean_variance_26042016.npz
    OUTPUT:
        VER       - Estimated VER (ph/cm^-3/s)
    NOTES:
        ==> VER Prior is loaded for default path. This can be changed at will inside the function
    HISTORY:
        18-Apr-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:

    '''

    # Check whether or not the prior is there. If not raise exception
    if not hasattr(VER_reshape, "__len__"):
        if (VER_reshape==0):
            print 'PRIOR IS MISSING'
            raise Exception
    

    # Interpolate brightness to remove zeros and calculate the covariance matrix
    x = range(0,len(Bright))
    xp = np.where(Bright!=0)[0]
    Bright_interp= np.interp(x, xp, Bright[xp])
    
    
    # To calculate the variance of the noise we just take a single Brightness profile and we find the variance. More accurate estimation can be performed if we sum the 6 stripes before feeding them here. We also need to have the sensitivity of the instrument to calculate the number of counts. 
    #W_var = np.sqrt(Bright * 0.0873);
    #W_var = Bright * 0.0873
    W_var = Bright_interp * 0.0873 * 12
    # Create Diagonal matrix of the VER variances.
    #Rw = np.diag(W_var)
    W_var = W_var/((12*0.0873)**2)
    
    #plt.errorbar(x,Bright_interp, yerr=np.sqrt(W_var))
    
    #maxvar = max(W_var)
    #W_var[:] = maxvar
    
    # Create Diagonal matrix of the Brightness variances.
    Rw = np.diag(W_var)

     # Load VER priors
    #data = np.load(path)

    # VER_var = data['arr_3']
    #Rx_ver = VER_var[:len(Bright),:len(Bright)]
    
    #VER_mean = data['arr_0']
    #VER_reshape = data['arr_1']

    VER_reshape = VER_reshape[:len(Bright),:]
    VER_mean =  VER_mean[:len(Bright)]
    
    VER_var = np.cov(VER_reshape);
    Rx_ver = VER_var

    covd12=sqrtm(np.linalg.inv(Rw))
    covm12=sqrtm(np.linalg.inv(Rx_ver))
    #covd12=np.linalg.inv(sqrtm(Rw))
    #covm12=np.linalg.inv(sqrtm(Rx_ver))
    
    A = np.concatenate((covd12.dot(S), covm12), axis=0)
    rhs= np.concatenate((covd12.dot(Bright), covm12.dot(VER_mean)), axis=0)

    
    #VER = np.linalg.pinv(A).dot(rhs)
    VER,_,_,_ = np.linalg.lstsq(A,rhs)
    # Check if the values of the VER are negative. 

    
    # Close the data
    #data.close()
    
    # Calculations for the MAP Estimation
    #Rx_map = np.linalg.inv(np.linalg.inv(Rx_ver) +np.dot(np.transpose(S), np.dot(np.linalg.inv(Rw),S)))
    #Rx_map = np.linalg.inv(np.linalg.inv(Rx_ver) +np.dot(np.transpose(S), np.dot(np.linalg.inv(Rw),S)))


    #b = np.dot(np.dot(np.transpose(S),np.linalg.inv(Rw)),Bright)
    #b = np.dot(np.transpose(S),np.dot(np.linalg.inv(Rw),Bright))
    #b = np.dot(np.transpose(S),np.dot(np.linalg.inv(Rw),(Bright-np.dot(S,VER_mean))))

    #VER = np.dot(Rx_map,b)+VER_mean
    
    # Check if the values of the VER are negative. 
    for i in range(0,np.size(VER)):
            if VER[i]<0:
                VER[i] = 0
    
    return VER

# Calculate roughening matrices 1st and 2nd degree
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
        raise Exception('Invalid degree')

# Calculate regularized inverse rolution
def calc_solution(A,b,lamda,L):
    '''
    Calculates the solution using non-negative least square for the regularization problem
    INPUTS:
        A     - Distance matrix
        b     - Brightness Profile (Rayleigh)
        lamda - The selected regularization parameter
        L     - Roughening matrix
    OUTPUT:
        sol       - Estimated electron density (cm^-3)       
    NOTES:
    HISTORY:
        15-May-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:
    '''
    C = np.concatenate((A, lamda*L))
    if np.size(np.shape(b))==2:
        b = b[:,0]
    elif np.size(np.shape(b))>2:
        raise Exception('Invalid data vector')
        
    # Constuct vector d = [b;0]
    d_temp = np.zeros(len(L))
    d = np.concatenate((b, d_temp))

    #print "Condition Number of Distacne Matrix:", np.linalg.cond(A)
    #print "Condition Number of Augmented Matrix:", np.linalg.cond(C)
    sol,rnorm = optimize.nnls(C,d)
    #sol = guassian_elinimation(C,d)
    
    return sol

def find_hm_Nm_F2(NE,rbot):
    '''
    Calculates the NmF2 and hmF2 of a given electron density profile
    INPUTS:
        NE   - Electron density profile [cm^{-3}]
        rbot - altitude vector [km]
    OUTPUT:
        hmF2 -  altitude of the maximum intesity value of the Ne profile
        NmF2 -  peak intensity value of the altitude profile
    NOTES:
        NONE
    HISTORY:
        06-Jun-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    indexf2= NE.argmax()

    hmF2 = rbot[indexf2]
    NmF2 = NE[indexf2]
    
    return hmF2,NmF2

# Function to calculate the NmF2 and hmF2 difference.
def find_hm_Nm_F2_Diff(Ne_e,rbot_e,Ne,rbot,abs_flag=0.,interp = 1,int_increase = 5.,kind ='cubic'):
    '''
    Calculates the F2 peak for a single vector
    INPUTS:
        NE_e     - Estimated Electron density profile [cm^{-3}]
        rbot_e   - altitude vector for estimated electron density profile [km]
        Ne       - Original Electron density profile [cm^{-3}]
        rbot     - altitude vector for original electron density profile [km]
        abs_flag - flag that indicates if the return differences will be the absolute values [0:Regular Difference (default), 1:Absolute Difference]
        interp   - flag indicating if interpolation will be used around the hmF2 region in order to increase accuracy
        int_increase - Number indicating the amount of points to be created. => len(hmf2Region)* int_increase
        kind     - Goes as input to the interp1d function and declares the type of interpolation
    OUTPUT:
        Hmf2   -  mean altitude of the maximum intesity values of the Ne profile
        Nmf2   -  mean peak intensity value of the altitude profiles
        Hmf2_s -  standard deviation of altitude of the maximum intesity values of the Ne profile [In case of a single profile return 0]
        Nmf2_s -  standard deviation of peak intensity value of the altitude profiles [In case of a single profile return 0]
    NOTES:
    HISTORY:
        18-Aug-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        18-Sep-2015: Modified function to work for single or multiple Ne profiles
        29-Apr-2016: Added interpolation function
    '''

    if np.count_nonzero(Ne_e)==0:
        print "Estimated Electron Density zeros"
        raise Exception

    if interp==1:
        Ne,rbot = hmF2_region_interpolate(Ne,rbot,int_increase = int_increase,kind =kind) 
    
    hmF2o,Nmf2o = find_hm_Nm_F2(Ne,rbot)
   
    if Ne_e.ndim==1:
        Nmdift = np.zeros(np.size(Ne_e,0))
        Hmdift = np.zeros(np.size(Ne_e,0))
        
        if interp==1:
            Ne_e,rbot_e = hmF2_region_interpolate(Ne_e,rbot_e,int_increase = int_increase,kind =kind) 
    
        Hmf2_t,Nmf2_t = find_hm_Nm_F2(Ne_e,rbot_e)

        if abs_flag == 0:
            Nmdift = (Nmf2_t - Nmf2o)/(Nmf2o)
            Hmdift =Hmf2_t-hmF2o
        else:
            Nmdift = abs((Nmf2_t - Nmf2o)/(Nmf2o))
            Hmdift = abs(Hmf2_t- hmF2o)

        return Hmdift,Nmdift,0,0
    else:
        
        Nmdift = np.zeros(np.size(Ne_e,1))
        Hmdift = np.zeros(np.size(Ne_e,1))

        for i in range(0,np.size(Ne_e,1)):
            
            if interp==1:
                Ne_e_interp,rbot_e_interp  = hmF2_region_interpolate(Ne_e[:,i],rbot_e,int_increase = int_increase,kind =kind) 
            else:
                Ne_e_interp = Ne_e[:,i]
                rbot_e_interp = rbot_e
                
            Hmf2_t,Nmf2_t = find_hm_Nm_F2(Ne_e_interp,rbot_e_interp)

            if abs_flag == 0:
                Nmdift[i] = (Nmf2_t - Nmf2o)/(Nmf2o)
                Hmdift[i] = Hmf2_t-hmF2o
            else:
                Nmdift[i] = abs((Nmf2_t - Nmf2o)/(Nmf2o))
                Hmdift[i] = abs(Hmf2_t- hmF2o)

        Hmf2 = np.mean(Hmdift)
        Nmf2 = np.mean(Nmdift)

        Hmf2_s = np.std(Hmdift)
        Nmf2_s = np.std(Nmdift)
        
        return Hmf2,Nmf2, Hmf2_s,Nmf2_s    

def hmF2_region_interpolate(Ne,rbot,int_increase = 5.,kind ='cubic'):
    '''
    Interpolates the values around the hmF2 region in order to increase altitude vector interval which helps in minimizing the hmF2 error
    INPUTS:
        Ne       - Original Electron density profile [cm^{-3}]
        rbot     - altitude vector for original electron density profile [km]
        int_increase - Number indicating the amount of points to be created. => len(hmf2Region)* int_increase
        kind     - Goes as input to the interp1d function and declares the type of interpolation
    OUTPUT:
        NE_inter   -  Interpolated altitude profile for the estimated electron density
        xnew       -  Interpolated altitude vector for the true electron density
    NOTES:
    HISTORY:
        29-Apr-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    
    nm_true = np.argmax(Ne)
    if nm_true>5 and nm_true<len(Ne)-5:
        y_true = Ne[nm_true-5:nm_true+5]
        x_true = rbot[nm_true-5:nm_true+5]

        f2 = interp1d(x_true, y_true, kind='cubic')
        
        xnew = np.linspace(x_true[0], x_true[-1], num=len(x_true)*5, endpoint=True)
        
        NE_inter = f2(xnew)
        
        return NE_inter,xnew
    else:
        print 'Ne maximum ill defined - No interpolation performed'
        
        return Ne,rbot
            
# From this on the function focus on creating the solution and writing the NetCDF file
def FUV_Level_2_OutputProduct_Calculation(Bright,h,satlatlonalt,az,ze,O,cont =1, regu=1,regu_order = 2,method=1,ireg=1,S_mp = 1,Spherical =1,low_ta = 150.,VER_mean=0.,VER_reshape=0,dn = datetime(2012, 9, 26, 20, 0, 0) ):
    '''
    Within this function, given the LVL1 Input file the VER and Ne profiles for the tangent altitude point are
    calculated
    INPUTS:
        Bright      - Brightness Profile [Rayleigh]
        h           - Altitude profile [km]
        satlatlonalt- Satellite Coordinates [degrees,degrees,km]
        az          - Azimuth [degrees]
        ze          - Zenith [degrees]
        O           - Oxygen profile ( Assume known. We can also pull it from MSIS)
        cont        - Contribution, 1=> RR+MN, 2=>RR ( This will be removed propably) [int]
        regu        - Choose regularization, 0=> Gaussian Elimination, 1=> Tikhonov L-Curve, 2=> GSVD, 3==> Bayesian MAP estimation (NEEDS PRIOR) [int]
        regu_order  - Regularization Order [int] (0,1,2 Possible Values)
        method      - Flag indicating which method for the Maximum Curvature is going to be used [0:Maximum Curvature 1: Second order derivative] or Flag indicating Gaussian Ellimination or Pseudoinverse [0:Gaussian Ellimination 1: Pseudoinverse]
        ireg        - Flag indicating if we want to print the parameter chosen when executing [1:print 0:noprint]
        S_mp        - Flag indication if multiprocessing will be used for the calculation of the Distance matrix for the Non-Spherical Earth
        Spherical   - Flag indication if the assumption for the Earth is spherical [1:Non-Spherical, 0: Spherical] - is needed for the calculation of the distance matrix
        low_ta       - Lower tangent altitude that we need to consider for our measurements (km)[default = 150 km :limb, for sublimb we can put -500 to cover all zenith] 
        VER_mean     - Prior mean
        VER_reshape  - Prior variance matrix  
        path         - Path that contains the prior file. Default value is the IRI file  '/home/dimitris/public_html/Datafiles/All_files_simulations/Priors/Prior_mean_variance_26042016.npz' 
    OUTPUT:
        VER         - Volume Emission Rate tangent altitude profile
        Counts      - Electron Density tangent altitude profile
    NOTES:
        The distance matrix S is calculated assuming WGS84 coordinate system. 
        => Further work might be need to ensure Right dimentionality of bright h etc
        => Further Regularization methods will be introduced.
        => Wont work for multiple noise levels because S will be calculated for every iteration 
        => Thoughts on whether the NetCDF file will be created here or outside
    HISTORY:
        14-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        24-Feb-2016: Added Spherical earth flag
        04-Apr-2016: Added Bayesian Method
        18-Apr-2016: Added Lowest Tangent altitude flag and find_mid_cell dunction call.

    '''
    
    
    # Calculate the mid tangent altitude for each cell
    h,rbot,_,h_loc_bot,h_loc = find_mid_cell(satlatlonalt,ze,az,Spherical,low_ta)

    # We need to do that so we dont calculate extra length of certain paths for lower altitudes
    #az = az[0:len(h)]
    #ze = ze[0:len(h)]

    if Spherical==1:
        if S_mp ==1:
            S,rmid = dmat.Calculate_D_Matrix_WGS84_mp(satlatlonalt,az,ze)
        elif S_mp == 0:
            S = dmat.Calculate_D_Matrix_WGS84(satlatlonalt,az,ze)
        else:
            raise Exception('Not valid S multiprocessing choice')
    else:
        S,rmid,_,_ = dmat.create_cells_Matrix_spherical_symmetry(np.deg2rad(ze),satlatlonalt[2])

        
    # Concatinate the matrix to include only lowest tangent altitude input
    print 'Distance Matrix Calculated'
    Bright= Bright[0:len(h)]
    S = S[0:len(h),0:len(h)]
    
    if regu ==0: # Gaussian Elimination
        print 'Gaussian Elimination Chosen'
        VER = guassian_elinimation(S,Bright,method)
        #counts= calc_electron_density(VER,O,cont)
        counts= calc_electron_density(VER,satlatlonalt,h,dn,cont)
        print 'VER & Ne profiles calculated'
    elif regu == 1: # LCurve
        print 'Tiknonov - Lcurve Chosen'
        VER,_,_,_,_ = Tikhonov(S,Bright,regu_order,0,method,ireg)
        counts = calc_electron_density(VER,satlatlonalt,h,dn,cont)
        print 'VER & Ne profiles calculated'
    elif regu == 2: # GSVD
        print 'Tikhonov-GCV Chosen'
        VER = gcv(S,Bright,regu_order,0,ireg)
        counts = calc_electron_density(VER,satlatlonalt,h,dn,cont)
        print 'VER & Ne profiles calculated'
    elif regu ==3: # Bayesian
        print 'MAP Estimation Method Chose'
        VER = MAP_Estimation(S,Bright,VER_mean=VER_mean,VER_reshape=VER_reshape) 
        counts =calc_electron_density(VER,satlatlonalt,h,dn,cont)
        print 'VER & Ne profiles calculated'
              
    return VER,counts,h

def FUV_Level_2_OutputProduct_NetCDF(dn,satlatlonalt,az,ze,tanlatlonalt,Bright,VER,Ne,NmF2,hmF2,path='/home/dimitris/public_html/Datafiles/LVL2TEST/'):
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
    f = netcdf.netcdf_file(path+'lvl2_5_Output.nc', 'w')
    
    # Create the history of the file (created date)
    f.history = "File Created on: " + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
    
    f.createDimension('profile', len(Bright))
    f.createDimension('profileVG', len(az))
    f.createDimension('single', 1)
    f.createDimension('time', 6)
    # Maybe I need to add dimensions
    
    # Save the variables
    '''
    UT_Date_Time = f.createVariable('UT_DATE_TIME', 'i',('time',))
    UT_Date_Time[:] = [dn.year, dn.month, dn.day, dn.hour, dn.minute, dn.second]
    UT_Date_Time.units = 'Year-Month-Day-Hour-Mins-Secs'
    UT_Date_Time.long_name = 'Universal date and time'
    '''
    # This is for SCOTTs UT variable. No idea what it represents
    UT_Date_Time = f.createVariable('UT_DATE_TIME', 'i',('profile',))
    UT_Date_Time[:] = dn
    UT_Date_Time.units = 'TBD' #YYYY:MM:DD:HH:MM:SS:MMMM
    UT_Date_Time.long_name = 'Universal date and time'
    
    # Satellite coordinates
    satlat = f.createVariable('ICON_WGS84_LATITUDE','float',('single',))
    satlat[:] = satlatlonalt[0]
    satlat.units = 'degrees'
    satlat.long_name = 'Satellite latitude coordinate'
    
    satlon = f.createVariable('ICON_WGS84_LONGITUDE','float',('single',))
    satlon[:] = satlatlonalt[1]
    satlon.units = 'degrees'
    satlon.long_name = 'Satellite longitude coordinate'
    
    satalt = f.createVariable('ICON_WGS84_ALTITUDE','float',('single',))
    satalt[:] = satlatlonalt[2]
    satalt.units = 'km'
    satalt.long_name = 'Satellite altitude coordinate'
    
    # Azimuth and Zenith angles
    az_f = f.createVariable('ICON_WGS84_AZ','float',('profileVG',))
    az_f[:] = az
    az_f.units = 'degrees'
    az_f.long_name = 'Viewing geometry of observation, azimuth'
    
    ze_f = f.createVariable('ICON_WGS84_ZE','float',('profileVG',))
    ze_f[:]= ze
    ze_f.units = 'degrees'
    ze_f.long_name = 'Viewing geometry of observation, zenith'
    
    # Tangent point coordinates
    tanlat = f.createVariable('FUV_TANGENT_LATITUDE','float',('profile',))
    tanlat[:] = tanlatlonalt[:,0]
    tanlat.units = 'degrees'
    tanlat.long_name = 'Tangent points latitude coordinates'
    
    tanlon = f.createVariable('FUV_TANGENT_LONGITUDE','float',('profile',))
    tanlon[:] = tanlatlonalt[:,1]
    tanlon.units = 'degrees'
    tanlon.long_name = 'Tangent points longitude coordinates'
    
    tanalt = f.createVariable('FUV_TANGENT_ALTITUDE','float',('profile',))
    tanalt[:] = tanlatlonalt[:,2]
    tanalt.units = 'km'
    tanalt.long_name = 'Tangent points altitude coordinates'
    
    # Brightness profile    
    Bright_f = f.createVariable('FUV_TANGENT_BRIGHTNESS','float',('profile',))
    Bright_f[:] = Bright
    Bright_f.units = 'Rayleigh'
    Bright_f.long_name = 'FUV Tangent altitude Brightness profile'
    
    # VER profile    
    VER_f = f.createVariable('FUV_TANGENT_VER','float',('profile',))
    VER_f[:] = VER
    VER_f.units = 'ph/cm^3/sec'
    VER_f.long_name = 'FUV Tangent altitude Volume Emission Rate profile'
    
    # Ne profile    
    NE_f = f.createVariable('FUV_TANGENT_NE','float',('profile',))
    NE_f[:] = Ne
    NE_f.units = 'cm^(-3)'
    NE_f.long_name = 'FUV Tangent altitude Electron Density profile'

    # F2 peaks
    NmF2_f = f.createVariable('FUV_NMF2','float',('single',))
    NmF2_f[:] = NmF2
    NmF2_f.units = 'cm^(-3)'
    NmF2_f.long_name = 'NmF2 Value'
    
    hmF2_f = f.createVariable('FUV_HMF2','float',('single',))
    hmF2_f[:] = hmF2
    hmF2_f.units = 'km'
    hmF2_f.long_name = 'hmF2 Value'
    
    f.close()
    
def Get_lvl2_5_product(path_input='/home/dimitris/Data_Files/ICON_FUV_ray_UT_15sec_night.nc',path_output='/home/dimitris/public_html/Datafiles/LVL2TEST/'):
    '''
    Operational Code that reads Lvl1 file and creates the corresponding Lvl2.5
    INPUTS:
        path_input  - Input file path 
        path_output - Output file path
    OUTPUT:
        Creates a NetCDF file on the desired output_path
    NOTES:
        This versions uses paths as input and output. That can change in case needed. Also, the lvl1 file that is used
        is elementary. No actual LVL1 file has been given yet.
    HISTORY:
        24-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    # Open input NetCDF file
    data = netcdf.netcdf_file(path_input,mode='r')
    
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

    ver,Ne,h = FUV_Level_2_OutputProduct_Calculation(bright,h,satlatlonalt,az,ze,O,2,1,2,1,1,1,0)
    hmF2,NmF2 = find_hm_Nm_F2(Ne,h)
    NmF2 = NmF2
    
    FUV_Level_2_OutputProduct_NetCDF(dn,satlatlonalt,az,ze,tanlatlonalt,bright,ver,Ne,NmF2,hmF2,path_output)
    
    print 'LVL2.5 Processing Terminated. File created!'
    
def print_netcdf_decription(path):
    '''
    Reads the Lvl2.5 and prints the NetCDF file description
    INPUTS:
        path - Input file path 

    OUTPUT:
        Prints the details of the NetCDF file
    NOTES:
    HISTORY:
        24-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    g = netcdf.netcdf_file(path, 'r')
    variables = g.variables
    no_info_variables = []

    template = "{Variable:22}|{Units:19}|{Description:15}" # same, but named
    print template.format(Variable="VARIABLE NAME", Units="UNITS", Description="DESCRIPTION")
    print template.format(Variable= 22*"=",Units= 19*"=",Description=22*"="   )
    for var_name in variables.iterkeys():
        try:
            print template.format(Variable=var_name, Units="["+variables[var_name].units+"]", Description= variables[var_name].long_name)
            #print '%010s [%s]:%s' % (var_name, variables[var_name].units, variables[var_name].long_name)
        except Exception as e:
            no_info_variables.append(var_name)
            print '\n The following variables had no information:'
            print no_info_variables
    g.close()


def find_mid_cell(satlatlonalt,ze,az,spherical,limb = 130.):
    '''
    Calculates the top bottom and mid of the cell determined by the tangents altitudes
    INPUTS:
        satlatlonalt - vector containing the satellite coordinates [lat-lon-alt] (km)
        ze           - zenith angle of look direction (degrees)
        az           - azimuth angle of look direction (degrees) 
        shperical    - flag indicating if spherical [0] or WGS84 model[1] will be used
        limb         - altitude where the limb starts [default 130km]
    OUTPUT:
        rmid - mid point of the cell[km]
        rbot - bottom boundaries[km] (tangent altitudes)
        rtop - bottom boundaries[km] (tangent altitudes + top bound)
        h_loc_bot - WGS84 coordinate vector for non-spherical earth with bottom boundaries as altitudes
        h_loc_mid - WGS84 coordinate vector for non-spherical earth with mid cell point as altitudes
    NOTES:
    HISTORY:
        02-Sep-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    CALLS:
    '''
    
    az_v = np.deg2rad(az)
    ze_v = np.deg2rad(ze)
    
    if (spherical==0):
        h = np.zeros(np.size(ze))
        RE = 6371.
        h = ic.angle2tanht(ze_v, satlatlonalt[2], RE) 
        h = h[np.where(h>limb)]
        h_loc_bot = 0
        h_loc_mid = 0
    else:
        # Tangent_poing function returns coordinates [lat-lon-alt] for the tangent point
        h_coord = np.zeros((len(ze),3))
       
        for i in range(0,len(ze)):
            h_coord[i,:] = ic.tangent_point(satlatlonalt,az[i],ze[i])
            
        # Initialize matrix to make sure dimensions are correct
        h_loc_bot = np.zeros((np.size(np.where(h_coord[:,2]>=limb)),3))
        h_loc_bot[:,:] = h_coord[np.where(h_coord[:,2]>=limb),:]
        h = h_coord[np.where(h_coord[:,2]>=limb),2]
        h = h[0,:]
    
    rbot = h
    rtop = rbot.copy()

    rtop[1:] = rbot[:-1]

    #rtop[0] = Horbit -1
    rtop[0] = satlatlonalt[2]
    # Define midpt of each layer
    rmid = (rbot + rtop)/2  
    h_loc_mid = h_loc_bot
    
    if np.size(h_loc_mid)!=1:
        h_loc_mid[:,2] = rmid
    
    return rmid,rbot,rtop, h_loc_bot, h_loc_mid

