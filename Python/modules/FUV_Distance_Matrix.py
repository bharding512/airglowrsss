import numpy as np
from datetime import datetime, timedelta
from scipy import stats
import multiprocessing
from multiprocessing import Pool
import time
import ICON as ic
reload(ic)
import math
from sys import exit

def create_cells_Matrix_spherical_symmetry(theta,Horbit,RE=6371.):
    '''
    Returns the boundaries for the created cells of the atmosphere assuming sperical symmetry of earth.
    INPUTS:
        theta - zenith agnle vector, rad
        Horbit - satellite altitude, km
        RE - radius of Earth, km
    OUTPUTS:
        S - Weight Matrix
        rtop - Upper Bound for each cell
        rmid - Mid Point for each cell
        rbot - Bottom Bound for each cell
    Comments:
        --  rbot is the tangent altitude for each zenith angle. 
        --  We want to have the same number of cells as the number 
            of zenith angle measurements in order to get an invertible matrix.
    '''
    
    Ntheta = np.size(theta)
    
    #for angle in theta:
    rbot= ic.angle2tanht(theta, Horbit, RE)
    # Define top of each layer.
    
    rtop = rbot.copy()

    rtop[1:] = rbot[:-1]

    rtop[0] = Horbit-1
    #rtop[0] = rtop[1] + (rtop[1] - rtop[2])
    # Define midpt of each layer
    rmid = (rbot + rtop)/2

    # Build observation matrix
    S = np.zeros((Ntheta, len(rmid)))
    for i in range(Ntheta):
        for j in range(len(rmid)):
            th = theta[i]
            rb = rbot[j]
            rt = rtop[j]
            sb2 = -math.sin(th)**2  + ((RE+rb)/(RE+Horbit))**2
            st2 = -math.sin(th)**2  + ((RE+rt)/(RE+Horbit))**2   
            if sb2 < 0: # there is no intersection of LOS with altitude rb. Set term to 0.
                # Note: this might be due to numerical rounding for tangent altitude. 
                # Do the same thing either way.
                sb2 = 0.
            if st2 < 0: # there is no intersection of LOS with altitude rt. Set term to 0.
                st2 = 0.
            path_len_km = 2*(RE+Horbit) * ( math.sqrt(st2) - math.sqrt(sb2) )
            path_len_cm = path_len_km*1e5
            # Result should be a conversion from VER to Rayleighs. path_length_cm * 10^-6 to match paper
            S[i,j] = path_len_cm * 1e-6
            
    return S,rmid,rbot,rtop


def Calculate_D_Matrix_WGS84_mp_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return ic.distance_to_shell(*a_b)

def Calculate_D_Matrix_WGS84_mp(satlatlonalt,az,ze,azze_def = 1):
    '''
    Returns the Distance matrix calculated given the satellite coordinates and viewing geometry for the WGS84 model.
    INPUTS:
        satlatlonalt -  vector containing the satellite coordinates [lat-lon-alt] (km)
        az -  azimuth vector (deg)
        ze -  zenith vector (deg)
        azze_def - falg indicating if default az and ze fov angles will be used from instrument parameters
    OUTPUT:
        S  - Distance matrix assuming WGS84 (km 10^{-1})
    NOTES:
    HISTORY:
        03-Jun-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        26-Aug-2015: Added the exception for the Pool
    '''
    try:
        rbot = np.zeros(np.size(ze,0))

        for i in range(0,np.size(ze,0)):
            _,_,rbot[i] =  ic.tangent_point(satlatlonalt,az[i],ze[i])

        rtop = rbot.copy()
        rtop[1:] = rbot[:-1]
        rtop[0] = satlatlonalt[2] - 1
        #rtop[0] = rtop[1] + (rtop[1] - rtop[2])
        rmid = (rbot + rtop)/2
        S = np.zeros((np.size(ze),np.size(rbot)))
        k = 0

        N = multiprocessing.cpu_count()

        # Create the pool.  Be nice.  Don't use all the cores!
        pool = Pool(processes=16)
        t0 = time.time()
        for i in range(0,np.size(ze)):

            job_args = [(satlatlonalt, az[i], ze[i],rtop[j]) for j in range(0,len(rbot))]
            job_args2 = [(satlatlonalt, az[i], ze[i],rbot[j]) for j in range(0,len(rbot))]
            job_args3 = [(satlatlonalt, az[i], ze[i],rtop[j],'second') for j in range(0,len(rbot))]

            ub = pool.map(Calculate_D_Matrix_WGS84_mp_star,job_args)
            lb = pool.map(Calculate_D_Matrix_WGS84_mp_star,job_args2)
            ub1 = np.array(ub)
            '''
            if np.sum(np.isnan(lb)) == len(lb):
                lb2 = pool.map(Calculate_D_Matrix_WGS84_mp_star,job_args3)
                lb21 = array(lb2)
                S[i,np.where(np.isnan(ub)==True)and(np.where(np.isnan(lb)==True))]=0
                diff = (ub1[np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb2)==False))] - lb21[np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb2)==False))])
                S[i,np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb2)==False))]=abs(diff)
            else:
                lb1 = array(lb)
                S[i,np.where(np.isnan(ub)==True)and(np.where(np.isnan(lb)==True))]=0
                diff = (ub1[np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb)==False))] - lb1[np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb)==False))])
                S[i,np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb)==False))]=2*abs(diff)
            '''
            lb2 = pool.map(Calculate_D_Matrix_WGS84_mp_star,job_args3)
            lb21 = np.array(lb2)
            S[i,np.where(np.isnan(ub)==True)and(np.where(np.isnan(lb)==True))]=0
            diff2 = (ub1[np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb2)==False))] - lb21[np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb2)==False))])
            lb1 = np.array(lb)
            diff = (ub1[np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb)==False))] - lb1[np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb)==False))])
            S[i,np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb)==False))]=2*abs(diff)
            diago = np.where(np.isnan(ub)==False)and(np.where(np.isnan(lb)==True))[0][0]
            S[i,diago]=abs(diff2[diago])
            
        t1 = time.time()
        #print t1-t0
        pool.close()
        pool.join()
        S = S*1e-1
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
    
    return S,rmid

def Calculate_D_Matrix_WGS84(satlatlonalt,az,ze,azze_def = 1):
    '''
    Returns the Distance matrix calculated given the satellite coordinates and viewing geometry for the WGS84 model.
    INPUTS:
        satlatlonalt -  vector containing the satellite coordinates [lat-lon-alt] (km)
        az -  azimuth vector (deg)
        ze -  zenith vector (deg)
        azze_def - falg indicating if default az and ze fov angles will be used from instrument parameters
    OUTPUT:
        S  - Distance matrix assuming WGS84 (km 10^{-1})
    NOTES:
    HISTORY:
        03-Jun-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    
    rbot = np.zeros(np.size(ze,0))
    
    for i in range(0,np.size(ze,0)):
        _,_,rbot[i] =  ic.tangent_point(satlatlonalt,az[i],ze[i])
    
    rtop = rbot.copy()
    rtop[1:] = rbot[:-1]
    rtop[0] = satlatlonalt[2] - 1
    rmid = (rbot + rtop)/2
    S = np.zeros((np.size(ze),np.size(rbot)))
    k = 0
    t0 = time.time()
    for i in range(0,np.size(ze)):
        for j in range(0,np.size(rbot)):
            ub = ic.distance_to_shell(satlatlonalt, az[i], ze[i], rtop[j])
            lb = ic.distance_to_shell(satlatlonalt, az[i], ze[i], rbot[j])
            if np.isnan(ub) and np.isnan(lb):
                S[i,j] = 0
            elif np.isnan(lb):
                lb = ic.distance_to_shell(satlatlonalt, az[i], ze[i], rtop[j],'second')
                S[i,j] = abs(ub - lb)
            else:
                S[i,j] = 2*abs(ub - lb)
    t1 = time.time()
    #print t1-t0
    S = S*1e-1
    
    return S

