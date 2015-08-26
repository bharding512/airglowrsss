# Common code for ICON users at Illinois. Geometrical transformations, airglow code, etc.

import numpy as np
from pyglow import pyglow
from datetime import datetime, timedelta
from scipy import stats
import multiprocessing
from multiprocessing import Pool
import time
import math
import ICON as ic
import Regularization as reg
from scipy.io import netcdf
from time import gmtime, strftime

# ICON FUV
def get_FUV_instrument_constants():
    '''
    Notes:
        12-Dec-2014: Last known parameters for ICON FUV   
    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu)
        01-Apr-2015: Last known parameters for ICON FUV - Rescell Dimensions  
    '''

    instrument = {    'npixelx': 6,            # number of pixels per slice of the interferogram in a rescell(horizontal) [512px total]
                      'npixely': 256,           # number of rescells in vertical (altitude) direction on CCD.
                'aperture_area': 0.006*0.032,   # [cm^2]
                   'coneangle1': 18,            # fov of instrument in horizontal direction, [deg]
                   'coneangle2': 24,            # fov of instrument in vertical direction, [deg]
                 'coneangle1_r': 0.314159265359,# fov of instrument in horizontal direction, [rad]
                 'coneangle2_r': 0.418879020479,# fov of instrument in vertical direction, [rad]
                    'exposure' : 12,            # exposure [sec]
                  'Sensitivity': 0.083,         # combined transmittance of all optics. [counts/res_cell/s/R]
                        'fov_l': 98,            # Field of View (Vertical) Lower Bound (degrees)
                        'fov_u': 123,           # Field of View (Vertical) Upper Bound (degrees)
                       'fovr_l': 1.71042266695, # Field of View (Vertical) Lower Bound (radians)
                       'fovr_u': 2.14675497995, # Field of View (Vertical) Upper Bound (radians)
                 'stripes_used': 1,             # Number of stripes used in the ccd the CCD
                   'Dark Noise': 4000,          # Dark noise [Not to be of our concern, only for reference]
                   'Read Noise': 60,            # Read noise [e-][Not to be of our concern, only for reference]
                  }

    return instrument

# ICON FUV
def calc_1356_nighttime(lat,lon,alt,dn,Ne_scaling = 1.):
    '''
    Return the radiative recombination (RR) and mutual neutralization (MN) for the
    135.6-nm emission according to a requested location and time using IRI and MSIS
    and the emission chemistry from Melendez-Alvira, 1999 
    (http://onlinelibrary.wiley.com/doi/10.1029/1999JA900136/abstract)
    INPUTS:
        lat -  latitude (deg)
        lon -  longitude (deg)
        alt -  altitude (km)
        dn  -  UT date and time to use (datetime)
        Ne_scaling - Scaling factor for Ne [default 1 => No scaling]
    OUTPUT:
        RR  - volume emission rate due to radiative recombination (1/cm^3/s)
        MN  - volume emission rate due to mutual neutralization (1/cm^3/s)
        Ne  - Electron Density Profile (1/cm^3/s)
        O   - Oxygen Profile (1/cm^3/s)
    NOTES:
    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu)
        13-Jul-2015: Add O as return value
    '''
    # Coefficients from Table 3 of Melendez-Alvira (1999)
    b1356 = 0.54    # yield parameter (unitless)
    a1356 = 7.3e-13 # radiative recombination rate (cm^3/s)
                    # consider that to be constant through altitude, normally it should change with temperature
    k1 = 1.3e-15    # radiative attachment rate (cm^3/s)
    k2 = 1e-7       # ion-ion neutralization rate (cm^3/s)
    k3 = 1.4e-10    # ion-atom neutralization rate (cm^3/2)
    
    # Run pyglow at the requested location and time.  IRI and MSIS calls are needed to generate
    # model outputs of needed constituents
    pt = pyglow.Point(dn, lat, lon, alt)
    pt.run_iri()
    pt.run_msis()
    
    # Pull necessary constitutents
    O_p = pt.ni['O+']  # O+ density (1/cm^3)
    O = pt.nn['O']     # O density (1/cm^3)
    Ne = pt.ne         # electron density (1/cm^3)
    Ne = Ne*Ne_scaling # scale Ne

    # Calcualte radiative recombination (equation 17) and mutual neutralization (equation 18)
    RR = a1356*Ne*O_p  # radiatvie recombination (1/cm^3/s)
    MN  = (b1356*(k1*k2)) *((Ne*O*O_p)/(k2*O_p + k3*O)) # mutual neutralization (1/cm^3/s)

    return RR,MN,Ne,O

def calculate_VER_1356_nighttime(satlat,satlon,satalt,dn,Ne_scaling = 1.):
    '''
    Step along a desired look direction from a given observing location and
    calculate the 135.6-nm intensity. Lines-of-sight are calculated using a spherical earth.
    INPUTS:
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        Ne_scaling - Scaling factor for Ne [default 1 => No scaling]
    OUTPUTS:
        VER         - A VER altitude profile for given lat/lon coords
        Ne          - A Electron density altitude profile for given lat/lon coords
    NOTES:

    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu) 
        08-Jul-2015: Added check for satalt scalar values (Dimitrios Iliou)
        13-Jul-2015: Add O to calc_1356_nighttime to have O tanget altitude profile
    '''
    VER = np.zeros(np.size(satalt))
    MN = np.zeros(np.size(satalt))
    VER_true = np.zeros(np.size(satalt))
    NE = np.zeros(np.size(satalt))
    O = np.zeros(np.size(satalt))
    
    if (np.size(satalt)==1):
        VER,MN,NE,O= calc_1356_nighttime(satlat,satlon,satalt,dn,Ne_scaling)
    else:
        for i in range(0,np.size(satalt)):
            VER[i],MN[i],NE[i],O[i]= calc_1356_nighttime(satlat,satlon,satalt[i],dn,Ne_scaling)

    VER_true = VER + MN

    return VER_true, VER,MN, NE, O

# ICON FUV
def calculate_pixel_1356_nighttime(ze,az,satlat,satlon,satalt,dn,cont=1,Ne_scaling=1., step=10.):
    '''
    Step along a desired look direction from a given observing location and
    calculate the 135.6-nm intensity. Lines-of-sight are calculated using a spherical earth.
    INPUTS:
        ze          - zenith angle of look direction (radians)
        az          - azimuth angle of look direction (0 is north, pi/4 is east) (radians)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
        Ne_scaling  - Scaling factor for Ne [default 1 => No scaling]
        step        - resolution of the integration along the line of sight (km, default = 10 km)
    OUTPUTS:
        Brightness  - the intensity of the integrated emission (R) or 0 if the contribution is not set correctly
    NOTES:
        A spherically symmetric atmosphere is currently assumed.  So, all calls to the calc_1356
        routine are done at the lat/lon of the satellite.  This needs to be generalized at some point.
    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu) and Brian J. Harding
       (bhardin2@illinois.edu) 
        03-Jun-2015: Changed Limb altitude from 90 to 150km (Dimitrios Iliou) 
        13-Jul-2015: Add _ to calc_1356_nighttime since O is not needed for the fwd model
    '''
    # Define constants and step sizes.  We work in a 2D coordinate system at this point
    RE = 6371.      # radius of earth (km)
    x = 0.          # This measures perpendicular to the satellite nadir
    y = RE + satalt # This measures parallel to satellite nadir
    dx = np.sin(ze) # unit step in x
    dy = np.cos(ze) # unit step in y
      
    IRR = 0.
    IMN = 0.
    Ne = 0
    Rayleigh = 0.   

    # We will trace along the ray path until we get back to the satellite altitude
    # NOTE: (1) Generalize this by calculating points and using a for loop.
    #       (2) Generalize this by allowing a ray path to be passed in.
    while (np.sqrt(x**2 + y**2) <= RE+satalt): # until you get far enough away from earth that it doesn't matter
             
        alt = np.sqrt(x**2 + y**2) - RE
        lat = np.nan # for now, this doesn't matter
        lon = np.nan # for now, this doesn't matter    
        
        # After this altitude the Raypath goes through the disk so no calculation is done for these altitudes
        #if alt < 150:    # 150 USUALLY 100
            #break
       
        # Calculate 1356 Emission for a single point
        VER,MN,Ne,_ = calc_1356_nighttime(satlat,satlon,alt,dn,Ne_scaling)
        
        IRR = IRR + VER * step*10**5
        IMN = IMN + MN * step*10**5
        
        # take the next step
        x += step*dx
        y += step*dy

    # WE dont divide with 4*pi
    IRR = (10**(-6)*IRR)
    IMN = 10**(-6) * IMN
    
    # Calculate Brightness depending on the contribution
    # Other contribution may be added here in the future
    # Contribution from Radiative Recombination and Mutual Nutralization
    if cont == 1:
        Rayleigh = IRR + IMN
    # Contribution only from Radiative Recombination
    elif cont ==2:
        Rayleigh = IRR # This is needed only when we need RR 
    else:
        print 'Invalid Choice of contribution. 1 for RR+MN or 2 for RR'
        return 0,0
    
    return Rayleigh

# ICON FUV WGS84
def calculate_pixel_1356_nighttime_WGS84(ze,az,satlat,satlon,satalt,dn,symmetry = 0,cont=1,Ne_scaling=1., step=10.):
    '''
    Step along a desired look direction from a given observing location and
    calculate the 135.6-nm intensity. Lines-of-sight are calculated using a non-spherical non-symmetric earth.
    INPUTS:
        ze          - zenith angle of look direction (deg)
        az          - azimuth angle of look direction (0 is north, pi/4 is east) (deg)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        symmetry    - flag indicating if symmetry [0] or non-symmetry[1] will be used
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
        Ne_scaling  - Scaling factor for Ne [default 1 => No scaling]
        step        - resolution of the integration along the line of sight (km, default = 10 km)
    OUTPUTS:
        Brightness  - the intensity of the integrated emission (R) or 0 if the contribution is not set correctly
    NOTES:
        
    HISTORY:
        08-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu) 
    '''
    
    # Initialize the integration variables to zeros
    IRR = 0.
    IMN = 0.
    Ne = 0
    Rayleigh = 0.   

    satlatlonalt = [satlat,satlon,satalt]
    #print '%f,%f'%(az,ze)
    # Call the function for a single az and ze angle. It returns all the stepping points along the line
    xyz, latlonalt = ic.project_line_of_sight(satlatlonalt, az,ze,step)
      
    # For all the points on the line   
    for j in range(0,np.size(latlonalt,1)):
        
        # Continue iterating till you reach the conjugate point of the satellite.
        # This is not necessery since for high altitudes the contributions are not significant.
        if (latlonalt[2,j]>latlonalt[2,0]):
            break
        
        # Calculate 1356 Emission for a single point        
        # Return Radiative Recombination, Mutual Neutralization and Ne value for each point on the line
        # The Ne can be scalled to test algorithm for various pertubations
        # The symmetry check changes the call so the lat lon coordinates will be the same as nadir direction or not. 
        if symmetry == 0:
            VER,MN,Ne,_ = calc_1356_nighttime(satlat,satlon,latlonalt[2,j],dn,Ne_scaling)
        elif symmetry == 1:
            VER,MN,Ne,_ = calc_1356_nighttime(latlonalt[0,j],latlonalt[1,j],latlonalt[2,j],dn,Ne_scaling)
        
        #print '%f,%f,%f' %(latlonalt[0,j],latlonalt[1,j],latlonalt[2,j])
        
        # Iterative summation (integration of raypath) 
        # step is multipled with 10**5 because we want the results in cm(^3)
        IRR = IRR + VER * step*10**5
        IMN = IMN + MN * step*10**5

    # WE dont divide with 4*pi
    # 10**(-6) is in the formula of calculating Brightness in Rayleigh. 
    IRR = 10**(-6) * IRR
    IMN = 10**(-6) * IMN
    
    # Calculate Brightness depending on the contribution
    # Other contribution may be added here in the future
    # Contribution from Radiative Recombination and Mutual Nutralization
    if cont == 1:
        Rayleigh = IRR + IMN
    # Contribution only from Radiative Recombination
    elif cont ==2:
        Rayleigh = IRR # This is needed only when we need RR 
    else:
        print 'Invalid Choice of contribution. 1 for RR+MN or 2 for RR'
        return 0,0

    return Rayleigh

# ICON FUV
def get_Photons_from_Brightness_1356_nighttime(ze,az,satlat,satlon,satalt,dn,symmetry =0,shperical=1,exposure=0.,TE=0.,cont=1,Ne_scaling = 1.,step = 10, stripes_used = 0):
    '''
    Calls 'calculate_pixel_1356_nighttime' which calculates the Brightness through integrated VER for a given zenith angle.
    INPUTS:
        ze          - zenith angle of look direction (radians)
        az          - azimuth angle of look direction (0 is north, pi/4 is east) (radians)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        symmetry    - flag indicating if symmetry [0] or non-symmetry[1] will be used
        shperical   - flag indicating if spherical [0] or WGS84 model[1] will be used
        exposure    - Exposure time of the FUV CCD, if input is zero the default parameters are loaded
        TE          - Total Efficiency of the FUC optical system, if input is zero the default parameters are loaded
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
        Ne_scaling - Scaling factor for Ne [default 1 => No scaling]
        step        - resolution of the integration along the line of sight (km, default = 10 km)
        stripes_used- number of stripes used of the CCD [default 1] [if zero on input loads default value]
    OUTPUTS:
        Brightness  - the intensity of the integrated emission (R) or 0 if the contribution is not set correctly
        photons     - the number of photons seen on the detecter for a given zenith angle [counts]
    NOTES:
        A spherically symmetric atmosphere is currently assumed.  So, all calls to the calc_1356
        routine are done at the lat/lon of the satellite.  This needs to be generalized at some point.

        Brightness and Photons are calculated for a specified number of pixel in a rescell: Default = 8 (One rescell)
        
        photon can be also calculated by multipling the radiance with the above parameters
        photons = r * solidangle * aperture * exposure * OE *0.03
    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu) 
        01-Apr-2015: pixels_per_rescell -> stripes used 
        08-Jul-2015: Changed the calculate_pixel_1356_nighttime to calculate_pixel_1356_nighttime_WGS84 to include non-spherical non-symmetric earth
        26-Jul-2015: Added the spherical and symmetry parameters on the function


    '''
    
    # Check if instrument parameters are zero from input to load default values
    params = get_FUV_instrument_constants()
    if exposure==0:
        exposure = params['exposure']
    if TE==0:
        TE =  params['Sensitivity']
    if stripes_used==0:
        stripes_used = params['stripes_usedl']

    # Get Brightness for a given zenith angle 
    if shperical == 0:
        Brightness = calculate_pixel_1356_nighttime(ze,az,satlat,satlon,satalt,dn,cont,Ne_scaling,step)
    elif shperical == 1:
        Brightness = calculate_pixel_1356_nighttime_WGS84(ze,az,satlat,satlon,satalt,dn,symmetry,cont,Ne_scaling,step)
    
    # Number of pixels in rescell = 8
    r = TE * exposure * stripes_used
  
    photons = r*Brightness
      
    return Brightness,photons

#Multiprocessing
def get_Photons_from_Brightness_1356_nighttime_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return get_Photons_from_Brightness_1356_nighttime(*a_b)

# ICON FUV
def get_Photons_from_Brightness_Profile_1356_nighttime(ze,az,satlat,satlon,satalt,dn,symmetry=0,shperical=1,exposure=0.,TE=0.,cont=1,Ne_scaling = 1.,step = 10,stripes_used = 0,proc=16):
    '''
    Calls 'get_Photons_from_Brightness_1356_nighttime' which calculates the Brightness and Electron density by going through integrated VER for all given zenith angles.
    INPUTS:
        ze          - zenith angle of look direction (radians)
        az          - azimuth angle of look direction (0 is north, pi/4 is east) (radians)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        symmetry    - flag indicating if symmetry [0] or non-symmetry[1] will be used
        shperical   - flag indicating if spherical [0] or WGS84 model[1] will be used
        exposure    - Exposure time of the FUV CCD, if input is zero the default parameters are loaded
        TE          - Total Efficiency of the FUC optical system, if input is zero the default parameters are loaded
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
        Ne_scaling - Scaling factor for Ne [default 1 => No scaling]
        step        - resolution of the integration along the line of sight (km, default = 10 km)
        stripes_used- number of stripes used of the CCD [default 1] [if zero on input loads default value]
        proc        - Number of cores used for multiprocessing
    OUTPUTS:
        Brightness  - the intensity of the integrated emission (R) or 0 if the contribution is not set correctly
        photons     - the number of photons seen on the detecter for a given zenith angle [counts]
    NOTES:
        A spherically symmetric atmosphere is currently assumed.  So, all calls to the calc_1356
        routine are done at the lat/lon of the satellite.  This needs to be generalized at some point.

        Brightness and Photons are calculated for a specified number of pixel in a rescell: Default = 8 (One rescell)

    HISTORY:
        12-Dec-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        01-Apr-2015: pixels_per_rescell -> stripes used 
        26-Jul-2015: Added the spherical and symmetry parameters on the function
    '''

    params = get_FUV_instrument_constants()
    if exposure==0:
        exposure = params['exposure']
    if TE==0:
        TE =  params['Sensitivity']
    if stripes_used==0:
        stripes_used = params['stripes_used']

    photons = np.zeros(np.size(ze)) # Number of Counts without Noise
    Rayl = np.zeros(np.size(ze))
    
    # i put index on azimuth
    job_args = [(ze[i],az[i] ,satlat,satlon,satalt, dn,symmetry,shperical,exposure,TE,cont,Ne_scaling,step,stripes_used) for i in range(0,len(ze))]
    N = multiprocessing.cpu_count()

    # Create the pool.  Be nice.  Don't use all the cores!
    pool = Pool(processes=16)
    
    t0 = time.time()
    results = pool.map(get_Photons_from_Brightness_1356_nighttime_star,job_args)
    for i in range(0,len(results)):
        Rayl[i] = results[i][0]
        photons[i] = results[i][1]

    t1 = time.time()
    pool.close()
    pool.join()
    #print t1-t0
    
    return Rayl,photons

# ICON FUV
def add_noise_to_photon_and_brightness(photons,exposure=0.,TE=0.,stripes_used = 0,reps=1):
    '''
    Returns Brightness and Electron Density profile(s) after adding shot noise
    INPUTS:
        photons     - number of photons (counts) [Noiseless Electron density profile]
        Exposure    - FUV Exposure time (s)   [if zero on input loads default value]
        TE          - Total Efficiency        [if zero on input loads default value]
        stripes_used- number of stripes used of the CCD [default 1] [if zero on input loads default value]
        reps        - Number of noise realizations to return [default = 1]
    OUTPUTS:
        Brightness- Noisy Brightess Profile(s) [Rayleighs]
        Shot_noise- Noisy Photon Profile(s) [count]
    Comments:
        noise - input noise, Poisson Dist (counts)
    HISTORY:
        12-Dec-2015: Written by Dimitrios Iliou (iliou2@illinois.edu) 
        01-Apr-2015: pixels_per_rescell -> stripes used 
    '''
    
     # Load default FUV parameters in case are not given in the input
    params = get_FUV_instrument_constants()
    if exposure==0:
        exposure = params['exposure']
    if TE==0:
        TE =  params['Sensitivity']
    if stripes_used==0:
        stripes_used = params['stripes_used']

    shot_noise = np.zeros((reps,np.size(photons))) 
    Rayl_ = np.zeros((reps,np.size(photons)))
    
    for rep in range(0,reps):
        for i in range(0,np.size(photons)):
            '''
            The variance of the photon noise, which is a measure of the expected difference 
            between the number of photons collected and the average number, is equal to the 
            square root of the average number of photons 
            '''
            if photons[i]==0:
                shot_noise[rep,i]
                print 'add_noise_to_photon_and_brighntess: zero value in signal'
            else:
                shot_noise[rep,i] = stats.poisson.rvs(photons[i],1)

        Rayl_[rep,:] = shot_noise[rep,:]/(TE*exposure*stripes_used)
    
    return Rayl_,shot_noise

def run_forward_modelling(satlatlonalt,date,symmetry = 0.,shperical=1, exp=12,reps=1000,sens=0,cont=2,Ne_sc=1.,step=10,stripes_used=0,proc=16):
    '''
    Top forward model modules. Creates brightness profile and adds noise returning the noisy Brighntess profile(s) to 
    be used for the inversion process. 
    INPUTS:
        satlatlonalt - vector containing the satellite coordinates [lat-lon-alt] (km)
        date         - datetime input
        symmetry     - flag indicating if symmetry [0] or non-symmetry[1] will be used
        shperical    - flag indicating if spherical [0] or WGS84 model[1] will be used
        exp          - exposude time in order to be used as proxy for SNR increase or decrease (s)
        reps         - number of noise realizations for the Brightness profile. Determines the size of the ouput vector.
        sens         - Sensitivity of the FUV optical system, if input is zero the default parameters are loaded
        cont         - contribution factors, cont=1 -> RR+MN, 2-> RR
        Ne_sc        - Scaling factor for Ne [default 1 => No scaling]
        step         - resolution of the integration along the line of sight (km, default = 10 km)
        stripes_used - number of stripes used of the CCD [default 1] [if zero on input loads default value]
        proc         - Number of cores used for multiprocessing [default 16 cores]        
    OUTPUT:
        NE           - Electron Density profile corresponding the the specified datetime and tangent altitude (cm^-3)
        Bright       - Brightness profile (Rayleigh)
        Bright_n     - Noisy Brightness profile (Rayleigh) - Size determined from the number of noise realizations
        h            - Tangent altitudes (lower cell boundaries) (km)
    NOTES:
        - Currently the forward model doesnt account for different azimuth angles. Azimuth is constant at 0 degrees 
          (S/C looking north)
        - Azimuth calculation to be added
        - Current ICON.py "calculate_pixel_1356" function assumes shperical symmetric earth instead of WGS84 needs to 
          be changed
        - NmF2 Difference between tangents points calculated for WGS84 and Spherical Symmetric Earth is 0.0016% 
          and hmF2 -1.12
        - VER and NE are calculated for fixed lat,lon-> Need to change to follow the WGS84 model
    HISTORY:
        03-Jun-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        08-Jul-2015: Changed everything to assume non-symmetric non-spherical earth. Choice for symmetry is given.
        22-Jul-2015: Changed VER in Symmetry to include satalt and satlot instead of zero-zero
        26-Jul-2015: Added the spherical parameter on the function
    '''
    
    # Satellite coordinates
    satlat = satlatlonalt[0]
    satlon = satlatlonalt[1]
    satalt = satlatlonalt[2]
    
    # Load default az and ze angles from the fov of the instrument
    az,ze = get_azze_default()
    az_v = np.deg2rad(az)
    ze_v = np.deg2rad(ze)

    '''
    Locate limb lower bound ~150km 
    '''
        
    check = 0
    
    # Locate limb lower bound ~150km 
    if (shperical==0):
        h = np.zeros(np.size(ze))
        RE = 6371.
        h = ic.angle2tanht(ze_v, satalt, RE) 
        h = h[np.where(h>150)]
        disk = len(h)
    else:
        # Tangent_poing function returns coordinates [lat-lon-lat] for the tangent point
        h_coord = np.zeros((len(ze),3))
        # Locate limb lower bound ~150km 
        for i in range(0,len(ze)):
            h_coord[i,:] = ic.tangent_point(satlatlonalt,az[i],ze[i])
        h_loc = h_coord[np.where(h_coord[:,2]>=150),:]
        h = h_coord[np.where(h_coord[:,2]>=150),2]
        h = h[0,:]
        disk = len(h)

    # Resize vectors for limb FOV
    ze = ze[0:disk]
    az = az[0:disk]
    
    ze_v = ze_v[0:disk]
    az_v = az_v[0:disk]
    
    
    # Initialize vectors
    photons = np.zeros(np.size(ze_v)) # Number of Counts without Noise
    Bright = np.zeros(np.size(ze_v))
    Bright_n = np.zeros(np.size(ze_v))
    Counts = np.zeros(np.size(ze_v)) # Number of Counts with Noise    

    '''
    Multicore Processing for Calculate Brighness
    '''
    ### ___ I NEED TO ADD THE SPHERICAL THING HERE
    
    # Takes s/c coordinates as input for poing location to calculate the brightness profile. This must change to tan. point
    if (shperical==0):
        Bright,photons = get_Photons_from_Brightness_Profile_1356_nighttime(ze_v,az_v,satlat,satlon,satalt,date,symmetry,shperical,exp,sens,cont,Ne_sc,step,stripes_used,proc)
    else:
        Bright,photons = get_Photons_from_Brightness_Profile_1356_nighttime(ze,az,satlat,satlon,satalt,date,symmetry,shperical,exp,sens,cont,Ne_sc,step,stripes_used,proc)
        #Bright,photons = ic.get_Photons_from_Brightness_Profile_1356_nighttime(ze_v,az_v,satlat,satlon,satalt,date,exp,sens,cont,Ne_sc,step,stripes_used,proc)
    #Bright,photons = ic.get_Photons_from_Brightness_Profile_1356_nighttime(ze_v,az_v,h[:,0],h[:,1],satalt,date,exp,0,2)
    
    '''
    Calculate Noisy Counts
    reps determines the number of noise profiles that the function will produce
    '''
    Bright_n= np.zeros((reps,np.size(ze_v)))
    #Bright_n,Counts = ic.add_noise_to_photon_and_brightness(photons,exp,sens,stripes_used,reps)
    Bright_n,Counts = add_noise_to_photon_and_brightness(photons,exp,sens,stripes_used,reps)
    Bright_n = np.transpose(Bright_n)

    '''
    Caclulate VER 
    '''    
    VER_true = np.zeros(np.size(h,0))
    MN = np.zeros(np.size(h,0))
    VER = np.zeros(np.size(h,0))
    NE = np.zeros(np.size(h,0))
    O = np.zeros(np.size(h,0))

    '''
    Calculates fixed NE for location 0,0 (S/C location)-> Needs to be changed for the WGS84 case where each tangent
    point has different lat and lon.
    '''
    if (shperical==0):
        VER_true,VER,MN,NE,O = calculate_VER_1356_nighttime(satlat,satlon,h,date)
    else:
        if symmetry == 0:
            _,_,_,NE,O = calculate_VER_1356_nighttime(satlat,satlon,h,date)
        else:
            for i in range(0,len(h)):
                _,_,_,NE[i],O[i] = calculate_VER_1356_nighttime(h_loc[0,i,0],h_loc[0,i,1],h_loc[0,i,2],date)

    return NE,Bright,Bright_n,h,O#,h_coord

def get_azze_default():
    '''
    Load the default fov parameters for FUV and calculates ze and az by evenly distributing the FOV over the size
    of the CCD
    INPUTS:
        NONE
    OUTPUT:
        az -  azimuth vector (deg)
        ze -  zenith vector (deg)
    NOTES:
        This function is used only if we want to load the default values when dealing with the forward modelling. When 
        data are give we will have to calculate the ze and az angle depending on that data.
    HISTORY:
        03-Jun-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''
    # Load instrument parameters
    params = get_FUV_instrument_constants()
    fov_l =  params['fov_l']
    fov_u =  params['fov_u']
    fovr_l =  params['fovr_l']
    fovr_u =  params['fovr_u']
    npixely = params['npixely']
    
    fov_az = params['coneangle1']
    npixelx = params['npixelx']

    # Create zenith angle vector by evenly distributing the FUV FOV angles 
    # For the azimuth angle we just assume that the satellite is facing north
    az = np.zeros(npixely) 
    az[:] = 0.
    ze = np.linspace(fov_l,fov_u,npixely)

    return az,ze

def decs(x, pos, sc = 1e-5):
    'The two args are the value and tick position'
    return '%1.f' % (x*sc)

def FUV_Level_2_OutputProduct_Calculation(Bright,h,satlatlonalt,az,ze,O,spherical=1,cont =1, regu=1,regu_order = 2,S_mp = 1):
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
        shperical   - Spherical Earth Assumption, 0=> Spherical, 1=>WGS84
        cont        - Contribution, 1=> RR+MN, 2=>RR ( This will be removed propably) [int]
        regu        - Choose regularization, 0=> Gaussian Elimination, 1=> Tikhonov L-Curve, 2=> GSVD [int]
        regu_order  - Regularization Order [int] (0,1,2 Possible Values)

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
    '''
    
    az = az[0:len(h)]
    ze = ze[0:len(h)]

    if (spherical==0):
        S,_,rbot,_ = ic.create_cells_Matrix_spherical_symmetry(np.deg2rad(ze),satlatlonalt[2])
    elif (spherical ==1):
        if S_mp ==1:
            S,rmid = Calculate_D_Matrix_WGS84_mp(satlatlonalt,az,ze)
        elif S_mp == 0:
            S = Calculate_D_Matrix_WGS84(satlatlonalt,az,ze)
        else:
            raise Exception('Not valid S multiprocessing choice')
    else:
        raise Exception('Not valid shperical earth choice')

    #S = S[0:len(h),0:len(h)]
    #print 'Distance Matrix Calculated'
    Bright= Bright[0:len(h)]
       
    if regu ==0: # Gaussian Elimination
        #print 'Gaussian Elimination Chosen'
        VER = reg.guassian_elinimation(S,Bright)
        if cont==1:
            counts= reg.calc_electron_density(VER,S,O)
        else:
            counts= reg.calc_electron_density(VER,S,0,2)
        #print 'VER & Ne profiles calculated'
    elif regu == 1: # LCurve
        #print 'Tiknonov - Lcurve Chosen'
        VER = reg.Tikhonov(S,Bright,regu_order,0)
        if cont==1:
            counts = reg.calc_electron_density(VER,S,O)
        else: 
            counts = reg.calc_electron_density(VER,S,0,2)
        #print 'VER & Ne profiles calculated'
    elif regu == 1: # GSVD
        #print 'Tikhonov-GSVD Chosen'
        VER = reg.gsvd(S,Bright,regu_order,0)
        if cont==1:
            counts = reg.calc_electron_density(VER,S,O)
        else: 
            counts = reg.calc_electron_density(VER,S,0,2)
        #print 'VER & Ne profiles calculated'
            
            
    return VER,counts,h

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

def find_hm_nm_f2_Diff(Ne_e,rbot_e,Ne,rbot,abs_flag=0.):
    
    '''
    Calculates the difference for the NmF2 and hmF2 of a given electron density profile and the estimated one.
    The estimated vector can contain many profiles and thus we get the average value for both hmF2 and NmF2.
    INPUTS:
        NE_e     - Estimated Electron density profile [cm^{-3}]
        rbot_e   - altitude vector for estimated electron density profile [km]
        Ne       - Original Electron density profile [cm^{-3}]
        rbot     - altitude vector for original electron density profile [km]
        abs_flag - flag that indicates if the return differences will be the absolute values [0:Regular Difference (default), 1:Absolute Difference]
    OUTPUT:
        Hmf2   -  mean altitude of the maximum intesity values of the Ne profile
        Nmf2   -  mean peak intensity value of the altitude profiles
        Hmf2_s -  standard deviation of altitude of the maximum intesity values of the Ne profile
        Nmf2_s -  standard deviation of peak intensity value of the altitude profiles
    NOTES:
        Till now it works for multiple Ne. Need to fix it for single comparison
    HISTORY:
        06-Jun-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''

    hmF2o,Nmf2o = find_hm_Nm_F2(Ne,rbot)
    
    Nmdift = np.zeros(np.size(Ne_e,1))
    Hmdift = np.zeros(np.size(Ne_e,1))
    
    for i in range(0,np.size(Ne_e,1)):
        indexf2= Ne_e[:,i].argmax()

        Hmf2_t = rbot_e[indexf2]
        Nmf2_t = Ne_e[indexf2,i]
        if abs_flag == 0:
            Nmdift[i] = (Nmf2o-Nmf2_t)/(Nmf2o)
            Hmdift[i] =hmF2o-Hmf2_t
        else:
            Nmdift[i] = abs((Nmf2o-Nmf2_t)/(Nmf2o))
            Hmdift[i] = abs(hmF2o-Hmf2_t)

    Hmf2 = np.mean(Hmdift)
    Nmf2 = np.mean(Nmdift)
    
    Hmf2_s = np.std(Hmdift)
    Nmf2_s = np.std(Nmdift)
    
    return Hmf2,Nmf2, Hmf2_s,Nmf2_s

def find_hm_nm_f2_aoe(Ne_e,rbot_e,Ne,rbot,abs_flag=0.):
    
    '''
    Calculates the difference for the NmF2 and hmF2 of a given electron density profile and the estimated one.
    INPUTS:
        NE_e     - Estimated Electron density profile [cm^{-3}]
        rbot_e   - altitude vector for estimated electron density profile [km]
        Ne       - Original Electron density profile [cm^{-3}]
        rbot     - altitude vector for original electron density profile [km]
        abs_flag - flag that indicates if the return differences will be the absolute values [0:Regular Difference (default), 1:Absolute Difference]
    OUTPUT:
        Nmdift   -  all of the maximum intesity values of each Ne profile
        Hmdift   -  all of  peak intensity value of the altitude profiles
    NOTES:
        NONE
    HISTORY:
        27-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''

    hmF2o,Nmf2o = find_hm_Nm_F2(Ne,rbot)
    
    Nmdift = np.zeros(np.size(Ne_e,1))
    Hmdift = np.zeros(np.size(Ne_e,1))
    
    for i in range(0,np.size(Ne_e,1)):
        indexf2= Ne_e[:,i].argmax()

        Hmf2_t = rbot_e[indexf2]
        Nmf2_t = Ne_e[indexf2,i]
        if abs_flag == 0:
            Nmdift[i] = (Nmf2o-Nmf2_t)/(Nmf2o)
            Hmdift[i] =hmF2o-Hmf2_t
        else:
            Nmdift[i] = abs((Nmf2o-Nmf2_t)/(Nmf2o))
            Hmdift[i] = abs(hmF2o-Hmf2_t)

    
    return Nmdift,Hmdift

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
    UT_Date_Time = f.createVariable('UT_DATE_TIME', 'i',('time',))
    UT_Date_Time[:] = [dn.year, dn.month, dn.day, dn.hour, dn.minute, dn.second]
    UT_Date_Time.units = 'Year-Month-Day-Hour-Mins-Secs'
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
    ze_f.units = 'km'
    ze_f.long_name = 'Viewing geometry of observation, zenith'
    
    # Tangent point coordinates
    tanlat = f.createVariable('FUV_TANGENT_LATITUDE','float',('profile',))
    tanlat[:] = tanlatlonalt[:,0]
    tanlat.units = 'degrees'
    tanlat.long_name = 'Tangent point latitude coordinate'
    
    tanlon = f.createVariable('FUV_TANGENT_LONGITUDE','float',('profile',))
    tanlon[:] = tanlatlonalt[:,1]
    tanlon.units = 'degrees'
    tanlon.long_name = 'Tangent point longitude coordinate'
    
    tanalt = f.createVariable('FUV_TANGENT_ALTITUDE','float',('profile',))
    tanalt[:] = tanlatlonalt[:,2]
    tanalt.units = 'km'
    tanalt.long_name = 'Tangent point altitude coordinate'
    
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

    bright = FUV_1356_IMAGE[limb,STRIPE][0][::-1]
    # h needs to be changed when non-symmetric earth is assumed. For now we have symmetry
    h = FUV_TANGENT_ALTITUDES[limb,STRIPE][0][::-1]
    # That might not work well now if non-symmetric earth is assumed from Scott
    O = FUV_TANGENT_O1[limb,STRIPE][0][::-1]

    ver,Ne,h = FUV_Level_2_OutputProduct_Calculation(bright,h,satlatlonalt,az,ze,O,cont =1, regu=1,regu_order = 2,S_mp = 1)
    hmF2,NmF2 = find_hm_Nm_F2(Ne,h)
    NmF2 = NmF2*100
    
    FUV_Level_2_OutputProduct_NetCDF(dn,satlatlonalt,az,ze,tanlatlonalt,bright,ver,Ne,NmF2,hmF2,path_output)
    
    print 'LVL2.5 Processing Terminated. File created!'

