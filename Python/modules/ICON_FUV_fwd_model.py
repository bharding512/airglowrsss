# Code containing the forward model for the ICON FUV instrument. Includes options for (non)spherical (non)symmetric earth-atmosphere.

# NOTE OF CAUTION: for the calculation of the Brightness multiple cores are used. Need to fix that option if running in a single core machine!

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
from sys import exit
from scipy.interpolate import interp1d,interp2d

import FUV as fuv
from scipy.io import netcdf
import FUV_Distance_Matrix as dmat



# ICON FUV
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
                  'Sensitivity': 0.0873,         # combined transmittance of all optics. [counts/res_cell/s/R]
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
def calc_1356_nighttime(lat,lon,alt,dn,Ne_scaling = 1.,testing=0):
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
        Testing - Flag indicate using a Gaussian Ne instead of IRI [0: Pyglow, 1: Gaussian]
    OUTPUT:
        RR  - volume emission rate due to radiative recombination (1/cm^3/s)
        MN  - volume emission rate due to mutual neutralization (1/cm^3/s) - 0 in case of testing
        Ne  - Electron Density Profile (1/cm^3/s) - 0 in case of testing
        O   - Oxygen Profile (1/cm^3/s) - 0 in case of testing
    NOTES:
    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu)
        13-Jul-2015: Add O as return value
        02-Sep-2015: Add testing input
        14-Sep-2015: Assumed Ne==Op or testing reasons.
    CALLS:
        Pyglow: IRI-MSIS
    '''
    # Coefficients from Table 3 of Melendez-Alvira (1999)
    b1356 = 0.54    # yield parameter (unitless)
    a1356 = 7.3e-13 # radiative recombination rate (cm^3/s)
                    # consider that to be constant through altitude, normally it should change with temperature
    k1 = 1.3e-15    # radiative attachment rate (cm^3/s)
    k2 = 1e-7       # ion-ion neutralization rate (cm^3/s)
    k3 = 1.4e-10    # ion-atom neutralization rate (cm^3/2)

    if testing==0:
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
        #RR = a1356*Ne*O_p  # radiatvie recombination (1/cm^3/s)
        #MN  = (b1356*(k1*k2)) *((Ne*O*O_p)/(k2*O_p + k3*O)) # mutual neutralization (1/cm^3/s)
        RR = a1356*Ne**2  # radiatvie recombination (1/cm^3/s)
        MN  = (b1356*(k1*k2)) *((Ne**2*O)/(k2*Ne + k3*O)) # mutual neutralization (1/cm^3/s)
    else:
        Ne = 1.6657129144319999*10**6*np.exp(-(alt-350)**2/100**2)
        # Calcualte radiative recombination (equation 17) and mutual neutralization (equation 18)
        RR = a1356*Ne**2  # radiatvie recombination (1/cm^3/s)
        MN = 0
        O = 0

    return RR,MN,Ne,O

def calculate_VER_1356_nighttime(satlat,satlon,satalt,dn,Ne_scaling = 1.,testing = 0,cont =2):
    '''
    Step along a desired look direction from a given observing location and
    calculate the 135.6-nm intensity. Lines-of-sight are calculated using a spherical earth.
    INPUTS:
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        Ne_scaling  - Scaling factor for Ne [default 1 => No scaling]
        Testing     - Flag indicate using a Gaussian Ne instead of IRI [0: Pyglow, 1: Gaussian]
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
    OUTPUTS:
        VER         - A VER altitude profile for given lat/lon coords
        RR          - A RR altitude profile for given lat/lon coords
        MN          - A MN altitude profile for given lat/lon coords
        Ne          - A Electron density altitude profile for given lat/lon coords
        O           - A O density altitude profile for given lat/lon coords
    NOTES:

    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu)
        08-Jul-2015: Added check for satalt scalar values (Dimitrios Iliou)
        13-Jul-2015: Add O to calc_1356_nighttime to have O tanget altitude profile
        02-Sep-2015: Add testing input - Gaussian Ne instead of IRI
        10-Sep-2015: Add the contribution parameter as an input.
    CALLS:
        - calc_1356_nighttime
    '''

    VER = np.zeros(np.size(satalt))
    MN = np.zeros(np.size(satalt))
    VER_true = np.zeros(np.size(satalt))
    NE = np.zeros(np.size(satalt))
    O = np.zeros(np.size(satalt))

    if (np.size(satalt)==1):
        VER,MN,NE,O= calc_1356_nighttime(satlat,satlon,satalt,dn,Ne_scaling,testing)
    else:
        for i in range(0,np.size(satalt)):
            VER[i],MN[i],NE[i],O[i]= calc_1356_nighttime(satlat,satlon,satalt[i],dn,Ne_scaling,testing)

    if cont == 1:
        VER_true = VER + MN
    else:
        VER_true = VER

    return VER_true, VER,MN, NE, O

# ICON FUV
def calculate_pixel_1356_nighttime(ze,az,satlat,satlon,satalt,dn,cont=1,Ne_scaling=1., step=10. , testing = 0):
    '''
    Step along a desired look direction from a given observing location and
    calculate the 135.6-nm intensity. Lines-of-sight are calculated using a spherical earth.
    INPUTS:
        ze          - zenith angle of look direction (deg)
        az          - azimuth angle of look direction (0 is north, 90 is east) (deg)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
        Ne_scaling  - Scaling factor for Ne [default 1 => No scaling]
        step        - resolution of the integration along the line of sight (km, default = 10 km)
        testing     - Flag indicate using a Gaussian Ne instead of IRI [0: Pyglow, 1: Gaussian]
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
        02-Sep-2015: Add testing input - Gaussian Ne instead of IRI
	12-Jul-2017: Changed ze/az input to deg
    CALLS:
        -calc_1356_nighttime
    '''
    # Define constants and step sizes.  We work in a 2D coordinate system at this point
    RE = 6371.      # radius of earth (km)
    x = 0.          # This measures perpendicular to the satellite nadir
    y = RE + satalt # This measures parallel to satellite nadir
    dx = np.sin(ze*np.pi/180.) # unit step in x
    dy = np.cos(ze*np.pi/180.) # unit step in y

    IRR = 0.
    IMN = 0.
    Ne = 0
    Rayleigh = 0.
    alti =[]
    # We will trace along the ray path until we get back to the satellite altitude
    # NOTE: (1) Generalize this by calculating points and using a for loop.
    #       (2) Generalize this by allowing a ray path to be passed in.
    while (np.sqrt(x**2 + y**2) <= RE+satalt): # until you get far enough away from earth that it doesn't matter

        alt = np.sqrt(x**2 + y**2) - RE
        lat = np.nan # for now, this doesn't matter
        lon = np.nan # for now, this doesn't matter
        alti.append(alt)
        # After this altitude the Raypath goes through the disk so no calculation is done for these altitudes
        #if alt < 150:    # 150 USUALLY 100
            #break

        # Calculate 1356 Emission for a single point
        VER,MN,Ne,_ = calc_1356_nighttime(satlat,satlon,alt,dn,Ne_scaling,testing)

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
def calculate_pixel_1356_nighttime_WGS84(ze,az,satlat,satlon,satalt,dn,symmetry = 0,cont=1,Ne_scaling=1., step=10.,testing = 0,total_distance = 5000.):
    '''
    Step along a desired look direction from a given observing location and
    calculate the 135.6-nm intensity. Lines-of-sight are calculated using a non-spherical non-symmetric earth.
    INPUTS:
        ze          - zenith angle of look direction (deg)
        az          - azimuth angle of look direction (0 is north, 90 is east) (deg)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        symmetry    - flag indicating if symmetry [0] or non-symmetry[1] will be used
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
        Ne_scaling  - Scaling factor for Ne [default 1 => No scaling]
        step        - resolution of the integration along the line of sight (km, default = 10 km)
        testing     - Flag indicate using a Gaussian Ne instead of IRI [0: Pyglow, 1: Gaussian]
        total_distance - length of the projected line for each raypath(km).
    OUTPUTS:
        Brightness  - the intensity of the integrated emission (R) or 0 if the contribution is not set correctly
    NOTES:

    HISTORY:
        08-Jul-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        02-Sep-2015: Add testing input - Gaussian Ne instead of IRI
        04-Sep-2015: Add total_distance input - goes to the WGS84 raypath calculations
    CALLS:
        -calc_1356_nighttime
        -ICON.project_line_of_sight
    '''

    # Initialize the integration variables to zeros
    IRR = 0.
    IMN = 0.
    Ne = 0
    Rayleigh = 0.

    satlatlonalt = [satlat,satlon,satalt]
    #print '%f,%f'%(az,ze)

    # Call the function for a single az and ze angle. It returns all the stepping points along the line
    xyz, latlonalt = ic.project_line_of_sight(satlatlonalt, az,ze,step,total_distance)

    # For all the points on the line
    for j in range(0,np.size(latlonalt,1)):

        # Continue iterating till you reach the conjugate point of the satellite.
        # This is not necessery since for high altitudes the contributions are not significant.
        #if (latlonalt[2,j]>latlonalt[2,0]):
            #break

        # Calculate 1356 Emission for a single point
        # Return Radiative Recombination, Mutual Neutralization and Ne value for each point on the line
        # The Ne can be scalled to test algorithm for various pertubations
        # The symmetry check changes the call so the lat lon coordinates will be the same as nadir direction or not.
        if symmetry == 0:
            VER,MN,Ne,_ = calc_1356_nighttime(satlat,satlon,latlonalt[2,j],dn,Ne_scaling,testing)
        elif symmetry == 1:
            VER,MN,Ne,_ = calc_1356_nighttime(latlonalt[0,j],latlonalt[1,j],latlonalt[2,j],dn,Ne_scaling,testing)

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
def get_Photons_from_Brightness_1356_nighttime(ze,az,satlat,satlon,satalt,dn,symmetry =0,shperical=1,exposure=0.,testing = 0,cont=1,TE=0.,Ne_scaling = 1.,step = 10, total_distance = 5000., stripes_used = 0):
    '''
    Calls 'calculate_pixel_1356_nighttime' which calculates the Brightness through integrated VER for a given zenith angle.
    INPUTS:
        ze          - zenith angle of look direction (deg)
        az          - azimuth angle of look direction (0 is north, 90 is east) (deg)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
        symmetry    - flag indicating if symmetry [0] or non-symmetry[1] will be used
        shperical   - flag indicating if spherical [0] or WGS84 model[1] will be used
        exposure    - Exposure time of the FUV CCD, if input is zero the default parameters are loaded
        testing     - Flag indicate using a Gaussian Ne instead of IRI [0: Pyglow, 1: Gaussian]
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
        TE          - Total Efficiency of the FUC optical system, if input is zero the default parameters are loaded
        Ne_scaling - Scaling factor for Ne [default 1 => No scaling]
        step        - resolution of the integration along the line of sight (km, default = 10 km)
        total_distance - length of the projected line for each raypath(km).
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
        02-Sep-2015: Add testing input - Gaussian Ne instead of IRI
        04-Sep-2015: Add total_distance input - goes to the WGS84 raypath calculations
	12-Jul-2017: Changed ze/az input to degrees
    CALLS:
        -get_FUV_instrument_constants
        -calculate_pixel_1356_nighttime
        -calculate_pixel_1356_nighttime_WGS84
    '''

    # Check if instrument parameters are zero from input to load default values
    params = get_FUV_instrument_constants()
    if exposure==0:
        exposure = params['exposure']
    if TE==0:
        TE =  params['Sensitivity']
    if stripes_used==0:
        stripes_used = params['stripes_used']

    # Get Brightness for a given zenith angle
    if shperical == 0:
        Brightness = calculate_pixel_1356_nighttime(ze,az,satlat,satlon,satalt,dn,cont,Ne_scaling,step,testing)
    elif shperical == 1:
        Brightness = calculate_pixel_1356_nighttime_WGS84(ze,az,satlat,satlon,satalt,dn,symmetry,cont,Ne_scaling,step,testing,total_distance)

    # Number of pixels in rescell = 8
    r = TE * exposure * stripes_used

    photons = r*Brightness

    return Brightness,photons

#Multiprocessing
def get_Photons_from_Brightness_1356_nighttime_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return get_Photons_from_Brightness_1356_nighttime(*a_b)

# ICON FUV
def get_Photons_from_Brightness_Profile_1356_nighttime(ze,az,satlat,satlon,satalt,dn,symmetry=0,shperical=1,exposure=0.,testing = 0,cont=1,TE=0.,Ne_scaling = 1.,step = 10,total_distance = 5000.,stripes_used = 0,proc=16):
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
        testing     - Flag indicate using a Gaussian Ne instead of IRI [0: Pyglow, 1: Gaussian]
        cont        - contribution factors, cont=1 -> RR+MN, 2-> RR
        TE          - Total Efficiency of the FUC optical system, if input is zero the default parameters are loaded
        Ne_scaling - Scaling factor for Ne [default 1 => No scaling]
        step        - resolution of the integration along the line of sight (km, default = 10 km)
        total_distance - length of the projected line for each raypath(km).
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
        26-Aug-2015: Exception for Pool added.
        02-Sep-2015: Add testing input - Gaussian Ne instead of IRI
        04-Sep-2015: Add total_distance input - goes to the WGS84 raypath calculations
    CALLS:
        -get_FUV_instrument_constants
        -Opens multiprocessing pool
        -get_Photons_from_Brightness_1356_nighttime_star
    '''
    try:
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
        job_args = [(ze[i],az[i] ,satlat,satlon,satalt, dn,symmetry,shperical,exposure,testing,cont,TE,Ne_scaling,step,total_distance,stripes_used) for i in range(0,len(ze))]
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

    return Rayl,photons

# ICON FUV
def add_noise_to_photon_and_brightness(photons,exposure=0., TE=0., stripes_used = 0, reps=1, ret_cov=False):
    '''
    Returns Brightness and Electron Density profile(s) after adding shot noise
    INPUTS:
        photons     - number of photons (counts) [Noiseless Electron density profile]
        Exposure    - FUV Exposure time (s)   [if zero on input loads default value]
        TE          - Total Efficiency        [if zero on input loads default value]
        stripes_used- number of stripes used of the CCD [default 1] [if zero on input loads default value]
        reps        - Number of noise realizations to return [default = 1]
        ret_cov     - If True, also return the covariance matrix of Brightness
    OUTPUTS:
        Brightness- Noisy Brightess Profile(s) [Rayleighs]
        Shot_noise- Noisy Photon Profile(s) [count]
        ret_cov   - (OPTIONAL) covariance matrix of Brightness. Only outputted if ret_cov == True
    Comments:
        noise - input noise, Poisson Dist (counts)
    HISTORY:
        12-Dec-2015: Written by Dimitrios Iliou (iliou2@illinois.edu)
        01-Apr-2015: pixels_per_rescell -> stripes used
        08-Oct-2015: changed the svs.poisson. removed the second argument
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
    sigma_Rayl = np.zeros(np.size(photons)) # 1-sigma uncertainty

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
                #shot_noise[rep,i] = stats.poisson.rvs(photons[i],1)
                shot_noise[rep,i] = stats.poisson.rvs(photons[i])

        Rayl_[rep,:] = shot_noise[rep,:]/(TE*exposure*stripes_used)

    for i in range(np.size(photons)):
        sigma_Rayl[i] = np.sqrt(photons[i]) / (TE*exposure*stripes_used)
    cov_Rayl = np.diag(sigma_Rayl**2)

    if ret_cov:
        return Rayl_, shot_noise, cov_Rayl
    else:
        return Rayl_,shot_noise

def run_forward_modelling(satlatlonalt,date,ze=0.,az=0.,symmetry = 0.,shperical=1, exp=12,reps=1000,sens=0,cont=2,testing = 0,Ne_sc=1.,step=10,total_distance = 6000.,stripes_used=0,proc=16,low_ta = 150.):
    '''
    Top forward model modules. Creates brightness profile and adds noise returning the noisy Brighntess profile(s) to
    be used for the inversion process.
    INPUTS:
        satlatlonalt - vector containing the satellite coordinates [lat-lon-alt] (km)
        date         - datetime input
        ze           - zenith angle of look direction (degrees)
        az           - azimuth angle of look direction (degrees)
        symmetry     - flag indicating if symmetry [0] or non-symmetry[1] will be used
        shperical    - flag indicating if spherical [0] or WGS84 model[1] will be used
        exp          - exposude time in order to be used as proxy for SNR increase or decrease (s)
        reps         - number of noise realizations for the Brightness profile. Determines the size of the ouput vector.
        sens         - Sensitivity of the FUV optical system, if input is zero the default parameters are loaded
        cont         - contribution factors, cont=1 -> RR+MN, 2-> RR
        testing     - Flag indicate using a Gaussian Ne instead of IRI [0: Pyglow, 1: Gaussian]
        Ne_sc        - Scaling factor for Ne [default 1 => No scaling]
        step         - resolution of the integration along the line of sight (km, default = 10 km)
        total_distance - length of the projected line for each raypath(km).
        stripes_used - number of stripes used of the CCD [default 1] [if zero on input loads default value]
        proc         - Number of cores used for multiprocessing [default 16 cores]
        low_ta       - Lower tangent altitude that we need to consider for our measurements (km)[default = 150 km :limb, for sublimb we can put -500 to cover all zenith]
    OUTPUT:
        NE           - Electron Density profile corresponding the the specified datetime and tangent altitude (cm^-3)
        Bright       - Brightness profile (Rayleigh)
        Bright_n     - Noisy Brightness profile (Rayleigh) - Size determined from the number of noise realizations
        h            - Tangent altitudes (mid cell points) (km)
        rbot         - bottom boundaries[km] (tangent altitudes) (km)
        O            - Oxygen profile corresponding the the specified datetime and tangent altitude (cm^-3)
        VER_true     - VER profile corresponding the the specified datetime and tangent altitude (ph/cm^-3)
        h_loc        - Location of mid altitude points [lat,lon,alt](degrees,degrees,km)
        h_loc_bot    - Location of bot tangent altitudes [lat,lon,alt](degrees,degrees,km)
        Sigma        - Covariance matrix from photons
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
        11-Aug-2015: Added VER_true,h_coord on the ouput.
        18-Aug-2015: Added Sigma variable on the output
        02-Sep-2015: Add testing input - Gaussian Ne instead of IRI
        03-Sep-2015: Added find_mid_cell function and fixed outputs to contain mid, bot vectors for both Spherical and Ellipsoid earth.
        04-Sep-2015: Add total_distance input - goes to the WGS84 raypath calculations
        25-Jan-2016: Add low_ta
    CALLS:
        -angle2tanht
        -ICON.tangent_point
        -get_Photons_from_Brightness_Profile_1356_nighttime
        -add_noise_to_photon_and_brightness
        -calculate_VER_1356_nighttime
        -find_mid_cell
    '''

    # Satellite coordinates
    satlat = satlatlonalt[0]
    satlon = satlatlonalt[1]
    satalt = satlatlonalt[2]

    if np.size(ze)==1 and np.size(az)==1:
        if ze==0 and az ==0:
            # Load default az and ze angles from the fov of the instrument
            az,ze = get_azze_default()
        elif ze==0:
            _,ze = get_azze_default()
            ze = ze[0:len(az)]
        elif az==0:
            az,_ = get_azze_default()
            az = az[0:len(ze)]
    elif np.size(az)==1:
        if az==0:
            az,_ = get_azze_default()
            az = az[0:len(ze)]
    elif np.size(ze)==1:
        if ze==0:
            _,ze = get_azze_default()
            ze = ze[0:len(az)]

    az_v = np.deg2rad(az)
    ze_v = np.deg2rad(ze)


    '''
    Locate limb lower bound ~150km
    '''

    check = 0
    '''
    # Locate limb lower bound ~130km
    limb = 130

    if (shperical==0):
        h = np.zeros(np.size(ze))
        RE = 6371.
        h = ic.angle2tanht(ze_v, satalt, RE)
        h = h[np.where(h>limb)]
        disk = len(h)
    else:
        # Tangent_poing function returns coordinates [lat-lon-lat] for the tangent point
        h_coord = np.zeros((len(ze),3))
        # Locate limb lower bound ~150km
        for i in range(0,len(ze)):
            h_coord[i,:] = ic.tangent_point(satlatlonalt,az[i],ze[i])
        h_loc = h_coord[np.where(h_coord[:,2]>=limb),:]
        h = h_coord[np.where(h_coord[:,2]>=limb),2]
        h = h[0,:]
        disk = len(h)


    # Resize vectors for limb FOV
    ze = ze[0:disk]
    az = az[0:disk]

    ze_v = ze_v[0:disk]
    az_v = az_v[0:disk]

    '''
    try:
        h,rbot,_,h_loc_bot,h_loc = find_mid_cell(satlatlonalt,ze,az,shperical,low_ta)
        #h,rbot,_,h_loc_bot,h_loc = find_mid_cell(satlatlonalt,ze,az,shperical,150.)

        ze = ze[0:len(h)]
        az = az[0:len(h)]

        ze_v = ze_v[0:len(h)]
        az_v = az_v[0:len(h)]

        # Initialize vectors
        photons = np.zeros(np.size(ze_v)) # Number of Counts without Noise
        Bright = np.zeros(np.size(ze_v))
        #Bright_n = np.zeros(np.size(ze_v))
        Counts = np.zeros(np.size(ze_v)) # Number of Counts with Noise

        '''
        Multicore Processing for Calculate Brighness
        '''
        ### ___ I NEED TO ADD THE SPHERICAL THING HERE

        # Takes s/c coordinates as input for poing location to calculate the brightness profile. This must change to tan. point
        if (shperical==0):
            Bright,photons = get_Photons_from_Brightness_Profile_1356_nighttime(ze_v,az_v,satlat,satlon,satalt,date,symmetry,shperical,exp,testing,cont,sens,Ne_sc,step,total_distance,stripes_used,proc)
        else:
            Bright,photons = get_Photons_from_Brightness_Profile_1356_nighttime(ze,az,satlat,satlon,satalt,date,symmetry,shperical,exp,testing,cont,sens,Ne_sc,step,total_distance,stripes_used,proc)

        '''
        Calculate Noisy Counts
        reps determines the number of noise profiles that the function will produce
        '''
        Bright_n= np.zeros((reps,np.size(ze_v)))
        Bright_n,Counts = add_noise_to_photon_and_brightness(photons,exp,sens,stripes_used,reps)
        Bright_n = np.transpose(Bright_n)

        Sigma = np.diag(1/np.sqrt(photons))

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
            VER_true,VER,MN,NE,O = calculate_VER_1356_nighttime(satlat,satlon,h,date,Ne_sc,testing,cont)
            h_coord = 0
        else:
            if symmetry == 0:
                VER_true,_,_,NE,O = calculate_VER_1356_nighttime(satlat,satlon,h,date,Ne_sc,testing,cont)
            else:
                for i in range(0,len(h)):
                    VER_true[i],_,_,NE[i],O[i] = calculate_VER_1356_nighttime(h_loc[i,0],h_loc[i,1],h_loc[i,2],date,Ne_sc,testing,cont)


    except Exception:
        print "Something Happened :("
        print type(inst)
        print inst
        exit(1)
    return NE,Bright,Bright_n,h,rbot,O,VER_true,h_loc,h_loc_bot,Sigma


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
    CALLS:
     -get_FUV_instrument_constants
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
    return '%2.1f' % (x*sc)

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



## TIEGCM Forward MODEL
def TIEGCM_Altitude_interpolate(Ne,rbot,xnew):
    '''
    Given a TIEGCM altitude profile we make cubic intepolation and linear extrapolation to calculate the values for the altitude vector resulting from the tangent altitudes
    INPUTS:
        Ne       - Original Electron density profile [cm^{-3}]
        rbot     - altitude vector for original electron density profile [km]
        xnew     - Altitude vector that we want to intepolate the values to
    OUTPUT:
        NE_inter   -  Interpolated altitude profile for the estimated electron density
        xnew       -  Interpolated altitude vector for the true electron density
    NOTES:
    HISTORY:
        16-May-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''

    # TIEGCM vectors
    y_true = Ne
    x_true = rbot

    # Create two interpolators, one linear for the extrapolation and one cubic for the intepolation
    # Cubic is more preferable since we want smoothness
    f_linear = interp1d(x_true, y_true, kind='linear',fill_value='extrapolate')
    f_cubic = interp1d(x_true,y_true,kind ='cubic')

    # Calculate the intepolations
    linear = f_linear(xnew)
    #cubic = f_cubic(xnew[np.where(xnew<max(x_true))])

    ## HERE I NEED TO ADD LINEAR INTERPOLATION BENEATH THE TANGENT POINT. MAYBE EXTRAPOLATE BELOW THE MIN VALUE OF TIEGCM
    cubic_indices = np.intersect1d(np.where(xnew>min(x_true)),np.where(xnew<max(x_true)))
    cubic = f_cubic(xnew[cubic_indices])

    # Replace the values in the linear vector with the corresponding values of the cubic one to get the correct dimensions for the altitude vector wanted
    # I need to find the indices that correspond to the the cubic and
    #linear[len(linear)-len(cubic):] = cubic
    linear[cubic_indices] = cubic

    linear[np.where(linear<=0)]=0

    return linear,xnew

def fix_lat_lon_interp_indices(x,value):

    xnew = np.zeros(len(x)+1)

    les = np.where(x<value)[0][-1]
    mor = np.where(x>value)[0][0]

    xnew[:les+1] = x[:les+1]
    xnew[les+1] = value
    xnew[les+2:] = x[mor:]

    return xnew

def TIEGCM_LatLon_2Dinterpolate(Ne,alt,lat_vec,lon_vec,lat,lon):
    '''
    Given a TIEGCM altitude profile we make cubic intepolation and linear extrapolation to calculate the values for the altitude vector resulting from the tangent altitudes
    INPUTS:
        Ne       - Original Electron density profile [cm^{-3}] TIEGCM
        alt      - Altitude profile TIEGCM [cm]
        lat_vec  - Latitude value vector
        lon_vec  - Longitude value vector
        lat      - latitude coordinate that we want to interpolate
        lon      - longitude coordinate that we want to interpolate
    OUTPUT:
        NE_interp   -  Interpolated altitude profile for the estimated electron density
        alt_interp  -  Interpolated altitude vector for the true electron density
        lat_new     - New latitude vector that includes the location for which we interpolated
        lon_new     - New longitude vector that includes the location for which we interpolated
    NOTES:
    HISTORY:
        16-May-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)
    '''

    if lat in lat_vec:
        lat_new = lat_vec
    else:
        lat_new = fix_lat_lon_interp_indeces(lat_vec,lat)
    if lon in lon_vec:
        lon_new = lon_vec
    else:
        lon_new = fix_lat_lon_interp_indeces(lon_vec,lon)

    Ne_interp = np.zeros((np.size(Ne,0),len(lat_new),len(lon_new)))

    alt_interp = np.zeros((np.size(alt,0),len(lat_new),len(lon_new)))

    for i in range(0,np.size(Ne,0)):
        f = interp2d( lon_vec,lat_vec,Ne[i,:,:], kind='cubic')
        Ne_interp[i,:,:] = f(lon_new,lat_new)

        f_alt = interp2d( lon_vec,lat_vec,alt[i,:,:], kind='cubic')
        alt_interp[i,:,:] = f_alt(lon_new,lat_new)

    return Ne_interp,alt_interp,lat_new,lon_new


def find_nearest(array,value):
    '''
    Function that finds the closest value in the array.
    This will be used to find the interval in which we want to interpolate
    in when we call TIEGCM
    INPUTS:
        array  - array that contains the values that we have
        value  - value that we want to find its closest neighbour
    OUTPUTS:
        inx    - index position of the matrix
    HISTORY:
        16-May-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)

    '''
    idx = (np.abs(array-value)).argmin()
    return idx


def TIEGCM_Brightness_Calculation(satlatlonalt,ze,az,hour = 0,spherical=0,limb = -500.,path = '/home/dimitris/public_html/Datafiles/TIEGMC/tiegcm_icon_merg2.0_DpertCbgrd.s_081.nc'):

    '''
    Foward model for  simulated Brightness measurements using TIEGCM files
    INPUTS:
        satlatlonalt    - Satellite coordinates [lat,lon,alt]
        ze              - zenith viewing geometry to calculate tangent altitudes
        az              - azimuth viewing geometry to calculate tangent altitudes
        hour            - Hour index for the TIEGCM model
        spherical       - Spherical earth assumption [0->Spherical, 1-> Non-Spherical]
        limb            - Lower tangent altitude
        path            - Path containing the TIEGCM file
    OUTPUTS:
        Brightness_tiegcm - Simulated Brightness
        ver_tiegcm        - Simulated VER
        Ne_tiegcm         - Electron density extracted from TIEGCM
        tan_alt           - Calculated tangent altitudes
    HISTORY:
        18-May-2016: Written by Dimitrios Iliou (iliou2@illinois.edu)
    NOTES:


    '''

    # Read the NETCDF file containing the TIEGCM model

    data = netcdf.netcdf_file(path, 'r')

    Ne = data.variables['NE'][:]
    Ne = Ne[:,::-1,:,:]
    lon = data.variables['lon'][:]
    lat = data.variables['lat'][:]
    day = data.variables['day'][:] # 81
    year = data.variables['year'][:] # 2009
    lev = data.variables['lev'][:]
    ilev = data.variables['ilev'][:]
    mtime = data.variables['mtime'][:]
    time_ = data.variables['time'][:]
    Z = data.variables['Z'][:]
    Z = Z[:,::-1,:,:] * 1e-5

    f107d = data.variables['f107d'][0]

    data.close()

    # Assumes that the start date is first month and first date of the year declared on the TIEGCM file
    start_date = datetime(year[0],1,1,0,0,0)
    dn =  start_date + timedelta(minutes=time_[hour])

    a1356 = 7.3e-13

    h_temp,_,_,_,_ = fuv.find_mid_cell(satlatlonalt = satlatlonalt,ze = ze,az =az,spherical = spherical,limb=limb)

    if spherical == 0:
        S,_,_,_ = dmat.create_cells_Matrix_spherical_symmetry(np.deg2rad(ze),satlatlonalt[2])
    else:
        S,rmid = dmat.Calculate_D_Matrix_WGS84_mp(satlatlonalt,az,ze)

    Ne_interp,Z_interp,lat_new,lon_new = TIEGCM_LatLon_2Dinterpolate(Ne[hour,:,:,:],Z[hour,:,:,:],lat,lon,satlatlonalt[0],satlatlonalt[1])

    lat_index = np.where(lat_new==satlatlonalt[0])[0][0]
    lon_index = np.where(lon_new==satlatlonalt[1])[0][0]

    Ne_tiegcm,tan_alt = TIEGCM_Altitude_interpolate(Ne_interp[:,lat_index,lon_index],Z_interp[:,lat_index,lon_index],h_temp)
    ver_tiegcm = a1356*Ne_tiegcm**2
    Brightness_tiegcm= S.dot(ver_tiegcm)

    return Brightness_tiegcm,ver_tiegcm,Ne_tiegcm,tan_alt
