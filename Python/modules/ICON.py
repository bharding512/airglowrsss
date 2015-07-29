# Common code for ICON users at Illinois. Geometrical transformations, airglow code, etc.

import numpy as np
from pyglow import pyglow
from datetime import datetime, timedelta
from scipy import stats
import multiprocessing
from multiprocessing import Pool
import time
import math

# ICON FUV
def get_FUV_instrument_constants():
    '''
    Notes:
        12-Dec-2014: Last known parameters for ICON FUV   
    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu)
        01-Apr-2015: Last known parameters for ICON FUV - Rescell Dimensions  
    '''

    instrument = {    'npixelx': 64,            # number of pixels per slice of the interferogram in a rescell(horizontal) [512px total]
                      'npixely': 256,           # number of rescells in vertical (altitude) direction on CCD.
                'aperture_area': 0.006*0.032,   # [cm^2]
                   'coneangle1': 18,            # fov of instrument in horizontal direction, [deg]
                   'coneangle2': 24,            # fov of instrument in vertical direction, [deg]
                 'coneangle1_r': 0.314159265359,# fov of instrument in horizontal direction, [rad]
                 'coneangle2_r': 0.418879020479,# fov of instrument in vertical direction, [rad]
                    'exposure' : 12,            # exposure [sec]
             'Total_efficiency': 0.083,         # combined transmittance of all optics. [counts/res_cell/s/R]
                        'fov_l': 98,            # Field of View (Vertical) Lower Bound (degrees)
                        'fov_u': 123,           # Field of View (Vertical) Upper Bound (degrees)
                       'fovr_l': 1.71042266695, # Field of View (Vertical) Lower Bound (radians)
                       'fovr_u': 2.14675497995, # Field of View (Vertical) Upper Bound (radians)
                 'stripes_used': 6,             # Number of stripes used in the ccd the CCD
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
    NOTES:
    HISTORY:
        12-Dec-2014: Written by Dimitrios Iliou (iliou2@illinois.edu)
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

    return RR,MN,Ne

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
    '''
    VER = np.zeros(np.size(satalt))
    MN = np.zeros(np.size(satalt))
    VER_true = np.zeros(np.size(satalt))
    NE = np.zeros(np.size(satalt))
    
    for i in range(0,np.size(satalt)):
        VER[i],MN[i],NE[i]= calc_1356_nighttime(satlat,satlon,satalt[i],dn,Ne_scaling)

    VER_true = VER + MN

    return VER_true, VER,MN, NE

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
        if alt < 90:    # USUALLY 100
            break
       
        # Calculate 1356 Emission for a single point
        VER,MN,Ne = calc_1356_nighttime(satlat,satlon,alt,dn,Ne_scaling)
        
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

# ICON FUV
def get_Photons_from_Brightness_1356_nighttime(ze,az,satlat,satlon,satalt,dn,exposure=0.,TE=0.,cont=1,Ne_scaling = 1.,step = 10, stripes_used = 0):
    '''
    Calls 'calculate_pixel_1356_nighttime' which calculates the Brightness through integrated VER for a given zenith angle.
    INPUTS:
        ze          - zenith angle of look direction (radians)
        az          - azimuth angle of look direction (0 is north, pi/4 is east) (radians)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
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
    '''
    
    # Check if instrument parameters are zero from input to load default values
    params = get_FUV_instrument_constants()
    if exposure==0:
        exposure = params['exposure']
    if TE==0:
        TE =  params['Total_efficiency']
    if stripes_used==0:
        stripes_used = params['stripes_usedl']

    # Get Brightness for a given zenith angle 
    Brightness = calculate_pixel_1356_nighttime(ze,az,satlat,satlon,satalt,dn,cont,Ne_scaling,step)
    
    # Number of pixels in rescell = 8
    r = TE * exposure * stripes_used
  
    photons = r*Brightness
      
    return Brightness,photons

#Multiprocessing
def get_Photons_from_Brightness_1356_nighttime_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return get_Photons_from_Brightness_1356_nighttime(*a_b)

# ICON FUV
def get_Photons_from_Brightness_Profile_1356_nighttime(ze,az,satlat,satlon,satalt,dn,exposure=0.,TE=0.,cont=1,Ne_scaling = 1.,step = 10,stripes_used = 0,proc=16):
    '''
    Calls 'get_Photons_from_Brightness_1356_nighttime' which calculates the Brightness and Electron density by going through integrated VER for all given zenith angles.
    INPUTS:
        ze          - zenith angle of look direction (radians)
        az          - azimuth angle of look direction (0 is north, pi/4 is east) (radians)
        satlat      - latitude of satellite (deg)
        satlon      - longitude of satellite (deg)
        satalt      - altitude of satellite (km)
        dn          - UT date and time to be run (datetime)
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
    '''

    params = get_FUV_instrument_constants()
    if exposure==0:
        exposure = params['exposure']
    if TE==0:
        TE =  params['Total_efficiency']
    if stripes_used==0:
        stripes_used = params['stripes_used']

    photons = np.zeros(np.size(ze)) # Number of Counts without Noise
    Rayl = np.zeros(np.size(ze))
    
    job_args = [(ze[i],az ,satlat,satlon,satalt, dn,exposure,TE,cont,Ne_scaling,step,stripes_used) for i in range(0,len(ze))]
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
    print t1-t0
    
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
        TE =  params['Total_efficiency']
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
                print 'zero value in signal'
            else:
                shot_noise[rep,i] = stats.poisson.rvs(photons[i],1)

        Rayl_[rep,:] = shot_noise[rep,:]/(TE*exposure*stripes_used)
    

        Rayl_[rep,:] = shot_noise[rep,:]/(TE*exposure*stripes_used)
    
    return Rayl_,shot_noise

  
def angle2tanht(theta, H, RE=6371.):
    '''
    Return the tangent altitude of the observation with zenith angle theta.
    INPUTS:
        theta - zenith angle, rad ( pi/2 < theta < pi )
        H - satellite altitude, km
        RE - radius of Earth, km
    OUTPUT:
        h - tangent altitude, km

    '''
    if np.any( theta < np.pi/2 ) or np.any( theta > np.pi ):
        raise Exception('Angle must be between pi/2 and pi, exclusive.')
    h = (H+RE)*np.sin(theta) - RE
    return h 


def calc_earth_radius_WGS84(lat):
    '''
    Step along a desired look direction from a given observing location and
    calculate the 135.6-nm intensity. Lines-of-sight are calculated using a spherical earth.
    INPUTS:
        lat  - latitude of a single point in WGS84 coordinates [deg]
    OUTPUTS:
        Earths Radius - Resulting earths radius for a specific latitude point [km]
    NOTES:
        
    HISTORY:
        24-Mar-2015: Written by Dimitrios Iliou (iliou2@illinois.edu) 
    '''

    # The Earth's equatorial radius a, or semi-major axis, is the distance from its center to the equator
    a = 6378.1370
    #The Earth's polar radius b, or semi-minor axis, is the distance from its center to the North and South Poles
    b = 6356.7523
    # Convert the input from degrees to radians
    lat_r = np.deg2rad(lat)
    
    RE = np.sqrt((a**4*np.cos(lat_r)**2 + b**4*np.sin(lat_r)**2)/(a**2*np.cos(lat_r)**2 + b**2*np.sin(lat_r)**2))

    return RE   

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
    rbot= angle2tanht(theta, Horbit, RE)
    # Define top of each layer.
    
    rtop = rbot.copy()

    rtop[1:] = rbot[:-1]

    rtop[0] = Horbit -1
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

def wgs84constants():
    # http://en.wikipedia.org/wiki/World_Geodetic_System (Nov 10, 2014)
    a = 6378.137 # semi-major axis of earth [km]
    b = 6356.75231424518 # semi-minor axis of earth [km]
    e = np.sqrt(1-b**2/a**2) # eccentricity of earth
    return a,b,e

def ecef_to_wgs84(ecef_xyz):
    '''
    Convert from earth-centered earth-fixed (ECEF)
    coordinates (x,y,z) to WGS-84 latitude, longitude, and altitude.
    INPUTS:
       ecef_xyz - a length-3 array containing the X, Y, and Z locations in ECEF
                  coordinates in kilometers.

    OUTPUTS:
       latlonalt - a length-3 array containing the WGS-84 coordinates:
                   [latitude (degrees), longitude (degrees), altitude (km)]
                   Altitude is defined as the height above the reference
                   ellipsoid.

    HISTORY:
       11-Jun-2006: Initial MATLAB template created by Jonathan J. Makela
       (jmakela@uiuc.edu)
       17-July-2006: Algorithm implemented by Dwayne P. Hagerman
       (dhagerm2@uiuc.ed)
       10-Nov-2014: Translated to Python by Brian J. Harding
       (bhardin2@illinois.edu)
       19-Jan-2015: Changed from iterative to closed-form implementation (BJH)
    '''
    
    
    a,b,e = wgs84constants()

    x = 1.0*ecef_xyz[0]
    y = 1.0*ecef_xyz[1]
    z = 1.0*ecef_xyz[2]
   
    lon = np.arctan2(y,x)
    # longitude is defined in [0,360) or [0,2*pi)
    if lon < 0.0:
        lon += 2*np.pi
     
    ep = np.sqrt((a**2-b**2)/b**2)
    p = np.sqrt(x**2 + y**2)
    th = np.arctan2(z*a,(p*b))

    lat = np.arctan2((z + ep**2 * b * np.sin(th)**3),(p - e**2 * a * np.cos(th)**3))

    N = a/np.sqrt(1 - e**2*np.sin(lat)**2)
    alt = p/np.cos(lat) - N
    lla = np.array([lat*180/np.pi, lon*180/np.pi, alt])

    return lla



def wgs84_to_ecef(latlonalt):
    '''
    Convert from WGS84 coordinates [latitude, longitude, altitude] to 
    earth-centered earth-fixed coordinates (ECEF) [x,y,z]
    
    INPUTS:
       latlonalt - a length-3 array containing the WGS-84 coordinates:
                   [latitude (degrees), longitude (degrees), altitude (km)]    

    OUTPUTS:
       ecef_xyz - a length-3 array containing the X, Y, and Z locations in ECEF
                  coordinates in kilometers.
                  
    HISTORY:
       11-Jun-2006: Initial MATLAB template created by Jonathan J. Makela
       (jmakela@uiuc.edu)
       17-July-2006: Algorithm implemented by Dwayne P. Hagerman
       (dhagerm2@uiuc.ed)
       10-Nov-2014: Translated to Python by Brian J. Harding
       (bhardin2@illinois.edu)    
    '''
    
    a,b,e = wgs84constants()
    
    lat = latlonalt[0]*np.pi/180.
    lon = latlonalt[1]*np.pi/180.
    alt = latlonalt[2]*1.0
    
    x = a * np.cos(lon) / np.sqrt(1 + (1-e**2) * np.tan(lat)**2) + alt*np.cos(lon)*np.cos(lat)
    y = a * np.sin(lon) / np.sqrt(1 + (1-e**2) * np.tan(lat)**2) + alt*np.sin(lon)*np.cos(lat)
    z = a * (1-e**2) * np.sin(lat) / np.sqrt(1 - e**2 * np.sin(lat)**2) + alt*np.sin(lat)
    
    return np.array([x,y,z])
  
  
  
  
    
def ven_to_ecef(latlonalt, ven):
    '''
    Convert a direction vector in local vertical-east-north (VEN) coordinates, defined at the 
    location given in WGS84 coordinates [latitude, longitude, altitude],
    to a unit vector in earth-centered earth-fixed (ECEF) coordinates [x,y,z]
    
    INPUTS:
        latlonalt - a length-3 array containing the WGS-84 coordinates of the location:
                   [latitude (degrees), longitude (degrees), altitude (km)]
        ven - the local direction vector [vertical, east, north] which will be converted to 
             ECEF coordinates. This will be converted to a unit vector
       
    OUTPUTS:
        xyz - the unit direction vector in ECEF coordinates [x,y,z]
        
    HISTORY:
       10-Nov-2014: Written by Brian J. Harding (bhardin2@illinois.edu)
    '''
    
    
    step = 1e-6 # step size for determining transformation from VEN to ECEF 
                # (units of degrees lat/lon or km altitude)
                
    # Convert to float values, in case they are integers
    latlonalt = np.array([float(x) for x in latlonalt])

    # Convert reference location to ECEF
    ecef_xyz = wgs84_to_ecef(latlonalt)

    # Find ecef vectors for single steps in vertical, east, and north
    # (1) vertical step in ECEF coordinates
    latlonalt2 = latlonalt.copy()
    latlonalt2[2] += step
    v_ecef = wgs84_to_ecef(latlonalt2) - ecef_xyz
    v_ecef = v_ecef/np.linalg.norm(v_ecef)
    # (2) east step in ECEF coordinates
    latlonalt2 = latlonalt.copy()
    latlonalt2[1] += step
    e_ecef = wgs84_to_ecef(latlonalt2) - ecef_xyz
    e_ecef = e_ecef/np.linalg.norm(e_ecef)
    # (3) north step in ECEF coordinates
    latlonalt2 = latlonalt.copy()
    latlonalt2[0] += step
    n_ecef = wgs84_to_ecef(latlonalt2) - ecef_xyz
    n_ecef = n_ecef/np.linalg.norm(n_ecef)

    # Create matrix that will transform a vector in VEN to ECEF
    M = np.zeros((3,3))
    M[:,0] = v_ecef
    M[:,1] = e_ecef
    M[:,2] = n_ecef

    # Transform
    s = np.linalg.norm(ven)
    ven = ven / s
    xyz = M.dot(ven)

    return xyz




def ecef_to_ven(latlonalt, ecef):
    '''
    Convert a direction vector in earth-centered earth-fixed (ECEF) coordinates [dx,dy,dz], 
    defined at the location given in WGS84 coordinates [latitude, longitude, altitude],
    to a unit vector in local vertical-east-north (VEN) coordinates.
    
    INPUTS:
        latlonalt - a length-3 array containing the WGS-84 coordinates of the location:
                   [latitude (degrees), longitude (degrees), altitude (km)]
        ecef - the direction vector [dx,dy,dz] which will be converted to 
               VEN coordinates. This will be converted to a unit vector.
       
    OUTPUTS:
        ven - the unit direction vector in [vertical, east, north] coordinates.
        
    HISTORY:
       06-Jan-2015: Written by Brian J. Harding (bhardin2@illinois.edu)
    '''
    step = 1e-6 # step size for determining transformation from ECEF to VEN 
                # (units of degrees lat/lon or km altitude)
                
    # Convert to float values, in case they are integers
    latlonalt = np.array([float(x) for x in latlonalt])

    # Convert reference location to ECEF
    ecef_xyz = wgs84_to_ecef(latlonalt)
    
    # Find ecef vectors for single steps in vertical, east, and north
    # (1) vertical step in ECEF coordinates
    latlonalt2 = latlonalt.copy()
    latlonalt2[2] += step
    v_ecef = wgs84_to_ecef(latlonalt2) - ecef_xyz
    v_ecef = v_ecef/np.linalg.norm(v_ecef)
    # (2) east step in ECEF coordinates
    latlonalt2 = latlonalt.copy()
    latlonalt2[1] += step
    e_ecef = wgs84_to_ecef(latlonalt2) - ecef_xyz
    e_ecef = e_ecef/np.linalg.norm(e_ecef)
    # (3) north step in ECEF coordinates
    latlonalt2 = latlonalt.copy()
    latlonalt2[0] += step
    n_ecef = wgs84_to_ecef(latlonalt2) - ecef_xyz
    n_ecef = n_ecef/np.linalg.norm(n_ecef)

    # Create matrix that will transform a vector in VEN to ECEF
    M = np.zeros((3,3))
    M[:,0] = v_ecef
    M[:,1] = e_ecef
    M[:,2] = n_ecef
    
    # Transform
    s = np.linalg.norm(ecef)
    ecef = ecef / s
    ven = np.linalg.solve(M,ecef)
    
    return ven




def ecef_to_azze(latlonalt, ecef):
    '''
    Convert a direction vector in earth-centered earth-fixed (ECEF) coordinates [dx,dy,dz], 
    defined at the location given in WGS84 coordinates [latitude, longitude, altitude],
    to the azimuth and zenith angles of the direction vector. This function is similar
    to ecef_to_ven.
    
    INPUTS:
        latlonalt - a length-3 array containing the WGS-84 coordinates of the location:
                   [latitude (degrees), longitude (degrees), altitude (km)]
        ecef - the direction vector [dx,dy,dz] which will be converted to 
               azimuth and zenith angles. This will be converted to a unit vector.
       
    OUTPUTS:
        az,ze - the azimuth (degrees East from North) and zenith (degrees down
                from Vertical) angles of the direction vector.
        
    HISTORY:
       24-Feb-2015: Written by Brian J. Harding (bhardin2@illinois.edu)
    '''
    
    # First convert to VEN
    ven = ecef_to_ven(latlonalt, ecef)
    
    # Then convert VEN to az, ze
    ze = np.arccos(ven[0])*180.0/np.pi
    az = np.arctan2(ven[1],ven[2])*180.0/np.pi
    
    return az,ze
    
    
    
  
def azze_to_ecef(latlonalt, az, ze):
    '''
    Convert a line of sight in (az,ze) coordinates, defined at the 
    location given in WGS84 coordinates [latitude, longitude, altitude],
    to a unit vector in earth-centered earth-fixed (ECEF) coordinates [x,y,z]. This
    function is very similiar to ven_to_ecef.
    
    INPUTS:
        latlonalt - a length-3 array containing the WGS-84 coordinates of the location:
                   [latitude (degrees), longitude (degrees), altitude (km)]
        az - direction of line of sight. degrees east of north.
        ze - direction of line of sight. degrees down from zenith
       
    OUTPUTS:
        xyz - the unit direction vector in ECEF coordinates [x,y,z]
        
    HISTORY:
       15-May-2014: Written by Brian J. Harding (bhardin2@illinois.edu)
    '''
    zer = ze*np.pi/180.
    azr = az*np.pi/180.    
    
    # First convert az,ze to VEN
    ven = np.zeros(3)
    ven[0] = np.cos(zer)
    ven[1] = np.sin(zer)*np.sin(azr)
    ven[2] = np.sin(zer)*np.cos(azr)
    # Then run ven_to_ecef
    return ven_to_ecef(latlonalt, ven)
    
    
    


def project_line_of_sight(satlatlonalt, az, ze, step_size = 1., total_distance = 4000.):
    '''
    Starting at the satellite, step along a line of sight, and return an array of 
    equally-spaced points at intervals along this line of sight, in WGS84 latitude,
    longitude, and altitude coordinates.
    
    INPUTS:
        satlatlonalt - array [latitude (deg), longitude (deg), altitude (km)] of the satellite
        az - direction of line of sight. degrees east of north.
        ze - direction of line of sight. degrees down from zenith
        
    OPTIONAL INPUTS:
        step_size - spacing between points in the returned array (km).
        total_distance - length of the projected line (km).
        
    OUTPUTS:
        xyz - array of size 3xN, where N is floor(total_distance/step_size). Contains
              the ECEF coordinates of every point along the line of sight, in step_size
              intervals.
        latlonalt - array of size 3xN, where each column contains the latitude, longitude
                    and altitude corresponding to the column of xyz, in WGS84 coordinates
    
    HISTORY:
        10-Nov-2014: Written by Brian J. Harding (bhardin2@illinois.edu)
    '''
    
    # Convert to radians
    zer = ze*np.pi/180.
    azr = az*np.pi/180.

    # Create unit vector describing the look direction in Vertical-East-North (VEN) coordinates
    lookven = np.array([np.cos(zer), np.sin(zer)*np.sin(azr), np.sin(zer)*np.cos(azr)])
    # Convert satellite location to ecef
    satxyz = wgs84_to_ecef(satlatlonalt)
    # Convert look direction to ecef. This is a unit vector.
    lookxyz = ven_to_ecef(satlatlonalt, lookven)

    # Step along this line of sight
    step_sizes = np.arange(0,total_distance,step_size)
    N = len(step_sizes)
    xyz = np.zeros((3,N))
    latlonalt = np.zeros((3,N))
    for i in range(N):
        xyzi = satxyz + step_sizes[i]*lookxyz
        xyz[:,i] = xyzi
        latlonalt[:,i] = ecef_to_wgs84(xyzi)  
        
    return xyz, latlonalt
 
 
    
    
def tangent_point(satlatlonalt, az, ze, tol=1e-6):
    '''
    Find the location (lat, lon, alt) of the tangent point of a ray from the satellite.
    
    INPUTS:
        satlatlonalt - array [latitude (deg), longitude (deg), altitude (km)] of the satellite
        az - direction of line of sight. degrees east of north.
        ze - direction of line of sight. degrees down from zenith
        
    OPTIONAL INPUTS:
        tol - km, stopping criterion for the solver. Stop when the solution is not moving
              by more than a horizontal distance of tol.
        
    OUTPUTS:
        latlonalt - array [latitude (deg), longitude (deg), altitude (km)] of the tangent location
    
    HISTORY:
        10-Dec-2014: Written by Brian J. Harding (bhardin2@illinois.edu)
    '''
    maxiters = 1e4 # Just in case, so it doesn't hang forever

    # Convert to radians
    zer = ze*np.pi/180.
    azr = az*np.pi/180.

    # Create unit vector describing the look direction in Vertical-East-North (VEN) coordinates
    lookven = np.array([np.cos(zer), np.sin(zer)*np.sin(azr), np.sin(zer)*np.cos(azr)])
    # Convert satellite location to ecef
    satxyz = wgs84_to_ecef(satlatlonalt)
    # Convert look direction to ecef. This is a unit vector.
    lookxyz = ven_to_ecef(satlatlonalt, lookven)

    # Find the step size which minimizes the altitude. This problem is convex,
    # so it's easy.
    def altitude(step_size):
        xyzi = satxyz + step_size*lookxyz
        latlonalt = ecef_to_wgs84(xyzi)
        return latlonalt[2]

    # Or, find the step size for which the slope is zero.
    # Define function to calculate the slope numerically
    def slope(step_size):
        numerical_step = 1e-4 # km
        alt1 = altitude(step_size - numerical_step/2)
        alt2 = altitude(step_size + numerical_step/2)
        return (alt2-alt1)/numerical_step

    # Find two values which straddle the solution.
    s0 = 0. # start at the satellite
    s1 = 1. # iterate on the next value until we've straddled the solution
    f0 = slope(s0)
    f1 = slope(s1)

    # Throw an error if there doesn't appear to be a tangent altitude.
    # (i.e., if the slope is positive, which would happen if the line
    # of sight is looking upward instead of downward)
    if f0 >= 0.0:
        raise Exception('No Tangent Altitude: Altitude not decreasing along line of sight')

    M = 2 # multiplier for line search to find straddle points
    iters = 0
    while(np.sign(f0)==np.sign(f1)):
        iters += 1
        if iters > maxiters:
            raise Exception('Initial line search failed: Maximum iterations reached')
        s1 = s1 * M
        f1 = slope(s1)

    # Straddle points found. Use bisection to converge to the answer.
    iters = 0
    while(abs(s0-s1) > tol):

        iters += 1
        if iters > maxiters:
            raise Exception('Bisection method failed: Maximum iterations reached')

        sn = (s0+s1)/2
        fn = slope(sn)
        if np.sign(fn)==np.sign(f0):
            s0 = sn
            f0 = fn
        else:
            s1 = sn
            f1 = fn

    xyzi = satxyz + sn*lookxyz
    latlonalt = ecef_to_wgs84(xyzi)
    return latlonalt






def distance_to_shell(satlatlonalt, az, ze, shell_altitude, intersection='first', tol=1e-5):
    '''
    Along the line of sight, find the distance from the satellite to the
    first intersection with the shell at a given altitude. If no intersection
    exists, return np.nan. If the satellite is below this shell, return np.nan.
    
    INPUTS:
        satlatlonalt - array [latitude (deg), longitude (deg), altitude (km)] of the satellite
        az - direction of line of sight. degrees east of north.
        ze - direction of line of sight. degrees down from zenith
        shell_altitude - the altitude of the shell in km.
        
    OPTIONAL INPUTS:
        intersection - str, 'first' or 'second'. In general there are two intersections of the 
                       ray with the shell. This argument specifies which will be returned.
        tol - km, stopping criterion for the solver. Stop when the solution is not changing
              by more than  tol.
        
    OUTPUTS:
        d - distance in km. If no intersection exists, d is nan.
    
    HISTORY:
        10-Dec-2014: Written by Brian J. Harding (bhardin2@illinois.edu)
        31-Mar-2015: Included "intersection" parameter.
    '''
    maxiters = 1e4 # Just in case, so it doesn't hang forever

    # Convert to radians
    zer = ze*np.pi/180.
    azr = az*np.pi/180.

    # Create unit vector describing the look direction in Vertical-East-North (VEN) coordinates
    lookven = np.array([np.cos(zer), np.sin(zer)*np.sin(azr), np.sin(zer)*np.cos(azr)])
    # Convert satellite location to ecef
    satxyz = wgs84_to_ecef(satlatlonalt)
    # Convert look direction to ecef. This is a unit vector.
    lookxyz = ven_to_ecef(satlatlonalt, lookven)
    
    # First find the tangent altitude
    tanglatlonalt = tangent_point(satlatlonalt, az, ze)
    tangent_altitude = tanglatlonalt[2]
    
    # If there's no intersection, return nan.
    if shell_altitude <= tangent_altitude or shell_altitude > satlatlonalt[2]:
        return np.nan

    # Find the step size which results in an altitude equal to the shell altitude
    def my_func(step_size): # We want this function value to be zero
        xyzi = satxyz + step_size*lookxyz
        latlonalt = ecef_to_wgs84(xyzi)
        return latlonalt[2] - shell_altitude 
    
    # We need different initializations for the 'first' and 'second' intersections.
    if intersection=='first':
        # Find two values which straddle the solution.
        # (0) Satellite location (distance of zero)
        # (1) Tangent location (we need to calculate this distance)
        # (0)
        s0 = 0.
        f0 = my_func(s0)
        # (1)
        tangxyz = wgs84_to_ecef(tanglatlonalt)
        s1 = np.linalg.norm(tangxyz - satxyz) # distance from satellite to tangent point
        f1 = my_func(s1)
    elif intersection=='second':
        # Find two values which straddle the solution.
        # (0) Tangent location (we need to calculate this distance)
        # (1) A distance past the tangent location by a large amount
        # (0)
        tangxyz = wgs84_to_ecef(tanglatlonalt)
        s0 = np.linalg.norm(tangxyz - satxyz) # distance from satellite to tangent point
        f0 = my_func(s0)    
        # (1) Twice again as far should do it, but we'll do three times to be safe
        s1 = 3*s0
        f1 = my_func(s1)
    else:
        raise Exception('Unrecognized argument: intersection="%s". Try "first" or "second"' % intersection)
    
    # Check if straddle points are indeed straddle points
    if np.sign(f0)==np.sign(f1):
        raise Exception('Something went horribly wrong')

    # Straddle points found. Use bisection to converge to the answer.
    iters = 0
    while(abs(s0-s1) > tol):
        iters += 1
        if iters > maxiters:
            raise Exception('Bisection method failed: Maximum iterations reached')

        sn = (s0+s1)/2
        fn = my_func(sn)
        if np.sign(fn)==np.sign(f0):
            s0 = sn
            f0 = fn
        else:
            s1 = sn
            f1 = fn
    return sn




def distance_to_tangent_point(satlatlonalt, az, ze):
    '''
    Along the line of sight, find the distance from the satellite to the
    tangent point.
    
    INPUTS:
        satlatlonalt - array [latitude (deg), longitude (deg), altitude (km)] of the satellite
        az - direction of line of sight. degrees east of north.
        ze - direction of line of sight. degrees down from zenith
        
    OUTPUTS:
        d - distance in km.
    
    HISTORY:
        07-Jan-2014: Written by Brian J. Harding (bhardin2@illinois.edu)
    '''
    tanglatlonalt = tangent_point(satlatlonalt, az, ze)
    tangxyz = wgs84_to_ecef(tanglatlonalt)
    satxyz  = wgs84_to_ecef(satlatlonalt)
    d = np.linalg.norm(tangxyz - satxyz)
    return d



def azze_to_lla(satlatlonalt, az, ze, ht, tol=1e-6):
    '''
    Find the location (lat, lon, alt) where the requested look direction (az, ze) from the requested location
    intersects the specified altitude
    
    INPUTS:
        satlatlonalt - array [latitude (deg), longitude (deg), altitude (km)] of the satellite
        az - direction of line of sight. degrees east of north.
        ze - direction of line of sight. degrees down from zenith.
        ht - altitude to calculate the intersection point for.
        
    OPTIONAL INPUTS:
        tol - km, stopping criterion for the solver. Stop when the solution is not moving
              by more than a horizontal distance of tol.
        
    OUTPUTS:
        latlonalt - array [latitude (deg), longitude (deg), altitude (km)] of the tangent location
    
    HISTORY:
        22-Jan-2015: Written by Jonathan J. Makela (jmakela@illinois.edu) based on ICON.tangent_point
    '''
    maxiters = 1e4 # Just in case, so it doesn't hang forever

    # Convert to radians
    zer = ze*np.pi/180.
    azr = az*np.pi/180.

    # Create unit vector describing the look direction in Vertical-East-North (VEN) coordinates
    lookven = np.array([np.cos(zer), np.sin(zer)*np.sin(azr), np.sin(zer)*np.cos(azr)])
    # Convert satellite location to ecef
    satxyz = wgs84_to_ecef(satlatlonalt)
    # Convert look direction to ecef. This is a unit vector.
    lookxyz = ven_to_ecef(satlatlonalt, lookven)

    # Find the step size which minimizes the altitude. This problem is convex,
    # so it's easy.
    def altitude(step_size):
        xyzi = satxyz + step_size*lookxyz
        latlonalt = ecef_to_wgs84(xyzi)
        return latlonalt[2]
    
    # Or, find the step size for which the slope is zero.
    # Define function to calculate the slope numerically
    def slope(step_size):
        numerical_step = 1e-4 # km
        alt1 = altitude(step_size - numerical_step/2)
        alt2 = altitude(step_size + numerical_step/2)
        return (alt2-alt1)/numerical_step

    # Find two values which straddle the solution.
    s0 = 0. # start at the satellite
    s1 = 1. # iterate on the next value until we've straddled the solution
    f0 = slope(s0)
    f1 = slope(s1)
        
    if (f0 > 0) & (satlatlonalt[2] > ht):
        raise Exception('Ray path will not intersect desired altitude')
        
    if (f0 < 0) & (satlatlonalt[2] < ht):
        raise Exception('Ray path will not intersect desired altitude')

    M = 2 # multiplier for line search to find straddle points
    iters = 0
    while(np.sign(f0)==np.sign(ht-altitude(s1))):
        iters += 1
        if iters > maxiters:
            raise Exception('Initial line search failed: Maximum iterations reached')
        s0 = s1
        s1 = s1 * M

    # Straddle points found. Use bisection to converge to the answer.
    iters = 0
    while(abs(s0-s1) > tol):

        iters += 1
        if iters > maxiters:
            raise Exception('Bisection method failed: Maximum iterations reached')

        # The next step
        sn = (s0+s1)/2
                
        # Figure out what side of the desired altitude this new step is    
        if(f0 < 0):
            if altitude(sn) < ht:
                s1 = sn
            else:
                s0 = sn
        else:
            if altitude(sn) < ht:
                s0 = sn
            else:
                s1 = sn

    xyzi = satxyz + sn*lookxyz
    latlonalt = ecef_to_wgs84(xyzi)
    return latlonalt
