
# TODO: Check this into svn.

# TODO: Incorporate some of these functions into pyglow.Point? Probably not - it'll be nice
# to have everything here.

import numpy as np


def calculate_slice(ze,az,satlat,satlon,satalt,ver,v,t,instr_params, step=1.):
    '''
    Step through the atmosphere modeled by the functions A, V, T and accumulate
    a MIGHTI interferogram along the line-of-sight specified. Lines-of-sight are
    calculated using a spherical earth.
    INPUTS:
        ze      - radians, zenith angle of look direction ( >pi/2 )
        az      - radians, azimuth angle of look direction (0 is north, 90 is east)
        satlat  - degrees, latitude of satellite
        satlon  - degrees, longitude of satellite
        satalt  - km,      altitude of satellite
        ver     - function of (lat,lon,alt) which specifies the volume emission rate
                  of the observed emission, in phot/cm^3/s
        v       - function of (lat,lon,alt) which specifies the neutral velocity
        t       - function of (lat,lon,alt) which specifies the neutral temperature
        params  - dictionary containing MIGHTI instrument parameters (see "interferogram")
        step    - km, resolution of the integration along the line of sight
    OUTPUTS:
        interferogram  - array of length instr_params['npixel']. 
    
    TODO: Eliminate the spherical-symmetry assumption. (i.e., actually use satlat, satlon, az)
    TODO: Account for nonzero FOV in horizontal direction. ... How? Maybe better to just parameterize this.
    '''
    RE = 6371. # km, radius of earth
    x = 0.          # This measures perpendicular to the satellite nadir
    y = RE + satalt # This measures parallel to satellite nadir
    dx = np.sin(ze) # unit step in x
    dy = np.cos(ze) # unit step in y
    igram = np.zeros(instr_params['npixelx'])
    
    while (np.sqrt(x**2 + y**2) <= RE+satalt): # until you get far enough away from earth that it doesn't matter
        # Step along line-of-sight and accumulate fringe patterns
        # calculate altitude
        # TODO: step in real 3-D space, calculating lat,lon as well
        alt = np.sqrt(x**2 + y**2) - RE
        lat = np.nan # for now, this doesn't matter
        lon = np.nan # for now, this doesn't matter

        # create interferogram at this location
        instr_params['V'] = v(lat, lon, alt)
        instr_params['T'] = t(lat, lon, alt)
        # ver is in phot/cm^3/s
        # multiply by step in cm and by 1e-6 to put units in Mphot/cm^2/s, i.e., Rayleighs
        instr_params['I'] = ver(lat, lon, alt) * step * 1.0e5 / 1.0e6

        # add the resulting interferogram to the accumulating interferogram
        igram += interferogram(instr_params)
        
        # take the next step
        x += step*dx
        y += step*dy

    return igram
        

    
def add_noise(I, instr_params):
    signal = I.copy() # in electrons

    # Add dark signal
    # (Adding shot noise on top of this will account for "dark noise")
    signal += instr_params['darkcurrent']*instr_params['exptime'] # still "electrons" 

    # Add shot noise (doesn't matter if you do this before/after binning)
    # This accounts for dark noise since dark current is in the signal.
    shotnoise = np.sqrt(signal)*np.random.randn(*np.shape(signal))
    signal += shotnoise

    # On-chip binning is already implicitly done.

    # Read noise
    readnoise = instr_params['readnoise']*np.random.randn(*np.shape(signal))
    signal += readnoise
    
    return signal



def interferogram(params):
    '''
    Function to simulate the noiseless, unquantized, interferogram produced 
    by the MIGHTI instrument.
    Creates a monochromatic interferogram (=1+cos) -- the envelope of the
    modulated part is multiplied by the apodization functions that are due 
    to the linewidth and the non-perfect instrument visibility.
    Note that the output of this function is in units of electrons,
    not CCD counts.

    INPUTS: a dictionary with the following entries
        frq - the spatial frequency to simulate [fringes per grating width]
        mass - the atomic mass of the line [amu]
        sigma - wavenumber of center wavelength [m^-1]
        T - temperature of the line [K]
        V - line of sight wind velocity [m/s]
        amplitude - interferogram amplitude at zero path for perfect visibility
        npixel - number of pixels in one fringe
        startpath - optical path difference on left edge of CCD [m]
        endpath - optical path difference on right edge of CCD [m]
        and others... (see below)
    OUTPUTS:
        interferogram - array of interference pattern, in units of electrons on the CCD

    HISTORY:
      Written by Jonathan J. Makela based on the routine "interferogram"
      provided by Christoph Englert in MIGHTIPerformanceModel_v4.1.pro
      17 Jan 2014: Updated by Brian Harding to match the specifications in
                   MIGHTIPerformanceModel_v11.pro provided by Kenneth Marr.
   
    '''

    # Physical constants
    c = 299792458.          # The speed of light [m/s]
    k = 1.3806503e-23      # Boltzmann's constant [J/K]
    Mamu = 1.66053886e-27  # atomic mass unit [kg]
    
    # unpack the parameter dictionary
    startpath = params['startpath']              
    endpath = params['endpath']
    npixelx = params['nx']
    npixely = params['ny']
    frq = params['frq']
    mass = params['mass']
    lam = params['lam'] # wavelength [m]
    T = params['T'] # Neutral Temperature [K]
    V = params['V'] # Neutral velocity [m/s]
    I = params['I'] # Intensity of emission, in Rayleighs.
                          # I think this should be "Intensity," though
                          # I've never been clear about this.
                          # e.g., VER*path put into units of Mphot/cm2/sec
    gain = params['gain']
    aperture_area = params['aperture_area'] # m^2
    coneangle1 = params['coneangle1']
    coneangle2 = params['coneangle2']
    opteff     = params['opteff']
    exptime    = params['exptime']
    
    
    # Instrument visibility (assumption: independent of wavenumber and fringe
    # frequency!  Real instrument will have some frequency dependence)
    # Question: Where do these numbers (0.8, 0.05) come from?
    InstrumentVisibility = 0.8-0.05*(np.arange(0.,npixelx))/(npixelx-1)

    # Calculate the thermal width of the line [cm^-1]
    sigma = 1./lam
    sigmaD = sigma*np.sqrt(k*T/mass/Mamu/c**2)
    
    path = startpath + (np.arange(0.,npixelx))/(npixelx-1)*(endpath-startpath)
    phase = (np.arange(0.,npixelx))/(npixelx-1) * frq * 2.* np.pi
    offs = np.mod((startpath/(endpath-startpath) * frq*2.*np.pi),2.*np.pi)
    phase = phase + offs
    phase = phase + 2.*np.pi*path*sigma*V/c        # add Doppler shift
    envelope = np.exp(-2.*np.pi**2*sigmaD**2.*path**2) # thermal envelope
    
    # Construct the interferogram with thermal envelope and instrument
    # visibility. No meaningful units.
    igram = ((np.cos(phase) * envelope) * InstrumentVisibility + 1.)
    
    # Account for etendue, emission brightness, and all other
    # instrumental effects that turn Rayleighs into electrons on the CCD.
    # Assume equi-angle lens projection in vertical direction.
    # (TODO: make sure this is right)
    solidangle = np.sin(coneangle1)*np.sin(coneangle2)/npixely
    photons = I*1.e10                         # phot/m^2/sec into 4*pi steradians
    photons = photons * solidangle/(4*np.pi)  # phot/m^2/sec into detector
    photons = photons * aperture_area         # phot/sec into detector
    photons = photons * exptime               # total photons entering MIGHTI
    electrons = photons * opteff              # total electrons on this row of CCD
    # split up electrons into pixels in this slice of the interferogram
    ccdslice = electrons * igram / np.sum(igram)
    
    return ccdslice


def get_emission_constants():

    emissions = {'green': {'frq': 40.2, # fringe frequency across array [cycles/pixel]
                           'mass': 16., # atomic mass of emitting specie [amu]
                           'lam': 557.7e-9, # center wavelength of emission [m]
                           }, 
                 'red':   {'frq': 85.7, # fringe frequency across array [total cycles]
                           'mass': 16.,# atomic mass of emitting specie [amu]
                           'lam': 630.0e-9, # center wavelength of emission [m]
                           }
                 }

    return emissions


def get_instrument_constants():

    instrument = {    #'nx': 512,  # number of superpixels per slice of the interferogram (horizontal)
                      'nx': 512,  # number of superpixels per slice of the interferogram (horizontal)
                      'ny': 1456/16,  # number of superpixels in vertical (altitude) direction on CCD.
                    'startpath': 4.52e-2, # optical path difference at start of interferometer [m]
                      'endpath': 5.50e-2, # optical path difference at end of interferometer [m]
                'aperture_area': 3.078e-2*3.078e-2, # cm^2
                   'coneangle1': 3.22, # fov of instrument in horizontal direction, [deg]
                   'coneangle2': 5.74, # fov of instrument in vertical direction, [deg]
                        'gain' : 1, # gain of CCD (counts per electron)
                       'opteff': 0.16, # combined transmittance of all optics.
                  'darkcurrent': 0.1*64, # electrons/sec/pixel (TODO: per superpixel or pixel?)
                    'readnoise': 7, # electrons rms, per 2.5 km superpixel
                  }

    return instrument


def tanht2angle(h, H, RE=6371.):
    '''
    Return the zenith angle of the look direction with tangent altitude h.
    INPUTS:
        h - tangent altitude, km
        H - satellite altitude, km (H > h)
        RE - radius of Earth, km
    OUTPUT:
        theta - zenith angle, rad
    This should work on scalar or vector input, which is why it's messy
    '''
    if hasattr(h,"__len__"):
        h = np.array(h)
        if any(H <= h):
            raise Exception('Tangent altitude must be below satellite altitude')
    elif H <= h:
        raise Exception('Tangent altitude must be below satellite altitude')
        
    theta = np.pi - np.arcsin( (h+RE)/(H+RE) )
    return theta


def angle2tanht(theta, H, RE=6371.):
    '''
    Return the tangent altitude of the observation with zenith angle theta.
    INPUTS:
        theta - zenith angle, rad ( pi/2 < theta < pi )
        H - satellite altitude, km
        RE - radius of Earth, km
    OUTPUT:
        h - tangent altitude, km
    This should work on scalar or vector input, which is why it's messy
    '''
    if hasattr(theta,"__len__"):
        theta = np.array(theta)
        if any( theta < np.pi/2 ) or any( theta > np.pi ):
            raise Exception('Angle must be between pi/2 and pi, exclusive.')
    elif ( theta < np.pi/2 ) or ( theta > np.pi ):
        raise Exception('Angle must be between pi/2 and pi, exclusive.')
    h = (H+RE)*np.sin(theta) - RE
    return h


def get_solar_zenith_angle(pt):
    '''
    Calculate the angle from zenith to the sun at the time and
    location of the pyglow.Point pt.
    INPUT:
        pt - pyglow.Point
    OUTPUT:
        sza - Solar zenith angle (radians)
    '''
    import ephem
    sun = ephem.Sun()
    obs = ephem.Observer()
    obs.lon = str(pt.lon)
    obs.lat = str(pt.lat)
    obs.date = pt.dn.strftime('%Y/%m/%d %H:%M:%S')
    obs.pressure = 0. # ignore refraction. This makes a negligible difference.
    obs.elevation = 1000*pt.alt # This makes a negligible difference.
    sun.compute(obs)
    sza = np.pi/2 - float(sun.alt) # radians
    return sza


def get_redline_airglow(pt):
    '''
    Return the red-line airglow for the pyglow.Point pt in photons/cm^3/s.
    Daytime climatological model from Zhang and Shepherd 2005 (Proc of SPIE).
    Nighttime model from Link and Cogger 1989 using MSIS concentrations.
    TODO: smoothly interpolate between twilight and nighttime somehow.
    Valid below 60 degrees latitude.
    '''
    
    # Calculate SZA to determine if it is daytime, twilight, or nighttime.
    sza = get_solar_zenith_angle(pt)
    szadeg = sza * 180/np.pi

    if szadeg < 87: ################## daytime
        # Model from Zhang and Shepherd 2005 doi:10.1117/12.627150
        cossza = np.cos(sza)
        F107 = pt.f107
        Vp = (1.791*F107 + 76.4)*cossza**(1/np.e) + 35
        I = (24.8*F107 + 608)*cossza**(1/np.e) + 390
        hp = (0.0006*F107 - 0.232)*Vp + (0.180*F107 + 237)
        W = 8.03*cossza + 0.079*F107 + 34
        V_6300 = Vp*np.exp(-(pt.alt - hp)**2/(2*W**2))
        
    elif szadeg < 104.5: ############# twilight
        # Model from Zhang and Shepherd 2005 doi:10.1117/12.627150
        cossza = np.cos(sza)
        F107 = pt.f107
        Vp = (4.703*F107 + 437)*(cossza + 0.25)**1.8 + (0.095*F107 + 9.4)
        I = (84.2*F107 + 2497)*(cossza + 0.25)**1.8 + (1.174*F107 + 91)
        hp = (0.0006*F107 - 0.232)*Vp + (0.180*F107 + 237)
        W = 8.03*cossza + 0.079*F107 + 34
        V_6300 = Vp*np.exp(-(pt.alt - hp)**2/(2*W**2))
        
    else: ############################ nighttime
        # let's see if IRI and MSIS have been executed
        # if not, run the appropriate models:
        if np.isnan(pt.ne):
            pt.run_iri()
        if np.isnan(pt.nn['O2']):
            pt.run_msis()

        Ne = pt.ne;       # electron density [cm^-3]
        Tn = pt.Tn_msis;  # neutral temperature [K]
        Ti = pt.Ti;       # ion temperature [K]
        Te = pt.Te;       # electron temperature [K]
        O2 = pt.nn['O2']; # O2 density [cm^-3]
        N2 = pt.nn['N2']; # N2 density [cm^-3]

        te = Te/300;
        ti = Ti/300;

        # These coefs are from Link and Cogger, JGR 93(A9), 988309892, 1988
        K1_6300 = 3.23e-12*np.exp(3.72/ti - 1.87/ti**2);
        K2_6300 = 2.78e-13*np.exp(2.07/ti - 0.61/ti**2);
        K3_6300 = 2.0e-11*np.exp(111.8/Tn);
        K4_6300 = 2.9e-11*np.exp(67.5/Tn);
        K5_6300 = 1.6e-12*Te**0.91;
        b6300 = 1.1;
        a1D = 7.45e-3;   # Corrected value form Link and Cogger, JGR, 94(A2), 1989
        a6300 = 5.63e-3; # Corrected value form Link and Cogger, JGR, 94(A2), 1989

        # Calculate O+ assuming mixture of ions (also from Link and Cogger, 1988)
        a1 = 1.95e-7*te**-0.7;
        a2 = 4.00e-7*te**-0.9;
        Oplus = Ne/(1.+K2_6300*N2/a2/Ne + K1_6300*O2/a1/Ne);

        AGNumerator = a6300/a1D*b6300*K1_6300*Oplus*O2;
        AGDenominator = 1.+(K3_6300*N2+K4_6300*O2+K5_6300*Ne)/a1D;
        V_6300 = AGNumerator / AGDenominator;
        
    return V_6300



def get_greenline_airglow(pt):
    '''
    Return the green-line airglow for the pyglow.Point pt in photons/cm^3/s.
    Daytime climatological model from Zhang and Shepherd 2005 (JGR).
    Nighttime model from Vargas et. al 2007 using MSIS concentrations.
    TODO: smoothly interpolate between twilight and nighttime somehow.
    Valid below 60 degrees latitude.
    '''
    
    # Calculate SZA to determine if it is daytime, twilight, or nighttime.
    sza = get_solar_zenith_angle(pt)
    szadeg = sza * 180/np.pi

    if szadeg < 80: ################## daytime
        # Model from Zhang and Shepherd 2005 doi:10.1029/2004JA010887.
        # Two peaks, in E and F layer
        cossza = np.cos(sza)
        F107 = pt.f107
        # F
        Vf = (3.56*F107 + 290)*cossza**(1.5) + 20
        If = (28.19*F107 + 1375)*cossza**(1.5) + 244
        hf = (-0.023*F107-17.6)*np.log(Vf) + (0.34*F107 + 230)
        Wf = -0.0147*Vf + (0.037*F107+26)
        bf = (pt.alt - hf)/Wf
        V_5577F = Vf*np.exp(1 - bf - np.exp(-bf))
        # E
        Ve = (4.5*F107 + 236)*cossza**1.2 + 136
        Ie = (8.09*F107 + 503)*cossza**1.2 + (1.22*F107 + 158)
        he = (0.023*F107 + 95)
        We = 8 # "Varies between 6 and 10 km" apparently in a random way
        be = (pt.alt - he)/We
        V_5577E = Ve*np.exp(1 - be - np.exp(-be))
        
        V_5577 = V_5577E + V_5577F
        
        
    
    else: ############################ nighttime
        # let's see if IRI and MSIS have been executed
        # if not, run the appropriate models:
        if np.isnan(pt.ne):
            pt.run_iri()
        if np.isnan(pt.nn['O2']):
            pt.run_msis()

        Ne = pt.ne;       # electron density [cm^-3]
        Tn = pt.Tn_msis;  # neutral temperature [K]
        Ti = pt.Ti;       # ion temperature [K]
        Te = pt.Te;       # electron temperature [K]
        O = pt.nn['O']    # O density [cm^-3]
        O2 = pt.nn['O2']; # O2 density [cm^-3]
        N2 = pt.nn['N2']; # N2 density [cm^-3]

        te = Te/300;
        ti = Ti/300;
        
        # Coefs and model from Vargas et al 2007 doi:10.1029/2006JD007642,
        # but are ultimately from Mcdade et al 1986 and Murtagh et al 1990.
        
        k1 = 4.7e-33*(300/Tn)**2
        k5 = 4.0e-12*np.exp(-865/Tn)
        A5 = 1.18
        A6 = 1.35
        
        V_5577 = (A5*k1*O**3*(O2 + N2)) / ((A6 + k5*O2)*(15*O2 + 211*O))
        
        
    return V_5577

    
    
    
    
    
    
    
