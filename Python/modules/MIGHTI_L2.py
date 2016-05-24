# A module for functions used for the conversion of MIGHTI Level 1 files to Level 2.1 and 2.2 files

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import ICON
import bisect # for quicker interpolation
from scipy import integrate



def get_emission_constants():
    '''
    The physical emission constants needed for the MIGHTI L2.1 analysis.
    INPUTS:
        none
    OUTPUTS:
        emissions --TYPE:dict, dictionary of red and green emission
                         parameters. (see below)
    '''

    emissions = {'green': {'mass': 16., # atomic mass of emitting specie [amu]
                           'lam': 557.7e-9, # center wavelength of emission [m]
                           }, 
                 'red':   {'mass': 16.,# atomic mass of emitting specie [amu]
                           'lam': 630.0e-9, # center wavelength of emission [m]
                           }
                 }

    return emissions




def get_instrument_constants(emission_color):
    '''
    The instrument constants needed for the MIGHTI L2.1 analysis.
    INPUTS:
        emission_color --TYPE:str, 'green' or 'red'
    OUTPUTS:
        instr_params   --TYPE:dict, dictionary of instrument parameters. (see below)
    '''
    emission = get_emission_constants()[emission_color]
    
    start_path = 4.62e-2 # optical path difference at start of interferometer [m]
    end_path   = 5.50e-2 # optical path difference at end of interferometer [m]
    
    # Set phase offset and zero phase depending on the color.
    # Zero phase:    The phase measurement which corresponds to a wind of zero.
    #                This needs to be found empirically.  ### JJM in rad?
    # Phase offset:  This constant is added to all phases so that there is no chance of a pi/-pi crossover.
    #                In practice, we'll have to see how the interferometer turns out before deciding on this.
    #                Changing this will require an equal change in zero_phase (already implemented below). ### JJM in rad?
    # TODO: How will this be done in practice, and what's the best way to hand off from NRL to Illinois?
    if emission_color == 'red':
        phase_offset = 2.309731
        zero_phase = 128.648466754 # updated 2015-11-24 using my sim
                                                       # (MIGHTI_Zero_wind_issues.ipynb)
    elif emission_color == 'green':
        phase_offset = 1.920212
        zero_phase = 54.66333601873041 # updated 2015-11-24 using my sim
                                                      # (MIGHTI_Zero_wind_issues.ipynb)
    else:
        raise Exception('emission_color = %s not understood' % emission_color)

    # Calculate phase-to-wind conversion factor
    # dphi = 2*pi*OPD*sigma*v/c (Eqn 1 of Englert et al 2007 Appl Opt.)
    # Use constant phase-to-wind factor for now, but we may eventually need to use
    # column-dependent values in some cases, such as satellite velocity removal.   
    c         = 299792458.0 # m/s, speed of light
    sigma     = 1.0/emission['lam']
    meanopd   = (start_path + end_path) / 2     ### JJM 2->2.?
    phase_to_wind_factor = c / (2*np.pi*sigma) * 1/meanopd   ### JJM 2->2., 1->1.?  This appears to be v/dphi, yes?  Why is meanopd treated outside of the second term?

    instr_params = { 'start_path': start_path,
                       'end_path': end_path,
                        'Nignore': 20, # The number of columns at the beginning and end of the interferogram
                                       # to ignore due to phase distortion from the filtering.
                   'phase_offset': phase_offset,
                     'zero_phase': zero_phase,
                        'min_amp': 0.0, # TODO. Is this the best way to implement this?
           'phase_to_wind_factor': phase_to_wind_factor,
                    
                 }

    return instr_params




def unwrap(x):
    '''
    Unwrap a monotonically increasing phase signal to remove -2*pi jumps.
    INPUTS:
        x     -- TYPE:array, UNITS:rad. Signal that has -2*pi jumps to remove
    OUTPUTS:
        xnew  -- TYPE:array, UNITS:rad. Copy of x with -2*pi jumps removed
    '''
    dx = np.diff(x)
    xnew = np.zeros(len(x))
    xnew[0] = x[0]
    idx = dx < 0
    dx[idx] = dx[idx] + 2*np.pi   ### JJM 2->2.?
    xnew[1:] = xnew[0] + np.cumsum(dx)
    return xnew




def circular_mean(angle0,angle1):
    '''
    Find the mean angle, taking into account -180/180 crossover. For example,
    circular_mean(10,50) is 30, but circular_mean(-170,160) is 175.
    INPUTS:
        angle0  -- TYPE:float or array, UNITS:deg. An angle in degrees.
        angle1  -- TYPE:float or array, UNITS:deg. An angle in degrees.
    OUTPUTS:
        angle   -- TYPE:float or array, UNITS:deg. The circular mean of the two
                   input angles.
    '''
    return 180/np.pi*np.angle((np.exp(1j*angle0*np.pi/180.) + np.exp(1j*angle1*np.pi/180.))/2) ### JJM 180->180., 2->2., use np.deg2rad instead of manually converting?




def tang_alt_to_ze(tang_alt, sat_alt, RE):
    '''
    Return the zenith angle(s) of the look direction(s) with given tangent 
    altitude(s). Uses spherical Earth approximation.
    INPUTS:
        tang_alt -- TYPE:array or float, UNITS:km. Tangent altitude(s) of the ray(s)
        sat_alt  -- TYPE:float,          UNITS:km. Satellite altitude (sat_alt > tang_alt)
        RE       -- TYPE:float,          UNITS:km. Effective radius of Earth
    OUTPUT:
        ze       -- TYPE:array or float, UNITS:deg. Zenith angle(s) of the ray(s)
    '''
    if hasattr(tang_alt,"__len__"):
        tang_alt = np.array(tang_alt)
        if any(sat_alt <= tang_alt):
            raise Exception('Tangent altitude must be below satellite altitude')
    elif sat_alt <= tang_alt:
        raise Exception('Tangent altitude must be below satellite altitude')
        
    ze = 180 - 180/np.pi*np.arcsin( (tang_alt+RE)/(sat_alt+RE) )  ### JJM 180->180.?
    return ze




def ze_to_tang_alt(ze, sat_alt, RE):
    '''
    Return the tangent altitude(s) of the look direction(s) with given zenith 
    angle(s). Uses spherical Earth approximation.
    INPUTS:
        ze       -- TYPE:array or float, UNITS:deg. Zenith angle(s) of the ray(s)
        sat_alt  -- TYPE:float,          UNITS:km.  Satellite altitude
        RE       -- TYPE:float,          UNITS:km.  Effective radius of Earth
    OUTPUT:
        tang_alt -- TYPE:array or float, UNITS:km.  Tangent altitude(s) of the ray(s)
    '''
    if hasattr(ze,"__len__"):
        ze = np.array(ze)
        if any( ze < 90 ) or any( ze > 180 ):
            raise Exception('Angle must be between 90 and 180, exclusive.')
    elif ( ze < 90 ) or ( ze > 180):
        raise Exception('Angle must be between 90 and 180, exclusive.')
    tang_alt = (sat_alt+RE)*np.sin(ze*np.pi/180.) - RE    ### JJM use np.deg2rad?
    return tang_alt




def remove_satellite_velocity(I, sat_velocity, sat_velocity_vector, mighti_ecef_vectors, phase_to_wind_factor,
                              show_plot = False):
    '''
    Modify the interferogram to remove the effect of satellite velocity upon the phase. 
    INPUTS:
        I                   -- TYPE:array(ny,nx), UNITS:arb.  The MIGHTI interferogram. 
        sat_velocity        -- TYPE:float,        UNITS:m/s.  ICON velocity.
        sat_velocity_vector -- TYPE:array(3,),    UNITS:none. The unit ECEF vector describing
                               the direction of ICON's velocity vector.
        mighti_vectors      -- TYPE:array(ny,3),  UNITS:none. Each row is a unit 3-vector in
                               ECEF coordinates defining the look direction of each measurement.
        phase_to_wind_factor-- TYPE:float,        UNITS:m/s/rad. The factor to multiply by a phase
                               change to get a wind change.
        show_plot           -- TYPE:bool,         UNITS:none. If True, display a plot
    OUTPUTS:
        I                   -- TYPE:array(ny,nx), UNITS:arb.  The MIGHTI interferogram, corrected
                               for the effects of satellite motion on the phase.
    '''
    ny,nx = np.shape(I)
    I2 = I.copy() # make a copy so that the input isn't overwritten
    
    # Create a column vector with projected satellite velocity. 
    # Remember, we are ignoring horizontal extent for now.  ### Is this in plans for improvements?  Should it be a TODO?
    proj_sat_vel = np.zeros(ny)
    for i in range(ny):
        look_ecef = mighti_ecef_vectors[i,:] # look direction of this pixel in ECEF
        proj_sat_vel[i] = sat_velocity * np.dot(sat_velocity_vector, look_ecef)

    ### JJM The comment below is not terribly clear to me.  The v/dphi factor seems to be passed in, so giving its equation here confused me.
    # Transform velocity to phase 
    # dphi = 2*pi*OPD*sigma*v/c (Eqn 1 of Englert et al 2007 Appl Opt.)
    # Use constant phase-to-wind factor for now, but eventually may need to use
    # actual OPD if we include horizontal extent of
    # the velocity correction. Is average satellite velocity enough?
    icon_vel_phase = proj_sat_vel/phase_to_wind_factor

    # Subtract phase from the interferogram (Change this when horizontal extent is included)  ### JJM What would this change look like?
    corr = np.exp(-1j*icon_vel_phase)
    for jj in range(nx):
        I2[:,jj] = I2[:,jj]*corr
        
    # Show a plot of satellite velocity, if requested.
    if show_plot:
        plt.plot(proj_sat_vel, range(ny), 'k')
        plt.xlabel('Projected Satellite Velocity [m/s]')
        plt.ylabel('Row Number')
        
    return I2




def find_phase_discontinuity(I, min_amp, Nignore):
    '''
    Raise a warning if a phase discontinuity is expected. A phase discontinuity 
    could cause large non-physical jumps in the reported wind.
    INPUTS:
        I          -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, MIGHTI interferogram.
        min_amp    -- TYPE:float,        UNITS:arb.  Rows of the interferogram with 
                      an amplitude less than this will not be analyzed.
        Nignore    -- TYPE:int,          UNITS:pixel. The number of columns at the
                      beginning and end of the interferogram to ignore due to phase
                      distortion from the filtering.
    OUTPUTS:
        None
    TODO:
        For now, warnings will be printed to stdout. In the future we may want to change that.

    '''
    
    ny = np.shape(I)[0]
    init_phase = np.zeros((ny))
    for j in range(ny):
        irow = I[j,:]
        phase = np.angle(irow)
        phaseu = unwrap(phase[Nignore:-Nignore])
        init_phase[j] = phaseu[0]
        ampl = abs(irow)
        if max(ampl) < min_amp:
            init_phase[j] = nan
    if any(abs(init_phase) > 0.8*np.pi):
        print 'WARNING: phase near pi/-pi crossover. Consider changing phase_offset'
    if any(np.diff(init_phase) > 0.8*2*np.pi):
        print 'WARNING: phase pi/-pi crossover detected. Recommend changing phase_offset'






def create_observation_matrix(tang_alt, icon_alt, top_layer='thin'):
    '''
    Define the matrix D whose inversion is known as "onion-peeling." The forward model is:
        I = D * Ip
    where I is the measured interferogram, D is the observation matrix (alternatively known
    as the distance matrix or path length matrix), and Ip is the onion-peeled interferogram.
    The matrix is created by approximating the line of sight integration as a summation using
    the trapezoidal rule, and collecting the weights for each node.
    
    INPUTS:
        tang_alt   -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
        icon_alt   -- TYPE:float,        UNITS:km.   Altitude of the satellite.
    OPTIONAL INPUTS:
        top_layer  -- TYPE:str,                      'thin': model top layer as a thin shell of constant emission (default)
                                                     'exp':  model top layer as an exponential falloff in altitude
    OUTPUTS:
        D          -- TYPE:array(ny,ny), UNITS:km.   Observation matrix. Also called the "path matrix"
                                                     or "distance matrix"
    '''
    
    H = 26. # km, assumed scale height of VER falloff with altitude (used if top_layer=='exp')
            # This was found by fitting many profiles for which there was significant
            # emission above 300 km. Profiles were generated from Zhang/Shepherd model and
            # from photochemical model fed by IRI/MSIS. (See MIGHTI SSR paper for details on
            # airglow models).
    
    def q(x,rm,r): ### JJM What does this function do?
        return 0.5*x*np.sqrt(rm**2 + x**2) + 0.5*rm**2 * np.log(2*(np.sqrt(rm**2 + x**2)+x)) - r*x ### JJM 2->2. in log?
    
    M = len(tang_alt)   ### JJM What is M?

    RE = 6371. # km, assume the Earth is locally spherical with an effective radius RE.
               # (The estimated winds are barely sensitive to the choice of RE. This
               #  approximation introduces an error < 1mm/s)
    
    D = np.zeros((M,M))
    for m in range(M):
        rm   = RE + tang_alt[m]
        # Loop over regions
        for k in range(m,M-1):
            # Region k is between nodes (i.e., tangent altitudes) k and k+1
            rk   = RE + tang_alt[k]
            rkp1 = RE + tang_alt[k+1]
            # Note that to use formulas from notes, I have to swap rk and rkp1  ### JJM What notes?
            # Contribution to node k
            wkkp1 = 2/(rk-rkp1) * ( q(np.sqrt(rk**2  -rm**2),rm,rk)   - q(np.sqrt(rkp1**2-rm**2),rm,rk  ) )   ### JJM 2->2.?
            wkk   = 2/(rk-rkp1) * ( q(np.sqrt(rkp1**2-rm**2),rm,rkp1) - q(np.sqrt(rk**2  -rm**2),rm,rkp1)  )  ### JJM 2->2.?

            D[m,k] += wkk
            D[m,k+1] += wkkp1
            
        # Handle contributions from above 300km differently, depending on top_layer='thin' or 'exp':
        if top_layer == 'thin': # Use assumption that airglow goes to zero just above top altitude
            # Calculate contribution to top node from above top tangent altitude
            # (Assuming thin shell: zero above a certain altitude)
            rk   = RE + tang_alt[M-1]
            rkp1 = RE + tang_alt[M-1] + (tang_alt[M-1]-tang_alt[M-2])
            wkk = 2/(rk-rkp1) * ( q(np.sqrt(rkp1**2-rm**2),rm,rkp1) - q(np.sqrt(rk**2  -rm**2),rm,rkp1)  )  ### JJM 2->2.?
            D[m,M-1] += wkk            
            
        elif top_layer == 'exp': # Use exponential falloff model
            rt = tang_alt[m] + RE 
            r0 = tang_alt[-1] + RE
            
            def func(x, rt):  ### JJM Describe function
                return np.exp(-1/H*(np.sqrt(x**2 + rt**2) - r0))   ### JJM 1->1.?
            
            x0 = np.sqrt(r0**2- rt**2)
            D[m,M-1] += 2*integrate.quad(func, x0, np.inf, args=(rt))[0]   ### JJM 2->2.?
    
    return D






def create_local_projection_matrix(tang_alt, icon_alt):
    '''
    Define the matrix B whose entries give the factor by which a horizontal wind
    would be projected onto the line of sight. This has the same shape as the
    observation matrix (i.e., distance matrix). At the tangent point, this factor
    is 1.0. Far from the tangent point, this factor is smaller. If this effect is 
    accounted for, it makes a small change in the winds (less than 5 m/s).
    INPUTS:
        tang_alt   -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
        icon_alt   -- TYPE:float,        UNITS:km.   Altitude of the satellite.
    OUTPUTS:
        B          -- TYPE:array(ny,ny), UNITS:km.   Local projection matrix. B[i,j] = cos(angle 
                                                     between ray i and the tangent of shell j
                                                     at the point where they intersect)
    '''
    
    # Assume the Earth is locally spherical with an effective radius RE.
    # (The estimated winds are barely sensitive to the choice of RE. This
    #  approximation introduces an error < 1mm/s)
    RE = 6371.
    theta = tang_alt_to_ze(tang_alt, icon_alt, RE)
    
    ny = len(tang_alt)
    
    # Calculate local-horizontal projection factors
    B = np.nan*np.zeros((ny,ny)) # matrix to hold cosine correction factors
    for i in range(ny):
        for j in range(i,ny): # only calculate upper triangular part
            th = theta[i]
            r = tang_alt[j]
            B[i,j] = (RE+icon_alt)/(RE+r) * np.sin(th*np.pi/180)
    return B
     

    
    
    
def extract_phase_from_row(row, zero_phase, phase_offset, Nignore):
    '''
    Given a 1-D interference pattern (i.e., a row of the intererogram), 
    analyze it to get a single phase value, which represents the wind.
    INPUTS:
        row          -- TYPE:array(nx), UNITS:arb.   A row of the complex-valued, MIGHTI interferogram.
        zero_phase   -- TYPE:float,     UNITS:rad.   The phase angle which is equivalent 
                                                     to a wind value of zero.
        phase_offset -- TYPE:float,     UNITS:rad.   An offset to avoid 2pi ambiguities.
        Nignore      -- TYPE:int,       UNITS:pixel. The number of columns at the
                       beginning and end of the interferogram to ignore due to phase
                       distortion from the filtering.
    OUTPUTS:
        phase        -- TYPE:float,     UNITS:rad. 
    '''
    
    row_phase = np.angle(row)

    # average phase and then take delta. Need unwrapping for this.
    phaseu = unwrap(row_phase[Nignore:-Nignore]-phase_offset) + phase_offset
    meanphase = np.mean(phaseu)
    phase = meanphase - zero_phase
    return phase

    
    
    
    
def perform_inversion(I, tang_alt, icon_alt, account_for_local_projection=False, zero_phase=None, phase_offset=None, Nignore=None, top_layer='thin'):
    '''
    Perform the onion-peeling inversion on the interferogram to return
    a new interferogram, whose rows refer to specific altitudes. In effect,
    this function undoes the integration along the line of sight. There is an
    option to use a simple inversion, or a complicated inversion that accounts for
    the fact that the line of sight is not always parallel to the ground.
    INPUTS:
        I           -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, MIGHTI interferogram.
        tang_alt    -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
        icon_alt    -- TYPE:float,        UNITS:km.   Altitude of the satellite.
    OPTIONAL INPUTS:
        account_for_local_projection -- TYPE:bool.   If False, a simple inversion is used. (default)
                                        If True, the inversion accounts for the fact that the ray is not 
                                        perfectly tangent to each shell at each point along the ray. If True,
                                        the following variables are needed:
        zero_phase  -- TYPE:float,      UNITS:rad.   The phase angle which is equivalent 
                                                     to a wind value of zero.
        phase_offset-- TYPE:float,      UNITS:rad.   An offset to avoid 2pi ambiguities.
        Nignore     -- TYPE:int,        UNITS:pixel. The number of columns at the
                                        beginning and end of the interferogram to ignore due to phase
                                        distortion from the filtering.
        top_layer   -- TYPE:str,        'thin': model top layer as a thin shell of constant emission (default)
                                        'exp':  model top layer as an exponential falloff in altitude
    OUTPUTS:
        Ip          -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, onion-peeled interferogram.
    '''
    
    ny,nx = np.shape(I)
    
    # Create the path matrix
    D = create_observation_matrix(tang_alt, icon_alt, top_layer)
    
    # The inversion will proceed in different ways depending on whether
    # we will try to account for the local horizontal projection.
    if not account_for_local_projection:
        
        # This is implemented with a simple linear inversion
        Ip = np.linalg.solve(D,I)
        
    else:
        
        B = create_local_projection_matrix(tang_alt, icon_alt)
        Ip = np.zeros((ny,nx), dtype=complex) # onion-peeled interferogram
        phase = np.zeros(ny) # phases at each altitude
        
        for i in range(ny)[::-1]: # onion-peel from the top altitude down
            dii = D[i,i] # path length
            Li = I[i,:] # we will peel off the other layers from this row
            # Loop over layers above this one
            for j in range(i+1,ny):
                dij = D[i,j]
                # Calculate the normalized jth row without the wind component
                Ij = Ip[j,:]*np.exp(-1j*phase[j])
                # Calculate it with the projected wind component
                Ij_proj = Ij*np.exp(1j*phase[j]*B[i,j])
                # Remove this contribution from the current layer
                Li = Li - dij*Ij_proj
            # final normalization by this layer's path length
            Li = Li/dii
            Ip[i,:] = Li
            # Analyze the layer to get the phase, and store it.
            # Note that the zero-phase/zero-wind is needed for this step.
            phase[i] = extract_phase_from_row(Li, zero_phase, phase_offset, Nignore)
            
    return Ip




def extract_wind(Ip, zero_phase, phase_offset, min_amp, Nignore, phase_to_wind_factor):
    '''
    Analyze the onion-peeled interferogram to extract wind and apparent intensity.
    For now, we include a placeholder for the uncertainty.
    INPUTS:
        Ip                  -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, onion-peeled 
                               interferogram.
        zero_phase          -- TYPE:float,        UNITS:rad.  The phase angle which is equivalent 
                               to a wind value of zero.
        phase_offset        -- TYPE:float,        UNITS:rad.   An offset to avoid 2pi ambiguities.
        min_amp             -- TYPE:float,        UNITS:arb.  Rows of the interferogram with 
                               an amplitude less than this will not be analyzed.
        Nignore             -- TYPE:int,          UNITS:pixel. The number of columns at the
                               beginning and end of the interferogram to ignore due to phase
                               distortion from the filtering.
        phase_to_wind_factor-- TYPE:float,        UNITS:m/s/rad. The factor to multiply by a phase
                               change to get a wind change.
    OUTPUTS:
        v                   -- TYPE:array(ny)     UNITS:m/s.  The estimated line of sight wind.
        ve                  -- TYPE:array(ny)     UNITS:m/s.  Estimated uncertainty in v.
        a                   -- TYPE:array(ny)     UNITS:arb.  The estimated apparent intensity. Not
                               to be confused with actual intensity, since temperature plays a role.
        ae                  -- TYPE:array(ny)     UNITS:arb.  Estimated uncertainty in a.
        p                   -- TYPE:array(ny)     UNITS:rad.  TEMPORARY. line of sight phase 
    '''
    ny,nx = np.shape(Ip)
    
    p = np.nan*np.zeros(ny) # analyzed phase of each shell
    a = np.nan*np.zeros(ny) # analyzed intensity of each shell
    for j in range(ny):
        irow = Ip[j,:]
        ampl  = abs(irow)
        if np.nanmax(ampl) > min_amp:
            dphase = extract_phase_from_row(irow, zero_phase, phase_offset, Nignore)
            p[j] = dphase
            a[j] = ampl[Nignore:-Nignore].mean()

    # Convert phase to velocity
    v = phase_to_wind_factor * p 
    
    ve = np.nan*np.zeros(np.shape(v))
    ae = np.nan*np.zeros(np.shape(a))

    return v, ve, a, ae, p




def fix_longitudes(lons, lon_target):
    '''
    Unwrap the list of longitudes to avoid 360-deg jumps. The list will
    be fixed so that it contains a value within 180 deg of lon_target and
    is otherwise continuous.
    INPUTS:
        lons       -- TYPE:array, UNITS:deg. An ordered list of longitudes to be unwrapped.
        lon_target -- TYPE:float, UNITS:deg. See above.
    OUTPUTS:
        lons_new   -- TYPE:array, UNITS:deg. An ordered list of longitudes with jumps removed.
    '''
    lons_new = np.array(lons).copy()
    
    # Find the index with value closest to lon_target (mod 360)
    diff_vec = np.mod(lons_new - lon_target + 180, 360) - 180  ### JJM int->float?
    k = np.argmin(abs(diff_vec))
    # Change the entire array up or down by 360 (or a multiple) if necessary, keying off of target_lon.
    n = round((lons_new[k] - lon_target)/360.)
    lons_new = lons_new - n*360 ### JJM 360->360.?
        
    # Define function to remove jumps
    def fix_jump(jump, val):
        n = round(jump/360.)
        return val - n*360 ### JJM 360->360.?
    # Traverse right, removing jumps > +/- 180
    for i in range(k+1,len(lons_new)):
        jump = lons_new[i] - lons_new[i-1]
        lons_new[i] = fix_jump(jump, lons_new[i])
    # Traverse left, removing jumps > +/- 180
    for i in range(k-1,-1,-1):
        jump = lons_new[i] - lons_new[i+1]
        lons_new[i] = fix_jump(jump, lons_new[i])   

    return lons_new




def attribute_measurement_location(tang_lat, tang_lon, tang_alt):
    '''
    Determine the geographical location to which the measurement will be attributed. The 
    current implementation of the inversion, which uses trapezoidal integration, means
    that we should simply return the tangent locations.
    
    NOTE: If the implementation of the following functions are changed, this function may need
    to change: create_observation_matrix, perform_inversion, extract_wind.
    INPUTS:
        tang_lat    -- TYPE:array(ny), UNITS:deg.   Tangent latitudes.
        tang_lon    -- TYPE:array(ny), UNITS:deg.   Tangent longitudes.
        tang_alt    -- TYPE:array(ny), UNITS:km.    Tangent altitudes.
    OUTPUTS:
        lat         -- TYPE:array(ny), UNITS:deg.   Measurement latitudes.
        lon         -- TYPE:array(ny), UNITS:deg.   Measurement longitudes.
        alt         -- TYPE:array(ny), UNITS:km.    Measurement altitudes.
    '''
    
    # The following two functions are saved for posterity, in case
    # we switch back to Riemann integration.
    def shift_up_by_half(vec):
        '''
        Shift the input vector up by half the resolution. Extrapolate for the top entry.
        '''
        bottom = vec
        top = bottom.copy()
        top[:-1] = top[1:]
        top[-1] = top[-1] + (top[-2]-bottom[-2])
        return (0.5*top + 0.5*bottom)
    
    def shift_up_by_half_angle(vec):
        '''
        Shift the input vector up by half the resolution. Extrapolate for the top entry.
        Use circular mean instead of arithmetic mean. This is intended for longitude
        calculations.
        '''
        # First, unwrap angles so there are no 360-deg jumps
        vec_new = fix_longitudes(vec, vec[0])   ### JJM Is this variable used anywhere?
        bottom = vec
        top = bottom.copy()
        top[:-1] = top[1:]
        top[-1] = top[-1] + (top[-2]-bottom[-2])
        mid = np.zeros(len(bottom))
        for i in range(len(mid)):
            mid[i] = circular_mean(top[i],bottom[i])
        # Un-unwrap angles, so that they are all in (-180,180)
        mid = np.mod(mid+180, 360) - 180  # JJM int->float?
        return mid
        
    lat = tang_lat
    lon = tang_lon
    alt = tang_alt
    
    return lat, lon, alt





def los_az_angle(sat_latlonalt, lat, lon, alt):
    '''
    Calculate the azimuth angle of the line of sight, evaluated at the 
    measurement location (lat, lon, alt). Assumes WGS84 Earth.
    INPUTS:
        sat_latlonalt -- TYPE:array(3),  UNITS:(deg,deg,km). Satellite location in WGS84.
        lat           -- TYPE:array(ny), UNITS:deg.          Measurement latitudes.
        lon           -- TYPE:array(ny), UNITS:deg.          Measurement longitudes.
        alt           -- TYPE:array(ny), UNITS:km.           Measurement altitudes.
    OUTPUTS:
        az            -- TYPE:array(ny), UNITS:deg.          Azimuth angle of line of sight
                         from the satellite to the measurement location, evaluated at the 
                         measurement location. Degrees East of North.
    '''
    ny = len(lat)
    local_az = np.zeros(ny)
    sat_xyz = ICON.wgs84_to_ecef(sat_latlonalt)
    for i in range(ny):
        meas_latlonalt = np.array([lat[i], lon[i], alt[i]]) # where the measurement is attributed to
        meas_xyz = ICON.wgs84_to_ecef(meas_latlonalt)
        look_xyz = meas_xyz - sat_xyz # look direction
        loc_az, loc_ze = ICON.ecef_to_azze(meas_latlonalt, look_xyz) # look direction in az, ze at measurement point.
        local_az[i] = loc_az  
    return local_az





def remove_Earth_rotation(v_inertial, az, lat, lon, alt):
    '''
    Transform wind measurement from inertial coordinates to a reference
    frame rotating with the Earth. This can be thought of as "removing 
    Earth rotation from the line-of-sight measurement."
    INPUTS:
        v_inertial    -- TYPE:array(ny), UNITS:m/s.   Line-of-sight velocity in inertial
                         coordinates.
        az            -- TYPE:array(ny), UNITS:deg.   Azimuth angle of line of sight
                         from the satellite to the measurement location, evaluated at the 
                         measurement location. Degrees East of North. See los_az_angle() above.
        lat           -- TYPE:array(ny), UNITS:deg.   Measurement latitudes.
        lon           -- TYPE:array(ny), UNITS:deg.   Measurement longitudes.
        alt           -- TYPE:array(ny), UNITS:km.    Measurement altitudes.
    OUTPUTS:
        v             -- TYPE:array(ny), UNITS:m/s.   Line-of-sight velocity in Earth-fixed
                         coordinates.
    '''
    ny = len(v_inertial)
    corot_contribution = np.zeros(ny)
    for i in range(ny):
        meas_latlonalt = np.array([lat[i], lon[i], alt[i]]) # where the measurement is attributed to
        meas_xyz = ICON.wgs84_to_ecef(meas_latlonalt)
        rho = np.sqrt(meas_xyz[0]**2 + meas_xyz[1]**2)
        sidereal_day_length = 23*60*60 + 56*60 + 4 # sidereal day is 23 hrs 56 min 4 sec ### JJM int->float?
        corot_vel = 2*np.pi*rho/sidereal_day_length*1e3 ### JJM 2->2.
        # Compute component along LoS
        corot_contribution[i] = corot_vel * np.sin(np.pi/180*az[i])
    v = v_inertial - corot_contribution
    return v






def interpolate_linear(x,y,x0,extrapolation='hold'):
    '''
    Linear interpolation of the function y = f(x) to the location x0.
    x and y are vectors comprising samples of this function. This function is
    5 times faster than scipy.interpolate.interp1d, and allows for
    zero-order-hold extrapolation. If you are interpolating to many points, 
    then scipy.interpolate.interp1d is likely faster.

    INPUTS:
        x    -- TYPE:array(n), UNITS:arb. Independent variable of samples of function.
        y    -- TYPE:array(n), UNITS:arb. Dependent variable of samples of function.
        x0   -- TYPE:float,    UNITS:arb. Independent variable of interpolation point.
    OPTIONAL INPUTS:
        extrapolation -- TYPE:str,        'hold': extrapolate by using values at end points
                                          'none': do not extrapolate. Points will be np.nan
    OUTPUTS:
        y0   -- TYPE:float,    UNITS:arb. Interpolated value.
    
    '''
    j0 = bisect.bisect(x,x0) - 1 # index to the left
    j1 = j0 + 1 # index to the right
    # Handle extrapolations
    if j0 == -1:
        if extrapolation=='hold':
            y0 = y[0]
        elif extrapolation == 'none':
            y0 = np.nan
        else: 
            raise Exception('"%s" not understood' % extrapolation)
    elif j1 == len(x):
        if extrapolation=='hold':
            y0 = y[-1]
        elif extrapolation == 'none':
            y0 = np.nan
        else: 
            raise Exception('"%s" not understood' % extrapolation)
    else: # linear interpolation
        y0 = (y[j1]-y[j0])/(x[j1]-x[j0])*(x0-x[j0]) + y[j0]    
    return y0




