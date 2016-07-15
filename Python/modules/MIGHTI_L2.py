# A module for functions used for the conversion of MIGHTI Level 1 files to Level 2.1 and 2.2 files

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import ICON
import bisect
from scipy import integrate
from datetime import datetime, timedelta
import netCDF4




def get_instrument_constants(emission_color, start_path, end_path):
    '''
    The instrument constants needed for the MIGHTI L2.1 analysis.
    INPUTS:
        emission_color --TYPE:str, 'green' or 'red'
        start_path     --TYPE:float, UNITS:m. optical path difference at left edge of interferogram
        end_path       --TYPE:float, UNITS:m. optical path difference at right edge of interferogram
    OUTPUTS:
        instr_params   --TYPE:dict, dictionary of instrument parameters. (see below)
    TODO:
        - What's the best way to store, maintain, and update these values, which may
          be time-dependent? (e.g., zero phase)
    '''
    
    
    # Set phase offset and zero phase depending on the color.
    # Zero phase:    radians. The phase measurement which corresponds to a wind of zero.
    #                In practice, this needs to be found empirically, and will change with time. 
    #                Most likely, this value will be drawn from a calibration database.
    #                For analyzing simulated data, we can determine exactly what it should be.
    # Phase offset:  radians. This constant is added to all phases so that there is no chance 
    #                of a pi/-pi crossover. In practice, we'll have to see how the interferometer 
    #                turns out before deciding on this. (TODO)

    if emission_color == 'red':
        lam          = 630.0304e-9   # center wavelength of emission [m] (Osterbrock et al. 1996)
        phase_offset = -1.139683     # updated 2016-06-29 using my instrument model
        zero_phase   = 128.648464318 # updated 2016-05-26 using my instrument model
                                                       # (MIGHTI_Zero_wind_issues.ipynb)
        bin_size     = 13            # How many rows to bin together in post-processing
        
    elif emission_color == 'green':
        lam          = 557.7338e-9   # center wavelength of emission [m]
        phase_offset = -1.630175     # updated 2016-06-29 using my instrument model
        zero_phase = 54.6633360579   # updated 2016-05-26 using my instrument model
                                                      # (MIGHTI_Zero_wind_issues.ipynb)
        bin_size     = 2             # How many rows to bin together in post-processing
    else: 
        raise Exception('emission_color = %s not understood' % emission_color)

    # Calculate phase-to-wind conversion factor
    # dphi = 2*pi*OPD*sigma*v/c (Eq 1 of Englert et al 2007 Appl Opt., 
    #                           and Eq 2 of Harding et al. 2016 SSR)
    # Use constant phase-to-wind factor for now, but we may eventually need to use
    # column-dependent values in some cases, such as satellite velocity removal.   
    c         = 299792458.0 # m/s, speed of light
    sigma     = 1.0/lam
    meanopd   = (start_path + end_path) / 2.
    phase_to_wind_factor = c / (2.*np.pi*sigma*meanopd)

    instr_params = {    'Nignore': 20, # The number of columns at the beginning and end of the interferogram
                                       # to ignore due to phase distortion from the filtering.
                   'phase_offset': phase_offset,
                     'zero_phase': zero_phase,
                        'min_amp': 0.0, # TODO. What's the best way to implement this? For now, do nothing.
           'phase_to_wind_factor': phase_to_wind_factor,
                       'bin_size': bin_size,
                 }

    return instr_params




def unwrap(x):
    '''
    Unwrap a monotonically increasing phase signal to remove -2*pi jumps.
    This is very similar to np.unwrap, but only unwraps negative jumps. 
    INPUTS:
        x     -- TYPE:array, UNITS:rad. Signal that has -2*pi jumps to remove
    OUTPUTS:
        xnew  -- TYPE:array, UNITS:rad. Copy of x with -2*pi jumps removed
    '''
    dx = np.diff(x)
    xnew = np.zeros(len(x))
    xnew[0] = x[0]
    idx = dx < 0
    dx[idx] = dx[idx] + 2.*np.pi
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
    return np.rad2deg(np.angle((np.exp(1j*np.deg2rad(angle0)) + np.exp(1j*np.deg2rad(angle1)))/2.))




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
        
    ze = 180. - np.rad2deg(np.arcsin( (tang_alt+RE)/(sat_alt+RE) ))
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
    tang_alt = (sat_alt+RE)*np.sin(np.deg2rad(ze)) - RE  
    return tang_alt




def remove_satellite_velocity(I, sat_velocity, sat_velocity_vector, mighti_ecef_vectors, phase_to_wind_factor):
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
    OUTPUTS:
        I                   -- TYPE:array(ny,nx), UNITS:arb.  The MIGHTI interferogram, corrected
                               for the effects of satellite motion on the phase.
    '''
    # TODO: Perhaps, remove on a column-by-column basis to account for the varying 
    # velocity projection (and maybe also the varying phase_to_wind_factor) with x.
    
    ny,nx = np.shape(I)
    I2 = I.copy() # make a copy so that the input isn't overwritten
    
    # Create a column vector with projected satellite velocity. 
    # Remember, we are ignoring horizontal extent for now. 
    proj_sat_vel = np.zeros(ny)
    for i in range(ny):
        look_ecef = mighti_ecef_vectors[i,:] # look direction of this pixel in ECEF
        proj_sat_vel[i] = sat_velocity * np.dot(sat_velocity_vector, look_ecef)


    # Convert the projected satellite velocity to a phase
    icon_vel_phase = proj_sat_vel/phase_to_wind_factor

    # Subtract phase from the interferogram
    corr = np.exp(-1j*icon_vel_phase)
    for jj in range(nx):
        I2[:,jj] = I2[:,jj]*corr
        
    return I2




def bin_array(b, y, lon = False):
    '''
    Downsample y by binning it, improving statistics. Every b
    elements of y will be averaged together to create a new array, y_b, 
    of length ny_b = ceil(len(y)/b). Binning starts at the end of the array, 
    so the first element of y_b may not represent exactly b samples of y.
    INPUTS:
        b    -- TYPE:int,          The number of rows to bin together
        y    -- TYPE:array(ny),    The array to be binned
    OPTIONAL INPUTS:
        lon  -- TYPE:bool,         If True, 360-deg discontinuities will
                                   be removed before averaging (e.g., for
                                   longitude binning).
    OUTPUTS:
        y_b  -- TYPE:array(ny_b),  The binned array
    '''
    ny = len(y)
    ny_b = int(np.ceil(1.0*ny/b))
    y_b = np.zeros(ny_b, dtype=y.dtype)
    for i in range(0,ny_b): # bin from the end to the beginning.
        i_new   = ny_b-i-1
        i_start = ny-(i+1)*b
        i_stop  = ny-i*b
        
        # grab the samples to be binned
        if np.mod(ny,b)!=0 and i_new==0: # special case in case ny is not divisible by b
            y_samps = y[:i_stop]
        else: # grab 
            y_samps = y[i_start:i_stop]

        if lon:
            y_samps = fix_longitudes(y_samps, 180.)
        y_b[i_new] = np.mean(y_samps)
        
    return y_b
    
    
    
    
def bin_image(b, I):
    '''
    Downsample the interferogram in altitude to improve statistics while
    degrading vertical resolution. Every b rows will be averaged together. 
    Binning starts at high altitudes, so the lower rows of I_b may not represent 
    exactly b rows of I.
    INPUTS:
        b           -- TYPE:int,                        The number of rows to bin together
        I           -- TYPE:array(ny,nx),   UNITS:arb.  The MIGHTI interferogram
    OUTPUTS:
        I_b         -- TYPE:array(ny_b,nx), UNITS:arb.  The binned MIGHTI interferogram
    '''
    
    ny,nx = np.shape(I)
    # Initial call to bin_array to see what the size of the new image will be
    tmp = bin_array(b, I[:,0])
    ny_b = len(tmp)
    
    # Bin the interfogram column by column
    I_b = np.zeros((ny_b,nx),dtype=I.dtype)
    for i in range(nx):
        I_b[:,i] = bin_array(b,I[:,i])
    return I_b


def create_observation_matrix(tang_alt, icon_alt, top_layer='exp', integration_order=0):
    '''
    Define the matrix D whose inversion is known as "onion-peeling." The forward model is:
        I = D * Ip
    where I is the measured interferogram, D is the observation matrix, and Ip is the 
    onion-peeled interferogram. If integration_order is 1, the observation matrix is 
    created by assuming the spectrum (and thus the interferogram) is a piecewise linear 
    function of altitude, treating the values of the interferogram at the tangent locations
    as the unknowns, and writing the measurements as a linear function of the unknowns.
    If integration_order is 0, the same recipe is followed, except the spectrum is 
    assumed to be a piecewise constant function of altitude, and the unknowns are the 
    values of the interferogram at the midpoint between two tangent altitudes.
    
    INPUTS:
        tang_alt   -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
        icon_alt   -- TYPE:float,        UNITS:km.   Altitude of the satellite.
    OPTIONAL INPUTS:
        top_layer  -- TYPE:str,          'thin': assume VER goes to zero above top layer
                                         'exp':  assume VER falls off exponentially in altitude (default)
        integration_order -- TYPE:int,   0: Use Riemann-sum rule for discretizing line-of-sight integral (default)
                                         1: Use trapezoidal rule for discretizing line-of-sight integral
    OUTPUTS:
        D          -- TYPE:array(ny,ny), UNITS:km.   Observation matrix. Also called the "path matrix"
                                                     or "distance matrix"
    '''
    
    
    H = 26. # km, assumed scale height of VER falloff with altitude (used if top_layer=='exp')
            # This was found by fitting many profiles for which there was significant
            # emission above 300 km. Profiles were generated from Zhang/Shepherd model and
            # from photochemical model fed by IRI/MSIS. (See MIGHTI SSR paper for details on
            # airglow models).
    
    def q(x,rm,r): 
        # antiderivative of (sqrt(x**2 + rm**2) - r)   w.r.t. x
        return 0.5*x*np.sqrt(rm**2 + x**2) + 0.5*rm**2 * np.log(2.*(np.sqrt(rm**2 + x**2)+x)) - r*x
    
    M = len(tang_alt)   # Number of rows of interferogram

    RE = 6371. # km, assume the Earth is locally spherical with an effective radius RE.
               # (The estimated winds are barely sensitive to the choice of RE. This
               #  approximation introduces an error < 1mm/s)
               
    D = np.zeros((M,M))
    
    #################### Zero-order integration #######################
    # Assume airglow is constant within thin altitude shells. This is
    # analogous to Riemann sum integration
    if integration_order == 0:
    
        theta = np.deg2rad(tang_alt_to_ze(tang_alt, icon_alt, RE))
        
        # Define grid. Bottom of each layer is defined by tangent height of observation.
        rbottom = tang_alt
        # Define top of each layer.
        rtop = rbottom.copy()
        rtop[:-1] = rbottom[1:]
        rtop[-1] = rbottom[-1] + (rtop[1]-rbottom[1])
        # Define midpt of each layer
        rmid = (rbottom + rtop)/2

        # Build observation matrix
        for m in range(M):
            for k in range(M):
                th = theta[m]
                rb = rbottom[k]
                rt = rtop[k]
                sb = np.cos(th)**2 - 1 + ((RE+rb)/(RE+icon_alt))**2
                st = np.cos(th)**2 - 1 + ((RE+rt)/(RE+icon_alt))**2
                if sb < 0: # there is no intersection of LOS with altitude rb. Set term to 0.
                    # Note: this might be due to numerical rounding for tangent altitude. 
                    # Do the same thing either way.
                    sb = 0.
                if st < 0: # there is no intersection of LOS with altitude rt. Set term to 0.
                    st = 0.
                D[m,k] = 2*(RE+icon_alt) * ( np.sqrt(st) - np.sqrt(sb) )
           
            if top_layer == 'exp': # Use exponential falloff model
                rt = tang_alt[m] + RE
                r0 = tang_alt[-1] + RE
                
                def func(x, rt):
                    # The extrapolation function to be numerically integrated. (Eq 6 in Harding et al. 2016 SSR)
                    return np.exp(-1./H*(np.sqrt(x**2 + rt**2) - r0))
                
                x0 = np.sqrt(r0**2- rt**2)
                D[m,M-1] = 2.*integrate.quad(func, x0, np.inf, args=(rt))[0]
                
    #################### First-order integration #######################
    # Assume airglow varies linearly within thin altitude shells. This is
    # analogous to trapezoidal rule integration
    elif integration_order == 1:
        for m in range(M):
            rm   = RE + tang_alt[m]
            # Loop over regions
            for k in range(m,M-1):
                # Region k is between nodes (i.e., tangent altitudes) k and k+1
                rk   = RE + tang_alt[k]
                rkp1 = RE + tang_alt[k+1]
                # Compile the contribution from this region to the nodes below and above, using the
                # analytical evaluation of the Abel integral.
                wkkp1 = 2./(rk-rkp1) * ( q(np.sqrt(rk**2  -rm**2),rm,rk)   - q(np.sqrt(rkp1**2-rm**2),rm,rk  ) )
                wkk   = 2./(rk-rkp1) * ( q(np.sqrt(rkp1**2-rm**2),rm,rkp1) - q(np.sqrt(rk**2  -rm**2),rm,rkp1)  )

                D[m,k] += wkk
                D[m,k+1] += wkkp1
                
            # Handle contributions from above 300km differently, depending on top_layer='thin' or 'exp':
            if top_layer == 'thin': # Use assumption that airglow goes to zero just above top altitude
                # Calculate contribution to top node from above top tangent altitude
                rk   = RE + tang_alt[M-1]
                rkp1 = RE + tang_alt[M-1] + (tang_alt[M-1]-tang_alt[M-2])
                wkk = 2./(rk-rkp1) * ( q(np.sqrt(rkp1**2-rm**2),rm,rkp1) - q(np.sqrt(rk**2  -rm**2),rm,rkp1)  )
                D[m,M-1] += wkk
                
            elif top_layer == 'exp': # Use exponential falloff model
                rt = tang_alt[m] + RE
                r0 = tang_alt[-1] + RE
                
                def func(x, rt):
                    # The extrapolation function to be numerically integrated. (Eq 6 in Harding et al. 2016 SSR)
                    return np.exp(-1./H*(np.sqrt(x**2 + rt**2) - r0))
                
                x0 = np.sqrt(r0**2- rt**2)
                D[m,M-1] += 2.*integrate.quad(func, x0, np.inf, args=(rt))[0]
                
    else:
        raise Exception('"integration_order == %i" not supported. Use 0 or 1.' % integration_order)
    
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
            B[i,j] = (RE+icon_alt)/(RE+r) * np.sin(np.deg2rad(th))
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

    # Average phase and then take delta. Need unwrapping for this.
    # First, apply phase offset (keeping phase between -pi and pi)
    row_phase_off = np.mod(row_phase[Nignore:-Nignore] - phase_offset + np.pi, 2*np.pi) - np.pi
    # Second, unwrap, and undo phase offset
    phaseu = unwrap(row_phase_off) + phase_offset
    meanphase = np.mean(phaseu)
    phase = meanphase - zero_phase
    return phase

    
    
    
    
def perform_inversion(I, tang_alt, icon_alt, top_layer='exp', integration_order=0, 
                      account_for_local_projection=True, zero_phase=None, 
                      phase_offset=None, Nignore=None):
    '''
    Perform the onion-peeling inversion on the interferogram to return
    a new interferogram, whose rows refer to specific altitudes. In effect,
    this function undoes the integration along the line of sight.
    INPUTS:
        I           -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, MIGHTI interferogram.
        tang_alt    -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
        icon_alt    -- TYPE:float,        UNITS:km.   Altitude of the satellite.
    OPTIONAL INPUTS:
        top_layer   -- TYPE:str,        'thin': assume VER goes to zero above top layer
                                        'exp':  assume VER falls off exponentially in altitude (default)
        integration_order -- TYPE:int,   0: Use Riemann-sum rule for discretizing line-of-sight integral (default)
                                         1: Use trapezoidal rule for discretizing line-of-sight integral
        account_for_local_projection -- TYPE:bool.   If False, a simple inversion is used.
                                        If True, the inversion accounts for the fact that the ray is not 
                                        perfectly tangent to each shell at each point along the ray. 
                                        (default True)
                                        If True, the following variables are needed:
        zero_phase  -- TYPE:float,      UNITS:rad.   The phase angle which is equivalent 
                                                     to a wind value of zero.
        phase_offset-- TYPE:float,      UNITS:rad.   An offset to avoid 2pi ambiguities.
        Nignore     -- TYPE:int,        UNITS:pixel. The number of columns at the
                                        beginning and end of the interferogram to ignore due to phase
                                        distortion from the filtering.
    OUTPUTS:
        Ip          -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, onion-peeled interferogram.
    '''
    
    ny,nx = np.shape(I)
    
    # Create the path matrix
    D = create_observation_matrix(tang_alt, icon_alt, top_layer=top_layer, integration_order=integration_order)
    
    # The inversion will proceed in different ways depending on whether
    # we will try to account for the local horizontal projection.
    if not account_for_local_projection:
        
        # This is implemented with a simple linear inversion
        Ip = np.linalg.solve(D,I)
        
    else:
        # The problem becomes nonlinear, but still solvable in closed form.
        # This code implements Eq (9) in the MIGHTI L2 Space Science Reviews
        # paper (Harding et al. 2016).
        
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
        a                   -- TYPE:array(ny)     UNITS:arb.  The estimated apparent VER. Not
                               to be confused with actual VER, since temperature plays a role.
        ae                  -- TYPE:array(ny)     UNITS:arb.  Estimated uncertainty in a.
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
            a[j] = ampl[Nignore:-Nignore].sum()

    # Convert phase to velocity
    v = phase_to_wind_factor * p 
    
    # TODO: propagate uncertainties
    ve = np.nan*np.zeros(np.shape(v))
    ae = np.nan*np.zeros(np.shape(a))

    return v, ve, a, ae





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
    diff_vec = np.mod(lons_new - lon_target + 180., 360.) - 180. 
    k = np.argmin(abs(diff_vec))
    # Change the entire array up or down by 360 (or a multiple) if necessary, keying off of target_lon.
    n = round((lons_new[k] - lon_target)/360.)
    lons_new = lons_new - n*360.
        
    # Define function to remove jumps
    def fix_jump(jump, val):
        n = round(jump/360.)
        return val - n*360. 
    # Traverse right, removing jumps > +/- 180
    for i in range(k+1,len(lons_new)):
        jump = lons_new[i] - lons_new[i-1]
        lons_new[i] = fix_jump(jump, lons_new[i])
    # Traverse left, removing jumps > +/- 180
    for i in range(k-1,-1,-1):
        jump = lons_new[i] - lons_new[i+1]
        lons_new[i] = fix_jump(jump, lons_new[i])   

    return lons_new





def attribute_measurement_location(tang_lat, tang_lon, tang_alt, integration_order=0):
    '''
    Determine the geographical location to which the measurement will be attributed. The 
    current implementation of the inversion, which uses trapezoidal integration, means
    that we should simply return the tangent locations.
    
    NOTE: If the implementation of any of the following functions are changed, this 
    function may need to change: create_observation_matrix, perform_inversion, extract_wind.
    INPUTS:
        tang_lat    -- TYPE:array(ny), UNITS:deg.   Tangent latitudes.
        tang_lon    -- TYPE:array(ny), UNITS:deg.   Tangent longitudes.
        tang_alt    -- TYPE:array(ny), UNITS:km.    Tangent altitudes.
    OPTIONAL INPUTS:
        integration_order -- TYPE:int,   0: Use Riemann-sum rule for discretizing line-of-sight integral (default)
                                         1: Use trapezoidal rule for discretizing line-of-sight integral
    OUTPUTS:
        lat         -- TYPE:array(ny), UNITS:deg.   Measurement latitudes.
        lon         -- TYPE:array(ny), UNITS:deg.   Measurement longitudes.
        alt         -- TYPE:array(ny), UNITS:km.    Measurement altitudes.
    '''
    
    def shift_up_by_half(vec):
        """
        Shift the input vector up by half the resolution. Extrapolate for the top entry.
        """
        bottom = vec
        top = bottom.copy()
        top[:-1] = top[1:]
        top[-1] = top[-1] + (top[-2] - bottom[-2])
        return 0.5 * top + 0.5 * bottom

    def shift_up_by_half_angle(vec):
        """
        Shift the input vector up by half the resolution. Extrapolate for the top entry.
        Use circular mean instead of arithmetic mean. This is intended for longitude
        calculations.
        """
        vec_new = fix_longitudes(vec, vec[0])
        bottom = vec
        top = bottom.copy()
        top[:-1] = top[1:]
        top[-1] = top[-1] + (top[-2] - bottom[-2])
        mid = np.zeros(len(bottom))
        for i in range(len(mid)):
            mid[i] = circular_mean(top[i], bottom[i])

        mid = np.mod(mid + 180, 360) - 180
        return mid
    
    if integration_order == 1:
        lat = tang_lat
        lon = tang_lon
        alt = tang_alt
    else:
        lat = shift_up_by_half(tang_lat)
        lon = shift_up_by_half_angle(tang_lon)
        alt = shift_up_by_half(tang_alt)
        
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
        sidereal_day_length = 23.*60.*60. + 56.*60. + 4. # sidereal day is 23 hrs 56 min 4 sec 
        corot_vel = 2.*np.pi*rho/sidereal_day_length*1e3
        # Compute component along LoS
        corot_contribution[i] = corot_vel * np.sin(np.deg2rad(az[i]))
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
        extrapolation -- TYPE:str,        'hold': extrapolate by using values at end points (default)
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





def level1_to_dict(L1_fn, emission_color):
    '''
    Read a level 1 file and translate it into a dictionary that the 
    level 2.1 processing can use.
    
    INPUTS:
        L1_fn   -- TYPE:str.  The full path and filename of the level 1 file
        emission_color --TYPE:str, 'green' or 'red'
        
    OUTPUTS:
        L1_dict -- TYPE:dict. A dictionary containing information needed for
                              the level 2.1 processing. The keys are:
                                     L1_fn
                                     I_amp
                                     I_phase
                                     I_amp_uncertainty
                                     I_phase_uncertainty
                                     tang_alt_start
                                     tang_alt_end
                                     tang_lat_start
                                     tang_lat_end
                                     tang_lon_start
                                     tang_lon_end
                                     emission_color
                                     icon_alt_start
                                     icon_alt_end
                                     icon_lat_start
                                     icon_lat_end
                                     icon_lon_start
                                     icon_lon_end
                                     mighti_ecef_vectors_start
                                     mighti_ecef_vectors_end
                                     icon_ecef_ram_vector_start
                                     icon_ecef_ram_vector_end
                                     icon_velocity_start
                                     icon_velocity_end
                                     source_files
                                     time
                                     exp_time
                                     interferometer_start_path
                                     interferometer_end_path
    TODO:
        - lots of stuff, mostly just waiting on finalization of L1 file.
        - interferometer start path and end path
        - time in better format
        - exposure time different than end minus start?
    '''
    
    f = netCDF4.Dataset(L1_fn)
    L1_dict = {}
    L1_dict['L1_fn']                       = L1_fn
    L1_dict['I_amp']                       = f['%s_ENVELOPE_NOISY' % emission_color.upper()][:]
    L1_dict['I_phase']                     = f['%s_PHASE_NOISY' % emission_color.upper()][:]
    L1_dict['I_amp_uncertainty']           = f['%s_ENVELOPE_UNCERTAINTY' % emission_color.upper()][:]
    L1_dict['I_phase_uncertainty']         = f['%s_PHASE_UNCERTAINTY' % emission_color.upper()][:]
    L1_dict['tang_alt_start']              = f['TANGENT_ALTITUDES_START'][:]
    L1_dict['tang_alt_end']                = f['TANGENT_ALTITUDES_END'][:]
    L1_dict['tang_lat_start']              = f['TANGENT_LATITUDES_START'][:]
    L1_dict['tang_lat_end']                = f['TANGENT_LATITUDES_END'][:]
    L1_dict['tang_lon_start']              = f['TANGENT_LONGITUDES_START'][:]
    L1_dict['tang_lon_end']                = f['TANGENT_LONGITUDES_END'][:]
    L1_dict['emission_color']              = emission_color
    L1_dict['icon_alt_start']              = f['ICON_ALTITUDE_START'][:].item()
    L1_dict['icon_alt_end']                = f['ICON_ALTITUDE_END'][:].item()
    L1_dict['icon_lat_start']              = f['ICON_LATITUDE_START'][:].item()
    L1_dict['icon_lat_end']                = f['ICON_LATITUDE_END'][:].item()
    L1_dict['icon_lon_start']              = f['ICON_LONGITUDE_START'][:].item()
    L1_dict['icon_lon_end']                = f['ICON_LONGITUDE_END'][:].item()
    L1_dict['mighti_ecef_vectors_start']   = f['MIGHTI_ECEF_VECTORS_START'][:]
    L1_dict['mighti_ecef_vectors_end']     = f['MIGHTI_ECEF_VECTORS_END'][:]
    L1_dict['icon_ecef_ram_vector_start']  = f['RAM_ECEF_VECTOR_START'][:]
    L1_dict['icon_ecef_ram_vector_end']    = f['RAM_ECEF_VECTOR_END'][:]
    L1_dict['icon_velocity_start']         = f['ICON_VELOCITY_START'][:].item()
    L1_dict['icon_velocity_end']           = f['ICON_VELOCITY_END'][:].item()
    L1_dict['source_files']                = [] # TODO
    tsec                                   = f['IMAGE_TIME_START'][:].item()
    L1_dict['time']                        = datetime(2015,1,1) + timedelta(seconds=tsec)
    L1_dict['exp_time']                    = f['IMAGE_TIME_END'][:].item() - f['IMAGE_TIME_START'][:].item()
    L1_dict['interferometer_start_path']   = 4.62e-2
    L1_dict['interferometer_end_path']     = 5.50e-2
    
    
    f.close()
    
    return L1_dict





def level1_uiuc_to_dict(L1_uiuc_fn):
    '''
    Read a level 1 file in special "UIUC" format, and translate it into a 
    dictionary that the level 2.1 processing can use. This function is 
    intended to be temporary, and used for testing only, not in operations.
    The L1_uiuc file is a .npz file.
    
    INPUTS:
        L1_uiuc_fn   -- TYPE:str.  The full path and filename of the level 1 file
    OUTPUTS:
        L1_dict      -- TYPE:dict. A dictionary containing information needed for
                                   the level 2.1 processing. The keys are:
                                     L1_fn
                                     I_amp
                                     I_phase
                                     I_amp_uncertainty
                                     I_phase_uncertainty
                                     tang_alt_start
                                     tang_alt_end
                                     tang_lat_start
                                     tang_lat_end
                                     tang_lon_start
                                     tang_lon_end
                                     emission_color
                                     icon_alt_start
                                     icon_alt_end
                                     icon_lat_start
                                     icon_lat_end
                                     icon_lon_start
                                     icon_lon_end
                                     mighti_ecef_vectors_start
                                     mighti_ecef_vectors_end
                                     icon_ecef_ram_vector_start
                                     icon_ecef_ram_vector_end
                                     icon_velocity_start
                                     icon_velocity_end
                                     source_files
                                     time
                                     exp_time
                                     interferometer_start_path
                                     interferometer_end_path
    '''
    
    L1_dict = {}
    L1_dict['L1_fn'] = L1_uiuc_fn
    
    npzfile = np.load(L1_uiuc_fn)
    
    ####### Hack for backwards compatibility #######
    # If interferogram is given as real/imaginary, then change
    # it to amp/phase
    if 'Ir' in npzfile.keys():
        Ir = npzfile['Ir']
        Ii = npzfile['Ii']
        I = Ir + 1j*Ii
        L1_dict['I_amp'] = abs(I)
        L1_dict['I_phase'] = np.angle(I)
    
    for key in npzfile.keys():
        if key in ['emission_color','icon_alt_start','icon_alt_end',
                   'icon_lat_start','icon_lat_end','icon_lon_start','icon_lon_end',
                   'icon_velocity_start','icon_velocity_end','time','exp_time',
                   'interferometer_start_path','interferometer_end_path',]:
            L1_dict[key] = npzfile[key].item()
        else:
            L1_dict[key] = npzfile[key]
    npzfile.close()
    
    return L1_dict
    





def level1_dict_to_level21(L1_dict, L21_fn, zero_phase_addition = 0.0, top_layer = 'exp', 
                           integration_order = 0, account_for_local_projection = True, 
                           bin_size = None, bin_method = 'before' ):
    '''
    High-level function to convert a level 1 dictionary (which was generated from
    a level 1 file) into a level 2.1 file.
    
    INPUTS:
        L1_dict             -- TYPE:dict.  A dictionary containing variables needed for
                                           the level 2.1 processing:
                                             L1_fn
                                             I_amp
                                             I_phase
                                             I_amp_uncertainty
                                             I_phase_uncertainty
                                             tang_alt_start
                                             tang_alt_end
                                             tang_lat_start
                                             tang_lat_end
                                             tang_lon_start
                                             tang_lon_end
                                             emission_color
                                             icon_alt_start
                                             icon_alt_end
                                             icon_lat_start
                                             icon_lat_end
                                             icon_lon_start
                                             icon_lon_end
                                             mighti_ecef_vectors_start
                                             mighti_ecef_vectors_end
                                             icon_ecef_ram_vector_start
                                             icon_ecef_ram_vector_end
                                             icon_velocity_start
                                             icon_velocity_end
                                             source_files
                                             time
                                             exp_time
                                             interferometer_start_path
                                             interferometer_end_path  
                                             
        L21_fn              -- TYPE:str.  The path and filename where the level 2.1
                                          file will be saved.                             
    OPTIONAL INPUTS:
        zero_phase_addition -- TYPE:float, UNITS:rad. A value which we will be added
                                                      to the zero phase used in the 
                                                      inversion.
        top_layer           -- TYPE:str, 'thin': assume VER goes to zero above top layer
                                         'exp':  assume VER falls off exponentially in altitude (default)
        integration_order   -- TYPE:int, 0: Use Riemann-sum rule for discretizing line-of-sight integral (default)
                                         1: Use trapezoidal rule for discretizing line-of-sight integral
        account_for_local_projection -- TYPE:bool.   If False, a simple inversion is used.
                                        If True, the inversion accounts for the fact that the ray is not 
                                        perfectly tangent to each shell at each point along the ray. 
                                        (default True)
        bin_size            -- TYPE:int or None,  The number of rows of the interferogram to bin together to 
                                                  improve statistics at the cost of altitude resolution. If 
                                                  None, use the default (color-dependent) value specified
                                                  in get_instrument_constants().
        bin_method          -- TYPE:str,  'before': perform binning on the L1 interferogram (default)
                                          'after' : perform binning on the onion-peeled interferogram
                                          Using 'before' has better precision but worse bias. It is recommended
                                          to use 'before' because if better bias is desired (at the expense of
                                          precision), then the recommended course is to not bin at all.
                                        
                                        
    OUTPUTS:
        flag                -- TYPE:int,              Equals 0 on success.
            
    TODO:
        - save NetCDF file.
        - error handling (return 0 or error code)
    
    '''
    if bin_method not in ['before','after']:
        raise Exception('bin_method="%s" not understood. Use "before" or "after".' % bin_method)

    ###  Load parameters from input dictionary
    Iraw = L1_dict['I_amp']*np.exp(1j*L1_dict['I_phase'])
    emission_color = L1_dict['emission_color']
    source_files = L1_dict['source_files']
    time = L1_dict['time']
    exp_time = L1_dict['exp_time']
    L1_fn = L1_dict['L1_fn']
    # Load parameters which are averaged from start to end of exposure.
    icon_alt = (L1_dict['icon_alt_start'] + L1_dict['icon_alt_end'])/2
    icon_lat = (L1_dict['icon_lat_start'] + L1_dict['icon_lat_end'])/2
    icon_lon = circular_mean(L1_dict['icon_lon_start'], L1_dict['icon_lon_end'])
    mighti_ecef_vectors = (L1_dict['mighti_ecef_vectors_start'] + L1_dict['mighti_ecef_vectors_end'])/2
    tang_alt = (L1_dict['tang_alt_start'] + L1_dict['tang_alt_end'])/2
    tang_lat = (L1_dict['tang_lat_start'] + L1_dict['tang_lat_end'])/2
    tang_lon = circular_mean(L1_dict['tang_lon_start'], L1_dict['tang_lon_end'])
    icon_ecef_ram_vector = (L1_dict['icon_ecef_ram_vector_start'] + L1_dict['icon_ecef_ram_vector_end'])/2
    icon_velocity = (L1_dict['icon_velocity_start'] + L1_dict['icon_velocity_end'])/2
    start_path = L1_dict['interferometer_start_path']
    end_path   = L1_dict['interferometer_end_path']

    
    ### Load instrument constants 
    instrument = get_instrument_constants(emission_color, start_path, end_path)
    if bin_size is None:
        bin_size = instrument['bin_size']

    #### Remove Satellite Velocity
    I = remove_satellite_velocity(Iraw, icon_velocity, icon_ecef_ram_vector, mighti_ecef_vectors, 
                                  instrument['phase_to_wind_factor'])
                         
    #### Bin data
    if bin_method == 'before':
        I        = bin_image(bin_size, I)
        tang_lat = bin_array(bin_size, tang_lat)
        tang_lon = bin_array(bin_size, tang_lon, lon=True)
        tang_alt = bin_array(bin_size, tang_alt)
        mighti_ecef_vectors = bin_image(bin_size, mighti_ecef_vectors)

                                                 
    #### Determine geographical locations of inverted wind
    lat, lon, alt = attribute_measurement_location(tang_lat, tang_lon, tang_alt,
                                                   integration_order=integration_order)
    
    #### Onion-peel interferogram
    Ip = perform_inversion(I, tang_alt, icon_alt, top_layer=top_layer, integration_order=integration_order,
                           account_for_local_projection=account_for_local_projection,
                           zero_phase=instrument['zero_phase'], phase_offset=instrument['phase_offset'],
                           Nignore=instrument['Nignore'])
                           
    #### Bin data
    if bin_method == 'after':
        Ip  = bin_image(bin_size, Ip)
        lat = bin_array(bin_size, lat)
        lon = bin_array(bin_size, lon, lon=True)
        alt = bin_array(bin_size, alt)
        mighti_ecef_vectors = bin_image(bin_size, mighti_ecef_vectors)


    #### Extract wind
    v_inertial, ve_inertial, a, ae = extract_wind(Ip, instrument['zero_phase'] + zero_phase_addition,
                                                      instrument['phase_offset'],
                                                      instrument['min_amp'], 
                                                      instrument['Nignore'],
                                                      instrument['phase_to_wind_factor'])
        

    #### Calculate azimuth angles at measurement locations
    icon_latlonalt = np.array([icon_lat, icon_lon, icon_alt])
    az = los_az_angle(icon_latlonalt, lat, lon, alt)

    #### Transform from inertial to rotating coordinate frame
    v = remove_Earth_rotation(v_inertial, az, lat, lon, alt)
    ve = ve_inertial.copy() # No uncertainty added in this process
    

    np.savez(L21_fn,
             los_wind                     = v,
             los_wind_error               = ve,
             lat                          = lat,
             lon                          = lon,
             alt                          = alt,
             time                         = time,
             exp_time                     = exp_time,
             az                           = az,
             emission_color               = emission_color,
             resolution_along_track       = np.nan,  # TODO
             resolution_cross_track       = np.nan, # TODO
             resolution_alt               = np.nan, # TODO
             icon_alt                     = icon_alt,
             icon_lat                     = icon_lat,
             icon_lon                     = icon_lon,
             fringe_amplitude             = a,
             fringe_amplitude_error       = ae,
             mighti_ecef_vectors          = mighti_ecef_vectors,
             icon_velocity_ecef_vector    = icon_velocity * icon_ecef_ram_vector,
             file_creation_time           = datetime.now(),
             source_files                 = np.concatenate((source_files,[L1_fn])),
             )
    
    return 0

