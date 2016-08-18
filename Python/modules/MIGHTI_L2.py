# A module for functions used for the conversion of MIGHTI Level 1 files to Level 2.1 and 2.2 files

import numpy as np
import ICON
import bisect
from scipy import integrate
from datetime import datetime, timedelta
import netCDF4
import getpass # for determining who is running the script
import glob


############################################################################################################
##########################################       Level 2.1       ###########################################
############################################################################################################



def get_instrument_constants(emission_color, start_path, end_path):
    '''
    The instrument constants needed for the MIGHTI L2.1 analysis.
    INPUTS:
        emission_color -- TYPE:str, 'green' or 'red'
        start_path     -- TYPE:float, UNITS:m. optical path difference at left edge of interferogram
        end_path       -- TYPE:float, UNITS:m. optical path difference at right edge of interferogram
    OUTPUTS:
        instr_params   -- TYPE:dict, dictionary of instrument parameters. (see below)
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
    Find the mean angle, taking into account 0/360 crossover. For example,
    circular_mean(10,50) is 30, but circular_mean(350,20) is 5.
    INPUTS:
        angle0  -- TYPE:float or array, UNITS:deg. An angle in degrees.
        angle1  -- TYPE:float or array, UNITS:deg. An angle in degrees.
    OUTPUTS:
        angle   -- TYPE:float or array, UNITS:deg. The circular mean of the two
                   input angles.
    '''
    x = np.rad2deg(np.angle((np.exp(1j*np.deg2rad(angle0)) + np.exp(1j*np.deg2rad(angle1)))/2.))
    x = np.mod(x,360.)
    return x




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
        if any( ze < 90. ) or any( ze > 180. ):
            raise Exception('Angle must be between 90 and 180, exclusive.')
    elif ( ze < 90. ) or ( ze > 180.):
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
    
    

def bin_uncertainty(b, ye):
    '''
    Determine the uncertainty of a binned array from the uncertainty
    of the un-binned array. Specifically:
    If the array y has uncertainty given by array ye, then the array
        y_b = bin_array(b, y)
    has uncertainty given by the array
        ye_b = bin_uncertainty(b, ye)
    INPUTS:
        b    -- TYPE:int,          The number of rows to bin together
        ye   -- TYPE:array(ny),    The uncertainty of the pre-binned data
    OUTPUTS:
        ye_b -- TYPE:array(ny_b), The uncertainty of the binned data
    '''
    ny = len(ye)
    ny_b = int(np.ceil(1.0*ny/b))
    ye_b = np.zeros(ny_b, dtype=ye.dtype)
    for i in range(0,ny_b): # bin from the end to the beginning.
        i_new   = ny_b-i-1
        i_start = ny-(i+1)*b
        i_stop  = ny-i*b
        
        # grab the samples to be binned
        if np.mod(ny,b)!=0 and i_new==0: # special case in case ny is not divisible by b
            ye_samps = ye[:i_stop]
        else: # grab 
            ye_samps = ye[i_start:i_stop]

        ye_b[i_new] = 1.0/len(ye_samps) * np.sqrt(np.sum(ye_samps**2))
        
    return ye_b
    
    
    
    
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

    
    
    
    
def perform_inversion(I, tang_alt, icon_alt, I_phase_uncertainty, I_amp_uncertainty, 
                      zero_phase, phase_offset, Nignore,
                      top_layer='exp', integration_order=0, 
                      account_for_local_projection=True):
    '''
    Perform the onion-peeling inversion on the interferogram to return
    a new interferogram, whose rows refer to specific altitudes. In effect,
    this function undoes the integration along the line of sight.
    INPUTS:
        I           -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, MIGHTI interferogram.
        tang_alt    -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
        icon_alt    -- TYPE:float,        UNITS:km.   Altitude of the satellite.
        I_phase_uncertainty -- TYPE:array(ny), UNITS:rad. Uncertainty in the unwrapped, mean phase of each row of I.
                                                          This is provided in L1 file.
        I_amp_uncertainty   -- TYPE:array(ny), UNITS:arb. Uncertainty in the summed amplitude of each row of I.
                                                          This is provided in L1 file.
        zero_phase  -- TYPE:float,        UNITS:rad.   The phase angle which is equivalent 
                                                       to a wind value of zero.
        phase_offset-- TYPE:float,        UNITS:rad.   An offset to avoid 2pi ambiguities.
        Nignore     -- TYPE:int,          UNITS:pixel. The number of columns at the
                                          beginning and end of the interferogram to ignore due to phase
                                          distortion from the filtering.
    OPTIONAL INPUTS:
        top_layer   -- TYPE:str,          'thin': assume VER goes to zero above top layer
                                          'exp':  assume VER falls off exponentially in altitude (default)
        integration_order -- TYPE:int,     0: Use Riemann-sum rule for discretizing line-of-sight integral (default)
                                           1: Use trapezoidal rule for discretizing line-of-sight integral
        account_for_local_projection   -- TYPE:bool.   If False, a simple inversion is used.
                                          If True, the inversion accounts for the fact that the ray is not 
                                          perfectly tangent to each shell at each point along the ray. 
                                          (default True)

    OUTPUTS:
        Ip          -- TYPE:array(ny,nx),    UNITS:arb. The complex-valued, onion-peeled interferogram.
        phase       -- TYPE:array(ny),       UNITS:rad. The unwrapped, mean phase of each row of Ip.
        amp         -- TYPE:array(ny),       UNITS:arb. The amplitude of each row of Ip.
        phase_uncertainty -- TYPE:array(ny), UNITS:rad. The uncertainty of phase
        amp_uncertainty   -- TYPE:array(ny), UNITS:rad. The uncertainty of amp
    '''
    
    ny,nx = np.shape(I)
    
    # Create the path matrix
    D = create_observation_matrix(tang_alt, icon_alt, top_layer=top_layer, integration_order=integration_order)
    
    
    
    ######### Onion-peeling inversion and amp/phase extraction #########
    # The inversion will proceed in different ways depending on whether
    # we will try to account for the local horizontal projection.
    phase = np.zeros(ny) # phases at each altitude
    if not account_for_local_projection:
        
        # This is implemented with a simple linear inversion
        Ip = np.linalg.solve(D,I)
        for i in range(ny):
            phase[i] = extract_phase_from_row(Ip[i,:], zero_phase, phase_offset, Nignore)
        
    else:
        # The problem becomes nonlinear, but still solvable in closed form.
        # This code implements Eq (9) in the MIGHTI L2 Space Science Reviews
        # paper (Harding et al. 2016).
        
        B = create_local_projection_matrix(tang_alt, icon_alt)
        Ip = np.zeros((ny,nx), dtype=complex) # onion-peeled interferogram

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
            
    amp = np.sum(abs(Ip),axis=1)        
            



    ######### Uncertainty propagation #########
    # Uncertainties can be propagated using simple linear inversion formula
    # (i.e., as if account_for_local_projection=False) to a very good approximation
    # (less than 1% error).
    
    ### Step 0: Characterize L1 and L2.1 interferograms with a single amp/phase per row
    ph_L1 = np.zeros(ny)
    for i in range(ny):
        ph_L1[i] = extract_phase_from_row(I[i,:], zero_phase, phase_offset, Nignore)
    A_L1 = np.sum(abs(I),axis=1)
    ph_L2 = phase.copy() # this was calculated above
    A_L2 = amp.copy() # this was calculated above
    # If amp is exactly zero (unlikely in practice), then replace it with a small number
    # so that uncertainties can be calculated.
    A_L2[A_L2==0.0] = 1e-6
        
    ### Step 1: Transform amp/phase uncertainties to real/imag uncertainties
    # Each row will have a 2x2 covariance matrix describing the real and imaginary parts
    cov_real_imag_L1 = np.zeros((ny,2,2))
    for m in range(ny):
        # Jacobian of transformation from ampl/phase to real/imag.
        J = np.array([[np.cos(ph_L1[m]), -A_L1[m]*np.sin(ph_L1[m])],
                      [np.sin(ph_L1[m]),  A_L1[m]*np.cos(ph_L1[m])]])
        cov_amp_phase = np.diag([I_amp_uncertainty[m], I_phase_uncertainty[m]])**2 # assuming uncorrelated
        cov_real_imag_L1[m,:,:] = J.dot(cov_amp_phase).dot(J.T) # Error propagation

    ### Step 2: Propagate uncertainties through the path length inversion
    # Treat real and imaginary parts separately.
    # Build covariance matrix of vector of real parts and of imag parts.
    cov_real_L1 = np.diag(cov_real_imag_L1[:,0,0]) # assume rows are uncorrelated
    cov_imag_L1 = np.diag(cov_real_imag_L1[:,1,1]) # assume rows are uncorrelated
    # Standard uncertainty propagation for matrix multiplication
    Dinv = np.linalg.inv(D)
    cov_real_L2 = Dinv.dot(cov_real_L1).dot(Dinv.T)
    cov_imag_L2 = Dinv.dot(cov_imag_L1).dot(Dinv.T)
    sigma_real_L2 = np.sqrt(np.diag(cov_real_L2))
    sigma_imag_L2 = np.sqrt(np.diag(cov_imag_L2))
    
    ### Step 3: Transform back to amp/phase #########
    # Each row will have a 2x2 covariance matrix describing the amplitude and phase
    cov_amp_phase_L2 = np.zeros((ny,2,2))
    for m in range(ny):
        # Jacobian of transformation from ampl/phase to real/imag
        J = np.array([[np.cos(ph_L2[m]), -A_L2[m]*np.sin(ph_L2[m])],
                      [np.sin(ph_L2[m]),  A_L2[m]*np.cos(ph_L2[m])]])
        # Jacobian of transformation from real/imag to ampl/phase
        Jinv = np.linalg.inv(J)
        cov_real_imag = np.diag([sigma_real_L2[m], sigma_imag_L2[m]])**2 # assume uncorrelated
        cov_amp_phase_L2[m,:,:] = Jinv.dot(cov_real_imag).dot(Jinv.T)
    # Extract amplitude and phase uncertainties
    amp_uncertainty = np.sqrt(cov_amp_phase_L2[:,0,0])
    phase_uncertainty = np.sqrt(cov_amp_phase_L2[:,1,1])
            
    return Ip, phase, amp, phase_uncertainty, amp_uncertainty




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
    function may need to change: create_observation_matrix, perform_inversion.
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






def interpolate_linear(x, y, x0, extrapolation='hold', prop_err = False, yerr = None):
    '''
    Linear interpolation of the function y = f(x) to the location x0.
    x and y are vectors comprising samples of this function. There is also
    an option to propagate errors to the interpolated value. This function is
    5 times faster than scipy.interpolate.interp1d, and allows for
    zero-order-hold extrapolation. If you are interpolating to many points, 
    then scipy.interpolate.interp1d is likely faster. 

    INPUTS:
        x     -- TYPE:array(n), UNITS:arb. Independent variable of samples of function.
        y     -- TYPE:array(n), UNITS:arb. Dependent variable of samples of function.
        x0    -- TYPE:float,    UNITS:arb. Independent variable of interpolation point.
    OPTIONAL INPUTS:
        extrapolation -- TYPE:str,        'hold': extrapolate by using values at end points (default)
                                          'none': do not extrapolate. Points will be np.nan
        prop_err      -- TYPE:bool,       True:  propagate errors from original to interpolated
                                                 value, and return an extra output; yerr must
                                                 be specified as an input. 
                                          False: do not propagate errors, and return only one
                                                 output (default).
        yerr          -- TYPE:array(n), UNITS:arb. Error in y, to be propagated to interpolated value.
    OUTPUTS:
        y0    -- TYPE:float,    UNITS:arb. Interpolated value.
    OPTIONAL OUTPUT (if prop_err = True):
        y0err -- TYPE:float,    UNTIS:arb. Propagated error of y0.
    '''
    
    if prop_err and yerr is None:
        raise Exception('If prop_err=True, then yerr must be specified')
    
    j0 = bisect.bisect(x,x0) - 1 # index to the left
    j1 = j0 + 1 # index to the right
    y0err = np.nan
    # Handle extrapolations
    if j0 == -1:
        if extrapolation=='hold':
            y0 = y[0]
            if prop_err:
                y0err = yerr[0]
        elif extrapolation == 'none':
            y0 = np.nan
        else: 
            raise Exception('"%s" not understood' % extrapolation)
    elif j1 == len(x):
        if extrapolation=='hold':
            y0 = y[-1]
            if prop_err:
                y0err = yerr[-1]
        elif extrapolation == 'none':
            y0 = np.nan
        else: 
            raise Exception('"%s" not understood' % extrapolation)
    else: # linear interpolation
        w1 = (x0-x[j0]) / (x[j1]-x[j0]) # weight of y[j1]
        w0 = 1.0-w1 # weight of y[j0]
        y0 = w0*y[j0] + w1*y[j1]
        if prop_err:
            # What is the best way to interpolate errors? 
            # Statistically correct way, but yields counterintuitive results, such as
            # a higher error near the sample points than between them:
            #y0err = np.sqrt(w0**2*yerr[j0]**2 + w1**2*yerr[j1]**2)
            # Simple way: just interpolate errors
            y0err = w0*yerr[j0] + w1*yerr[j1]
    if prop_err:
        return y0, y0err
    else:
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
                              the level 2.1 processing. See documentation for 
                              level1_dict_to_level21(...) for required keys.
    TODO:
        - lots of stuff, mostly just waiting on finalization of L1 file.
        - interferometer start path and end path
        - time in better format
        - exposure time different than end minus start?
        - will ICON velocity be specified in m/s eventually?
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
    L1_dict['icon_velocity_start']         = f['ICON_VELOCITY_START'][:].item()*1e3 # convert to m/s
    L1_dict['icon_velocity_end']           = f['ICON_VELOCITY_END'][:].item()*1e3   # convert to m/s
    L1_dict['source_files']                = [] # TODO
    tsec_start                             = f['IMAGE_TIME_START'][:].item()
    tsec_end                               = f['IMAGE_TIME_END'][:].item()
    L1_dict['time_start']                  = datetime(2015,1,1) + timedelta(seconds=tsec_start)
    L1_dict['time_end']                    = datetime(2015,1,1) + timedelta(seconds=tsec_end)
    L1_dict['exp_time']                    = f['IMAGE_TIME_END'][:].item() - f['IMAGE_TIME_START'][:].item()
    L1_dict['interferometer_start_path']   = 4.62e-2
    L1_dict['interferometer_end_path']     = 5.50e-2
    
    # Dummy placeholder code for reading global attributes
    nc_attrs = f.ncattrs()
    for nc_attr in nc_attrs:
        a = f.getncattr(nc_attr)
    
    
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
                                   the level 2.1 processing. See documentation for 
                                   level1_dict_to_level21(...) for required keys.
    '''
    
    L1_dict = {}
    L1_dict['L1_fn'] = L1_uiuc_fn
    
    npzfile = np.load(L1_uiuc_fn)
    for key in npzfile.keys():
        if key in ['emission_color','icon_alt_start','icon_alt_end',
                   'icon_lat_start','icon_lat_end','icon_lon_start','icon_lon_end',
                   'icon_velocity_start','icon_velocity_end','time','exp_time',
                   'interferometer_start_path','interferometer_end_path',]:
            L1_dict[key] = npzfile[key].item()
        else:
            L1_dict[key] = npzfile[key]
    npzfile.close()
    
    ####### Hacks for backwards compatibility #######
    # If interferogram is given as real/imaginary, then change
    # it to amp/phase
    if 'Ir' in L1_dict.keys():
        Ir = L1_dict['Ir']
        Ii = L1_dict['Ii']
        I = Ir + 1j*Ii
        L1_dict['I_amp'] = abs(I)
        L1_dict['I_phase'] = np.angle(I)
    # If uncertainties aren't in the file, add placeholders (0.01% and 2 mrad)
    if 'I_amp_uncertainty' not in L1_dict.keys():
        L1_dict['I_amp_uncertainty'] = 0.0001 * np.sum(L1_dict['I_amp'],axis=1)
    if 'I_phase_uncertainty' not in L1_dict.keys():
        L1_dict['I_phase_uncertainty'] = 0.002 * np.ones(np.shape(L1_dict['I_phase'])[0])
    # Change time to time_start and time_end
    if 'time' in L1_dict:
        L1_dict['time_start'] = L1_dict['time']
        L1_dict['time_end']   = L1_dict['time'] + timedelta(seconds=L1_dict['exp_time'])
        del L1_dict['time']
    # Ensure longitudes are defined [0,360]
    for key in ['icon_lon_start','icon_lon_end','tang_lon_start','tang_lon_end']:
        L1_dict[key] = np.mod(L1_dict[key],360.)
        
    
    return L1_dict
    





def level1_dict_to_level21_dict(L1_dict, zero_phase_addition = 0.0, top_layer = 'exp', 
                                integration_order = 0, account_for_local_projection = True, 
                                bin_size = None):
    '''
    High-level function to run the Level 2.1 processing. It takes a dictionary (containing
    input variables extracted from a Level 1 file) and outputs a dictionary (containing 
    output variables, which can be written to a file with another function).
    
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
                                             time_start
                                             time_end
                                             exp_time
                                             interferometer_start_path
                                             interferometer_end_path  
                                                                  
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
                                                                          
    OUTPUTS:
        L21_dict            -- TYPE:dict. A dictionary containing output variables of the Level 2.1 processing:
                                             los_wind 
                                             los_wind_error
                                             lat 
                                             lon  
                                             alt    
                                             time_start 
                                             time_end
                                             exp_time
                                             az  
                                             emission_color 
                                             resolution_along_track 
                                             resolution_cross_track
                                             resolution_alt   
                                             icon_alt    
                                             icon_lat  
                                             icon_lon 
                                             fringe_amplitude 
                                             fringe_amplitude_error  
                                             mighti_ecef_vectors  
                                             icon_velocity_ecef_vector 
                                             file_creation_time
                                             source_files 
                                             bin_size 
                                             top_layer
                                             integration_order
                                             zero_phase
            
    TODO:
        - error handling
        - Nignore will be done in L1, not L2.1, eventually.
    
    '''

    ###  Load parameters from input dictionary
    Iraw = L1_dict['I_amp']*np.exp(1j*L1_dict['I_phase'])
    I_amp_uncertainty = L1_dict['I_amp_uncertainty']
    I_phase_uncertainty = L1_dict['I_phase_uncertainty']
    emission_color = L1_dict['emission_color']
    source_files = L1_dict['source_files']
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
    I        = bin_image(bin_size, I)
    tang_lat = bin_array(bin_size, tang_lat)
    tang_lon = bin_array(bin_size, tang_lon, lon=True)
    tang_alt = bin_array(bin_size, tang_alt)
    mighti_ecef_vectors = bin_image(bin_size, mighti_ecef_vectors)
    I_amp_uncertainty   = bin_uncertainty(bin_size, I_amp_uncertainty)
    I_phase_uncertainty = bin_uncertainty(bin_size, I_phase_uncertainty)
    
    
    
    #### Determine geographical locations of inverted wind
    lat, lon, alt = attribute_measurement_location(tang_lat, tang_lon, tang_alt,
                                                   integration_order=integration_order)
    
    
    #### Onion-peel interferogram
    zero_phase = instrument['zero_phase'] + zero_phase_addition
    Ip, phase, amp, phase_uncertainty, amp_uncertainty = perform_inversion(I, tang_alt, icon_alt, 
                           I_phase_uncertainty, I_amp_uncertainty,
                           zero_phase, instrument['phase_offset'], instrument['Nignore'],
                           top_layer=top_layer, integration_order=integration_order,
                           account_for_local_projection=account_for_local_projection)


    #### Transform from phase to wind
    v_inertial             = instrument['phase_to_wind_factor'] * phase
    v_inertial_uncertainty = instrument['phase_to_wind_factor'] * phase_uncertainty
        

    #### Calculate azimuth angles at measurement locations
    icon_latlonalt = np.array([icon_lat, icon_lon, icon_alt])
    az = los_az_angle(icon_latlonalt, lat, lon, alt)

    #### Transform from inertial to rotating coordinate frame
    v = remove_Earth_rotation(v_inertial, az, lat, lon, alt)
    v_uncertainty = v_inertial_uncertainty.copy() # No appreciable uncertainty added in this process
    
    # Make a L2.1 dictionary
    L21_dict = {
             'los_wind'                     : v,
             'los_wind_error'               : v_uncertainty,
             'lat'                          : lat,
             'lon'                          : lon,
             'alt'                          : alt,
             'time_start'                   : L1_dict['time_start'],
             'time_end'                     : L1_dict['time_end'],
             'exp_time'                     : exp_time,
             'az'                           : az,
             'emission_color'               : emission_color,
             'resolution_along_track'       : np.nan, # TODO
             'resolution_cross_track'       : np.nan, # TODO
             'resolution_alt'               : np.nan, # TODO
             'icon_alt'                     : icon_alt,
             'icon_lat'                     : icon_lat,
             'icon_lon'                     : icon_lon,
             'fringe_amplitude'             : amp,
             'fringe_amplitude_error'       : amp_uncertainty,
             'mighti_ecef_vectors'          : mighti_ecef_vectors,
             'icon_velocity_ecef_vector'    : icon_velocity * icon_ecef_ram_vector,
             'file_creation_time'           : datetime.now(),
             'source_files'                 : np.concatenate((source_files,[L1_fn])),
             'bin_size'                     : bin_size,
             'top_layer'                    : top_layer,
             'integration_order'            : integration_order,
             'zero_phase'                   : zero_phase,
    }
        
    return L21_dict
   
   
   
   
   
def _create_variable(ncfile, name, value, format_nc='f8', format_fortran='F', dimensions=(), zlib=True, complevel=6, 
                    shuffle=True,  depend_0=None, depend_1=None, depend_2=None, chunk_sizes=None, desc='', 
                    display_type='scalar', field_name='', fill_value=None,label_axis='', bin_location=0.5, 
                    time_base='FIXED: 1970 (POSIX)', time_scale='UTC', units='', valid_min=None, valid_max=None, 
                    notes='', var_type='data'):
    '''
    A helper function to write a variable to a netCDF file.
    INPUTS:
        See above. Notes:
            - fill_value = None --> default fill values will be used, if they exist. See netCDF4.default_fillvals
            - display_type: for now, 'scalar', 'altitude_profile', or 'image' will be used
            - var_type: one of 'data', 'support_data', 'metadata', 'ignore_data'
            - format_fortran: Used by ISTP. See http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
            - except as specified above, if a variable attribute is left as the default None, it will not be written to the file
    OUTPUT:
        The netCDF4._netCDF4.Variable object that was created and written to
    
    '''    
    
    # Rudimentary error-checking:
    valid_var_types = ['data','support_data','metadata','ignore_data']
    if var_type not in valid_var_types:
        raise Exception('var_type="%s" is not valid. Try one of: %s' % (var_type, valid_var_types) )
    if len(desc) > 80:
        raise Exception('"desc" is too long (%i chars). Shorten to 80 characters:\n"%s"' % (len(desc),desc))
    if len(field_name) > 30:
        raise Exception('field_name="%s" is too long (%i chars). Shorten to 30 characters:\n"%s"' % (len(field_name),field_name))
    if len(label_axis) > 10:
        raise Exception('label_axis="%s" is too long (%i chars). Shorten to 10 characters.' % (len(label_axis),label_axis))
    
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
    var.Display_Type       = display_type 
    var.FieldNam           = field_name
    # Note: t_var._FillValue not expliclity needed since that is set by the createVariable function argument "fill_value"
    #var._FillValue         = fill_value
    if fill_value is not None:
        var.FillVal        = var._FillValue
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

   
   
    
    
def save_nc_level21(path, L21_dict):
    '''
    Take the output of the Level 2.1 processing and save it as a NetCDF4 file in the official format.
    NetCDF4 file conventions taken from "Science Operations Center Data Product Conventions" Rev 0.4,
    authored by Tori Fae., including the new filenaming requirement that will presumably be in Rev 0.5.
    
    INPUTS:
        path        -- TYPE:str.  The directory the file will be saved in, including trailing "/"
                                  (e.g., '/home/user/')
        L21_dict    -- TYPE:dict. A dictionary containing output variables of the Level 2.1 processing.
                                  See documentation for level1_dict_to_level21_dict(...) for details.
    OUTPUTS:
        L21_fn      -- TYPE:str.  The full path to the saved file.
    TO-DO:
      - figure out how to make np.nan the fill value
      - Are attributes _FillValue and FillVal needed for variables that are strings? python-NetCDF4 doesn't seem to support this
      - Maybe: Fill in more notes for each variable
      - How will sensor and the 4 versions be specified?
      - Have a second set of eyes look at descriptions of each variable
      - Should dimensions be labeled the same as variables? Altitude/Vector/Epoch. Should Depend_0 point to vars or dims?
    '''


    data_version = 1 # TODO: how will this be calculated? It goes into global attr Data_Version
    software_version = 0.01 # TODO: how will this be determined and stored? It goes into global attr Software_Version
    version_for_filename = 1 # TODO: How is this related to the two above? Rec'd major software version number
    revision_for_filename = 0 # TODO: How to determine this?
    
    # TODO: How will sensor be determined? Will it be in L1 file?
    source_files = ' '.join(L21_dict['source_files'])
    if '_A_' in source_files:
        sensor = 'A'
    elif '_B_' in source_files:
        sensor = 'B'
    else:
        raise Exception('Cannot Determine Sensor from "%s"'%source_files) 

    #################### Compile variables to write in file ######################
    ### Timing:
    t_start = L21_dict['time_start']
    t_end   = L21_dict['time_end']
    t_mid   = t_start + timedelta(seconds=(t_end - t_start).total_seconds()/2) # middle of exposure
    t_start_msec = (t_start - datetime(1970,1,1)).total_seconds()*1e3 # milliseconds since epoch
    t_end_msec   = (t_end   - datetime(1970,1,1)).total_seconds()*1e3
    t_mid_msec   = (t_mid   - datetime(1970,1,1)).total_seconds()*1e3
    t_start_msec = np.int64(np.round(t_start_msec)) # cast to signed 64 bit integer
    t_end_msec   = np.int64(np.round(t_end_msec)) 
    t_mid_msec   = np.int64(np.round(t_mid_msec))
    t_file  = datetime.now()   # time this file was created  
    ### Who's running this process
    user_name = getpass.getuser()
    ### Parent files
    parents = '' # This will go in global attr Parents
    for source_fn in L21_dict['source_files']:
        s = source_fn.split('/')[-1].split('.')
        pre = '.'.join(s[:-1])
        post = s[-1].upper()
        parents += '%s > %s, ' % (post, pre)
    if parents: parents = parents[:-2] # trim trailing comma


    ######################### Open file for writing ##############################
    L21_fn = 'ICON_L2_1_MIGHTI-%s_%s_v%02i_r%02i.nc' % (sensor,t_mid.strftime('%Y-%m-%d_%H.%M.%S'),
                                                        version_for_filename, revision_for_filename)
    L21_full_fn = '%s%s'%(path, L21_fn)
    ncfile = netCDF4.Dataset(L21_full_fn,mode='w',format='NETCDF4') 

    try:
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
                                                "  * Neutral Wind profiles  Dr. Jon Makela\n",
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
        ncfile.Calibration_File =               '((TODO: zero wind cal file, and pass through L1 cal files?))'
        ncfile.Conventions =                    'SPDF ISTP/IACG Modified for NetCDF'
        ncfile.Data_Level =                     'L2.1'
        ncfile.Data_Type =                      'DP21 > Data Product 2.1: Line-of-sight Wind Profile'
        ncfile.Data_Version =                   np.uint16(data_version)
        ncfile.Date_End =                       t_mid.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC' # single measurement: use midpoint
        ncfile.Date_Start =                     t_mid.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC' # single measurement: use midpoint
        ncfile.Description =                    'ICON MIGHTI Line-of-sight Winds (DP 2.1)'
        ncfile.Descriptor =                     'MIGHTI-%s > Michelson Interferometer for Global High-resolution ' % sensor+\
                                                'Thermospheric Imaging, Sensor %s' % sensor
        ncfile.Discipline =                     'Space Physics > Ionospheric Science'
        ncfile.File =                           L21_fn
        ncfile.File_Date =                      t_file.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC'
        ncfile.Generated_By =                   'ICON SDC > ICON UIUC MIGHTI L2.1 Processor v%.2f, B. J. Harding' % software_version
        ncfile.Generation_Date =                t_file.strftime('%Y%m%d')
        ncfile.History =                        'Version %i, %s, %s, ' % (data_version, user_name, t_file.strftime('%Y-%m-%dT%H:%M:%S')) +\
                                                'MIGHTI L2.1 Processor v%.2f ' % software_version +\
                                                '((TODO: Prepend previous history of this file, if it exists))'
        ncfile.HTTP_LINK =                      'http://icon.ssl.berkeley.edu/Instruments/MIGHTI'
        ncfile.Instrument =                     'MIGHTI-%s' % sensor
        ncfile.Instrument_Type =                'Imagers (space)'
        ncfile.Link_Text =                      'MIGHTI Line-of-sight Wind Profile (DP 2.1)'
        ncfile.Link_Title =                     'ICON MIGHTI'
        ncfile.Logical_File_ID =                L21_fn[:-3]
        ncfile.Logical_Source =                 'ICON_L2_1_MIGHTI-%s_' % (sensor,)
        ncfile.Logical_Source_Description =     'MIGHTI Sensor %s - Line-of-sight Wind Profile'
        ncfile.Mission_Group =                  'Ionospheric Investigations'
        ncfile.MODS =                           ncfile.History
        ncfile.Parents =                        parents
        ncfile.PI_Affiliation =                 'UC Berkeley > SSL'
        ncfile.PI_Name =                        'T. J. Immel'
        ncfile.Project =                        'NASA > ICON'
        ncfile.Rules_of_Use =                   'Public Data for Scientific Use'
        ncfile.Software_Version =               'ICON SDC > ICON UIUC MIGHTI L2.1 Processor v%.2f' % software_version
        ncfile.Source_Name =                    'ICON > Ionospheric Connection Explorer'
        ncfile.Spacecraft_ID =                  'NASA > ICON - 493'
        ncfile.Text =                           'ICON explores the boundary between Earth and space  the ionosphere  ' +\
                                                'to understand the physical connection between our world and the immediate '+\
                                                'space environment around us. Visit \'http://icon.ssl.berkeley.edu\' for more details.'
        ncfile.Text_Supplement =                '((TODO: Insert reference to Harding et al [2016] paper))'
        ncfile.Time_Resolution =                '%.1f seconds' % L21_dict['exp_time']
        ncfile.Title =                          'ICON MIGHTI Line-of-sight Wind Profile (DP 2.1)'


        ################################## Dimensions ########################################
        n = len(L21_dict['alt'])
        ncfile.createDimension('Epoch',0)
        ncfile.createDimension('Altitude', n)
        ncfile.createDimension('Vector',3)


        ################################## Variables #########################################



        ######### Timing Variables #########

        # Time midpoint (the official required "Epoch" variable)
        var = _create_variable(ncfile, 'Epoch', t_mid_msec, 
                              dimensions=(),
                              format_nc='i8', format_fortran='I', desc='Sample time, midpoint of exposure. Number of msec since Jan 1, 1970.', 
                              display_type='scalar', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
                              units='ms', valid_min=0, valid_max=1000*365*86400e3, var_type='support_data', chunk_sizes=1,
                              notes='')

        # Time start
        var = _create_variable(ncfile, 'Epoch_Start', t_start_msec, 
                              dimensions=(),
                              format_nc='i8', format_fortran='I', desc='Sample time, start of exposure. Number of msec since Jan 1, 1970.', 
                              display_type='scalar', field_name='Time Start', fill_value=None, label_axis='Time Start', bin_location=0.0,
                              units='ms', valid_min=0, valid_max=1000*365*86400e3, var_type='support_data', chunk_sizes=1,
                              notes='')

        # Time end
        var = _create_variable(ncfile, 'Epoch_End', t_end_msec, 
                              dimensions=(),
                              format_nc='i8', format_fortran='I', desc='Sample time, end of exposure. Number of msec since Jan 1, 1970.', 
                              display_type='scalar', field_name='Time End', fill_value=None, label_axis='Time End', bin_location=1.0,
                              units='ms', valid_min=0, valid_max=1000*365*86400e3, var_type='support_data', chunk_sizes=1,
                              notes='')


        ######### Data Location and Direction Variables #########

        # Altitude
        val = L21_dict['alt']*1e3 # convert to meters
        var_alt = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Altitude'%sensor, val, 
                              dimensions=('Altitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='WGS84 altitude of each wind sample', 
                              display_type='altitude_profile', field_name='Altitude', fill_value=None, label_axis='Altitude', bin_location=0.5,
                              units='m', valid_min=0, valid_max=1e10, var_type='support_data', chunk_sizes=[1],
                              notes='')

        # Latitude
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Latitude'%sensor, L21_dict['lat'], 
                              dimensions=('Altitude'), depend_0 = var_alt.name,
                              format_nc='f8', format_fortran='F', desc='WGS84 latitude of each wind sample', 
                              display_type='altitude_profile', field_name='Latitude', fill_value=None, label_axis='Latitude', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='support_data', chunk_sizes=[1],
                              notes='')

        # Longitude
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Longitude'%sensor, L21_dict['lon'], 
                              dimensions=('Altitude'), depend_0 = var_alt.name,
                              format_nc='f8', format_fortran='F', desc='WGS84 longitude of each wind sample', 
                              display_type='altitude_profile', field_name='Longitude', fill_value=None, label_axis='Longitude', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[1],
                              notes='')

        # Azimuth angle of line of sight
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Line-of-sight_Azimuth'%sensor, L21_dict['az'], 
                              dimensions=('Altitude'), depend_0 = var_alt.name,
                              format_nc='f8', format_fortran='F', desc='Azimuth angle of the line of sight at the tangent point. Deg East of North.', 
                              display_type='altitude_profile', field_name='Line-of-sight Azimuth', fill_value=None, label_axis='Azimuth', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[1],
                              notes='')


        ######### Data Variables #########

        # Line-of-sight wind profile
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Line-of-sight_Wind'%sensor, L21_dict['los_wind'], 
                              dimensions=('Altitude'), depend_0 = var_alt.name,
                              format_nc='f8', format_fortran='F', desc='Line-of-sight horizontal wind profile. A positive wind is towards MIGHTI.', 
                              display_type='altitude_profile', field_name='Line-of-sight Wind', fill_value=None, label_axis='LoS Wind', bin_location=0.5,
                              units='m/s', valid_min=-1e10, valid_max=1e10, var_type='data', chunk_sizes=[1],
                              notes='')

        # Line-of-sight wind error profile
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Line-of-sight_Wind_Error'%sensor, L21_dict['los_wind_error'], 
                              dimensions=('Altitude'), depend_0 = var_alt.name,
                              format_nc='f8', format_fortran='F', desc='Line-of-sight Horizontal Wind Error Profile', 
                              display_type='altitude_profile', field_name='Line-of-sight Wind Error', fill_value=None, label_axis='Wind Error', bin_location=0.5,
                              units='m/s', valid_min=-1e10, valid_max=1e10, var_type='data', chunk_sizes=[1],
                              notes='')

        # Fringe amplitude profile (TODO: will this be replaced by a VER data product?)
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Fringe_Amplitude'%sensor, L21_dict['fringe_amplitude'], 
                              dimensions=('Altitude'), depend_0 = var_alt.name,
                              format_nc='f8', format_fortran='F', desc='Fringe Amplitude Profile', 
                              display_type='altitude_profile', field_name='Fringe Amplitude', fill_value=None, label_axis='Fringe Amp', bin_location=0.5,
                              units='arb', valid_min=-1e30, valid_max=1e30, var_type='data', chunk_sizes=[1],
                              notes='An approximate volume emission rate (VER) profile in arbitrary units. Technically this a profile of the visibility '+
                                    'of the fringes, which has a dependence on temperature and background emission.')

        # Fringe amplitude error profile (TODO: will this be replaced by a VER data product?)
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Fringe_Amplitude_Error'%sensor, L21_dict['fringe_amplitude_error'], 
                              dimensions=('Altitude'), depend_0 = var_alt.name,
                              format_nc='f8', format_fortran='F', desc='Fringe Amplitude Error Profile', 
                              display_type='altitude_profile', field_name='Fringe Amplitude Error', fill_value=None, label_axis='Amp Err', bin_location=0.5,
                              units='arb', valid_min=0, valid_max=1e30, var_type='data', chunk_sizes=[1],
                              notes='')


        ######### Other Metadata Variables #########

        # Emission color
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Emission_Color'%sensor, L21_dict['emission_color'].capitalize(), 
                              dimensions=(),
                              format_nc=str, format_fortran='A', desc='Emission used for wind estimate: "Red" or "Green"', 
                              display_type='scalar', field_name='Emission Color', fill_value=None, label_axis='Color', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=1,
                              notes='')

        # Zero wind phase
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Zero_Wind_Phase'%sensor, L21_dict['zero_phase'], 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The fringe phase corresponding to a wind of zero', 
                              display_type='scalar', field_name='Zero Wind Phase', fill_value=None, label_axis='Zero Phase', bin_location=0.5,
                              units='rad', valid_min=-1e30, valid_max=1e30, var_type='metadata', chunk_sizes=1,
                              notes='')

        # ICON velocity vector
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Spacecraft_Velocity_Vector'%sensor, L21_dict['icon_velocity_ecef_vector'], 
                              dimensions=('Vector'),# depend_0 = 'Vector',
                              format_nc='f8', format_fortran='F', desc='ICON\'s velocity vector in Earth-Centered, Earth-fixed coordinates', 
                              display_type='scalar', field_name='ICON Velocity Vector', fill_value=None, label_axis='S/C Vel', bin_location=0.5,
                              units='m/s', valid_min=-1e10, valid_max=1e10, var_type='metadata', chunk_sizes=[1],
                              notes='')

        # ICON latitude
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Spacecraft_Latitude'%sensor, L21_dict['icon_lat'], 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The WGS84 latitude of ICON', 
                              display_type='scalar', field_name='Spacecraft Latitude', fill_value=None, label_axis='S/C Lat', bin_location=0.5,
                              units='deg', valid_min=--90., valid_max=90., var_type='metadata', chunk_sizes=1,
                              notes='')

        # ICON longitude
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Spacecraft_Longitude'%sensor, L21_dict['icon_lon'], 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The WGS84 longitude of ICON', 
                              display_type='scalar', field_name='Spacecraft Longitude', fill_value=None, label_axis='S/C Lon', bin_location=0.5,
                              units='deg', valid_min=-0., valid_max=360., var_type='metadata', chunk_sizes=1,
                              notes='')

        # ICON altitude
        val = L21_dict['icon_alt']*1e3 # convert to m
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Spacecraft_Altitude'%sensor, val, 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The WGS84 altitude of ICON', 
                              display_type='scalar', field_name='Spacecraft Altitude', fill_value=None, label_axis='S/C Alt', bin_location=0.5,
                              units='m', valid_min=0., valid_max=1e10, var_type='metadata', chunk_sizes=1,
                              notes='')

        # Along-track resolution
        val = L21_dict['resolution_along_track']*1e3 # convert to meters
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Resolution_Along_Track'%sensor, val, 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The horizontal resolution in the spacecraft velocity direction', 
                              display_type='scalar', field_name='Along-Track Resolution', fill_value=None, label_axis='Hor Res AT', bin_location=0.5,
                              units='m', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
                              notes='')

        # Cross-track resolution
        val = L21_dict['resolution_cross_track']*1e3 # convert to meters
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Resolution_Cross_Track'%sensor, val, 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The horizontal resolution perpendicular to the spacecraft velocity direction', 
                              display_type='scalar', field_name='Cross-Track Resolution', fill_value=None, label_axis='Hor Res CT', bin_location=0.5,
                              units='m', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
                              notes='')

        # Altitude resolution
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Resolution_Altitude'%sensor, L21_dict['resolution_alt'], 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The vertical resolution', 
                              display_type='scalar', field_name='Vertical Resolution', fill_value=None, label_axis='Vert Res', bin_location=0.5,
                              units='m', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
                              notes='')

        # MIGHTI ECEF vectors
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Line-of-sight_Vector'%sensor, L21_dict['mighti_ecef_vectors'], 
                              dimensions=('Altitude','Vector'), depend_0 = var_alt.name, # depend_1 = 'Vector',
                              format_nc='f8', format_fortran='F', desc='The look direction of each MIGHTI line of sight, as a vector in ECEF', 
                              display_type='altitude_profile', field_name='Line-of-sight Vector', fill_value=None, label_axis='LoS Vec', bin_location=0.5,
                              units='', valid_min=-1., valid_max=1., var_type='metadata', chunk_sizes=[1,3],
                              notes='')


        # Bin Size
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Bin_Size'%sensor, L21_dict['bin_size'], 
                              dimensions=(),
                              format_nc='i1', format_fortran='I', desc='How many raw samples were binned vertically for each reported sample', 
                              display_type='scalar', field_name='Bin Size', fill_value=None, label_axis='Bin Size', bin_location=0.5,
                              units='', valid_min=0, valid_max=100000, var_type='metadata', chunk_sizes=1,
                              notes='')


        # Integration order
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Integration_Order'%sensor, L21_dict['integration_order'], 
                              dimensions=(),
                              format_nc='i1', format_fortran='I', desc='Order used to discretize the integral for inversion: 0=Riemann, 1=Trapezoidal', 
                              display_type='scalar', field_name='Order', fill_value=None, label_axis='Order', bin_location=0.5,
                              units='', valid_min=0, valid_max=10, var_type='metadata', chunk_sizes=1,
                              notes='')

        # How the top layer was handled in the inversion
        var = _create_variable(ncfile, 'ICON_L2_1_MIGHTI-%s_Top_Layer_Model'%sensor, L21_dict['top_layer'].capitalize(), 
                              dimensions=(),
                              format_nc=str, format_fortran='A', desc='How the top altitudinal layer is handled in the inversion: "Exp" or "Thin"', 
                              display_type='scalar', field_name='Top Layer', fill_value=None, label_axis='Top Layer', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=1,
                              notes='')    
            
            
        ncfile.close()
        
    except: # make sure the file is closed
        ncfile.close()
        raise
    
    return L21_full_fn
    
    
    
    
# def level1_to_level21()
# TODO: implement this high-level function
    
    
    
    
def level21_to_dict(L21_fns):
    ''' 
    Load a series of Level 2.1 files and return relevant variables in a dictionary.
    INPUTS:
        L21_fns  -- TYPE:list of str.  The paths to the Level 2.1 files to be loaded.
    OUTPUTS:
        L21_dict -- TYPE:dict. A dictionary containing the following variables. Most 
                               are provided as arrays of shape (ny,nt), where ny is the number
                               of altitude samples and nt is the number of time samples.
                    lat             -- TYPE:array(ny,nt), UNITS:deg. Sample latitudes.
                    lon             -- TYPE:array(ny,nt), UNITS:deg. Sample longitudes.
                    alt             -- TYPE:array(ny,nt), UNITS:km.  Sample altitudes.
                    los_wind        -- TYPE:array(ny,nt), UNITS:m/s. Line-of-sight wind component towards MIGHTI.
                    los_wind_error  -- TYPE:array(ny,nt), UNITS:m/s. Error in los_wind variable.
                    local_az        -- TYPE:array(ny,nt), UNITS:deg. Azimuth angle of vector pointing from 
                                       MIGHTI towards the sample location, at the sample location (deg E of N).
                    amp             -- TYPE:array(ny,nt), UNITS:arb. Fringe amplitude at sample locations
                    time            -- TYPE:array(nt).               Array of datetime objects, one per file.
                    exp_time        -- TYPE:array(nt),    UNITS:sec. Exposure time of each sample.
                    emission_color  -- TYPE:str,                     'red' or 'green'.
                    source_files    -- TYPE:list of str,             A copy of the input.
        
    '''
    # Open the first file to see how many altitude bins there are
    d = netCDF4.Dataset(L21_fns[0],mode='r')
    stub = 'ICON_L2_1_%s' % d.Instrument # The prefix of each variable
    ny = len(d.variables['%s_Altitude'%stub])
    nt = len(L21_fns)
    d.close()
    
    lat = np.zeros((ny,nt))
    lon = np.zeros((ny,nt))
    alt = np.zeros((ny,nt))
    los_wind = np.zeros((ny,nt))
    los_wind_error = np.zeros((ny,nt))
    local_az = np.zeros((ny,nt))   
    amp = np.zeros((ny,nt))
    time = np.zeros(nt, dtype=object)
    exp_time = np.zeros(nt)

    for i in range(nt):
        fn = L21_fns[i]
        d = netCDF4.Dataset(fn,mode='r')
        stub = 'ICON_L2_1_%s' % d.Instrument # The prefix of each variable
        emission_color = d.variables['%s_Emission_Color' % stub][...].lower()      

        lat[:,i] = d.variables['%s_Latitude' % stub][...]
        lon[:,i] = d.variables['%s_Longitude' % stub][...]
        alt[:,i] = 1e-3 * d.variables['%s_Altitude' % stub][...] # m to km
        los_wind[:,i] = d.variables['%s_Line-of-sight_Wind' % stub][...]
        los_wind_error[:,i] = d.variables['%s_Line-of-sight_Wind_Error' % stub][...]
        local_az[:,i] = d.variables['%s_Line-of-sight_Azimuth' % stub][...]
        amp[:,i] = d.variables['%s_Fringe_Amplitude' % stub][...]
        time_msec =         d.variables['Epoch'][...].item()
        time[i] = datetime(1970,1,1) + timedelta(seconds = 1e-3*time_msec)
        exp_time[i] = float(d.Time_Resolution.split(' ')[0])
        d.close()    
    
    L21_dict = {}
    L21_dict['lat'] = lat
    L21_dict['lon'] = lon
    L21_dict['alt'] = alt
    L21_dict['los_wind'] = los_wind
    L21_dict['los_wind_error'] = los_wind_error
    L21_dict['local_az'] = local_az
    L21_dict['amp'] = amp
    L21_dict['time'] = time
    L21_dict['exp_time'] = exp_time
    L21_dict['emission_color'] = emission_color
    L21_dict['source_files'] = L21_fns
    
    return L21_dict




    
    
    
############################################################################################################
##########################################       Level 2.2       ###########################################
############################################################################################################

    
   

    
def level21_dict_to_level22_dict(L21_A_dict, L21_B_dict):
    '''
    Given Level 2.1 data from MIGHTI A and MIGHTI B, process it with the Level 2.2 algorithm. 
    This entails interpolating line-of-sight wind measurements from MIGHTI A and B to a 
    common grid, and rotating to a geographic coordinate system to derive vector horizontal winds. 
    INPUTS:
        L21_A_dict  -- TYPE:dict.  The dictionary corresponding to three orbits of MIGHTI A measurements 
                                   for a single emission color. See "level21_to_dict" for required keys.
        L21_B_dict  -- TYPE:dict.  The dictionary corresponding to three orbits of MIGHTI B measurements
                                   for a single emission color, which is the same as for the A measurements.
                                   See "level21_to_dict" for required keys.
    OUTPUTS:
        L22_dict    -- TYPE:dict.  A dictionary containing the following variables. Most are given as 
                                   arrays with shape (ny,nx) where ny is the number of altitude grid
                                   points, and nx is the number of horizontal grid points.
                    lat             -- TYPE:array(ny,nx),    UNITS:deg.  Latitude of each point on the grid.
                    lon             -- TYPE:array(ny,nx),    UNITS:deg.  Longitude of each point on the grid.
                    alt             -- TYPE:array(ny,nx),    UNITS:km.   Altitude of each point on the grid.
                    u               -- TYPE:array(ny,nx),    UNITS:m/s.  Estimated zonal wind (positive eastward).
                    v               -- TYPE:array(ny,nx),    UNITS:m/s.  Estimated meridional wind (positive northward).
                    u_error         -- TYPE:array(ny,nx),    UNITS:m/s.  Uncertainty in u.
                    v_error         -- TYPE:array(ny,nx),    UNITS:m/s.  Uncertainty in v.
                    error_flags     -- TYPE:array(ny,nx,ne), UNITS:none. The error flags (either 0 or 1) for each point
                                                                         in the grid. Each point has ne=8 error flags, 
                                                                         which are set to 1 under the following 
                                                                         circumstances:
                                                                             0  = missing A file
                                                                             1  = missing B file
                                                                             2  = A did not sample this altitude
                                                                             3  = B did not sample this altitude
                                                                             4  = A sample exists but equals np.nan
                                                                             5  = B sample exists but equals np.nan
                                                                             6  = spherical asymmetry: A&B VER 
                                                                                  estimates disagree
                                                                             7  = unknown error
                    time            -- TYPE:array(ny,nt),    UNITS:none. The average between the time of the MIGHTI A 
                                                                         and B measurements that contribute to this 
                                                                         grid point, given as a datetime object. 
                    time_delta      -- TYPE:array(ny,nt),    UNITS:s.    The difference between the time of the MIGHTI
                                                                         A and B measurements that contribute to this 
                                                                         grid point.
                    amp_rel_diff    -- TYPE:array(ny,nt),    UNITS:none. The difference between the fringe amplitude 
                                                                         of the MIGHTI A and B measurements that 
                                                                         contribute to this grid point, divided by the
                                                                         mean. When this is high, it indicates that 
                                                                         spherical asymmetry may be a problem.
                    emission_color  -- TYPE:str,                         'red' or 'green'.        
                    source_files    -- TYPE:list of str,                 All the files used to create the data product,
                                                                         including the full paths.
                    
    TODO:
        - only save the parent files that were actually used.
        - explicitly pass in which files are previous, current, and next orbit
        - Will an orbit be -180 to 180 deg or 0 to 360? Who decides this? 
        - propagate uncertainties (easy).                                             
    '''

    
    lat_A      =       L21_A_dict['lat']
    lon_A      =       L21_A_dict['lon']
    alt_A      =       L21_A_dict['alt']
    los_wind_A =       L21_A_dict['los_wind']
    los_wind_A_err =   L21_A_dict['los_wind_error']
    local_az_A =       L21_A_dict['local_az']
    amp_A      =       L21_A_dict['amp']
    time_A     =       L21_A_dict['time']
    exp_time_A =       L21_A_dict['exp_time']
    emission_color =   L21_A_dict['emission_color']
    N_alts_A, N_times_A = np.shape(lon_A)

    lat_B      =       L21_B_dict['lat']
    lon_B      =       L21_B_dict['lon']
    alt_B      =       L21_B_dict['alt']
    los_wind_B =       L21_B_dict['los_wind']
    los_wind_B_err =   L21_B_dict['los_wind_error']
    local_az_B =       L21_B_dict['local_az']
    amp_B      =       L21_B_dict['amp']
    time_B     =       L21_B_dict['time']
    exp_time_B =       L21_B_dict['exp_time']
    emission_color =   L21_B_dict['emission_color']
    N_alts_B, N_times_B = np.shape(lon_B)
    

    ################### Parameters ###############################
    kcurr_A = N_times_A/2 # an index of the A files that is known to be in the "current" orbit
                          # (i.e., not "previous" or "next")
    kcurr_B = N_times_B/2 # same for B
    VER_REL_DIFF_THRESH = 0.19 # Relative difference in VER measured by A and B, beyond
                               # which the spherical-asymmetry flag will be raised.
                               # (Note that we really mean "fringe amplitude" which is
                               # not exactly proportional to VER due to the temperature
                               # dependence.)


    ############### Unwrap sample longitudes #####################
    # Eventually, this section will change when we finalize how the L2.1
    # files will be specified by Tori. A couple options are:
    #     1) Tori explicitly denotes which files are "previous", "current", and "next" orbit
    #     2) I look inside the file to determine the orbit number (assuming it is in there)
    # For now we define:
    #   - "current"  orbit is -180 to +180 deg
    #   - "previous" orbit is -540 to -180 deg
    #   - "next"     orbit is +180 to +540 deg
    # The inversion will be performed in the domain -180 to +180 deg.
    # This part of the code may need to be adjusted for use on real data.
    # It is currently defined like this so that Scott's simulation (which spans one
    # orbit from -180 to 180 deg). 
    for i in range(N_alts_A):
        for k in range(kcurr_A,N_times_A-1):
            if (lon_A[i,k+1] - lon_A[i,k]) < 0.0:
                lon_A[i,k+1:] += 360.
        for k in range(kcurr_A,0,-1):
            if (lon_A[i,k] - lon_A[i,k-1]) < 0.0:
                lon_A[i,:k] -= 360.
    for i in range(N_alts_B):
        for k in range(kcurr_B,N_times_B-1):
            if (lon_B[i,k+1] - lon_B[i,k]) < 0.0:
                lon_B[i,k+1:] += 360.
        for k in range(kcurr_B,0,-1):
            if (lon_B[i,k] - lon_B[i,k-1]) < 0.0:
                lon_B[i,:k] -= 360.

    ################# Define reconstruction grid: lon/alt ###################
    # Determine how finely to define longitude grid based on minimum L2.1 sampling rate
    lon_res_A = (np.diff(lon_A,axis=1)).min()
    lon_res_B = (np.diff(lon_B,axis=1)).min()
    lon_res = min(lon_res_A, lon_res_B)
    # Define longitude grid
    lon_vec = np.arange(-180.,180.,lon_res)
    # Define altitude grid similarly (predefine 90 to 300 km)
    alt_res_A = abs(np.diff(alt_A,axis=0)).min()
    alt_res_B = abs(np.diff(alt_B,axis=0)).min()
    alt_res = min(alt_res_A, alt_res_B)
    alt_vec = np.arange(90., 300., alt_res)
    # Define 2D reconstruction grid based on lon and alt
    lon,alt = np.meshgrid(lon_vec, alt_vec)
    N_alts, N_lons = np.shape(lon)

    ############### Interpolate values to reconstruction grid ##################
    # Use bilinear interpolation. This is somewhat complicated because
    # the sample grid is not exactly regular, but it's close, and 
    # we are approximating it as such. We're implementing our own
    # bilinear interpolation so we can control extrapolation in 
    # longitude and altitude as desired. Bilinear interpolation is 
    # used because it is insensitive to the units used to define
    # the sample grid.
    # This proceeds in 4 steps:
    # 1) Setup
    # 2) Error flagging
    # 3) Interpolation
    # 4) Inverting (i.e., rotating LoS winds to cardinal)
            
    # Output variables, which will be defined on the reconstruction grid
    U = np.nan*np.zeros(np.shape(lon))                # zonal wind
    V = np.nan*np.zeros(np.shape(lon))                # meridional wind
    U_err = np.nan*np.zeros(np.shape(lon))            # zonal wind uncertainty
    V_err = np.nan*np.zeros(np.shape(lon))            # meridional wind uncertainty
    lat = np.nan*np.zeros(np.shape(lon))              # latitude
    time = np.empty(np.shape(lon), dtype=object)      # time ascribed to L2.2 data point (as datetime objects)
    time_delta = np.nan*np.zeros(np.shape(lon))       # difference between A and B times used (seconds)
    ver_A = np.nan*np.zeros(np.shape(lon))            # fringe amplitude from A (related to VER)
    ver_B = np.nan*np.zeros(np.shape(lon))            # fringe amplitude from B (related to VER)
    ver_rel_diff = np.nan*np.zeros(np.shape(lon))     # relative difference in the above two
    error_flags = np.zeros((N_alts, N_lons, 8))       # Error flags, one set per grid point, defined as follows:
                                                      # 0  = missing A file
                                                      # 1  = missing B file
                                                      # 2  = A did not sample this altitude
                                                      # 3  = B did not sample this altitude
                                                      # 4  = A sample exists but equals np.nan
                                                      # 5  = B sample exists but equals np.nan
                                                      # 6  = Spherical asymmetry: A&B VER estimates disagree
                                                      # 7  = unknown error
    # Loop over the reconstruction altitudes
    for i in range(N_alts):   
        alt_pt = alt_vec[i]
        # Create a list of longitudes, one per A and B file, which have been
        # interpolated to this altitude.
        lon_list_A = np.zeros(N_times_A)
        lon_list_B = np.zeros(N_times_B)
        for k in range(N_times_A):
            lon_list_A[k] = interpolate_linear(alt_A[:,k], lon_A[:,k], alt_pt)
        for k in range(N_times_B):
            lon_list_B[k] = interpolate_linear(alt_B[:,k], lon_B[:,k], alt_pt)
        
        # Loop over the reconstruction longitudes
        for k in range(N_lons):
            lon_pt = lon_vec[k]
            
            # Find the file to the left and right in longitude. 
            kA0 = bisect.bisect(lon_list_A, lon_pt) - 1
            kA1 = kA0 + 1
            kB0 = bisect.bisect(lon_list_B, lon_pt) - 1
            kB1 = kB0 + 1
            
            
            
            ##################### Error Flagging ##########################
            # Never extrapolate in longitude. This error should not normally happen, and 
            # probably indicates an entire missing orbit or an extended calibration routine.
            # Mark as missing, and continue.
            if kA0 < 0 or kA1 >= N_times_A or kB0 < 0 or kB1 >= N_times_B:
                if kA0 < 0 or kA1 >= N_times_A:
                    error_flags[i,k,0] = 1
                if kB0 < 0 or kB1 >= N_times_B:
                    error_flags[i,k,1] = 1
                continue
            
            # Determine if there are "missing" files by checking the time between the straddling
            # files we just found and comparing to the exposure time of the files.
            # If so, throw error flag and continue.
            # Note that it's the exposure time of the first file that matters here.
            missing_A = (time_A[kA1] - time_A[kA0]).total_seconds() > 1.5*exp_time_A[kA0]
            missing_B = (time_B[kB1] - time_B[kB0]).total_seconds() > 1.5*exp_time_B[kB0]
            if missing_A or missing_B:
                if missing_A:
                    error_flags[i,k,0] = 1
                if missing_B:
                    error_flags[i,k,1] = 1
                continue
                    
            # If the desired altitude is outside the range of altitudes sampled by the 
            # instruments, throw error flag and continue.
            # For this, allow some wiggle room to handle case where MIGHTI A samples
            # at 90.01 km but we wanted 90.00 km.
            altmin_A = max(min(alt_A[:,kA0]), min(alt_A[:,kA1])) - alt_res
            altmin_B = max(min(alt_B[:,kB0]), min(alt_B[:,kB1])) - alt_res
            altmax_A = min(max(alt_A[:,kA0]), max(alt_A[:,kA1])) + alt_res
            altmax_B = min(max(alt_B[:,kB0]), max(alt_B[:,kB1])) + alt_res
            if alt_pt > min(altmax_A, altmax_B) or alt_pt < max(altmin_A, altmin_B):
                if alt_pt > altmax_A or alt_pt < altmin_A:
                    error_flags[i,k,2] = 1
                if alt_pt > altmax_B or alt_pt < altmin_B:
                    error_flags[i,k,3] = 1
                continue
            
            
            
            ######################## Interpolating ############################
            # If it passed all the error checks, perform bilinear interpolation (altitude, then longitude).
            # Variables to interpolate to this point:
            #   - los_wind (A and B)
            #   - az       (A and B)
            #   - lat      (A and B, to be averaged)
            #   - time     (A and B, to be averaged and subtracted)
            #   - ver      (A and B, to be compared)
            
            def bilinear_interp(lon_AB, alt_AB, val, prop_err = False, valerr = None):
                '''
                Helper function that will bilinearly interpolate the Nx2 array "val", sampled
                at the Nx2 array of points described by lon_AB and alt_AB, to the point 
                currently under consideration (i.e., lon_pt, alt_pt).
                
                Optional input to propagate error from original to interpolated value. If
                prop_err = True, valerr must be specified (as the error of each value of
                val), and the interpolated error will be provided an additional output.
                '''
                
                if prop_err and valerr is None:
                    raise Exception('If prop_err = True, then valerr must be specified')
                
                if not prop_err:
                    # Do interpolate of value to the desired altitude, for each longitude
                    val_0 = interpolate_linear(alt_AB[:,0], val[:,0], alt_pt)
                    val_1 = interpolate_linear(alt_AB[:,1], val[:,1], alt_pt)
                    # Interpolate the longitude coordinate to the desired altitude
                    lon_0 = interpolate_linear(alt_AB[:,0], lon_AB[:,0], alt_pt)
                    lon_1 = interpolate_linear(alt_AB[:,1], lon_AB[:,1], alt_pt)
                    # Do interpolation to the desired longitude
                    val_pt = interpolate_linear([lon_0, lon_1], [val_0, val_1], lon_pt,
                                                          extrapolation='none')
                    return val_pt
                
                else: # prop_err is True
                    # Do interpolate of value to the desired altitude, for each longitude
                    val_0, val_0_err = interpolate_linear(alt_AB[:,0], val[:,0], alt_pt, 
                                                          prop_err = True, yerr = valerr[:,0])
                    val_1, val_1_err = interpolate_linear(alt_AB[:,1], val[:,1], alt_pt, 
                                                          prop_err = True, yerr = valerr[:,1])
                    # Interpolate the longitude coordinate to the desired altitude
                    lon_0 = interpolate_linear(alt_AB[:,0], lon_AB[:,0], alt_pt)
                    lon_1 = interpolate_linear(alt_AB[:,1], lon_AB[:,1], alt_pt)
                    # Do interpolation to the desired longitude
                    val_pt, val_pt_err = interpolate_linear([lon_0, lon_1], [val_0, val_1], lon_pt,
                                                            extrapolation='none', 
                                                            prop_err = True, yerr = [val_0_err, val_1_err])
                    return val_pt, val_pt_err
                    
            
            los_wind_A_pt, los_wind_A_pt_err = \
                            bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], los_wind_A[:,kA0:kA1+1],
                                            prop_err = True, valerr = los_wind_A_err[:,kA0:kA1+1])
            los_wind_B_pt, los_wind_B_pt_err = \
                            bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], los_wind_B[:,kB0:kB1+1],
                                            prop_err = True, valerr = los_wind_B_err[:,kB0:kB1+1])
            local_az_A_pt = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], local_az_A[:,kA0:kA1+1])
            local_az_B_pt = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], local_az_B[:,kB0:kB1+1])
            lat_A_pt      = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], lat_A     [:,kA0:kA1+1])
            lat_B_pt      = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], lat_B     [:,kB0:kB1+1])
            ver_A_pt      = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], amp_A     [:,kA0:kA1+1])
            ver_B_pt      = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], amp_B     [:,kB0:kB1+1])
            
            # Interpolate time, which is more complicated because it's a datetime object
            t_A_0 = time_A[kA0]
            t_A_1 = time_A[kA1]
            t_B_0 = time_B[kB0]
            t_B_1 = time_B[kB1]
            # Interpolate the longitude coordinate to the desired altitude
            lon_A_0 = interpolate_linear(alt_A[:,kA0], lon_A[:,kA0], alt_pt)
            lon_A_1 = interpolate_linear(alt_A[:,kA1], lon_A[:,kA1], alt_pt)
            lon_B_0 = interpolate_linear(alt_B[:,kB0], lon_B[:,kB0], alt_pt)
            lon_B_1 = interpolate_linear(alt_B[:,kB1], lon_B[:,kB1], alt_pt)
            # Interpolate time gap, and add to time
            tgap_A = (t_A_1 - t_A_0).total_seconds()
            tgap_B = (t_B_1 - t_B_0).total_seconds()
            toff_A = interpolate_linear([lon_A_0, lon_A_1], [0, tgap_A], lon_pt, extrapolation='none')
            toff_B = interpolate_linear([lon_B_0, lon_B_1], [0, tgap_B], lon_pt, extrapolation='none')
            t_A = t_A_0 + timedelta(seconds=(toff_A))
            t_B = t_B_0 + timedelta(seconds=(toff_B))
            
            
            ############################ Inversion #############################
            # Coordinate transformation of winds from lines of sight to cardinal directions
            # Construct LoS winds in vector y
            y = np.array([-los_wind_A_pt, -los_wind_B_pt])
            # Coordinate transform (the heart of the L2.2 processing)
            azA_rad = np.deg2rad(local_az_A_pt)
            azB_rad = np.deg2rad(local_az_B_pt)
            A = np.array([[np.sin(azA_rad), np.cos(azA_rad)],
                          [np.sin(azB_rad), np.cos(azB_rad)]])
            invA = np.linalg.inv(A) # explicitly compute inverse
            x = invA.dot(y)
            u = x[0]
            v = x[1]
            # propagate uncertainties
            Sig_y = np.array([[los_wind_A_pt_err**2, 0.0],
                              [0.0, los_wind_B_pt_err**2]]) # covariance matrix of y
            Sig_x = invA.dot(Sig_y.dot(invA.T)) # standard linear error propagation
            u_err = np.sqrt(Sig_x[0,0])
            v_err = np.sqrt(Sig_x[1,1])

            
            ###################### Final error flagging #######################
            # Check spherical symmetry
            ver_rel_diff_pt = abs(ver_A_pt - ver_B_pt)/np.mean([ver_A_pt,ver_B_pt])
            if ver_rel_diff_pt > VER_REL_DIFF_THRESH:
                error_flags[i,k,6] = 1
        
            # Check if L2.1 data were nan
            if np.isnan(los_wind_A_pt):
                error_flags[i,k,4] = 1
            if np.isnan(los_wind_B_pt):
                error_flags[i,k,5] = 1
            if np.isnan(u) or np.isnan(v) and all(error_flags[i,k,:] == 0): # Unknown error
                error_flags[i,k,7] = 1
                
                
            # Fill in all the relevant variables at this grid point
            U[i,k] = u
            V[i,k] = v
            U_err[i,k] = u_err
            V_err[i,k] = v_err
            lat[i,k] = (lat_A_pt + lat_B_pt)/2
            time[i,k] = t_A + timedelta(seconds=(t_B-t_A).total_seconds()/2)
            time_delta[i,k] = (t_B-t_A).total_seconds()
            ver_A[i,k] = ver_A_pt
            ver_B[i,k] = ver_B_pt
            ver_rel_diff[i,k] = ver_rel_diff_pt
                
             
    # Create dictionary to be returned
    L22_dict = {}
    L22_dict['lat'] = lat
    L22_dict['lon'] = lon
    L22_dict['alt'] = alt
    L22_dict['u'] = U
    L22_dict['v'] = V
    L22_dict['u_error'] = U_err
    L22_dict['v_error'] = V_err
    L22_dict['error_flags'] = error_flags
    L22_dict['time'] = time
    L22_dict['time_delta'] = time_delta
    L22_dict['ver_rel_diff'] = ver_rel_diff
    L22_dict['emission_color'] = emission_color
    L22_dict['source_files'] = np.concatenate((L21_A_dict['source_files'], L21_B_dict['source_files']))
    
    return L22_dict





def save_nc_level22(path, L22_dict):
    '''
    Take the output of the Level 2.2 processing and save it as a NetCDF4 file in the official format.
    NetCDF4 file conventions taken from "Science Operations Center Data Product Conventions" Rev 0.4,
    authored by Tori Fae., including the new filenaming requirement that will presumably be in Rev 0.5.
    
    INPUTS:
        path        -- TYPE:str.  The directory the file will be saved in, including trailing "/"
                                  (e.g., '/home/user/')
        L22_dict    -- TYPE:dict. A dictionary containing output variables of the Level 2.2 processing.
                                  See documentation for level21_dict_to_level22_dict(...) for details.
    OUTPUTS:
        L22_fn      -- TYPE:str.  The full path to the saved file.
    TO-DO:
      - figure out how to make np.nan the fill value
      - Are attributes _FillValue and FillVal needed for variables that are strings? python-NetCDF4 doesn't seem to support this
      - Maybe: Fill in more notes for each variable
      - How will the 4 versions be specified?
      - Have a second set of eyes look at descriptions of each variable
      - Should dimensions be labeled the same as variables? Altitude/Vector/Epoch. Should Depend_0 point to vars or dims?
      - Parent files: Do we really want this? There are hundreds.
      - Time Resolution attribute? 30-60 seconds? Varied? N/A?
      - Add VER/fringe-amplitude estimate?
      - So far I haven't used Depend_0 or Depend_1 because there isn't an array for altitude or longitude. It's not a
        regular grid, in general. I'm not sure how NetCDF4/ISTP wants to handle that. 
    '''


    data_version = 1 # TODO: how will this be calculated? It goes into global attr Data_Version
    software_version = 0.01 # TODO: how will this be determined and stored? It goes into global attr Software_Version
    version_for_filename = 1 # TODO: How is this related to the two above? Rec'd major software version number
    revision_for_filename = 0 # TODO: How to determine this?
   
    
    #################### Compile variables to write in file ######################
    ### Timing:
    t_all = filter(None, L22_dict['time'].flatten()) # Extract all non-None grid times as a 1-D array
    t_start = min(t_all)
    t_end   = max(t_all)
    t_mid   = t_start + timedelta(seconds=(t_end - t_start).total_seconds()/2) # midpoint time
    t_start_msec = (t_start - datetime(1970,1,1)).total_seconds()*1e3 # milliseconds since epoch
    t_end_msec   = (t_end   - datetime(1970,1,1)).total_seconds()*1e3
    t_mid_msec   = (t_mid   - datetime(1970,1,1)).total_seconds()*1e3
    t_start_msec = np.int64(np.round(t_start_msec)) # cast to signed 64 bit integer
    t_end_msec   = np.int64(np.round(t_end_msec)) 
    t_mid_msec   = np.int64(np.round(t_mid_msec))
    t_file  = datetime.now()   # time this file was created  
    ### Who's running this process
    user_name = getpass.getuser()
    ### Parent files
    parents = '' # This will go in global attr Parents
    # TODO: parent files (there will be a lot...)
    for source_fn in L22_dict['source_files']:
        s = source_fn.split('/')[-1].split('.')
        pre = '.'.join(s[:-1])
        post = s[-1].upper()
        parents += '%s > %s, ' % (post, pre)
    if parents: parents = parents[:-2] # trim trailing comma


    ######################### Open file for writing ##############################
    L22_fn = 'ICON_L2_2_MIGHTI_%s_v%02i_r%02i.nc' % (t_mid.strftime('%Y-%m-%d_%H.%M.%S'),
                                                        version_for_filename, revision_for_filename)
    L22_full_fn = '%s%s'%(path, L22_fn)
    ncfile = netCDF4.Dataset(L22_full_fn,mode='w',format='NETCDF4') 
    
    try: # always close file if an error occurs
    
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
                                                "  * Neutral Wind profiles  Dr. Jon Makela\n",
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
        ncfile.Calibration_File =               '((TODO: zero wind cal file, and pass through L1 cal files?))'
        ncfile.Conventions =                    'SPDF ISTP/IACG Modified for NetCDF'
        ncfile.Data_Level =                     'L2.2'
        ncfile.Data_Type =                      'DP22 > Data Product 2.2: Cardinal Vector Winds'
        ncfile.Data_Version =                   np.uint16(data_version)
        ncfile.Date_End =                       t_end.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC' 
        ncfile.Date_Start =                     t_start.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC' 
        ncfile.Description =                    'ICON MIGHTI Cardinal Vector Winds (DP 2.2)'
        ncfile.Descriptor =                     'MIGHTI > Michelson Interferometer for Global High-resolution Thermospheric Imaging' 
        ncfile.Discipline =                     'Space Physics > Ionospheric Science'
        ncfile.File =                           L22_fn
        ncfile.File_Date =                      t_file.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC'
        ncfile.Generated_By =                   'ICON SDC > ICON UIUC MIGHTI L2.2 Processor v%.2f, B. J. Harding' % software_version
        ncfile.Generation_Date =                t_file.strftime('%Y%m%d')
        ncfile.History =                        'Version %i, %s, %s, ' % (data_version, user_name, t_file.strftime('%Y-%m-%dT%H:%M:%S')) +\
                                                'MIGHTI L2.2 Processor v%.2f ' % software_version +\
                                                '((TODO: Prepend previous history of this file, if it exists))'
        ncfile.HTTP_LINK =                      'http://icon.ssl.berkeley.edu/Instruments/MIGHTI'
        ncfile.Instrument =                     'MIGHTI'
        ncfile.Instrument_Type =                'Imagers (space)'
        ncfile.Link_Text =                      'MIGHTI Cardinal Vector Winds (DP 2.2)'
        ncfile.Link_Title =                     'ICON MIGHTI'
        ncfile.Logical_File_ID =                L22_fn[:-3]
        ncfile.Logical_Source =                 'ICON_L2_2_MIGHTI'
        ncfile.Logical_Source_Description =     'MIGHTI - Cardinal Vector Winds'
        ncfile.Mission_Group =                  'Ionospheric Investigations'
        ncfile.MODS =                           ncfile.History
        ncfile.Parents =                        parents
        ncfile.PI_Affiliation =                 'UC Berkeley > SSL'
        ncfile.PI_Name =                        'T. J. Immel'
        ncfile.Project =                        'NASA > ICON'
        ncfile.Rules_of_Use =                   'Public Data for Scientific Use'
        ncfile.Software_Version =               'ICON SDC > ICON UIUC MIGHTI L2.2 Processor v%.2f' % software_version
        ncfile.Source_Name =                    'ICON > Ionospheric Connection Explorer'
        ncfile.Spacecraft_ID =                  'NASA > ICON - 493'
        ncfile.Text =                           'ICON explores the boundary between Earth and space  the ionosphere  ' +\
                                                'to understand the physical connection between our world and the immediate '+\
                                                'space environment around us. Visit \'http://icon.ssl.berkeley.edu\' for more details.'
        ncfile.Text_Supplement =                '((TODO: Insert reference to Harding et al [2016] paper))'
        ncfile.Time_Resolution =                '30 or 60 seconds'
        ncfile.Title =                          'ICON MIGHTI Cardinal Vector Winds (DP 2.2)'


        ################################## Dimensions ########################################
        ny,nx = np.shape(L22_dict['alt'])
        ncfile.createDimension('Epoch',0)
        ncfile.createDimension('Altitude', ny)
        ncfile.createDimension('Longitude', nx)
        ncfile.createDimension('N_flags', np.shape(L22_dict['error_flags'])[2])
        


        ################################## Variables #########################################

        ######### Timing Variables #########

        # Time midpoint (the official required "Epoch" variable)
        # This is a little confusing since time is a dependent variable in our case, and the ISTP
        # format seems to want it to be the primary independent variable.
        t_msec = np.zeros((ny,nx),dtype=np.int64)
        t_fillval = np.int64(-1)
        for i in range(ny):
            for j in range(nx):
                if L22_dict['time'][i,j] is None:
                    t_msec[i,j] = t_fillval
                else:
                    t_msec[i,j] = np.int64(np.round((L22_dict['time'][i,j] - datetime(1970,1,1)).total_seconds()*1e3))
        var = _create_variable(ncfile, 'Epoch', t_msec, 
                              dimensions=('Altitude', 'Longitude'),
                              format_nc='i8', format_fortran='I', desc='Sample time, midpoint of A and B measurements. Number of msec since Jan 1, 1970.', 
                              display_type='scalar', field_name='Time', fill_value=t_fillval, label_axis='Time', bin_location=0.5,
                              units='ms', valid_min=0, valid_max=1000*365*86400e3, var_type='support_data', chunk_sizes=[1,1],
                              notes='')


        
        ######### Data Location and Direction Variables #########

        # Altitude
        val = L22_dict['alt']*1e3 # convert to meters
        var_alt = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Altitude', val, 
                              dimensions=('Altitude','Longitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='WGS84 altitude of each wind sample', 
                              display_type='image', field_name='Altitude', fill_value=None, label_axis='', bin_location=0.5,
                              units='m', valid_min=0, valid_max=1e10, var_type='support_data', chunk_sizes=[1,1],
                              notes='')

        
        # Longitude
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Longitude', L22_dict['lon'], 
                              dimensions=('Altitude','Longitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='WGS84 longitude of each wind sample', 
                              display_type='image', field_name='Longitude', fill_value=None, label_axis='', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[1,1],
                              notes='')
                              
        # Latitude
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Latitude', L22_dict['lat'], 
                              dimensions=('Altitude','Longitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='WGS84 latitude of each wind sample', 
                              display_type='image', field_name='Latitude', fill_value=None, label_axis='', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='support_data', chunk_sizes=[1,1],
                              notes='')
                              

        
        ######### Data Variables #########
        
        # Zonal Wind
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Zonal_Wind', L22_dict['u'], 
                              dimensions=('Altitude','Longitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Zonal component of the horizontal wind for an orbit. Positive Eastward.', 
                              display_type='image', field_name='Zonal Wind', fill_value=None, label_axis='', bin_location=0.5,
                              units='m/s', valid_min=-1e10, valid_max=1e10, var_type='data', chunk_sizes=[1,1],
                              notes='')
        
        # Meridional Wind
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Meridional_Wind', L22_dict['v'], 
                              dimensions=('Altitude','Longitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Meridional component of the horizontal wind for an orbit. Positive Northward.', 
                              display_type='image', field_name='Meridional Wind', fill_value=None, label_axis='', bin_location=0.5,
                              units='m/s', valid_min=-1e10, valid_max=1e10, var_type='data', chunk_sizes=[1,1],
                              notes='')    
                              
                              
        # Zonal Wind Error
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Zonal_Wind_Error', L22_dict['u_error'], 
                              dimensions=('Altitude','Longitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Error in the zonal wind estimate.', 
                              display_type='image', field_name='Zonal Wind Error', fill_value=None, label_axis='', bin_location=0.5,
                              units='m/s', valid_min=-1e10, valid_max=1e10, var_type='data', chunk_sizes=[1,1],
                              notes='')
        
        # Meridional Wind Error
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Meridional_Wind_Error', L22_dict['v_error'], 
                              dimensions=('Altitude','Longitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Error in the meridional wind estimate.', 
                              display_type='image', field_name='Meridional Wind Error', fill_value=None, label_axis='', bin_location=0.5,
                              units='m/s', valid_min=-1e10, valid_max=1e10, var_type='data', chunk_sizes=[1,1],
                              notes='')    


        ######### Other Metadata Variables #########
        
        # Emission color
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Emission_Color', L22_dict['emission_color'].capitalize(), 
                              dimensions=(),
                              format_nc=str, format_fortran='A', desc='Emission used for wind estimate: "Red" or "Green"', 
                              display_type='scalar', field_name='Emission Color', fill_value=None, label_axis='Color', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=1,
                              notes='')

        # Fringe amplitude relative difference
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Fringe_Amplitude_Relative_Difference', L22_dict['ver_rel_diff'], 
                              dimensions=('Altitude','Longitude'), # depend_0 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Difference in MIGHTI A&B\'s fringe amplitude estimates, divided by the mean', 
                              display_type='image', field_name='Fringe Amplitude Difference', fill_value=None, label_axis='', bin_location=0.5,
                              units='', valid_min=-1e10, valid_max=1e10, var_type='metadata', chunk_sizes=[1,1],
                              notes='This is the quantity used to determine if spherical asymmetry flag is raised')    
        # Error flags
        
        var = _create_variable(ncfile, 'ICON_L2_2_MIGHTI_Error_Flag', L22_dict['error_flags'], 
                              dimensions=('Altitude','Longitude','N_flags'), # depend_0 = 'Altitude',
                              format_nc='b', format_fortran='I', desc='Error flags. See Var_Notes attribute for description.', 
                              display_type='image', field_name='Error Flags', fill_value=None, label_axis='', bin_location=0.5,
                              units='', valid_min=0, valid_max=1, var_type='metadata', chunk_sizes=[1,1,1],
                              notes='Eight error flags for each grid point, each either 0 or 1:\n' +\
                                    '    0 = missing MIGHTI A file'+\
                                    '    1 = missing MIGHTI B file'+\
                                    '    2 = A did not sample this altitude'+\
                                    '    3 = B did not sample this altitude'+\
                                    '    4 = A sample exists but is NaN'+\
                                    '    5 = B sample exists but is NaN'+\
                                    '    6 = Spherical asymmetry: A&B VER estimates disagree'+\
                                    '    7 = Unknown Error')
        
        ncfile.close()
        
    except: # Make sure the file is closed
        ncfile.close()
        raise
            
    return L22_full_fn




def level21_to_level22(L21_path, L22_path):
    '''
    High-level function to apply the Level-2.1-to-Level-2.2 algorithm to a series of Level 2.1 files
    (in the L21_path directory) and create a Level 2.2 file (in the L22_path directory)
    
    INPUTS:
        L21_path    -- TYPE:str.  The directory containing the L2.1 files to be processed, including
                                  trailing "/" (e.g., '/home/user/'). This should contain a series of
                                  MIGHTI A and B files, over three orbits (previous, current, and next).
                                  
        L22_path    -- TYPE:str.  The directory the L2.2 file will be saved in, including trailing "/"
                                  (e.g., '/home/user/')    
    OUTPUTS:
        L22_fn      -- TYPE:str.  The full path to the saved L2.2 file.
    '''
    

    ##### Load L2.1 files into dictionaries
    Afns = glob.glob('%s*_MIGHTI-A_*.nc' % L21_path)
    Bfns = glob.glob('%s*_MIGHTI-B_*.nc' % L21_path)

    # Sort the files by time (same as alphanumeric sorting)
    Afns.sort()
    Bfns.sort()

    level21_A = level21_to_dict(Afns)
    level21_B = level21_to_dict(Bfns)

    ##### Run L2.2 processing to create L2.2 dictionary
    L22_dict = level21_dict_to_level22_dict(level21_A, level21_B)

    ##### Save L2.2 data to file
    L22_fn = save_nc_level22(L22_path, L22_dict)

    return L22_fn















    
    
