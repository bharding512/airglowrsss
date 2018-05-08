# A module for the conversion of MIGHTI Level 1 files to Level 2.1 and 2.2 files.
# Level 1 files - Calibrated MIGHTI interferograms
# Level 2.1 files - Line-of-sight wind profiles (this is where the onion-peeling inversion happens)
# Level 2.2 files - Vector wind profiles (this is where the A/B matchup happens)
# Altitudes and distances are expressed in km everywhere in the code, except when it's about to be
# saved in a netCDF file in meters.


####################################### VERSION CONTROL ############################################
# These need to be manually changed, when necessary.
# NOTE: When the major version is updated, you should change the History global attribute
# in both the L2.1 and L2.2 netcdf files, to describe the change (if that's still the convention)
software_version_major = 1 # Should only be incremented on major changes
software_version_minor = 8 # [0-99], increment on ALL published changes, resetting when the major version changes
__version__ = '%i.%02i' % (software_version_major, software_version_minor) # e.g., 2.03
####################################################################################################


############################## GLOBAL PROCESSING PARAMETERS ########################################
# Unless overridden by the user, the following parameters will be used:
global_params = {}

global_params['red'] = {
    'sigma'             : 1.0/630.0304e-9, # reciprocal of center wavelength of emission [m^-1] 
                                           # (from John Harlander)  
    'bin_size'          : 1,               # The number of rows of the interferogram to bin together to 
                                           # improve statistics at the cost of altitude resolution.   
    'account_for_local_projection': True,  # Whether to account for the fact that the line of sight is not
                                           # quite horizontal everywhere along the line of sight
    'integration_order' : 0,               # 0: Use Riemann-sum rule for discretizing line-of-sight integral
                                           # 1: Use trapezoidal rule for discretizing line-of-sight integral
    'top_layer'         : 'exp',           # 'thin': assume VER goes to zero above top layer
                                           # 'exp':  assume VER falls off exponentially in altitude with H=26km
    'sph_asym_thresh'   : 0.19,            # The relative difference in VER estimates from A&B,
                                           # beyond which the spherical asymmetry flag will be
                                           # raised in L2.2.
    'chi2_thresh'       : 0.5,             # [rad^2] Minimum value of chi^2, below which the line fit to extract phase 
                                           # should not be trusted. As defined, it is normalized by the number of pixels
                                           # in a row. If this threshold is being reached too often, then consider binning.
}

global_params['green'] = {
    'sigma'             : 1.0/557.7339e-9,
    'bin_size'          : 1,
    'account_for_local_projection': True,
    'integration_order' : 0,
    'top_layer'         : 'exp',
    'sph_asym_thresh'   : 0.19,
    'chi2_thresh'       : 0.5,
}


#####################################################################################################


import numpy as np
import ICON
import bisect
from scipy import integrate
from datetime import datetime, timedelta
import netCDF4
import getpass # for determining who is running the script
import glob
import traceback # for printing detailed error traces
import sys
import copy

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.colors import LogNorm

# Custom colormap for plots, and other aesthetics
from matplotlib.colors import LinearSegmentedColormap
cdict = {'red':   ((0.0, 0.7, 0.7),
                   (1.0, 0.7, 0.7)),
         'green': ((0.0, 0.7, 0.7),
                   (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 1.0, 1.0),
                   (1.0, 0.0, 0.0))}
cmap_custom = LinearSegmentedColormap('Custom', cdict)
matplotlib.rcParams['xtick.labelsize'] = 'small'
matplotlib.rcParams['ytick.labelsize'] = 'small'

# Ignore errors having to do with NaN. These clog up the log file.
np.seterr(invalid='ignore')


############################################################################################################
##########################################       Level 2.1       ###########################################
############################################################################################################




def phase_to_wind_factor(sigma_opd):
    '''
    Return the value f that satisfies w = f*p, where w is a wind change and p is a phase change.
    dphi = 2*pi*OPD*sigma*v/c (Eq 1 of Englert et al. 2007 Appl Opt., and Eq 2 of Harding et al. 2017 SSR)
                               
    INPUTS:
    
      *  sigma_opd   -- TYPE:float, UNITS:none.  sigma times opd: Optical path difference measured 
                                                 in wavelengths. If analyzing an entire row at once, 
                                                 the mean OPD of that row should be used.
                                                 
    OUTPUTS:
    
      *  f     -- TYPE:float, UNITS:m/s/rad. Phase to wind factor, described above.
      
    '''
    c      = 299792458.0 # m/s, speed of light
    return c / (2.*np.pi*sigma_opd)





def circular_mean(angle0,angle1):
    '''
    Find the mean angle, taking into account 0/360 crossover. For example,
    circular_mean(10,50) is 30, but circular_mean(350,20) is 5.
    
    INPUTS:
    
      *  angle0  -- TYPE:float or array, UNITS:deg. An angle in degrees.
      *  angle1  -- TYPE:float or array, UNITS:deg. An angle in degrees.
      
    OUTPUTS:
    
      *  angle   -- TYPE:float or array, UNITS:deg. The circular mean of the two input angles.
                   
    '''
    x = np.rad2deg(np.angle((np.exp(1j*np.deg2rad(angle0)) + np.exp(1j*np.deg2rad(angle1)))/2.))
    x = np.mod(x,360.)
    return x




def remove_satellite_velocity(I, sat_latlonalt, sat_velocity_vector, mighti_vectors, sigma_opd,):
    '''
    Modify the interferogram to remove the effect of satellite velocity upon the phase. 
    
    INPUTS:
    
      *  I                   -- TYPE:array(ny,nx),    UNITS:arb.  The MIGHTI interferogram. 
      *  sat_latlonalt       -- TYPE:array(3),        UNITS:(deg,deg,km). Satellite location in WGS84.
      *  sat_velocity_vector -- TYPE:array(3),        UNITS:m/s.  Satellite velocity vector in ECEF
      *  mighti_vectors      -- TYPE:array(ny,nx,3),  UNITS:none. mighti_vectors[i,j,:] is a unit 3-vector in ECEF
                                                                  coordinates defining the look direction of pixel (i,j).
      *  sigma_opd           -- TYPE:array(nx),       UNITS:none. The optical path difference (measured in wavelengths)
                                                                  for each column of the interferogram.
                                                                  
    OUTPUTS:
    
      *  I                   -- TYPE:array(ny,nx), UNITS:arb.  The MIGHTI interferogram, corrected
                                for the effects of satellite motion on the phase.
                                
    '''
    
    ny,nx = np.shape(I)
    
    # Loop over each pixel, calculating the look direction, projected satellite velocity
    # and the resulting phase shift in the interferogram
    sat_vel_phase = np.zeros((ny,nx))
    for j in range(nx): # Loop over columns
        phase_to_wind = phase_to_wind_factor(sigma_opd[j]) # use a different phase-to-wind factor for each column
        for i in range(ny): # loop over rows
            look_vector = mighti_vectors[i,j,:]
            proj_sat_vel = np.dot(sat_velocity_vector, look_vector) # positive apparent wind towards MIGHTI
            sat_vel_phase[i,j] = proj_sat_vel/phase_to_wind

    # Subtract phase from the interferogram
    I2 = I*np.exp(-1j*sat_vel_phase)
        
    return I2




def bin_array(b, y, lon = False):
    '''
    Downsample y by binning it, improving statistics. Every b
    elements of y will be averaged together to create a new array, y_b, 
    of length ny_b = ceil(len(y)/b). Binning starts at the end of the array, 
    so the first element of y_b may not represent exactly b samples of y.
    
    INPUTS:
    
      *  b    -- TYPE:int,          The number of rows to bin together
      *  y    -- TYPE:array(ny),    The array to be binned
      
    OPTIONAL INPUTS:
    
      *  lon  -- TYPE:bool,         If True, 360-deg discontinuities will
                                    be removed before averaging (e.g., for
                                    longitude binning).
                                    
    OUTPUTS:
    
      *  y_b  -- TYPE:array(ny_b),  The binned array
      
    '''
    # To save time, return quickly if b==1
    if b==1:
        return y
    
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
            
        if all(np.isnan(y_samps)): # If all the values are nan, return nan
            if y_b.dtype == complex:
                y_b[i_new] = np.nan + 1j*np.nan
            else:
                y_b[i_new] = np.nan
        else: # if there are non-nan values, take the nanmean
            y_b[i_new] = np.nanmean(y_samps)
        
    return y_b
    
    

def bin_uncertainty(b, ye):
    '''
    Determine the uncertainty of a binned array from the uncertainty
    of the un-binned array. Specifically:
    If the array y has uncertainty given by array ye, then the array
    ::
    
      y_b = bin_array(b, y)
        
    has uncertainty given by the array
    ::
    
      ye_b = bin_uncertainty(b, ye)
        
    INPUTS:
    
      *  b    -- TYPE:int,          The number of rows to bin together
      *  ye   -- TYPE:array(ny),    The uncertainty of the pre-binned data
      
    OUTPUTS:
    
      *  ye_b -- TYPE:array(ny_b), The uncertainty of the binned data
      
    '''
    # To save time, return quickly if b==1
    if b==1:
        return ye
    
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

        ye_b[i_new] = 1.0/sum(~np.isnan(ye_samps)) * np.sqrt(np.nansum(ye_samps**2))
        
    return ye_b
    
    
    
    
def bin_image(b, I):
    '''
    Downsample the interferogram in altitude to improve statistics while
    degrading vertical resolution. Every b rows will be averaged together. 
    Binning starts at high altitudes, so the lower rows of I_b may not represent 
    exactly b rows of I.
    
    INPUTS:
    
      *  b           -- TYPE:int,                        The number of rows to bin together
      *  I           -- TYPE:array(ny,nx),   UNITS:arb.  The MIGHTI interferogram
      
    OUTPUTS:
    
      *  I_b         -- TYPE:array(ny_b,nx), UNITS:arb.  The binned MIGHTI interferogram
      
    '''
    
    # To save time, return quickly if b==1
    if b==1:
        return I
    
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
    Define the matrix D whose inversion is known as "onion-peeling." 
    
    The forward model is:
    ::
    
        I = D * Ip
        
    where I is the measured interferogram, D is the observation matrix, and Ip is the 
    onion-peeled interferogram. If integration_order is 1, the observation matrix is 
    created by assuming the spectrum (and thus the interferogram) is a piecewise linear 
    function of altitude, treating the values of the interferogram at the tangent locations
    as the unknowns, and writing the measurements as a linear function of the unknowns.
    If integration_order is 0, the same recipe is followed, except the spectrum is 
    assumed to be a piecewise constant function of altitude, and the unknowns are the 
    values of the interferogram at the midpoint between two tangent altitudes.
    
    Setting integration_order==0 is better for precision.
    
    Setting integration_order==1 is better for accuracy.
    
    INPUTS:
    
      *  tang_alt   -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
      *  icon_alt   -- TYPE:float,        UNITS:km.   Altitude of the satellite.
      
    OPTIONAL INPUTS:
    
      *  top_layer  -- TYPE:str,          'thin': assume VER goes to zero above top layer
                                          'exp':  assume VER falls off exponentially in altitude (default)
      *  integration_order -- TYPE:int,
      
                      * 0: Use Riemann-sum rule for discretizing line-of-sight integral (default).
                      * 1: Use trapezoidal rule for discretizing line-of-sight integral.
                      
                                          
    OUTPUTS:
    
      *  D          -- TYPE:array(ny,ny), UNITS:km.   Observation matrix. Also called the "path matrix"
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
    
        theta = np.deg2rad(ICON.tang_alt_to_ze(tang_alt, icon_alt, RE))
        
        # Define grid. Bottom of each layer is defined by tangent height of observation.
        rbottom = tang_alt
        # Define top of each layer.
        rtop = rbottom.copy()
        rtop[:-1] = rbottom[1:]
        rtop[-1] = rbottom[-1] + (rtop[1]-rbottom[1])
        # Define midpt of each layer
        rmid = (rbottom + rtop)/2

        
        # Build observation matrix and other required matrix variables
        # Thanks to JJM for the vectorized formulation which improved runtime.
        th = np.tile(theta,(len(theta),1)).T
        rb = np.tile(rbottom,(len(rbottom),1))
        rt = np.tile(rtop,(len(rtop),1))
        sb2 = -np.sin(th)**2 + ((RE+rb)/(RE+icon_alt))**2
        st2 = -np.sin(th)**2 + ((RE+rt)/(RE+icon_alt))**2

        sb2[sb2<0] = 0. # there is no intersection of LOS with altitude rb. Set term to 0.
        st2[st2<0] = 0. # there is no intersection of LOS with altitude rt. Set term to 0.
        D = 2*(RE + icon_alt) * (np.sqrt(st2) - np.sqrt(sb2))
        
        if top_layer == 'exp': # Use exponential falloff model
            for m in range(M):
                rt = tang_alt[m] + RE
                r0 = tang_alt[-1] + RE
                def func(x, rt):
                    # The extrapolation function to be numerically integrated. (Eq 6 in Harding et al. 2017 SSR)
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
                    # The extrapolation function to be numerically integrated. (Eq 6 in Harding et al. 2017 SSR)
                    return np.exp(-1./H*(np.sqrt(x**2 + rt**2) - r0))
                
                x0 = np.sqrt(r0**2- rt**2)
                D[m,M-1] += 2.*integrate.quad(func, x0, np.inf, args=(rt))[0]
                
    else:
        raise Exception('"integration_order == %i" not supported. Use 0 or 1.' % integration_order)
    
    return D






def create_local_projection_matrix(tang_alt, icon_alt):
    '''
    Define the matrix B whose entries give the factor byF which a horizontal wind
    would be projected onto the line of sight. This has the same shape as the
    observation matrix (i.e., distance matrix). At the tangent point, this factor
    is 1.0. Far from the tangent point, this factor is smaller. If this effect is 
    accounted for, it makes a small change in the winds (less than 5 m/s). 
    Currently implemented as a simple unweighted average of the phase across the row.
    
    INPUTS:
    
      *  tang_alt   -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
      *  icon_alt   -- TYPE:float,        UNITS:km.   Altitude of the satellite.
      
    OUTPUTS:
    
      *  B          -- TYPE:array(ny,ny), UNITS:km.   Local projection matrix. B[i,j] = cos(angle 
                                                      between ray i and the tangent of shell j
                                                      at the point where they intersect)
                                                      
    '''
    
    # Assume the Earth is locally spherical with an effective radius RE.
    # (The estimated winds are barely sensitive to the choice of RE. This
    #  approximation introduces an error < 1mm/s)
    RE = 6371.
    theta = ICON.tang_alt_to_ze(tang_alt, icon_alt, RE)
    
    ny = len(tang_alt)
    
    # Calculate local-horizontal projection factors
    B = np.nan*np.zeros((ny,ny)) # matrix to hold cosine correction factors
    for i in range(ny):
        for j in range(i,ny): # only calculate upper triangular part
            th = theta[i]
            r = tang_alt[j]
            B[i,j] = (RE+icon_alt)/(RE+r) * np.sin(np.deg2rad(th))
    return B
     

    
    
    
def analyze_row(row):
    '''
    Given a 1-D interference pattern (i.e., a row of the complex intererogram), 
    analyze it to get a scalar phase value, which represents the wind, and a
    single amplitude value, which is roughly proportional to the emission rate.
    Also return normalized chi^2 for the fit used to determine the phase.
    
    INPUTS:
    
      *  row               -- TYPE:array(nx), UNITS:arb.   A row of the complex-valued, MIGHTI interferogram.
      
    OUTPUTS:
    
      *  phase             -- TYPE:float,     UNITS:rad.   A scalar phase which represents the observed wind.
      *  amp               -- TYPE:float,     UNITS:arb.   A scalar amplitude which represents the brightness.
      *  chi2              -- TYPE:float,     UNITS:rad^2. Normalized chi^2 for the line fit to unwrapped phase.
                                                           (normalized chi^2 = average value of (data - fit)**2)
       
    '''
    if all(np.isnan(row)):
        return np.nan, np.nan, np.nan
    
    row_phase = np.angle(row)
    
    tot_phase = np.nanmean(row_phase)
    
    resid = row_phase - tot_phase
    chi2 = np.nanmean(resid**2)
    
    # Evaluate total amplitude
    # Note: even though some rows are longer than others (near the bottom of the CCD),
    # they are apodized and thus the missing pixels wouldn't contribute much anyway.
    amp = np.nansum(abs(row))
        
    return tot_phase, amp, chi2

    
    
    
    
def perform_inversion(I, tang_alt, icon_alt, I_phase_uncertainty, I_amp_uncertainty,
                      top_layer='exp', integration_order=0, account_for_local_projection=True, chi2_thresh=np.inf,
                      linear_amp = True):
    '''
    Perform the onion-peeling inversion on the interferogram to return
    a new interferogram, whose rows refer to specific altitudes. In effect,
    this function undoes the integration along the line of sight.
    
    INPUTS:
    
      *  I           -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, MIGHTI interferogram.
      *  tang_alt    -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
      *  icon_alt    -- TYPE:float,        UNITS:km.   Altitude of the satellite.
      *  I_phase_uncertainty -- TYPE:array(ny), UNITS:rad. Uncertainty in the unwrapped, mean phase of each row of I.
                                                           This is provided in L1 file.
      *  I_amp_uncertainty   -- TYPE:array(ny), UNITS:arb. Uncertainty in the summed amplitude of each row of I.
                                                           This is provided in L1 file.
      
    OPTIONAL INPUTS:
    
      *  top_layer   -- TYPE:str,          'thin': assume VER goes to zero above top layer
                                           'exp':  assume VER falls off exponentially in altitude (default)
      *  integration_order -- TYPE:int,
      
                      * 0: Use Riemann-sum rule for discretizing line-of-sight integral (default).
                      * 1: Use trapezoidal rule for discretizing line-of-sight integral.
                      
      *  account_for_local_projection   -- TYPE:bool.   If False, a simple inversion is used.
                                           If True, the inversion accounts for the fact that the ray is not 
                                           perfectly tangent to each shell at each point along the ray. 
                                           (default True)
      *  chi2_thresh       -- TYPE:float. UNITS:rad^2. The chi^2 value for the line fit, below which the 
                                                       data will be ignored, since it is mostly noise.
                                                       Normalized by the number of samples. (default np.inf)
      *  linear_amp        -- TYPE:bool.   If True, a linear inversion is used on fringe amplitude. If False, the 
                                           fringe amplitude is estimated during the onion-peeling. These are the
                                           same in the absence of noise, but may make a difference at low SNR,
                                           especially for computing uncertainties. (default True)
                                           
                                                           
    OUTPUTS:
    
      *  Ip                -- TYPE:array(ny,nx), UNITS:arb. The complex-valued, onion-peeled interferogram.
      *  phase             -- TYPE:array(ny),    UNITS:rad. The mean phase of each row of Ip.
      *  amp               -- TYPE:array(ny),    UNITS:arb. The amplitude of each row of Ip.
      *  phase_uncertainty -- TYPE:array(ny),    UNITS:rad. The uncertainty of phase.
      *  amp_uncertainty   -- TYPE:array(ny),    UNITS:arb. The uncertainty of amp.
      
    '''
    
    if top_layer not in ['exp','thin']:
        raise ValueError('Argument top_layer=\'%s\' not recognized. Use \'exp\' or \'thin\'.' % top_layer)
    if integration_order not in [0,1]:
        raise ValueError('Argument integration_order=\'%s\' not recognized. Use 0 or 1')
          
    
    ny,nx = np.shape(I)
    
    # Create the path matrix
    D = create_observation_matrix(tang_alt, icon_alt, top_layer=top_layer, integration_order=integration_order)
    
    # Create local horizontal projection matrix (and set it to unity if we are to ignore this effect)
    B = np.ones((ny,ny))
    if account_for_local_projection:
        B = create_local_projection_matrix(tang_alt, icon_alt)
        
    
    ######### Onion-peeling inversion and phase extraction #########
    phase = np.zeros(ny) # phases at each altitude
    amp   = np.zeros(ny) # fringe amplitude at each altitude
    chi2  = np.zeros(ny) # chi^2 for line fit at each altitude
    Ip = np.zeros((ny,nx), dtype=complex) # onion-peeled interferogram

    # This code implements Eq (9) in the MIGHTI L2 Space Science Reviews
    # paper (Harding et al. 2017).

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
            # If the subtracted layer is nan, subtract nothing
            if not all(np.isnan(Ij_proj)):
                Li = Li - dij*Ij_proj
        # final normalization by this layer's path length
        Li = Li/dii
        Ip[i,:] = Li
        # Analyze the layer to get the phase, and store it.
        p,a,c = analyze_row(Ip[i,:])
        phase[i] = p
        amp[i] = a
        chi2[i] = c
            
        # Initial quality control
        # (Note: comprehensive quality control is done in a different function)
        if c > chi2_thresh:
            Ip[i,:] = 0.0 # Zeroed out to not contaminate other rows
            amp[i] = 0.0 # This is overly conservative. We probably trust
                         # amplitudes more than phase. Note that this will
                         # be overwritten if linear_amp=True
            
            
    if linear_amp: # Replace the onion-peeled amplitude with the linear inversion
        amp_L1 = np.zeros(ny) # fringe amplitude at each row of L1 interferogram
        for i in range(ny):
            _, a, _ = analyze_row(I[i,:])
            amp_L1[i] = a
        # Be careful with nans. If left unchecked they will contaminate other
        # altitudes. For purposes of inverting amplitude/VER, treat nan 
        # rows as 0 signal, but replace with nan at the end.
        inan = np.isnan(amp_L1)
        amp_L1[inan] = 0.0
        amp = np.linalg.solve(D, amp_L1) # fringe amplitude at each altitude
        amp[inan] = np.nan

    ######### Uncertainty propagation #########
    # Uncertainties can be propagated using simple linear inversion formula
    # (i.e., as if account_for_local_projection=False) to a very good approximation
    # (less than 1% error).
    
    ### Step 0: Characterize L1 and L2.1 interferograms with a single amp/phase per row
    ph_L1 = np.zeros(ny)
    A_L1 = np.zeros(ny)
    for i in range(ny):
        p,a,_ = analyze_row(I[i,:])
        ph_L1[i] = p
        A_L1[i] = a
    ph_L2 = phase.copy() # this was calculated above
    A_L2 = amp.copy() # this was calculated above

    ### Step 0.5: Handle nans and bad rows
    # For purposes of computing uncertainty, nans are like signal that we know is zero
    # to good precision.
    ph_L1[np.isnan(ph_L1)] = 0.0
    ph_L2[np.isnan(ph_L2)] = 0.0
    A_L1[np.isnan(A_L1)] = 0.0
    A_L2[np.isnan(A_L2)] = 0.0
    # If amp is exactly zero, then replace it with a small number
    # so that uncertainties can be calculated.
    A_L1[A_L1==0.0] = 1e-6
    A_L2[A_L2==0.0] = 1e-6
    
    ### Step 1: Transform amp/phase uncertainties to real/imag uncertainties
    # Each row will have a 2x2 covariance matrix describing the real and imaginary parts
    cov_real_imag_L1 = np.zeros((ny,2,2))
    for m in range(ny):
        # Jacobian of transformation from ampl/phase to real/imag.
        J = np.array([[np.cos(ph_L1[m]), -A_L1[m]*np.sin(ph_L1[m])],
                      [np.sin(ph_L1[m]),  A_L1[m]*np.cos(ph_L1[m])]])
        cov_amp_phase = np.diag([I_amp_uncertainty[m], I_phase_uncertainty[m]])**2 # assuming uncorrelated
        cov_amp_phase[np.isnan(cov_amp_phase)] = 0.0 # so the few nan *columns* don't corrupt the entire image
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
    
    if linear_amp: # Use linear propagation formula to overwrite
        DD = Dinv.dot(Dinv.T)
        amp_uncertainty = I_amp_uncertainty * np.sqrt(np.diag(DD))
        
    # Places where the values are nan should have nan uncertainty
    # (It was temporarily set to zero above for the purposes of 
    # propagating uncertainty, but should now be corrected)
    amp_uncertainty[np.isnan(amp)] = np.nan
    phase_uncertainty[np.isnan(phase)] = np.nan
            
    return Ip, phase, amp, phase_uncertainty, amp_uncertainty, chi2





def fix_longitudes(lons, lon_target):
    '''
    Unwrap the list of longitudes to avoid 360-deg jumps. The list will
    be fixed so that it contains a value within 180 deg of lon_target and
    is otherwise continuous.
    
    INPUTS:
    
      *  lons       -- TYPE:array, UNITS:deg. An ordered list of longitudes to be unwrapped.
      *  lon_target -- TYPE:float, UNITS:deg. See above.
      
    OUTPUTS:
    
      *  lons_new   -- TYPE:array, UNITS:deg. An ordered list of longitudes with jumps removed.
      
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
    Determine the geographical location to which the measurement will be attributed. Depending
    on integration_order (see function create_observation_matrix), this will either return
    the tangent locations, or the midpoint between two adjacent tangent locations.

    INPUTS:
    
      *  tang_lat    -- TYPE:array(ny), UNITS:deg.   Tangent latitudes.
      *  tang_lon    -- TYPE:array(ny), UNITS:deg.   Tangent longitudes.
      *  tang_alt    -- TYPE:array(ny), UNITS:km.    Tangent altitudes.
      
    OPTIONAL INPUTS:
    
      *  integration_order -- TYPE:int   
      
                      * 0: Use Riemann-sum rule for discretizing line-of-sight integral (default).
                      * 1: Use trapezoidal rule for discretizing line-of-sight integral.
                      
    OUTPUTS:
    
      *  lat         -- TYPE:array(ny), UNITS:deg.   Measurement latitudes.
      *  lon         -- TYPE:array(ny), UNITS:deg.   Measurement longitudes.
      *  alt         -- TYPE:array(ny), UNITS:km.    Measurement altitudes.
      
    '''
    if integration_order not in [0,1]:
        raise ValueError('integration_order = "%s" not recognized. Use 0 or 1')
    
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
    
      *  sat_latlonalt -- TYPE:array(3),  UNITS:(deg,deg,km). Satellite location in WGS84.
      *  lat           -- TYPE:array(ny), UNITS:deg.          Measurement latitudes.
      *  lon           -- TYPE:array(ny), UNITS:deg.          Measurement longitudes.
      *  alt           -- TYPE:array(ny), UNITS:km.           Measurement altitudes.
      
    OUTPUTS:
    
      *  az            -- TYPE:array(ny), UNITS:deg.          Azimuth angle of line of sight
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





def interpolate_linear(x, y, x0, extrapolation='hold', prop_err = False, yerr = None):
    '''
    Linear interpolation of the function y = f(x) to the location x0.
    x and y are vectors comprising samples of this function. There is also
    an option to propagate errors to the interpolated value. This function is
    5 times faster than scipy.interpolate.interp1d, and allows for
    zero-order-hold extrapolation. If you are interpolating to many points, 
    then scipy.interpolate.interp1d is likely faster. 

    INPUTS:
    
      *  x     -- TYPE:array(n), UNITS:arb. Independent variable of samples of function.
      *  y     -- TYPE:array(n), UNITS:arb. Dependent variable of samples of function.
      *  x0    -- TYPE:float,    UNITS:arb. Independent variable of interpolation point.
      
    OPTIONAL INPUTS:
    
      *  extrapolation -- TYPE:str,        'hold': extrapolate by using values at end points (default)
                                           'none': do not extrapolate. Points will be np.nan
      *  prop_err      -- TYPE:bool,
      
                                      * True:  propagate errors from original to interpolated
                                               value, and return an extra output; yerr must
                                               be specified as an input. 
                                      * False: do not propagate errors, and return only one
                                               output (default).
                                               
      *  yerr          -- TYPE:array(n), UNITS:arb. Error in y, to be propagated to interpolated value.
      
    OUTPUTS:
    
      *  y0    -- TYPE:float,    UNITS:arb. Interpolated value.
      
    OPTIONAL OUTPUT (if prop_err = True):
    
      *  y0err -- TYPE:float,    UNTIS:arb. Propagated error of y0.
      
    '''
    
    if prop_err and yerr is None:
        raise Exception('If prop_err=True, then yerr must be specified')    
        
    # Special corner case: x0 is exactly on the last grid point
    if x0==x[-1]:
        if prop_err:
            return y[-1], yerr[-1]
        else:
            return y[-1]
    
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





def level1_to_dict(L1_fn, emission_color, startstop = True):
    '''
    Read a level 1 file and translate it into a dictionary that the 
    level 2.1 processing can use.
    
    INPUTS:
    
      *  L1_fn          -- TYPE:str.  The full path and filename of the level 1 file.
      *  emission_color -- TYPE:str, 'green' or 'red'.
      
    OPTIONAL INPUTS:
      *  startstop      -- TYPE:bool, If False, read the midpoint value and write it to
                                      both the "start" and "stop" variables. (default True)
        
    OUTPUTS:
    
      *  L1_dict -- TYPE:dict. A dictionary containing information needed for
                               the level 2.1 processing. See documentation for 
                               level1_dict_to_level21(...) for required keys.
                               

    '''
    
    f = netCDF4.Dataset(L1_fn)
    
    # Is this A or B? There's no variable that says it (yet?) so we have to infer it from the file name
    sensor = None
    if 'MIGHTI-A' in L1_fn:
        sensor = 'A'
    elif 'MIGHTI-B' in L1_fn:
        sensor = 'B'
    else:
        raise Exception('Cannot determine sensor (A or B) from %s' % L1_fn)
        
    istop = 2  # index of stop in start/mid/stop
    istart = 0 # index of start in start/mid/stop
    if not startstop:
        istop = 1
        istart = 1
    
    L1_dict = {}
    L1_dict['L1_fn']                       = L1_fn
    L1_dict['sensor']                      = sensor
    L1_dict['I_amp']                       = f['ICON_L1_MIGHTI_%s_%s_Envelope' % (sensor, emission_color.capitalize())][0,:,:]
    L1_dict['I_phase']                     = f['ICON_L1_MIGHTI_%s_%s_Phase' % (sensor, emission_color.capitalize())][0,:,:]
    L1_dict['I_amp_uncertainty']           = f['ICON_L1_MIGHTI_%s_%s_Envelope_Uncertainties' % (sensor, emission_color.capitalize())][0,:]
    L1_dict['I_phase_uncertainty']         = f['ICON_L1_MIGHTI_%s_%s_Phase_Uncertainties' % (sensor, emission_color.capitalize())][0,:]
    # For tangent locations, only use the center, not the full horizontal distribution
    ny,nx = np.shape(L1_dict['I_amp'])
    tang_lla                               = f['ICON_L1_MIGHTI_%s_%s_Tangent_LatLonAlt' % (sensor, emission_color.capitalize())][0,:,:,:] 
    L1_dict['tang_alt_start']              = tang_lla[istart,2,:]
    L1_dict['tang_alt_stop']               = tang_lla[istop ,2,:]
    L1_dict['tang_lat_start']              = tang_lla[istart,0,:]
    L1_dict['tang_lat_stop']               = tang_lla[istop ,0,:]
    L1_dict['tang_lon_start']              = tang_lla[istart,1,:]
    L1_dict['tang_lon_stop']               = tang_lla[istop ,1,:]
    L1_dict['emission_color']              = emission_color
    # In the L1 file, the ECEF vectors are stored in multidimensional array: (time, start/mid/stop, vector_xyz, vert, horz)
    tmp                                    = f['ICON_L1_MIGHTI_%s_%s_ECEF_Unit_Vectors'% (sensor, emission_color.capitalize())][0,:,:,:]
    L1_dict['mighti_ecef_vectors']         = np.transpose(tmp, (1,2,0)) # V x H x vector
    L1_dict['icon_velocity_vector']        = f['ICON_L1_MIGHTI_%s_SC_Velocity_ECEF'% sensor][0,1,:] # at middle of exposure
    L1_dict['source_files']                = [f.Parents]
    L1_dict['acknowledgement']             = f.Acknowledgement
    tsec_start                             = f['ICON_L1_MIGHTI_%s_Image_Times'% sensor][0,istart]*1e-3
    tsec_stop                              = f['ICON_L1_MIGHTI_%s_Image_Times'% sensor][0,istop ]*1e-3
    L1_dict['time_start']                  = datetime(1970,1,1) + timedelta(seconds=tsec_start)
    L1_dict['time_stop']                   = datetime(1970,1,1) + timedelta(seconds=tsec_stop)
    L1_dict['exp_time']                    = tsec_stop - tsec_start # Should this be replaced with the actual exposure time?
    L1_dict['optical_path_difference']     = f['ICON_L1_MIGHTI_%s_%s_Array_OPD' % (sensor, emission_color.capitalize())][0,:].astype(float)*1e-2 # convert to m
    icon_ecef                              = f['ICON_L1_MIGHTI_%s_SC_Position_ECEF'% sensor][:][0,:,:] # timetable x [x,y,z]
    icon_latlonalt = np.zeros((3,3))
    for i in range(3):
        icon_latlonalt[i,:] = ICON.ecef_to_wgs84(icon_ecef[i,:])
    L1_dict['icon_alt_start'] = icon_latlonalt[istart,2]
    L1_dict['icon_alt_stop']  = icon_latlonalt[istop ,2]
    L1_dict['icon_lat_start'] = icon_latlonalt[istart,0]
    L1_dict['icon_lat_stop']  = icon_latlonalt[istop ,0]
    L1_dict['icon_lon_start'] = icon_latlonalt[istart,1]
    L1_dict['icon_lon_stop']  = icon_latlonalt[istop ,1]
    col = f.getncattr('Reference_Pixel_%s' % emission_color.capitalize())
    if not hasattr(col,'__len__'): # Only a single value was given. Expand to all rows
        col = col*np.ones(ny)
    
    # Dummy placeholder code for reading global attributes, if that matters
    nc_attrs = f.ncattrs()
    
    
    f.close()
    
    return L1_dict



    


def level1_dict_to_level21_dict(L1_dict, sigma = None, top_layer = None, 
                                integration_order = None, account_for_local_projection = None, 
                                bin_size = None, chi2_thresh = None, linear_amp = True):
    '''
    High-level function to run the Level 2.1 processing. It takes a dictionary (containing
    input variables extracted from a Level 1 file) and outputs a dictionary (containing 
    output variables, which can be written to a file using save_nc_level21). Basic quality
    control is performed.
    
    INPUTS:
    
      *  L1_dict       -- TYPE:dict.  A dictionary containing variables needed for
                                      the level 2.1 processing:
                                             
                                      * L1_fn                      -- TYPE:str.                      
                                                                      Level 1 filename
                                      * sensor                     -- TYPE:str.
                                                                      Which sensor took the data: 'A' or 'B'
                                      * I_amp                      -- TYPE:array(ny,nx), UNITS:arb.  
                                                                      Magnitude of interferogram
                                      * I_phase                    -- TYPE:array(ny,nx), UNITS:rad.  
                                                                      Phase of interferogram
                                      * I_amp_uncertainty          -- TYPE:array(ny),    UNITS:arb.  
                                                                      Uncertainty in the sum of each row of I_amp
                                      * I_phase_uncertainty        -- TYPE:array(ny),    UNITS:rad. 
                                                                      Uncertainty in the mean phase of each row of interferogram
                                      * tang_alt_start             -- TYPE:array(ny),    UNITS:km.   
                                                                      Tangent altitudes at beginning of exposure
                                      * tang_alt_stop              -- TYPE:array(ny),    UNITS:km.   
                                                                      Tangent altitudes at end of exposure
                                      * tang_lat_start             -- TYPE:array(ny),    UNITS:deg.  
                                                                      Tangent latitudes at beginning of exposure
                                      * tang_lat_stop              -- TYPE:array(ny),    UNITS:deg.  
                                                                      Tangent latitudes at end of exposure
                                      * tang_lon_start             -- TYPE:array(ny),    UNITS:deg.  
                                                                      Tangent longitudes at beginning of exposure
                                      * tang_lon_stop              -- TYPE:array(ny),    UNITS:deg.  
                                                                      Tangent longitudes at end of exposure
                                      * emission_color             -- TYPE:str.                    
                                                                      'red' or 'green'
                                      * icon_alt_start             -- TYPE:float,        UNITS:km.   
                                                                      Spacecraft altitude at beginning of exposure
                                      * icon_alt_stop              -- TYPE:float,        UNITS:km.   
                                                                      Spacecraft altitude at end of exposure
                                      * icon_lat_start             -- TYPE:float,        UNITS:deg.  
                                                                      Spacecraft latitude at beginning of exposure
                                      * icon_lat_stop              -- TYPE:float,        UNITS:deg.  
                                                                      Spacecraft latitude at end of exposure
                                      * icon_lon_start             -- TYPE:float,        UNITS:deg.  
                                                                      Spacecraft longitude at beginning of exposure
                                      * icon_lon_stop              -- TYPE:float,        UNITS:deg.  
                                                                       Spacecraft longitude at end of exposure
                                      * mighti_ecef_vectors        -- TYPE:array(ny,nx,3).           
                                                                      Unit ECEF vector of line of sight of each pixel at middle of exposure
                                      * icon_velocity_vector       -- TYPE:array(3),     UNITS:m/s.
                                                                      Spacecraft velocity vector in ECEF coordinates at middle of exposure
                                      * source_files               -- TYPE:list of strs.             
                                                                      All files that were used to generate this L1 file
                                      * time_start                 -- TYPE:datetime (timezone naive).                
                                                                      Start of exposure in UTC
                                      * time_stop                  -- TYPE:datetime (timezone naive).                 
                                                                      End of exposure in UTC
                                      * exp_time                   -- TYPE:float.        UNITS:s.    
                                                                      Length of exposure
                                      * optical_path_difference    -- TYPE:array(nx).    UNITS:m.    
                                                                      Optical path difference for each column of interferogram
                                      * acknowledgements           -- TYPE:str.
                                                                      The value in the global attribute "Acknowledgement" which will be propagated
                                             
    OPTIONAL INPUTS - If None, defaults from MIGHTI_L2.global_params will be used 
    
      *  sigma               -- TYPE:float, UNITS:m^-1. The wavenumber of the emission (1/wavelength)
      *  top_layer           -- TYPE:str, 'thin': assume VER goes to zero above top layer
                                          'exp':  assume VER falls off exponentially in altitude
      *  integration_order   -- TYPE:int, 0: Use Riemann-sum rule for discretizing line-of-sight integral
                                          1: Use trapezoidal rule for discretizing line-of-sight integral
      *  account_for_local_projection -- TYPE:bool. If False, a simple inversion is used.
                                         If True, the inversion accounts for the fact that the ray is not 
                                         perfectly tangent to each shell at each point along the ray.
      *  bin_size            -- TYPE:int, The number of rows of the interferogram to bin together to 
                                          improve statistics at the cost of altitude resolution.
      *  chi2_thresh         -- TYPE:float, UNITS:rad^2. The chi^2 value for the line fit (used to analyze
                                                         the inteferogram), below which the data will be 
                                                         ignored, since it is probably just noise.
      *  linear_amp          -- TYPE:bool. If True, a linear inversion is used on fringe amplitude. If False, the 
                                           fringe amplitude is estimated during the onion-peeling. These are the
                                           same in the absence of noise, but different with noise, especially
                                           for computing uncertainties. (default True)
                                                                          
    OUTPUTS:
    
      *  L21_dict            -- TYPE:dict. A dictionary containing output variables of the Level 2.1 processing:

                                * los_wind                  -- TYPE:array(ny),   UNITS:m/s.   Line-of-sight wind profile 
                                * los_wind_error            -- TYPE:array(ny),   UNITS:m/s.   Uncertainty of los_wind (1-sigma)
                                * lat                       -- TYPE:array(ny),   UNITS:deg.   Latitude of each point in profile
                                * lon                       -- TYPE:array(ny),   UNITS:deg.   Longitude of each point in profile
                                * alt                       -- TYPE:array(ny),   UNITS:alt.   Altitude of each point in profile
                                * time_start                -- TYPE:datetime (timezone naive) Time at start of exposure in UTC
                                * time_stop                 -- TYPE:datetime (timezone naive) Time at end of exposure in UTC
                                * exp_time                  -- TYPE:float,       UNITS:s.     Exposure time
                                * az                        -- TYPE:array(ny),   UNITS:deg.   The azimuth angle of the line of sight
                                                                                              at the tangent point (deg East of North)
                                * emission_color            -- TYPE:str.                      'red' or 'green'
                                * sensor                    -- TYPE:str.                      'A' or 'B'
                                * resolution_along_track    -- TYPE:array(ny),   UNITS:km.    Horizontal resolution along the line of sight
                                * resolution_cross_track    -- TYPE:array(ny),   UNITS:km.    Horizontal resolution perpendicular to line of sight
                                * resolution_alt            -- TYPE:array(ny),   UNITS:km.    Vertical resolution
                                * icon_alt                  -- TYPE:float,       UNITS:km.    Spacecraft altitude
                                * icon_lat                  -- TYPE:float,       UNITS:deg.   Spacecraft latitude
                                * icon_lon                  -- TYPE:float,       UNITS:deg.   Spacecraft longitude [0,360]
                                * fringe_amplitude          -- TYPE:array(ny),   UNITS:arb.   The fringe contrast, a proxy for volume emission rate
                                * fringe_amplitude_error    -- TYPE:array(ny),   UNITS:arb.   Uncertainty in fringe_amplitude (1-sigma)
                                * mighti_ecef_vectors       -- TYPE:array(ny,3).              ECEF unit vector for each line of sight
                                * icon_velocity_ecef_vector -- TYPE:array(3).    UNITS:m/s.   ECEF vector of spacecraft velocity
                                * file_creation_time        -- TYPE:datetime (timezone naive) Time this processing was run in UTC
                                * source_files              -- TYPE:list of str.              All science files that went into creating this file
                                * bin_size                  -- TYPE:int.                      Bin size used in the processing
                                * top_layer                 -- TYPE:str.                      How the top layer was handled: 'thin' or 'exp'
                                * integration_order         -- TYPE:int.                      Order of integration used in inversion: 0 or 1
                                * I                         -- TYPE:array(ny,nx) UNITS:arb.   The complex-valued, onion-peeled interferogram
                                * chi2                      -- TYPE:array(ny)    UNITS:rad^2. The normalized chi^2 from the phase
                                                                                              extraction from each row of the interferogram.
                                * quality                   -- TYPE:array(ny)    UNITS:arb.   A number from 0 (Bad) to 1 (Good) for each altitude.
                                * acknowledgement           -- TYPE:str.                      A copy of the Acknowledgement attribute in the L1 file.
    
    '''
    
    #### Parse input parameters and load defaults
    emission_color = L1_dict['emission_color']
    params = global_params[emission_color]
    if sigma is None:
        sigma = params['sigma']
    if top_layer is None:
        top_layer = params['top_layer']
    if integration_order is None:
        integration_order = params['integration_order']
    if account_for_local_projection is None:
        account_for_local_projection = params['account_for_local_projection']
    if bin_size is None:
        bin_size = params['bin_size']
    bin_size = int(bin_size)
    if chi2_thresh is None:
        chi2_thresh = params['chi2_thresh']

    ####  Load parameters from input dictionary
    Iraw = L1_dict['I_amp']*np.exp(1j*L1_dict['I_phase'])
    I_amp_uncertainty = L1_dict['I_amp_uncertainty']
    I_phase_uncertainty = L1_dict['I_phase_uncertainty']
    source_files = L1_dict['source_files']
    exp_time = L1_dict['exp_time']
    L1_fn = L1_dict['L1_fn']
    opd = L1_dict['optical_path_difference']
    sigma_opd = sigma * opd # Optical path difference, in units of wavelengths
    sensor = L1_dict['sensor']
    mighti_ecef_vectors = L1_dict['mighti_ecef_vectors']
    icon_velocity_vector = L1_dict['icon_velocity_vector']
    
    # Load parameters which are averaged from start to stop of exposure.
    icon_alt = (L1_dict['icon_alt_start'] + L1_dict['icon_alt_stop'])/2
    icon_lat = (L1_dict['icon_lat_start'] + L1_dict['icon_lat_stop'])/2
    icon_lon = circular_mean(L1_dict['icon_lon_start'], L1_dict['icon_lon_stop'])
    tang_alt = (L1_dict['tang_alt_start'] + L1_dict['tang_alt_stop'])/2
    tang_lat = (L1_dict['tang_lat_start'] + L1_dict['tang_lat_stop'])/2
    tang_lon = circular_mean(L1_dict['tang_lon_start'], L1_dict['tang_lon_stop'])
    
    #### Remove Satellite Velocity
    icon_latlonalt = np.array([icon_lat, icon_lon, icon_alt])
    I = remove_satellite_velocity(Iraw, icon_latlonalt, icon_velocity_vector, mighti_ecef_vectors, sigma_opd)
                         
    #### Bin data: average nearby rows together
    I        = bin_image(bin_size, I)
    tang_lat = bin_array(bin_size, tang_lat)
    tang_lon = bin_array(bin_size, tang_lon, lon=True)
    tang_alt = bin_array(bin_size, tang_alt)
    ny, nx = np.shape(I)
    mighti_ecef_vectors_center = mighti_ecef_vectors[:,nx/2,:] # For reporting in output file, determine ecef vector at center of row
    mighti_ecef_vectors_center = bin_image(bin_size, mighti_ecef_vectors_center) # bin each component separately
    I_amp_uncertainty   = bin_uncertainty(bin_size, I_amp_uncertainty)
    I_phase_uncertainty = bin_uncertainty(bin_size, I_phase_uncertainty)
    
    #### Determine geographical locations of inverted wind
    lat, lon, alt = attribute_measurement_location(tang_lat, tang_lon, tang_alt,
                                                   integration_order=integration_order)
    
    #### Onion-peel interferogram
    Ip, phase, amp, phase_uncertainty, amp_uncertainty, chi2 = perform_inversion(I, tang_alt, icon_alt, 
                           I_phase_uncertainty, I_amp_uncertainty,
                           top_layer=top_layer, integration_order=integration_order,
                           account_for_local_projection=account_for_local_projection,
                           chi2_thresh=chi2_thresh, linear_amp=linear_amp)
    
    #### Transform from phase to wind
    f = phase_to_wind_factor(np.mean(sigma_opd)) # Use average OPD to analyze entire row
    v             = f * phase
    v_uncertainty = f * phase_uncertainty
    
    #### Quality control (TODO: this ought to be its own function)
    ## Chi^2 filtering: mask phase, but not amplitude
    idx_chi2 = chi2 > chi2_thresh
#    if sum(idx_chi2) > 0.0:
#        print 'Removing %i wind samples (chi^2 > %.3f) with the following characteristics:' % (sum(idx_chi2), chi2_thresh)
#        print '\t chi^2   wind [m/s]   wind_error [m/s]'
#        for i in np.where(idx_chi2)[0]:
#            print'\t %6.3f    %7.1f       %7.1f' % (chi2[i], v[i], v_uncertainty[i])
#        print '\n'
    v[idx_chi2] = np.nan
    v_uncertainty[idx_chi2] = np.nan
  
    #### Calculate azimuth angles at measurement locations
    az = los_az_angle(icon_latlonalt, lat, lon, alt)
    
    # Make a L2.1 dictionary
    L21_dict = {
             'los_wind'                     : v,
             'los_wind_error'               : v_uncertainty,
             'lat'                          : lat,
             'lon'                          : lon,
             'alt'                          : alt,
             'time_start'                   : L1_dict['time_start'],
             'time_stop'                    : L1_dict['time_stop'],
             'exp_time'                     : exp_time,
             'az'                           : az,
             'emission_color'               : emission_color,
             'sensor'                       : sensor,
             'resolution_along_track'       : np.nan, # TODO
             'resolution_cross_track'       : np.nan, # TODO
             'resolution_alt'               : np.nan, # TODO
             'icon_alt'                     : icon_alt,
             'icon_lat'                     : icon_lat,
             'icon_lon'                     : icon_lon,
             'fringe_amplitude'             : amp,
             'fringe_amplitude_error'       : amp_uncertainty,
             'mighti_ecef_vectors'          : mighti_ecef_vectors_center,
             'icon_velocity_ecef_vector'    : icon_velocity_vector,
             'file_creation_time'           : datetime.now(),
             'source_files'                 : [L1_fn],
             'bin_size'                     : bin_size,
             'top_layer'                    : top_layer,
             'integration_order'            : integration_order,
             'I'                            : Ip,
             'chi2'                         : chi2,
             'wind_quality'                 : np.nan, # TODO
             'acknowledgement'              : L1_dict['acknowledgement'],
    }
        
    return L21_dict
   
   



def _create_variable(ncfile, name, value, format_nc='f8', format_fortran='F', dimensions=(), zlib=True, complevel=6, 
                    shuffle=True,  depend_0=None, depend_1=None, depend_2=None, chunk_sizes=None, desc='', 
                    display_type='scalar', field_name='', fill_value=None,label_axis='', bin_location=0.5, 
                    time_base='FIXED: 1970 (POSIX)', time_scale='UTC', units='', valid_min=None, valid_max=None, 
                    notes='', var_type='data', monoton=None):
    '''
    A helper function to write a variable to a netCDF file.
    
    INPUTS:
    
      *  Self evident from the parameters above. Notes:
      
            * fill_value = None --> default fill values will be used, if they exist. See netCDF4.default_fillvals
            * display_type: for now, 'scalar', 'time_series', 'altitude_profile', or 'image' will be used
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
        raise Exception('field_name="%s" is too long (%i chars). Shorten to 30 characters:\n"%s"' % (field_name,len(field_name)))
    if len(label_axis) > 10:
        raise Exception('label_axis="%s" is too long (%i chars). Shorten to 10 characters.' % (label_axis,len(label_axis)))
    if monoton not in ['Increase', 'Decrease', None]:
        raise Exception('Input "monoton" must be either "Increase", "Decrease" or None.')
    
    # If fill value was not specified, use the default value, if it exists.
    # It will not exist for strings, for example, for which fill values
    # cannot be set. (TODO: is this right?)
    if fill_value is None and format_nc in netCDF4.default_fillvals.keys():
        fill_value = netCDF4.default_fillvals[format_nc]
    
    var = ncfile.createVariable(name, format_nc, dimensions=dimensions, zlib=zlib, complevel=complevel,
                               shuffle=shuffle, chunksizes=chunk_sizes, fill_value=fill_value)
    var.CatDesc            = desc
    var.Long_Name          = desc
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
    elif fill_value is None and format_nc == str: 
        # Special case for strings. Make sure to set FillVal even though _FillValue can't be set
        var.FillVal        = ''
        
    var.Format             = format_fortran
    var.LablAxis           = label_axis
    if monoton is not None:
        var.MonoTon        = monoton
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
    if valid_min is not None and valid_max is not None:
        var.Valid_Range    = [valid_min, valid_max]
    var.setncattr_string('Var_Notes', notes) # to allow for multi-strings
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

   
   
    
    
def save_nc_level21(path, L21_dict, data_revision=0):
    '''
    Take the output of the Level 2.1 processing and save it as a NetCDF4 file in the official format.
    NetCDF4 file conventions taken from "Science Operations Center Data Product Conventions" Rev 0.5.
    
    INPUTS:
    
      *  path        -- TYPE:str.  The directory the file will be saved in, including trailing "/"
                                   (e.g., '/home/user/')
      *  L21_dict    -- TYPE:dict. A dictionary containing output variables of the Level 2.1 processing.
                                   See documentation for level1_dict_to_level21_dict(...) for details.
                                   
    OPTIONAL INPUTS:
    
      *  data_revision       -- TYPE:int,  The minor version of the data [0-999]. The major version is set
                                           by the software's major version.
                                   
    OUTPUTS:
    
      *  L21_fn      -- TYPE:str.  The full path to the saved file.
      
    TO-DO:
    
      * How can we append to global attributes History and MODS when the processing is re-run?
      
    '''
    
    L21_dict = copy.deepcopy(L21_dict) # because netCDF4 seems to change the input when nans are involved
    
    data_version_major = software_version_major # enforced as per Data Product Conventions Document

    #################### Compile variables to write in file ######################
    ### Sensor:
    sensor = L21_dict['sensor']
    ### Timing:
    t_start = L21_dict['time_start']
    t_stop  = L21_dict['time_stop']
    t_mid   = t_start + timedelta(seconds=(t_stop - t_start).total_seconds()/2) # middle of exposure
    t_start_msec = (t_start - datetime(1970,1,1)).total_seconds()*1e3 # milliseconds since epoch
    t_stop_msec  = (t_stop  - datetime(1970,1,1)).total_seconds()*1e3
    t_mid_msec   = (t_mid   - datetime(1970,1,1)).total_seconds()*1e3
    t_start_msec = np.int64(np.round(t_start_msec)) # cast to signed 64 bit integer
    t_stop_msec  = np.int64(np.round(t_stop_msec)) 
    t_mid_msec   = np.int64(np.round(t_mid_msec))
    t_file  = datetime.now()   # time this file was created  
    ### Who's running this process
    user_name = getpass.getuser()
    ### Parent files
    parents = [] # This will go in global attr Parents
    for source_fn in L21_dict['source_files']:
        s = source_fn.split('/')[-1].split('.')
        pre = '.'.join(s[:-1])
        post = s[-1].upper()
        parents.append('%s > %s' % (post, pre))


    L21_fn = 'ICON_L2_MIGHTI-%s_Line-of-Sight-Wind-%s_%s_v%02ir%03i.NC' % (sensor,L21_dict['emission_color'].capitalize(),
                                                           t_mid.strftime('%Y-%m-%d_%H%M%S'),
                                                           data_version_major, data_revision)
    L21_full_fn = '%s%s'%(path, L21_fn)
    ncfile = netCDF4.Dataset(L21_full_fn,mode='w',format='NETCDF4') 

    try:
        ########################## Global Attributes #################################
        ncfile.setncattr_string('Acknowledgement',                L21_dict['acknowledgement'])
        ncfile.setncattr_string('ADID_Ref',                       'NASA Contract > NNG12FA45C')
        ncfile.setncattr_string('Calibration_File',               '')
        ncfile.setncattr_string('Conventions',                    'SPDF ISTP/IACG Modified for NetCDF')
        ncfile.setncattr_string('Data_Level',                     'L2.1')
        ncfile.setncattr_string('Data_Type',                      'DP21 > Data Product 2.1: Line-of-sight Wind Profile')
        ncfile.Data_Version_Major =                               np.uint16(data_version_major)
        ncfile.Data_Revision =                                    np.uint16(data_revision)
        ncfile.Data_Version =                                     data_version_major + 0.001 * data_revision
        ncfile.setncattr_string('Date_Stop',                      t_mid.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC') # single measurement: use midpoint
        ncfile.setncattr_string('Date_Start',                     t_mid.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC') # single measurement: use midpoint
        ncfile.setncattr_string('Description',                    'ICON MIGHTI Line-of-sight Winds (DP 2.1)')
        ncfile.setncattr_string('Descriptor',                     'MIGHTI-%s > Michelson Interferometer for Global High-resolution ' % sensor+\
                                                                  'Thermospheric Imaging, Sensor %s' % sensor)
        ncfile.setncattr_string('Discipline',                     'Space Physics > Ionospheric Science')
        ncfile.setncattr_string('File',                           L21_fn)
        ncfile.setncattr_string('File_Date',                      t_file.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Generated_By',                   'ICON SDC > ICON UIUC MIGHTI L2.1 Processor v%s, B. J. Harding' % __version__)
        ncfile.setncattr_string('Generation_Date',                t_file.strftime('%Y%m%d'))
        ncfile.setncattr_string('History',                        'v1.0: First release of MIGHTI L2.1/L2.2 software, B. J. Harding, 05 Mar 2018')
        ncfile.setncattr_string('HTTP_LINK',                      'http://icon.ssl.berkeley.edu/Instruments/MIGHTI')
        ncfile.setncattr_string('Instrument',                     'MIGHTI-%s' % sensor)
        ncfile.setncattr_string('Instrument_Type',                'Imagers (space)')
        ncfile.setncattr_string('Link_Text',                      'MIGHTI Line-of-sight Wind Profile (DP 2.1)')
        ncfile.setncattr_string('Link_Title',                     'ICON MIGHTI')
        ncfile.setncattr_string('Logical_File_ID',                L21_fn[:-3])
        ncfile.setncattr_string('Logical_Source',                 'ICON_L2_MIGHTI-%s_' % (sensor,))
        ncfile.setncattr_string('Logical_Source_Description',     'MIGHTI Sensor %s - Line-of-sight Wind Profile' % (sensor,))
        ncfile.setncattr_string('Mission_Group',                  'Ionospheric Investigations')
        ncfile.setncattr_string('MODS',                           ncfile.History)
        ncfile.setncattr_string('Parents',                        parents)
        ncfile.setncattr_string('PI_Affiliation',                 'UC Berkeley > SSL')
        ncfile.setncattr_string('PI_Name',                        'T. J. Immel')
        ncfile.setncattr_string('Project',                        'NASA > ICON')
        ncfile.setncattr_string('Rules_of_Use',                   'Public Data for Scientific Use')
        ncfile.setncattr_string('Software_Version',               'ICON SDC > ICON UIUC MIGHTI L2.1 Processor v%s' % __version__)
        ncfile.setncattr_string('Source_Name',                    'ICON > Ionospheric Connection Explorer')
        ncfile.setncattr_string('Spacecraft_ID',                  'NASA > ICON - 493')
        ncfile.setncattr_string('Text',                           'ICON explores the boundary between Earth and space - the ionosphere - '
                                                                  'to understand the physical connection between our world and the immediate '
                                                                  'space environment around us. Visit \'http://icon.ssl.berkeley.edu\' for more details.')
        ncfile.setncattr_string('Text_Supplement',                ["This data product contains an altitude profile of the line-of-sight winds for a single "
        "image taken by MIGHTI. There is one file for each sensor (A or B), for each color (red or green) and for each image (usually every 30 or "
        "60 seconds). The profile "
        "spans from ~90 km (for green) or ~150 km (for red) to ~300 km, though altitudes with low signal levels are masked out. This data product is "
        "generated from the Level 1 MIGHTI product, which comprises calibrated amplitudes and phases. The spacecraft velocity is removed from the "
        "interferogram phase, , then (optionally) the data are binned from their native altitude resolution (~2.5 km) to improve statistics. "
        "An onion-peeling inversion"
        "then an onion-peeling inversion is performed to remove the effect of the line-of-sight integration. After the inversion, each row (i.e., altitude) "
        "is analyzed to extract the phase, and thus the line-of-sight wind. Level 2.1 files from MIGHTI-A and MIGHTI-B are combined during the Level 2.2 "
        "processing (not discussed here). See Harding et al. [2017, doi:10.1007/s11214-017-0359-3] for more details of the inversion algorithm. One update "
        "to this paper is relevant: Zero wind removal is now performed prior to the creation of the Level 1 file, instead of during the L2.1 processing. "
        ])
        ncfile.setncattr_string('Time_Resolution',                '%.1f seconds' % L21_dict['exp_time'])
        ncfile.setncattr_string('Title',                          'ICON MIGHTI Line-of-sight Wind Profile (DP 2.1)')


        ################################## Dimensions ########################################
        prefix = 'ICON_L2_MIGHTI_%s_%s' % (sensor, L21_dict['emission_color'].capitalize()) # prefix of each variable,
                                                                                            # e.g., ICON_L2_MIGHTI_A_Red
        n = len(L21_dict['alt'])
        var_alt_name = '%s_Altitude'%prefix
        ncfile.createDimension('Epoch',1)
        ncfile.createDimension(var_alt_name, n)
        ncfile.createDimension('Vector',3)
        ncfile.createDimension('Start_Mid_Stop',3)


        ################################## Variables #########################################

   
        ######### Timing Variables #########

        # Time midpoint (the official required "Epoch" variable)
        var = _create_variable(ncfile, 'Epoch', t_mid_msec, 
                              dimensions=(),
                              format_nc='i8', format_fortran='I', desc='Sample time, midpoint of exposure. Number of msec since Jan 1, 1970.', 
                              display_type='scalar', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
                              units='ms', valid_min=0, valid_max=1000*365*86400*1000, var_type='support_data', chunk_sizes=1,
                              notes="This variable contains the time corresponding to the wind profile reported in this file. It is evaluated "
                                    "at the midpoint of the exposure time. It is in UTC and has units of milliseconds since "
                                    "Jan 1, 1970. A human-readable version of the time can be found in the global attributes "
                                    "Date_Start and Date_Stop."
                              )

        # Time start/mid/stop
        var = _create_variable(ncfile, '%s_Time'%prefix, np.array([t_start_msec, t_mid_msec, t_stop_msec]),
                              dimensions=('Start_Mid_Stop'),
                              format_nc='i8', format_fortran='I', desc='Sample time at start, mid, stop of exposure. Number of msec since Jan 1, 1970.', 
                              display_type='scalar', field_name='Time', fill_value=None, label_axis='Time', bin_location=[0.0,0.5,1.0],
                              units='ms', valid_min=0, valid_max=1000*365*86400*1000, var_type='support_data', chunk_sizes=[3],
                              notes="This variable is the same as Epoch, except it includes the start time and stop time. "
                                    "This variable has length 3: [start_time, mid_time, stop_time]."
                              )
        
        # Human readable version of midpoint time
        var = _create_variable(ncfile, '%s_UTC_Time'%prefix, t_mid.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3], 
                              dimensions=(),
                              format_nc=str, format_fortran='A', desc='Sample time, midpoint of exposure.', 
                              display_type='time_series', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='support_data', chunk_sizes=1,
                              notes="This variable is the same as Epoch but is formatted as a human-readable string.")
                               
            
        ######### Data Variables #########

        # Line-of-sight wind profile
        var = _create_variable(ncfile, '%s_Line_of_Sight_Wind'%prefix, L21_dict['los_wind'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Line-of-sight horizontal wind profile. A positive wind is towards MIGHTI.', 
                              display_type='altitude_profile', field_name='Line-of-sight Wind', fill_value=None, label_axis='LoS Wind', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[n],
                              notes="The wind is the primary data product in this file. This variable contains the projection of the horizontal wind "
                              "(at the tangent point) onto the line of sight. An entire altitude profile is observed simultaneously. An onion-peeling "
                              "inversion is used on the raw observations to remove the effects of the integration along the line of sight. The "
                              "line-of-sight wind is defined such that a positive value indicates motion towards the spacecraft. This direction is given "
                              "by the Line_of_Sight_Azimuth variable. It is assumed that the vertical wind is zero, but even large vertical winds "
                              "(~100 m/s) do not significantly affect this data product, since the line of sight is nearly horizontal everywhere. It "
                              "should be noted that while this measurement is ascribed to a particular latitude, longitude and altitude, it is actually "
                              "an average over many hundreds of kilometers horizontally, and 2.5-30 kilometers vertically (depending on the binning). "
                              "See Harding et al. [2017, doi:10.1007/s11214-017-0359-3] for a more complete discussion of the inversion algorithm."
                              )

        # Line-of-sight wind error profile
        var = _create_variable(ncfile, '%s_Line_of_Sight_Wind_Error'%prefix, L21_dict['los_wind_error'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Line-of-sight horizontal wind error profile', 
                              display_type='altitude_profile', field_name='Line-of-sight Wind Error', fill_value=None, label_axis='Wind Error', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[n],
                              notes="The statistical (1-sigma) error in the line-of-sight wind. This is usually dominated by shot noise, but "
                              "also includes the effects of dark and read noise, as well as calibrations errors (e.g., the zero wind calibration), "
                              "and spacecraft pointing error (which affects the uncertainty in removing the spacecraft velocity from the observed "
                              "velocity). Other systematic errors or biases may exist (e.g., the effect of gradients along the line of sight) "
                              "which are not included in this variable."
                              )
        
        # Quality code
        var = _create_variable(ncfile, '%s_Wind_Quality'%prefix, L21_dict['wind_quality'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='A quantification of the quality, from 0 (Bad) to 1 (Good)', 
                              display_type='altitude_profile', field_name='Wind Quality', fill_value=None, label_axis='Quality', bin_location=0.5,
                              units='arb', valid_min=0.0, valid_max=1.0, var_type='data', chunk_sizes=[n],
                              notes=["NOT YET IMPLEMENTED",
                                     "While the intent is that the variable Line_of_Sight_Wind_Error accurately characterizes the statistical "
                                     "error in the wind data, it is possible that systematic errors are present, or that the statistical error "
                                     "estimation is not accurate. If it is suspected that this is the case, the quality will be less than 1.0. If "
                                     "the data are definitely unusable, the the quality will be 0.0 and the sample will be masked. Users should "
                                     "exercise caution when the quality is less than 1.0."
                                     ])

        # Fringe amplitude profile (TODO: will this be replaced by a VER data product?)
        var = _create_variable(ncfile, '%s_Fringe_Amplitude'%prefix, L21_dict['fringe_amplitude'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Fringe amplitude profile', 
                              display_type='altitude_profile', field_name='Fringe Amplitude', fill_value=None, label_axis='Fringe Amp', bin_location=0.5,
                              units='arb', valid_min=-1e30, valid_max=1e30, var_type='data', chunk_sizes=[n],
                              notes="An approximate volume emission rate (VER) profile in arbitrary units. Technically this a profile of the "
                              "amplitude of the fringes, which has a dependence on thermospheric temperature and background emission. Thus, it does not "
                              "truly represent volume emission rate. However, it is a useful proxy. The units are arbitrary, but an attempt has "
                              "been made to cross-calibrate MIGHTI-A with MIGHTI-B. In contrast with the wind inversion, which is nonlinear due to the "
                              "phase extraction step, the amplitude inversion is purely linear. The Level 1 interferogram is analyzed to obtain a single "
                              "brightness value per observing angle, and this is inverted with the distance matrix to obtain a value of the amplitude "
                              "per altitude."
                              )

        # Fringe amplitude error profile (TODO: will this be replaced by a VER data product?)
        var = _create_variable(ncfile, '%s_Fringe_Amplitude_Error'%prefix, L21_dict['fringe_amplitude_error'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Fringe amplitude error profile', 
                              display_type='altitude_profile', field_name='Fringe Amplitude Error', fill_value=None, label_axis='Amp Err', bin_location=0.5,
                              units='arb', valid_min=0.0, valid_max=1e30, var_type='data', chunk_sizes=[n],
                              notes="The statistical (1-sigma) error in the fringe amplitude. As with the wind, systematic errors are not "
                              "included, but can arise from sources such as horizontal gradients and inaccurate calibration.")
        

        ######### Data Location and Direction Variables #########
        
        # Altitude
        val = L21_dict['alt']*1e3 # convert to meters
        var = _create_variable(ncfile, '%s_Altitude'%prefix, val, 
                              dimensions=(var_alt_name),  depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='WGS84 altitude of each wind sample', 
                              display_type='altitude_profile', field_name='Altitude', fill_value=None, label_axis='Altitude', bin_location=0.5,
                              units='m', valid_min=50e3, valid_max=1000e3, var_type='support_data', chunk_sizes=[n], monoton='Increase',
                              notes="The altitudes of each point in the wind profile, evaluated using the WGS84 ellipsoid. If the variable "
                                    "Integration_Order=0 (which is the default value), then these altitudes are one half sample above "
                                    "the tangent altitudes of each pixel's line of sight (consistent with the assumption implicit "
                                    "in the inversion that the wind and emission rate are constant within the layer between tangent "
                                    "altitudes). If Integration_Order=1, this variable contains the tangent altitudes."
                              )

        # Latitude
        var = _create_variable(ncfile, '%s_Latitude'%prefix, L21_dict['lat'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='WGS84 latitude of each wind sample', 
                              display_type='altitude_profile', field_name='Latitude', fill_value=None, label_axis='Latitude', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='support_data', chunk_sizes=[n],
                              notes="The latitudes of each point in the wind profile, evaluated using the WGS84 ellipsoid. The latitude only "
                              "varies by several degrees from the bottom of the profile to the top. It should be noted that while a single latitude "
                              "value (the tangent latitude) is given for each point, the observation is inherently a horizontal average "
                              "over many hundreds of kilometers."
                              )

        # Longitude
        var = _create_variable(ncfile, '%s_Longitude'%prefix, L21_dict['lon'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='WGS84 longitude of each wind sample', 
                              display_type='altitude_profile', field_name='Longitude', fill_value=None, label_axis='Longitude', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[n],
                              notes="The longitudes (0-360) of each point in the wind profile, evaluated using the WGS84 ellipsoid. The longitude only "
                              "varies by several degrees from the bottom of the profile to the top. It should be noted that while a single longitude "
                              "value (the tangent longitude) is given for each point, the observation is inherently a horizontal average "
                              "over many hundreds of kilometers."
                              )

        # Azimuth angle of line of sight
        var = _create_variable(ncfile, '%s_Line_of_Sight_Azimuth'%prefix, L21_dict['az'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Azimuth angle of the line of sight at the tangent point. Deg East of North.', 
                              display_type='altitude_profile', field_name='Line-of-sight Azimuth', fill_value=None, label_axis='Azimuth', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[n],
                              notes="Consider the vector pointing from the spacecraft to the tangent point (i.e., the line of sight). At the tangent "
                              "point, this vector is parallel to the ground. This variable contains the azimuth angle of this vector, evaluated at "
                              "the tangent point. It follows the typical geophysical convention of degrees East of North (North=0, East=90, South=180, "
                              "West=270). It can vary by a few degrees from the top of the profile to the bottom, so one value is reported per altitude. "
                              "MIGHTI-A and MIGHTI-B will have values approximately 90 degrees apart."
                              )



        ######### Other Metadata Variables #########

        # Chi^2 of line fit to phase
        var = _create_variable(ncfile, '%s_Chi2'%prefix, L21_dict['chi2'], 
                              dimensions=(var_alt_name), depend_0 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Chi^2 of the line fit to unwrapped phase used to determine wind', 
                              display_type='altitude_profile', field_name='Chi^2 of line fit to phase', fill_value=None, label_axis='Chi^2', bin_location=0.5,
                              units='rad^2', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=[n],
                              notes="To convert a complex interferogram to a wind value, a line is fit to the unwrapped phase vs optical path "
                              "difference. The chi^2 value of that fit (here defined as the square of the residual, averaged over optical path "
                              "difference) can be used as a quality-of-fit metric. It is used to find and mask values which have a unusably low "
                              "SNR."
                              )
        
        # ICON velocity vector
        var = _create_variable(ncfile, '%s_Spacecraft_Velocity_Vector'%prefix, L21_dict['icon_velocity_ecef_vector'], 
                              dimensions=('Vector'),
                              format_nc='f8', format_fortran='F', desc="ICON's velocity vector in Earth-centered, Earth-fixed coordinates", 
                              display_type='scalar', field_name='ICON Velocity Vector', fill_value=None, label_axis='S/C Vel', bin_location=0.5,
                              units='m/s', valid_min=-100e6, valid_max=100e6, var_type='metadata', chunk_sizes=[3],
                              notes="A length-3 vector [vx,vy,vz] of ICON's velocity in Earth-centered Earth-fixed (ECEF) coordinates at the "
                              "midpoint time of the observation. The effect of spacecraft velocity has already been removed from the Line_of_Sight_Wind "
                              "variable.")

        # ICON latitude
        var = _create_variable(ncfile, '%s_Spacecraft_Latitude'%prefix, L21_dict['icon_lat'], 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The WGS84 latitude of ICON', 
                              display_type='scalar', field_name='Spacecraft Latitude', fill_value=None, label_axis='S/C Lat', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='metadata', chunk_sizes=1,
                              notes='The latitude of ICON at the midpoint time of the observation, using the WGS84 ellipsoid.')

        # ICON longitude
        var = _create_variable(ncfile, '%s_Spacecraft_Longitude'%prefix, L21_dict['icon_lon'], 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The WGS84 longitude of ICON', 
                              display_type='scalar', field_name='Spacecraft Longitude', fill_value=None, label_axis='S/C Lon', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='metadata', chunk_sizes=1,
                              notes='The longitude (0-360) of ICON at the midpoint time of the observation, using the WGS84 ellipsoid.')

        # ICON altitude
        val = L21_dict['icon_alt']*1e3 # convert to m
        var = _create_variable(ncfile, '%s_Spacecraft_Altitude'%prefix, val, 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The WGS84 altitude of ICON', 
                              display_type='scalar', field_name='Spacecraft Altitude', fill_value=None, label_axis='S/C Alt', bin_location=0.5,
                              units='m', valid_min=100e3, valid_max=2000e3, var_type='metadata', chunk_sizes=1,
                              notes='The altitude of ICON at the midpoint time of the observation, using the WGS84 ellipsoid.')

        # Along-track resolution
        val = L21_dict['resolution_along_track']*1e3 # convert to meters
        var = _create_variable(ncfile, '%s_Resolution_Along_Track'%prefix, val, 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The horizontal resolution in the spacecraft velocity direction', 
                              display_type='scalar', field_name='Along-Track Resolution', fill_value=None, label_axis='Hor Res AT', bin_location=0.5,
                              units='m', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
                              notes=["NOT YET IMPLEMENTED",]
                                     )

        # Cross-track resolution
        val = L21_dict['resolution_cross_track']*1e3 # convert to meters
        var = _create_variable(ncfile, '%s_Resolution_Cross_Track'%prefix, val, 
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The horizontal resolution perpendicular to the spacecraft velocity direction', 
                              display_type='scalar', field_name='Cross-Track Resolution', fill_value=None, label_axis='Hor Res CT', bin_location=0.5,
                              units='m', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
                              notes=["NOT YET IMPLEMENTED",]
                                     )

        # Altitude resolution
        var = _create_variable(ncfile, '%s_Resolution_Altitude'%prefix, L21_dict['resolution_alt']*1e3, # km to meters
                              dimensions=(),
                              format_nc='f8', format_fortran='F', desc='The vertical resolution', 
                              display_type='scalar', field_name='Vertical Resolution', fill_value=None, label_axis='Vert Res', bin_location=0.5,
                              units='m', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
                              notes=["NOT YET IMPLEMENTED",]
                                     )

        # MIGHTI ECEF vectors
        var = _create_variable(ncfile, '%s_Line_of_Sight_Vector'%prefix, L21_dict['mighti_ecef_vectors'], 
                              dimensions=(var_alt_name,'Vector'), depend_0 = var_alt_name, # depend_1 = 'Vector',
                              format_nc='f8', format_fortran='F', desc='The look direction of each MIGHTI line of sight, as a vector in ECEF', 
                              display_type='altitude_profile', field_name='Line-of-sight Vector', fill_value=None, label_axis='LoS Vec', bin_location=0.5,
                              units='', valid_min=-1., valid_max=1., var_type='metadata', chunk_sizes=[n,3],
                              notes="The vector from the spacecraft to the tangent point (i.e., along MIGHTI's line of sight), as a unit "
                              "vector in Earth-centered Earth-fixed (ECEF) coordinates. A vector is provided for each tangent point. If this "
                              "vector is transformed to an azimuth and zenith angle at the tangent point, the zenith angle will be 90 deg, and "
                              "the azimuth angle will be the same as the Line_of_Sight_Azimuth variable."
                              )

        # Bin Size
        var = _create_variable(ncfile, '%s_Bin_Size'%prefix, L21_dict['bin_size'], 
                              dimensions=(),
                              format_nc='i1', format_fortran='I', desc='How many raw samples were binned vertically for each reported sample', 
                              display_type='scalar', field_name='Bin Size', fill_value=None, label_axis='Bin Size', bin_location=0.5,
                              units='', valid_min=np.int8(1), valid_max=np.int8(100), var_type='metadata', chunk_sizes=1,
                              notes="To improve statistics, adjacent rows of the interferogram can be averaged together before the inversion. "
                              "This improves precision at the cost of vertical resolution. If no binning is performed, this value will be 1, "
                              "corresponding to ~2.5 km resolution. A value of 2 corresponds to ~5 km resolution, etc."
                              )

        # Integration order
        var = _create_variable(ncfile, '%s_Integration_Order'%prefix, L21_dict['integration_order'], 
                              dimensions=(),
                              format_nc='i1', format_fortran='I', desc='Order used to discretize the integral for inversion: 0=Riemann, 1=Trapezoidal', 
                              display_type='scalar', field_name='Order', fill_value=None, label_axis='Order', bin_location=0.5,
                              units='', valid_min=np.int8(0), valid_max=np.int8(1), var_type='metadata', chunk_sizes=1,
                              notes="In formulating the inversion, an assumption must be made regarding the choice of basis functions, "
                              "which can be thought of as an assumption regarding the behavior of the wind and amplitude within each altitude "
                              "layer. The most basic assumption is that these quantities are constant within each altitude layer, which corresponds "
                              "to Integration_Order=0. However, if it is assumed that the variation within each layer is linear, "
                              "Integration_Order=1. This sacrifices precision to improve vertical resolution."
                              )

        # How the top layer was handled in the inversion
        var = _create_variable(ncfile, '%s_Top_Layer_Model'%prefix, L21_dict['top_layer'], 
                              dimensions=(),
                              format_nc=str, format_fortran='A', desc='How the top altitudinal layer is handled in the inversion: "exp" or "thin"', 
                              display_type='scalar', field_name='Top Layer', fill_value=None, label_axis='Top Layer', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=1,
                              notes="In formulating the inversion, an assumption must be made about the shape of the emission rate profile "
                              "above the top measured altitude, since this shape is not measured. It can be assumed to go to zero "
                              "(Top_Layer_Model=\"thin\") or assumed to fall off exponentially with a scale height of 26 km, a value extracted "
                              "from running a variety of airglow models (Top_Layer_Model=\"exp\"). Usually this choice will not affect the "
                              "inversion significantly. In cases where it does, the quality variable will be decreased."
                              )    
            
        ncfile.close()
        
    except: # make sure the file is closed
        ncfile.close()
        raise
    
    return L21_full_fn
    
    
    
    
    
    
def level1_to_level21_without_info_file(L1_fn, emission_color, L21_path , data_revision=0, sigma=None, top_layer=None, integration_order=None, 
                                        account_for_local_projection=None, bin_size=None, chi2_thresh=None):
    '''
    High-level function to apply the Level-1-to-Level-2.1 algorithm to a Level 1 file. This version
    of the function requires the user to input the arguments instead of specifying them with an 
    Information.TXT file, as will be done in the Science Data Center. A single Level 1 file can 
    become two Level 2.1 files (one per emission color).
    
    INPUTS:
    
      *  L1_fn               -- TYPE:str.  The full path to the Level 1 file to be processed
      *  emission_color      -- TYPE:str.  'red' or 'green'
      *  L21_path            -- TYPE:str.  The directory the Level 2.1 file will be saved in, including trailing "/"
                                           (e.g., '/home/user/')
      
    OPTIONAL INPUT:
    
      *  data_revision       -- TYPE:int,  The minor version of the data [0-999]. The major version is set
                                           by this software's major version. (default 0)

    MORE OPTIONAL INPUTS - If None, defaults from MIGHTI_L2.global_params will be used 
    
      *  sigma               -- TYPE:float, UNITS:m^-1. The wavenumber of the emission (1/wavelength)
      *  top_layer           -- TYPE:str, 'thin': assume VER goes to zero above top layer
                                          'exp':  assume VER falls off exponentially in altitude
      *  integration_order   -- TYPE:int, 0: Use Riemann-sum rule for discretizing line-of-sight integral
                                          1: Use trapezoidal rule for discretizing line-of-sight integral
      *  account_for_local_projection -- TYPE:bool. If False, a simple inversion is used.
                                         If True, the inversion accounts for the fact that the ray is not 
                                         perfectly tangent to each shell at each point along the ray.
      *  bin_size            -- TYPE:int, The number of rows of the interferogram to bin together to 
                                          improve statistics at the cost of altitude resolution.
      *  chi2_thresh         -- TYPE:float, UNITS:rad^2. The chi^2 value for the line fit (used to analyze
                                                         the inteferogram), below which the data will be 
                                                         ignored, since it is probably just noise.
                                           
    OUTPUTS:      
    
      *  L21_fn              -- TYPE:str.  The full path to the saved L2.1 file.

    '''
    
    # Parse inputs
    if emission_color not in ['red','green']:
        raise ValueError('Argument emission_color=\'%s\' not recognized. Use \'red\' or \'green\'.' % emission_color)
    # For other inputs, just pass them to the lower-level function. It will replace the Nones.
    
    # Read L1 file into a dictionary
    L1_dict = level1_to_dict(L1_fn, emission_color)
    
    # Perform L1 to L2.1 processing
    L21_dict = level1_dict_to_level21_dict(L1_dict, sigma, top_layer = top_layer, 
                                           integration_order = integration_order, 
                                           account_for_local_projection = account_for_local_projection, 
                                           bin_size = bin_size, chi2_thresh = chi2_thresh)
    
    
    # Save L2.1 file
    L21_fn = save_nc_level21(L21_path, L21_dict, data_revision)
    
    return L21_fn
    
    
    
    
    
def read_info_file(info_fn):
    '''
    Read the Information.TXT file that the Science Data Center provides, and return a dictionary
    of parameters and a list of files.
    
    INPUTS:
    
      * info_fn  -- TYPE:str.  Full path to an ASCII file in the following format:
      
                                        [PARAMETERS]
                                        Revision=001
                                        Directory=/path/to/wherever/
                                        <other parameters>

                                        [FILES]
                                        ICON_L1_MIGHTI-A_Science_2017-03-03_191803_v04r006.NC
                                        
                                        
    OUTPUTS:
    
      * info    -- TYPE:dict.         A dictionary of the parameters specified in the file. Keys
                                      and values are both strings.
      * files   -- TYPE:list of str.  A list of the files specified in the file.
      
    '''
    
    # Read the info file to extract parameters for the processing function
    info = {} # expected to have entries: 'Directory' and 'Revision'
    files = [] # files to process
    with open(info_fn, 'r') as f:
        # Read until the [Parameters] line
        line = f.readline()
        while line and '[parameters]' not in line.lower():
            line = f.readline()
        if not line:
            raise IOError('Information file format not understood: "[Parameters]" not found')
        # Read and save all the parameters, until a blank line
        line = f.readline().strip()
        while line:
            x = line.split('=')
            info[x[0]] = x[1]
            line = f.readline().strip()
        # Read until the [Parameters] line
        line = f.readline()
        while line and '[files]' not in line.lower():
            line = f.readline()
        if not line:
            raise IOError('Information file format not understood: "[Files]" not found')
        # Read until the end of the file, recording all the fns
        line = f.readline()
        while line:
            files.append(line.strip())
            line = f.readline()
            
    return info, files
                    
    
  

def level1_to_level21(info_fn):
    '''
    Highest-level function to apply the Level-1-to-Level-2.1 algorithm to a Level 1 file, with 
    input arguments specified via an information file. Many files may be specified in the information
    file; this routine loops over all files specified. For each file, the processing is run twice:
    once for 'red' and once for 'green'. The output L2.1 files will be saved to the Output folder in
    the directory specified in the information file. Summary plot(s) will also be created in the Output
    folder.
    
    INPUTS:
    
      * info_fn  -- TYPE:str.  Full path to an ASCII file in the following format:
      
                                        [PARAMETERS]
                                        Revision=001
                                        Directory=/path/to/wherever/
                                        <other parameters>

                                        [FILES]
                                        ICON_L1_MIGHTI-A_Science_2017-03-03_191803_v04r006.NC
                                        
    OUTPUTS:   
    
      *  ret     -- TYPE:str. '' if everything worked. If not, a human-readable error message for each file that failed
    
    '''    
    
    info, L1_fns = read_info_file(info_fn)
    
    # Parse the info
    # (0) Make sure there's a trailing "/" on the directory
    direc = info['Directory']
    if direc[-1] != '/':
        direc += '/'
    # (1) Add the directory to all the L1 files, which are in the Input folder
    L1_full_fns = []
    for L1_fn in L1_fns:
        L1_full_fns.append(direc + 'Input/' + L1_fn)
    # (2) Parse list of data revision numbers
    s = info['Revision'].split(',')
    data_revision = [int(x) for x in s]
    # if data_revision only has one entry, it should be applied to all input files
    if len(data_revision)==1:
        data_revision = [data_revision[0]]*len(L1_fns)
    assert len(L1_fns)==len(data_revision), "Length of revision list != Length of file list"
    
    tstart = datetime.now()
    # Loop and call the lower-level function which does all the real work.
    L21_fns = []
    failed_L1_fns = []
    failure_messages = []
    for L1_fn, rev in zip(L1_full_fns, data_revision):
        for emission_color in ['red','green']:
            try:
                L21_fn = level1_to_level21_without_info_file(L1_fn, emission_color, direc + 'Output/', data_revision = rev)
                L21_fns.append(L21_fn)
            except Exception as e:
                failed_L1_fns.append(L1_fn)
                failure_messages.append('Failed processing:\n\tL1_file = %s\n\tColor   = %s\n%s\n'%(L1_fn,emission_color,traceback.format_exc()))
    tfiles = datetime.now()
    #print 'Finished creating L2.1 files (%.2f min)' % ((tfiles-tstart).total_seconds()/60.)
    
    # Summary plots
    Afns = [fn for fn in L21_fns if 'MIGHTI-A' in fn]
    Bfns = [fn for fn in L21_fns if 'MIGHTI-B' in fn]
    try:
        if Afns:
            plot_level21(Afns, direc + 'Output/')
        if Bfns:
            plot_level21(Bfns, direc + 'Output/')
    except Exception as e:
        failure_messages.append('Failed creating summary plot:\n%s\n'%(traceback.format_exc()))
    tplots = datetime.now()
    #print 'Finished making L2.1 plots (%.2f min)' % ((tplots-tfiles).total_seconds()/60.)

    
    if not failure_messages: # Everything worked
        return ''
    
    else:
        s = '\n' + '\n'.join(failure_messages)
        return s
    
    


    
    
def level21_to_dict(L21_fns):
    ''' 
    Load a series of Level 2.1 files and return relevant variables in a dictionary. It is
    assumed that all files are from the same emission (red or green)
    
    INPUTS:
    
      *  L21_fns  -- TYPE:list of str.  The paths to the Level 2.1 files to be loaded.
      
    OUTPUTS:
    
      *  L21_dict -- TYPE:dict. A dictionary containing the following variables. Most 
                                are provided as arrays of shape (ny,nt), where ny is the number
                                of altitude samples and nt is the number of time samples.
                                
                  * lat             -- TYPE:array(ny,nt), UNITS:deg. Sample latitudes.
                  * lon             -- TYPE:array(ny,nt), UNITS:deg. Sample longitudes.
                  * alt             -- TYPE:array(ny,nt), UNITS:km.  Sample altitudes.
                  * los_wind        -- TYPE:array(ny,nt), UNITS:m/s. Line-of-sight wind component towards MIGHTI.
                  * los_wind_error  -- TYPE:array(ny,nt), UNITS:m/s. Error in los_wind variable.
                  * local_az        -- TYPE:array(ny,nt), UNITS:deg. Azimuth angle of vector pointing from 
                                       MIGHTI towards the sample location, at the sample location (deg E of N).
                  * amp             -- TYPE:array(ny,nt), UNITS:arb. Fringe amplitude at sample locations.
                  * amp_error       -- TYPE:array(ny,nt), UNITS:arb. Error in amp variable.
                  * time            -- TYPE:array(nt).               Array of datetime objects, one per file.
                  * icon_lat        -- TYPE:array(nt),    UNITS:deg. Spacecraft latitude.
                  * icon_lon        -- TYPE:array(nt),    UNITS:deg. Spacecraft longitude.
                  * icon_alt        -- TYPE:array(nt),    UNITS:km.  Spacecraft altitude
                  * exp_time        -- TYPE:array(nt),    UNITS:sec. Exposure time of each sample.
                  * emission_color  -- TYPE:str,                     'red' or 'green'.
                  * sensor          -- TYPE:list(nt),                each element is 'A' or 'B'
                  * source_files    -- TYPE:list of str,             A copy of the input.
                  * acknowledgement -- TYPE:str,                     The Acknowledgment attribute in the first file
                  
    '''
    
    if len(L21_fns)==0:
        raise ValueError('No files specified')
    
    # Initialize variables to empty list
    variables = ['lat',
                 'lon',
                 'alt',
                 'los_wind',
                 'los_wind_error',
                 'local_az',
                 'amp',
                 'amp_error',
                 'time',
                 'icon_lat',
                 'icon_lon',
                 'icon_alt',
                 'exp_time',
                 'sensor',]

    z = {}
    for v in variables:
        z[v] = []

    first_color = L21_fns[0].split('/')[-1].split('_')[3].split('-')[-1] # Red or Green
    for fn in L21_fns:
        try:
            d = netCDF4.Dataset(fn)
        except IOError as e:
            print 'Error reading file: %s' % fn
            raise
        sens  = d.Instrument[-1] # 'A' or 'B'
        color = fn.split('/')[-1].split('_')[3].split('-')[-1] # Red or Green
        prefix = 'ICON_L2_MIGHTI_%s_%s' % (sens, color.capitalize())    

        # Load all variables
        z['lat'].append(           d.variables['%s_Latitude' % prefix][...])
        z['lon'].append(           d.variables['%s_Longitude' % prefix][...])
        z['alt'].append(           d.variables['%s_Altitude' % prefix][...]*1e-3) # m to km
        z['los_wind'].append(      d.variables['%s_Line_of_Sight_Wind' % prefix][...])
        z['los_wind_error'].append(d.variables['%s_Line_of_Sight_Wind_Error' % prefix][...])
        z['local_az'].append(      d.variables['%s_Line_of_Sight_Azimuth' % prefix][...])
        z['amp'].append(           d.variables['%s_Fringe_Amplitude' % prefix][...])
        z['amp_error'].append(     d.variables['%s_Fringe_Amplitude_Error' % prefix][...])
        time_msec =                d.variables['Epoch'][...].item()
        z['time'].append(datetime(1970,1,1) + timedelta(seconds = 1e-3*time_msec))
        z['icon_lat'].append(      d.variables['%s_Spacecraft_Latitude' % prefix][...].item())
        z['icon_lon'].append(      d.variables['%s_Spacecraft_Longitude' % prefix][...].item())
        z['icon_alt'].append(      d.variables['%s_Spacecraft_Altitude' % prefix][...].item())
        z['exp_time'].append(float(d.Time_Resolution.split(' ')[0]))
        z['sensor'].append(sens)

        ack = d.Acknowledgement    
        assert color == first_color, "Input files do not have the same emission color"
        d.close()

    # Replace masked arrays with nan and other minor things
    for v in variables:
        tmp = [np.ma.filled(x, np.nan) for x in z[v]]
        z[v] = np.vstack(tmp).T
        # If the shape is 1 by something (i.e., only one value per time),
        # then just make it a 1D array instead of 2D
        if np.shape(z[v])[0] == 1:
            z[v] = z[v][0,:]

    # Add extra global information
    z['acknowledgement'] = ack
    z['emission_color'] = color.lower()
    z['source_files'] = L21_fns
    
    return z

    
    
    
    
############################################################################################################
##########################################       Level 2.2       ###########################################
############################################################################################################
   

    
def level21_dict_to_level22_dict(L21_A_dict, L21_B_dict, sph_asym_thresh = None, time_start = None, time_stop = None):
    '''
    Given Level 2.1 data from MIGHTI A and MIGHTI B, process it with the Level 2.2 algorithm. 
    This entails interpolating line-of-sight wind measurements from MIGHTI A and B to a 
    common grid, and rotating to a geographic coordinate system to derive vector horizontal winds. 
    The input files should span 24 hours (0 UT to 23:59:59 UT), with ~20 minutes of extra files
    on either side.
    
    INPUTS:
    
      *  L21_A_dict      -- TYPE:dict.  The dictionary corresponding to 24+ hours of MIGHTI A measurements 
                                        for a single emission color. See "level21_to_dict" for required keys.
      *  L21_B_dict      -- TYPE:dict.  The dictionary corresponding to 24+ hours of MIGHTI B measurements
                                        for a single emission color, which is the same as for the A measurements.
                                        See "level21_to_dict" for required keys.
                                        
    OPTIONAL INPUTS:
    
      *  sph_asym_thresh -- TYPE:float.    Relative difference in emission rate measured by A and B, beyond which
                                           the spherical-asymmetry flag will be raised. Technically, it should be
                                           "fringe amplitude" instead of "emission rate" due to the temperature
                                           dependence. Relative difference is defined as abs(A-B)/((A+B)/2). If 
                                           None (default), the default from MIGHTI_L2.global_params will be used.
      *  time_start      -- TYPE:datetime. A timezone-naive datetime in UT, specifying the beginning of the interval
                                           which the Level 2.2 data product should span. (Some L2.1 files from before
                                           the start time are needed to handle the crossover). If None (default), 
                                           the start time defaults to the first 0 UT after the first input file's time.
                                           Note that the time changes from the top to the bottom of a L2.2 profile. We
                                           use the top of the field of view as the reference. This means that it's
                                           possible for some time samples to occur before the start time or after the
                                           stop time.
      *  time_stop       -- TYPE:datetime. A timezone-naive datetime in UT, specifying the end of the interval
                                           which the Level 2.2 data product should span. (Some L2.1 files from after
                                           the start time are needed to handle the crossover). If None (default), 
                                           the stop time defaults to 24 hours after the start time.
                                        
    OUTPUTS:
    
      *  L22_dict        -- TYPE:dict.  A dictionary containing the following variables. Most are given as 
                                        arrays with shape (ny,nx) where ny is the number of altitude grid
                                        points, and nx is the number of horizontal (or time) grid points.
                                        
                    * lat             -- TYPE:array(ny,nx),    UNITS:deg.  Latitude of each point on the grid.
                    * lon             -- TYPE:array(ny,nx),    UNITS:deg.  Longitude of each point on the grid [0-360]
                    * lon_unwrapped   -- TYPE:array(ny,nx),    UNITS:deg.  Same as lon, but with 0/360 jumps removed
                    * alt             -- TYPE:array(nx),       UNITS:km.   Altitude of each row of the grid.
                    * u               -- TYPE:array(ny,nx),    UNITS:m/s.  Estimated zonal wind (positive eastward).
                    * v               -- TYPE:array(ny,nx),    UNITS:m/s.  Estimated meridional wind (positive northward).
                    * u_error         -- TYPE:array(ny,nx),    UNITS:m/s.  Uncertainty in u.
                    * v_error         -- TYPE:array(ny,nx),    UNITS:m/s.  Uncertainty in v.
                    * error_flags     -- TYPE:array(ny,nx,ne), UNITS:none. The error flags (either 0 or 1) for each point
                                                                           in the grid. Each point has a number of error flags, 
                                                                           which are set to 1 under the following 
                                                                           circumstances:
                                                                          
                                                                            * 0  = missing A file
                                                                            * 1  = missing B file
                                                                            * 2  = A signal too weak
                                                                            * 3  = B signal too weak
                                                                            * 4  = A did not sample this altitude
                                                                            * 5  = B did not sample this altitude
                                                                            * 6  = spherical asymmetry: A&B VER 
                                                                                   estimates disagree
                                                                            * 7  = unknown error
                                                                            
                   * epoch_full      -- TYPE:array(ny,nx),    UNITS:none. The average between the time of the MIGHTI A 
                                                                          and B measurements that contribute to this 
                                                                          grid point, given as a datetime object. 
                   * epoch           -- TYPE:array(nx),       UNITS:none. mean(epoch_full,axis=0), accounting for missing values, and
                                                                          given as a datetime object. This is useful because 
                                                                          while the time varies along a column, it doesn't vary 
                                                                          much (+/- 30 seconds or so).
                   * time_start      -- TYPE:datetime                     The start time for defining the reconstruction grid
                   * time_stop       -- TYPE:datetime                     The stop time for defining the reconstruction grid
                   * time_delta      -- TYPE:array(ny,nx),    UNITS:s.    The difference between the time of the MIGHTI
                                                                          A and B measurements that contribute to this 
                                                                          grid point.
                   * fringe_amp_A    -- TYPE:array(ny,nc),    UNITS:none. The fringe amplitude measured by MIGHTI A (roughly
                                                                          proportional to VER)
                   * fringe_amp_B    -- TYPE:array(ny,nc),    UNITS:none. The fringe amplitude measured by MIGHTI B (roughly
                                                                          proportional to VER)
                   * fringe_amp      -- TYPE:array(ny,nc),    UNITS:none. Mean of the above two quantitiies.
                   * fringe_amp_rel_diff -- TYPE:array(ny,nx),UNITS:none. The difference between the fringe amplitude 
                                                                          of the MIGHTI A and B measurements that 
                                                                          contribute to this grid point, divided by the
                                                                          mean. When this is high, it indicates that 
                                                                          spherical asymmetry may be a problem.
                   * emission_color  -- TYPE:str,                         'red' or 'green'.        
                   * source_files    -- TYPE:list of str,                 All the files used to create the data product,
                                                                          including the full paths.
                   * acknowledgement -- TYPE:str,                         The Acknowlegment copied from the L2.1 files
    '''

    ################################## Parse Inputs ####################################
    emission_color = L21_A_dict['emission_color']
    params = global_params[emission_color]
    if sph_asym_thresh is None:
        sph_asym_thresh = params['sph_asym_thresh']
    
    assert L21_A_dict['emission_color'] == L21_B_dict['emission_color'], "Files for A and B are for different emissions"
    
    
    lat_A      =       L21_A_dict['lat']
    lon_A_raw  =       L21_A_dict['lon']
    alt_A      =       L21_A_dict['alt']
    los_wind_A =       L21_A_dict['los_wind']
    los_wind_A_err =   L21_A_dict['los_wind_error']
    local_az_A =       L21_A_dict['local_az']
    amp_A      =       L21_A_dict['amp']
    time_A     =       L21_A_dict['time']
    exp_time_A =       L21_A_dict['exp_time']
    emission_color =   L21_A_dict['emission_color']
    N_alts_A, N_times_A = np.shape(lat_A)

    lat_B      =       L21_B_dict['lat']
    lon_B_raw  =       L21_B_dict['lon']
    alt_B      =       L21_B_dict['alt']
    los_wind_B =       L21_B_dict['los_wind']
    los_wind_B_err =   L21_B_dict['los_wind_error']
    local_az_B =       L21_B_dict['local_az']
    amp_B      =       L21_B_dict['amp']
    time_B     =       L21_B_dict['time']
    exp_time_B =       L21_B_dict['exp_time']
    emission_color =   L21_B_dict['emission_color']
    N_alts_B, N_times_B = np.shape(lat_B)
    
    if time_start is None:
        # Default to 0 UT on date of middle A image
        t_mid = time_A[N_times_A/2]
        time_start = datetime(t_mid.year, t_mid.month, t_mid.day)
        
    if time_stop is None:
        # Default to 24 hours after the start time
        time_stop = time_start + timedelta(hours=24)
    
    assert any(time_A < time_stop ), "No MIGHTI-A files found before time_stop=%s"%(time_stop)
    assert any(time_B < time_stop ), "No MIGHTI-B files found before time_stop=%s"%(time_stop)
    assert any(time_A > time_start), "No MIGHTI-A files found before time_start=%s"%(time_start)
    assert any(time_B > time_start), "No MIGHTI-B files found before time_start=%s"%(time_start)
    
    ####################### Define reconstruction grid: lon/alt ########################

    # Unwrap the sample longitude arrays to avoid 0/360 jumps
    lon_A = lon_A_raw.copy()
    lon_B = lon_B_raw.copy()
    lon_A[:,0] = fix_longitudes(lon_A[:,0], lon_A[-1,0]) # Use top first A longitude as target
    lon_B[:,0] = fix_longitudes(lon_B[:,0], lon_A[-1,0]) # Use top first A longitude as target
    for i in range(N_alts_A):
        lon_A[i,:] = fix_longitudes(lon_A[i,:], lon_A[i,0])
    for i in range(N_alts_B):
        lon_B[i,:] = fix_longitudes(lon_B[i,:], lon_B[i,0])

    # Determine start longitude
    # This should be the longitude of the tangent point at the top of the profile for the first
    # exposure after the start time, averaged between A and B.
    iA = bisect.bisect(time_A, time_start) # First A file after start time.
    iB = bisect.bisect(time_B, time_start) # First B file after start time.
    start_lon_A = lon_A[-1,iA]
    start_lon_B = lon_B[-1,iB]
    assert abs(start_lon_A - start_lon_B) < 90.,"A and B start longitudes are off by 360 deg."
    start_lon = np.mean([start_lon_A, start_lon_B])

    # Determine stop longitude
    # This should be the longitude of the tangent point at the top of the profile for the first
    # exposure before the stop time, averaged between A and B.
    iA = bisect.bisect(time_A, time_stop) - 1 # First A file before (or equal to) stop time
    iB = bisect.bisect(time_B, time_stop) - 1 # First B file before (or equal to) stop time
    stop_lon_A = lon_A[-1,iA]
    stop_lon_B = lon_B[-1,iB]
    assert abs(stop_lon_A - stop_lon_B) < 90.,"A and B start longitudes are off by 360 deg."
    stop_lon = np.mean([stop_lon_A, stop_lon_B])
    
    time_min = min(min(time_A), min(time_B))
    time_max = max(max(time_A), max(time_B))
    assert stop_lon > start_lon, 'No files found between time_start="%s" and time_stop="%s". Files provided span from "%s" to "%s"'%(time_start, \
                                                                                   time_stop, time_min, time_max)

    # Determine how finely to define longitude grid based on minimum L2.1 sampling rate
    lon_res_A = (np.diff(lon_A,axis=1)).min()
    lon_res_B = (np.diff(lon_B,axis=1)).min()
    lon_res = min(lon_res_A, lon_res_B)
    # Define longitude grid
    lon_vec = np.arange(start_lon, stop_lon, lon_res)

    # Define altitude grid based on the min and max in the L2.1 data (i.e., don't throw out data)
    # Define altitude resolution based upon the resolution of L2.1 data
    alt_min = min(alt_A.min(), alt_B.min())
    alt_max = max(alt_A.max(), alt_B.max())
    nalts = max(np.shape(alt_A)[0], np.shape(alt_B)[0])
    alt = np.linspace(alt_min, alt_max, nalts)
    alt_res = alt[1] - alt[0]

    # Create grid for longitude (since we might change how this is done in the future)
    lon,_ = np.meshgrid(lon_vec, alt)
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
    U = np.nan*np.zeros((N_alts, N_lons))                # zonal wind
    V = np.nan*np.zeros((N_alts, N_lons))                # meridional wind
    U_err = np.nan*np.zeros((N_alts, N_lons))            # zonal wind uncertainty
    V_err = np.nan*np.zeros((N_alts, N_lons))            # meridional wind uncertainty
    lat = np.nan*np.zeros((N_alts, N_lons))              # latitude
    time = np.empty((N_alts, N_lons), dtype=object)      # time ascribed to L2.2 data point (as datetime objects)
    time_delta = np.nan*np.zeros((N_alts, N_lons))       # difference between A and B times used (seconds)
    ver_A = np.nan*np.zeros((N_alts, N_lons))            # fringe amplitude from A (related to VER)
    ver_B = np.nan*np.zeros((N_alts, N_lons))            # fringe amplitude from B (related to VER)
    ver   = np.nan*np.zeros((N_alts, N_lons))            # fringe amplitude (mean of A and B)
    ver_rel_diff = np.nan*np.zeros((N_alts, N_lons))     # relative difference in A and B VER
    error_flags = np.zeros((N_alts, N_lons, 8), dtype=bool)      # Error flags, one set per grid point. See above for definition.
    
    # Loop over the reconstruction altitudes
    for i in range(N_alts):
            
        alt_pt = alt[i]
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
                    error_flags[i,k,4] = 1
                if alt_pt > altmax_B or alt_pt < altmin_B:
                    error_flags[i,k,5] = 1
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
                    # Do interpolation of value to the desired altitude, for each longitude
                    val_0, val_0_err = interpolate_linear(alt_AB[:,0], val[:,0], alt_pt, 
                                                          prop_err = True, yerr = valerr[:,0])
                    val_1, val_1_err = interpolate_linear(alt_AB[:,1], val[:,1], alt_pt, prop_err = True, yerr = valerr[:,1])
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
            
            # Check if the interpolated L2.1 data were np.nan, which would have to 
            # indicate weak signal (or is there another reason they might be masked?)
            if np.isnan(los_wind_A_pt):
                error_flags[i,k,2] = 1
            if np.isnan(los_wind_B_pt):
                error_flags[i,k,3] = 1
                
            # Check spherical symmetry
            ver_rel_diff_pt = abs(ver_A_pt - ver_B_pt)/np.mean([ver_A_pt,ver_B_pt])
            if ver_rel_diff_pt > sph_asym_thresh:
                error_flags[i,k,6] = 1
                # Should we nan the data too?
            
            # Unknown error - this should never happen
            if (np.isnan(u) or np.isnan(v)) and all(error_flags[i,k,:] == 0): # Unknown error
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
            ver[i,k]   = (ver_A_pt + ver_B_pt)/2.
            ver_rel_diff[i,k] = ver_rel_diff_pt
                
    # Create average time per column (i.e., per longitude)
    epoch = np.empty(N_lons, dtype=object)
    for k in range(N_lons):
        t_k = [ti for ti in time[:,k] if ti is not None]
        if len(t_k) > 0:
            # Take mean of that array
            dt_k = np.array([(ti - t_k[0]).total_seconds() for ti in t_k])
            epoch[k] = t_k[0] + timedelta(seconds=np.mean(dt_k))
    
    
    # Create dictionary to be returned
    L22_dict = {}
    L22_dict['lat'] = lat
    L22_dict['lon'] = np.mod(lon, 360.)
    L22_dict['lon_unwrapped'] = lon
    L22_dict['alt'] = alt
    L22_dict['u'] = U
    L22_dict['v'] = V
    L22_dict['u_error'] = U_err
    L22_dict['v_error'] = V_err
    L22_dict['error_flags'] = error_flags
    L22_dict['epoch_full'] = time # 2D
    L22_dict['epoch'] = epoch # 1D
    L22_dict['time_delta'] = time_delta
    L22_dict['time_start'] = time_start
    L22_dict['time_stop'] = time_stop
    L22_dict['fringe_amp_A'] = ver_A
    L22_dict['fringe_amp_B'] = ver_B
    L22_dict['fringe_amp'] = ver
    L22_dict['fringe_amp_rel_diff'] = ver_rel_diff
    L22_dict['emission_color'] = emission_color
    L22_dict['source_files'] = np.concatenate((L21_A_dict['source_files'], L21_B_dict['source_files']))
    L22_dict['acknowledgement'] = L21_A_dict['acknowledgement']
    L22_dict['wind_quality'] = np.nan*np.ones((N_alts,N_lons)) # TODO
    
    return L22_dict





def save_nc_level22(path, L22_dict, data_revision = 0):
    '''
    Take the output of the Level 2.2 processing and save it as a NetCDF4 file in the official format.
    NetCDF4 file conventions taken from "Science Operations Center Data Product Conventions" Rev 0.5.
    
    INPUTS:
    
      *  path        -- TYPE:str.  The directory the file will be saved in, including trailing "/"
                                   (e.g., '/home/user/')
      *  L22_dict    -- TYPE:dict. A dictionary containing output variables of the Level 2.2 processing.
                                   See documentation for level21_dict_to_level22_dict(...) for details.
                     
    OPTIONAL INPUTS:
    
      *  data_revision       -- TYPE:int,  The minor version of the data [0-999]. The major version is set
                                           by the software's major version.
                                           
    OUTPUTS:
    
      *  L22_fn      -- TYPE:str.  The full path to the saved file.

    '''
    
    L22_dict = copy.deepcopy(L22_dict) # because netCDF4 actually modifies its input
   
    data_version_major = software_version_major # enforced as per Data Product Conventions Document
    
    #################### Compile variables to write in file ######################
    ### Timing:
    t_all = filter(None, L22_dict['epoch_full'].flatten()) # Extract all non-None grid times as a 1-D array
    t_start = L22_dict['time_start']
    t_stop  = L22_dict['time_stop']
    t_mid   = t_start + timedelta(seconds=(t_stop - t_start).total_seconds()/2) # midpoint time
    t_start_msec = (t_start - datetime(1970,1,1)).total_seconds()*1e3 # milliseconds since epoch
    t_stop_msec  = (t_stop  - datetime(1970,1,1)).total_seconds()*1e3
    t_mid_msec   = (t_mid   - datetime(1970,1,1)).total_seconds()*1e3
    t_start_msec = np.int64(np.round(t_start_msec)) # cast to signed 64 bit integer
    t_stop_msec  = np.int64(np.round(t_stop_msec)) 
    t_mid_msec   = np.int64(np.round(t_mid_msec))
    t_file  = datetime.now()   # time this file was created  
    ### Who's running this process
    user_name = getpass.getuser()
    
    ### Parent files
    parents = [] # This will go in global attr Parents
    for source_fn in L22_dict['source_files']:
        s = source_fn.split('/')[-1].split('.')
        pre = '.'.join(s[:-1])
        post = s[-1].upper()
        parents.append('%s > %s' % (post, pre))

    ######################### Open file for writing ##############################
    L22_fn = 'ICON_L2_MIGHTI_Vector-Wind-%s_%s_v%02ir%03i.NC' % (L22_dict['emission_color'].capitalize(),
                                                        t_start.strftime('%Y-%m-%d'),
                                                        data_version_major, data_revision)
    L22_full_fn = '%s%s'%(path, L22_fn)
    ncfile = netCDF4.Dataset(L22_full_fn,mode='w',format='NETCDF4') 
    
    try: # always close file if an error occurs
    
        ########################## Global Attributes #################################
        ncfile.setncattr_string('Acknowledgement',                L22_dict['acknowledgement'])
        ncfile.setncattr_string('ADID_Ref',                       'NASA Contract > NNG12FA45C')
        ncfile.setncattr_string('Calibration_File',               '')
        ncfile.setncattr_string('Conventions',                    'SPDF ISTP/IACG Modified for NetCDF')
        ncfile.setncattr_string('Data_Level',                     'L2.2')
        ncfile.setncattr_string('Data_Type',                      'DP22 > Data Product 2.2: Cardinal Vector Winds')
        ncfile.Data_Version_Major =                               np.uint16(data_version_major)
        ncfile.Data_Revision =                                    np.uint16(data_revision)
        ncfile.Data_Version =                                     data_version_major + 0.001 * data_revision
        ncfile.setncattr_string('Date_Stop',                      t_stop.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Date_Start',                     t_start.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Description',                    'ICON MIGHTI Cardinal Vector Winds (DP 2.2)')
        ncfile.setncattr_string('Descriptor',                     'MIGHTI > Michelson Interferometer for Global High-resolution Thermospheric Imaging')
        ncfile.setncattr_string('Discipline',                     'Space Physics > Ionospheric Science')
        ncfile.setncattr_string('File',                           L22_fn)
        ncfile.setncattr_string('File_Date',                      t_file.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Generated_By',                   'ICON SDC > ICON UIUC MIGHTI L2.2 Processor v%s, B. J. Harding' % __version__)
        ncfile.setncattr_string('Generation_Date',                t_file.strftime('%Y%m%d'))
        ncfile.setncattr_string('History',                        'v1.0: First release of MIGHTI L2.1/L2.2 software, B. J. Harding, 05 Mar 2018')
        ncfile.setncattr_string('HTTP_LINK',                      'http://icon.ssl.berkeley.edu/Instruments/MIGHTI')
        ncfile.setncattr_string('Instrument',                     'MIGHTI')
        ncfile.setncattr_string('Instrument_Type',                'Imagers (space)')
        ncfile.setncattr_string('Link_Text',                      'MIGHTI Cardinal Vector Winds (DP 2.2)')
        ncfile.setncattr_string('Link_Title',                     'ICON MIGHTI')
        ncfile.setncattr_string('Logical_File_ID',                L22_fn[:-3])
        ncfile.setncattr_string('Logical_Source',                 'ICON_L2_MIGHTI_')
        ncfile.setncattr_string('Logical_Source_Description',     'MIGHTI - Cardinal Vector Winds')
        ncfile.setncattr_string('Mission_Group',                  'Ionospheric Investigations')
        ncfile.setncattr_string('MODS',                           ncfile.History)
        ncfile.setncattr_string('Parents',                        parents)
        ncfile.setncattr_string('PI_Affiliation',                 'UC Berkeley > SSL')
        ncfile.setncattr_string('PI_Name',                        'T. J. Immel')
        ncfile.setncattr_string('Project',                        'NASA > ICON')
        ncfile.setncattr_string('Rules_of_Use',                   'Public Data for Scientific Use')
        ncfile.setncattr_string('Software_Version',               'ICON SDC > ICON UIUC MIGHTI L2.2 Processor v%s' % __version__)
        ncfile.setncattr_string('Source_Name',                    'ICON > Ionospheric Connection Explorer')
        ncfile.setncattr_string('Spacecraft_ID',                  'NASA > ICON - 493')
        ncfile.setncattr_string('Text',                           'ICON explores the boundary between Earth and space - the ionosphere - ' +\
                                                                  'to understand the physical connection between our world and the immediate '+\
                                                                  'space environment around us. Visit \'http://icon.ssl.berkeley.edu\' for more details.')
        ncfile.setncattr_string('Text_Supplement',                ['This data product contains cardinal (i.e., zonal and meridional) thermospheric winds ' +
        'obtained by combining Level 2.1 (line-of-sight winds) from MIGHTI A and MIGHTI B. The cardinal winds are given as a function of time ' +
        '(spanning 24 hours) and altitude (spanning nominally 90-300 km). There is one file per emission color (red or green).',
         
        'Cardinal wind observations are enabled by the 90-degree offset between the two MIGHTI sensors. First, MIGHTI A measures a wind component along '+
        'its line of sight. Five to eight minutes later, depending on tangent point altitude, the spacecraft has moved to a position such that MIGHTI B '+
        'measures a nearly orthogonal wind component at approximately the same location. A coordinate rotation is performed on the two line-of-sight '+
        'components to obtain the northward and eastward components reported in this file. The assumption is that the thermospheric wind has not changed '+
        'during this interval. Because the Level 2.1 data are naturally on an irregular grid, '+
        'they are first interpolated to a regular grid of longitude and altitude before the coordinate rotation is performed. See Harding et al. [2017, '+
        'doi:10.1007/s11214-017-0359-3] for more details of the Level 2.2 algorithm.'])
        ncfile.setncattr_string('Time_Resolution',                '30 or 60 seconds')
        ncfile.setncattr_string('Title',                          'ICON MIGHTI Cardinal Vector Winds (DP 2.2)')
        


        ################################## Dimensions ########################################
        prefix = 'ICON_L2_MIGHTI_%s' % (L22_dict['emission_color'].capitalize()) # variable prefix
        var_alt_name = '%s_Altitude'%prefix
        
        ny,nx = np.shape(L22_dict['u'])
        ncfile.createDimension('Epoch',nx)
        ncfile.createDimension(var_alt_name, ny)
        ncfile.createDimension('N_Flags', np.shape(L22_dict['error_flags'])[2])
        


        ################################## Variables #########################################
        # Some variables are transposed to satisfy the ICON requirement that EPOCH is the first dimension
        
        ######### Timing Variables #########
        # Two timing variables are used: The simple one which is 1 dimensional, and the 
        # more correct 2 dimensional one, since the time can vary slightly along altitude,
        # but only by +/- 30-60 seconds.
        # Epoch = 1 dimensional
        # Epoch_Full = 2 dimensional
        # Each variable has two forms: "milliseconds since Jan 1, 1970" and a human-readable string
        
        epoch_msec = np.zeros(nx, dtype=np.int64)
        t_fillval = np.int64(-1)
        for j in range(nx):
            if L22_dict['epoch'][j] is None:
                epoch_msec[j] = t_fillval
            else:
                epoch_msec[j] = np.int64(np.round((L22_dict['epoch'][j] - datetime(1970,1,1)).total_seconds()*1e3))
        var = _create_variable(ncfile, 'Epoch', epoch_msec, 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i8', format_fortran='I', desc='Sample time, average of A and B measurements. Number of msec since Jan 1, 1970.', 
                              display_type='time_series', field_name='Time', fill_value=t_fillval, label_axis='Time', bin_location=0.5,
                              units='ms', valid_min=0, valid_max=1000*365*86400*1000, var_type='support_data', chunk_sizes=[nx],
                              notes="A one-dimensional array defining the time dimension of the two-dimensional data grid (the other dimension being "
                              "altitude). This is the average of the MIGHTI-A and "
                              "MIGHTI-B sample times, which differ by 5-8 minutes. The matchup between MIGHTI-A and B happens at slightly different "
                              "times at different altitudes, a complication which is ignored by this variable. The effect is small (plus or minus 30-60 "
                              "seconds), but in cases where it is important, it is recommended to use the alternative time variable "
                              "Epoch_Full, which is two dimensional and captures the variation with altitude."
                              )
        
        epoch_full_msec = np.zeros((ny,nx),dtype=np.int64)
        for i in range(ny):
            for j in range(nx):
                if L22_dict['epoch_full'][i,j] is None:
                    epoch_full_msec[i,j] = t_fillval
                else:
                    epoch_full_msec[i,j] = np.int64(np.round((L22_dict['epoch_full'][i,j] - datetime(1970,1,1)).total_seconds()*1e3))        
        var = _create_variable(ncfile, 'Epoch_Full', epoch_full_msec.T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='i8', format_fortran='I', desc='Sample time, midpoint of A and B measurements. Number of msec since Jan 1, 1970.', 
                              display_type='image', field_name='Time', fill_value=t_fillval, label_axis='Time', bin_location=0.5,
                              units='ms', valid_min=0, valid_max=1000*365*86400*1000, var_type='support_data', chunk_sizes=[nx,ny],
                              notes="See the notes for the variable Epoch. This variable is the same as Epoch but contains a second dimension, "
                              "which captures the small (30-60 second) variation of time with altitude. For most applications this is expected to "
                              "be negligible, and Epoch can be used instead of this variable. Also see the variable Time_Delta, which contains "
                              "the difference between the MIGHTI-A and MIGHTI-B times that contributed to each point. (This variable "
                              "contains the average)."
                              )
        
        # Human readable versions of those two variables
        utctime = []
        for t in L22_dict['epoch']:
            if t is None:
                utctime.append('')
            else:
                utctime.append(t.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3])
        utctime = np.array(utctime)
        var = _create_variable(ncfile, '%s_UTC_Time'%prefix, utctime, 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc=str, format_fortran='A', desc='Sample time, average of A and B measurements.', 
                              display_type='time_series', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='support_data', chunk_sizes=[nx],
                              notes="This variable is the same as Epoch but is formatted as a human-readable string. Missing "
                              "grid points are labeled with empty strings.")       
        
        utctime_full = np.empty((ny,nx), dtype='S23')
        for i in range(ny):
            for j in range(nx):
                t = L22_dict['epoch_full'][i,j]
                if t is None:
                    utctime_full[i,j] = ''
                else:
                    utctime_full[i,j] = t.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
        #utctime_full = np.array(utctime_full)
        var = _create_variable(ncfile, '%s_UTC_Time_Full'%prefix, utctime_full.T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc=str, format_fortran='A', desc='Sample time, average of A and B measurements.', 
                              display_type='image', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='support_data', chunk_sizes=[nx,ny],
                              notes="This variable is the same as Epoch_Full but is formatted as a human-readable string. Missing "
                              "grid points are labeled with empty strings.")       
        
        ######### Data Variables #########
        
        # Zonal Wind
        var = _create_variable(ncfile, '%s_Zonal_Wind'%prefix, L22_dict['u'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Zonal component of the horizontal wind. Positive Eastward.', 
                              display_type='image', field_name='Zonal Wind', fill_value=None, label_axis='', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="The zonal (positive eastward) and meridional (positive northward) winds are the primary "
                              "data product in this file. They are defined on a grid with dimensions of time and altitude, "
                              "spanning 24 hours and nominally 90-300 km (150-300 km for the red channel). The altitude, time, "
                              "latitude and longitude corresponding to each point in the grid are given by other variables in "
                              "this file. It should be noted that while each measurement is ascribed to a particular latitude, "
                              "longitude, altitude, and time, it is actually an average over many hundreds of kilometers "
                              "horizontally and 2.5-30 kilometers vertically (depending on the binning). It also assumes stationarity "
                              "over the 5-8 minutes between the MIGHTI-A and B measurements used for each point. See Harding et "
                              "al. [2017, doi:10.1007/s11214-017-0359-3] for a more complete discussion of the inversion algorithm."
                              )
        
        # Meridional Wind
        var = _create_variable(ncfile, '%s_Meridional_Wind'%prefix, L22_dict['v'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Meridional component of the horizontal wind. Positive Northward.', 
                              display_type='image', field_name='Meridional Wind', fill_value=None, label_axis='', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="The zonal (positive eastward) and meridional (positive northward) winds are the primary "
                              "data product in this file. They are defined on a grid with dimensions of time and altitude, "
                              "spanning 24 hours and nominally 90-300 km (150-300 km for the red channel). The altitude, time, "
                              "latitude and longitude corresponding to each point in the grid are given by other variables in "
                              "this file. It should be noted that while each measurement is ascribed to a particular latitude, "
                              "longitude, altitude, and time, it is actually an average over many hundreds of kilometers "
                              "horizontally and 2.5-30 kilometers vertically (depending on the binning). It also assumes stationarity "
                              "over the 5-8 minutes between the MIGHTI-A and B measurements used for each point. See Harding et "
                              "al. [2017, doi:10.1007/s11214-017-0359-3] for a more complete discussion of the inversion algorithm."
                              )    
                                      
        # Zonal Wind Error
        var = _create_variable(ncfile, '%s_Zonal_Wind_Error'%prefix, L22_dict['u_error'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Error in the zonal wind estimate.', 
                              display_type='image', field_name='Zonal Wind Error', fill_value=None, label_axis='', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="The statistical (1-sigma) error in the zonal wind, propagated from the error in the "
                              "L2.1 (line-of-sight wind) files. This is usually dominated by shot noise in the detectors, but "
                              "also includes the effects of dark and read noise, as well as calibrations errors  (e.g., the "
                              "zero wind calibration), and spacecraft pointing error (which affects the uncertainty in removing "
                              "the spacecraft velocity from the observed velocity). Other systematic errors or biases may exist "
                              "(e.g., the effect of gradients along the line of sight) which are not included in this variable."
                              )
        
        # Meridional Wind Error
        var = _create_variable(ncfile, '%s_Meridional_Wind_Error'%prefix, L22_dict['v_error'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Error in the meridional wind estimate.', 
                              display_type='image', field_name='Meridional Wind Error', fill_value=None, label_axis='', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="The statistical (1-sigma) error in the meridional wind, propagated from the error in the "
                              "L2.1 (line-of-sight wind) files. This is usually dominated by shot noise in the detectors, but "
                              "also includes the effects of dark and read noise, as well as calibrations errors  (e.g., the "
                              "zero wind calibration), and spacecraft pointing error (which affects the uncertainty in removing "
                              "the spacecraft velocity from the observed velocity). Other systematic errors or biases may exist "
                              "(e.g., the effect of gradients along the line of sight) which are not included in this variable."
                              ) 
        
        # Quality code
        var = _create_variable(ncfile, '%s_Wind_Quality'%prefix, L22_dict['wind_quality'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='A quantification of the quality, from 0 (Bad) to 1 (Good)', 
                              display_type='image', field_name='Wind Quality', fill_value=None, label_axis='Quality', bin_location=0.5,
                              units='arb', valid_min=0.0, valid_max=1.0, var_type='data', chunk_sizes=[nx,ny],
                              notes=["NOT YET IMPLEMENTED",
                                     "While the intent is that the Wind_Error variable accurately characterizes the statistical "
                                     "error in the wind data, it is possible that systematic errors are present, or that the statistical error "
                                     "estimation is not accurate. If it is suspected that this is the case, the quality will be less than 1.0. If "
                                     "the data are definitely unusable, the the quality will be 0.0 and the sample will be masked. Users should "
                                     "exercise caution when the quality is less than 1.0."
                                     ])
        
         # Fringe amplitude (combining MIGHTI-A and MIGHTI-B)
        var = _create_variable(ncfile, '%s_Fringe_Amplitude'%prefix, L22_dict['fringe_amp'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Fringe Amplitude', 
                              display_type='image', field_name='Fringe Amplitude A', fill_value=None, label_axis='', bin_location=0.5,
                              units='arb', valid_min=-1e30, valid_max=1e30, var_type='data', chunk_sizes=[nx,ny],
                              notes="An approximate volume emission rate (VER) profile in arbitrary units, estimated by combining "
                              "MIGHTI-A and MIGHTI-B data. Technically this is not the VER, but rather the amplitude of the fringes, "
                              "which has a dependence on thermospheric temperature and background emission. Thus, it does not truly "
                              "represent volume emission rate. However, it is a useful proxy. The units are arbitrary, as the "
                              "fringe amplitudes are not calibrated. See also variables Fringe_Amplitude_Relative_Difference, "
                              "Fringe_Amplitude_A, and Fringe_Amplitude_B."
                              )
        
        
        
        ######### Data Location Variables #########
        
        # Altitude
        val = L22_dict['alt']*1e3 # convert to meters
        var = _create_variable(ncfile, var_alt_name, val, 
                              dimensions=(var_alt_name), depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='WGS84 altitude of each wind sample', 
                              display_type='image', field_name='Altitude', fill_value=None, label_axis='', bin_location=0.5,
                              units='m', valid_min=50e3, valid_max=1000e3, var_type='support_data', chunk_sizes=[ny], monoton='Increase',
                              notes="A one-dimensional array defining the altitude dimension of the data grid (the other dimension "
                              "being time). For each Level 2.2 file, the altitude grid is defined based on the minimum and maximum altitudes, "
                              "(and highest resolution) in the Level 2.1 (line-of-sight wind) files. Altitude is defined using the WGS84 ellipsoid."
                              )

        
        # Longitude
        var = _create_variable(ncfile, '%s_Longitude'%prefix, L22_dict['lon'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='WGS84 longitude of each wind sample', 
                              display_type='image', field_name='Longitude', fill_value=None, label_axis='', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[nx,ny],
                              notes="A two-dimensional array defining the longitude (0-360 deg) of the two-dimensional data grid. In the "
                              "initial implementation, the longitude is constant with altitude, but this may change in the future to capture "
                              "the slight (few deg) variation with altitude. Longitude is defined using the WGS84 ellipsoid. It should be noted "
                              "that while a single longitude value is given for each point, the observation is "
                              "inherently a horizontal average over many hundreds of kilometers."
                              )
                              
        # Latitude
        var = _create_variable(ncfile, '%s_Latitude'%prefix, L22_dict['lat'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='WGS84 latitude of each wind sample', 
                              display_type='image', field_name='Latitude', fill_value=None, label_axis='', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='support_data', chunk_sizes=[nx,ny],
                              notes="A two-dimensional array defining the latitude of the two-dimensional data grid. The latitude varies "
                              "only slightly (a few deg) with altitude, but this variation is included. Latitude is defined using the WGS84 "
                              "ellipsoid. It should be noted that while a single longitude value is given for each point, the observation is "
                              "inherently a horizontal average over many hundreds of kilometers."
                              )
                              
            

        ######### Other Metadata Variables #########
        
        
        # Difference between the MIGHTI-A and MIGHTI-B times contributing to each point
        var = _create_variable(ncfile, '%s_Time_Delta'%prefix, L22_dict['time_delta'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Difference between MIGHTI-A and B times contributing to each point', 
                              display_type='image', field_name='Time Delta', fill_value=t_fillval, label_axis='Time Delta', bin_location=0.5,
                              units='s', valid_min=-30.*60, valid_max=30.*60, var_type='support_data', chunk_sizes=[nx,ny],
                              notes="To determine the cardinal wind at each point, a MIGHTI-A line-of-sight wind is combined "
                              "with a MIGHTI-B line-of-sight wind from several minutes later. This variable contains this time "
                              "difference for every point. During standard operations, this variable should be positive, but can "
                              "potentially become negative during conjugate operations or when ICON is observing to the south."
                              )
        
        # Fringe amplitude profile from MIGHTI-A
        var = _create_variable(ncfile, '%s_Fringe_Amplitude_A'%prefix, L22_dict['fringe_amp_A'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Fringe Amplitude from MIGHTI-A', 
                              display_type='image', field_name='Fringe Amplitude A', fill_value=None, label_axis='', bin_location=0.5,
                              units='', valid_min=-1e30, valid_max=1e30, var_type='metadata', chunk_sizes=[nx,ny],
                              notes="See Fringe_Amplitude. This variable contains the fringe amplitude measured "
                              "by MIGHTI-A, interpolated to the reconstruction grid. This is one of two variables "
                              "used to create Fringe_Amplitude. When A and B are significantly different, large "
                              "gradients are suspected, and the quality is reduced."
                              )
        
        # Fringe amplitude profile from MIGHTI-B
        var = _create_variable(ncfile, '%s_Fringe_Amplitude_B'%prefix, L22_dict['fringe_amp_B'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Fringe Amplitude from MIGHTI-B', 
                              display_type='image', field_name='Fringe Amplitude A', fill_value=None, label_axis='', bin_location=0.5,
                              units='', valid_min=-1e30, valid_max=1e30, var_type='metadata', chunk_sizes=[nx,ny],
                              notes="See Fringe_Amplitude. This variable contains the fringe amplitude measured "
                              "by MIGHTI-B, interpolated to the reconstruction grid. This is one of two variables "
                              "used to create Fringe_Amplitude. When A and B are significantly different, large "
                              "gradients are suspected, and the quality is reduced."
                              )
       

        # Fringe amplitude relative difference
        var = _create_variable(ncfile, '%s_Fringe_Amplitude_Relative_Difference'%prefix, L22_dict['fringe_amp_rel_diff'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc="Difference in MIGHTI A and B's amplitude estimates, divided by the mean", 
                              display_type='image', field_name='Fringe Amplitude Difference', fill_value=None, label_axis='', bin_location=0.5,
                              units='', valid_min=0.0, valid_max=1e10, var_type='metadata', chunk_sizes=[nx,ny],
                              notes="The absolute value of the difference between Fringe_Amplitude_A and Fringe_Amplitude_B, divided "
                              "by the average. Ideally, MIGHTI A and B should measure the same amplitude. When they do not, this is an "
                              "indication of potential violations of the spherical symmetry assumption inherent to the inversion. This "
                              "is the quantity used to determine if spherical asymmetry flag is raised."
                              )
        
        # Error flags      
        var = _create_variable(ncfile, '%s_Error_Flag'%prefix, np.transpose(L22_dict['error_flags'], (1,0,2)), 
                              dimensions=('Epoch', var_alt_name, 'N_Flags'), depend_0 = 'Epoch', depend_1 = var_alt_name, depend_2 = 'N_Flags',
                              format_nc='i8', format_fortran='I', desc='Error flags', 
                              display_type='image', field_name='Error Flags', fill_value=-1, label_axis='', bin_location=0.5,
                              units='', valid_min=0, valid_max=1, var_type='metadata', chunk_sizes=[nx,ny,8],
                              notes=["This variable provides information on why the Quality variable is reduced from 1.0. Eight error flags can "
                              "exist for each grid point, each either 0 or 1. More than one flag can be raised per point. This variable "
                              "is a three-dimensional array with dimensions time, altitude, and number of flags.",
                                    "    0 = missing MIGHTI A file",
                                    "    1 = missing MIGHTI B file",
                                    "    2 = A signal too weak",
                                    "    3 = B signal too weak",
                                    "    4 = A did not sample this altitude",
                                    "    5 = B did not sample this altitude",
                                    "    6 = Spherical asymmetry: A and B VER estimates disagree",
                                    "    7 = Unknown Error",
                                    ])
        
        ncfile.close()
        
    except: # Make sure the file is closed
        ncfile.close()
        raise
            
    return L22_full_fn




def level21_to_level22_without_info_file(L21_fns, L22_path, data_revision=0, sph_asym_thresh=None, time_start=None, time_stop=None):
    '''
    High-level function to apply the Level-2.1-to-Level-2.2 algorithm to a series of same-colored Level 2.1 files
    (in the L21_path directory) and create a single Level 2.2 file (in the L22_path directory). This version
    of the function requires the user to input parameters manually, instead of specifying an 
    Information.TXT file, like is done at the Science Data Center.
    
    INPUTS:
    
      *  L21_fns         -- TYPE:list of str.   A list of L2.1 files to be processed, including both MIGHTI-A and
                                                MIGHTI-B files, all from the same color.
      *  L22_path        -- TYPE:str.           The directory the L2.2 file will be saved in, including trailing "/"
                                                (e.g., '/home/user/')

    OPTIONAL INPUTS:
    
      *  data_revision   -- TYPE:int,      The minor version of the data [0-999]. The major version is set
                                           by the software's major version.
      *  sph_asym_thresh -- TYPE:float.    Relative difference in emission rate measured by A and B, beyond which
                                           the spherical-asymmetry flag will be raised. Technically, it should be
                                           "fringe amplitude" instead of "emission rate" due to the temperature
                                           dependence. Relative difference is defined as abs(A-B)/((A+B)/2). If 
                                           None (default), the default from MIGHTI_L2.global_params will be used.
      *  time_start      -- TYPE:datetime. A timezone-naive datetime in UT, specifying the beginning of the interval
                                           which the Level 2.2 data product should span. (Some L2.1 files from before
                                           the start time are needed to handle the crossover). If None (default), 
                                           the start time defaults to 0 UT on the date of the middle input file's time.
      *  time_stop       -- TYPE:datetime. A timezone-naive datetime in UT, specifying the end of the interval
                                           which the Level 2.2 data product should span. (Some L2.1 files from after
                                           the stop time are needed to handle the crossover). If None (default), 
                                           the stop time defaults to 24 hours after the start time. 
                                       
    OUTPUTS:
    
      *  L22_fn          -- TYPE:str.      The full path to the saved L2.2 file.
      
    '''
    
    assert len(L21_fns)>0, "No L2.1 files specified"
    

    ##### Load L2.1 files into dictionaries
    # Sort files by sensor
    Afns = [f for f in L21_fns if 'MIGHTI-A' in f]
    Bfns = [f for f in L21_fns if 'MIGHTI-B' in f]
    
    if len(Afns)==0:
        raise ValueError('No MIGHTI-A files found')
    if len(Bfns)==0:
        raise ValueError('No MIGHTI-B files found')

    # Sort the files by time (same as alphanumeric sorting)
    Afns.sort()
    Bfns.sort()

    # Read files (this is the bottleneck in terms of runtime)
    level21_A = level21_to_dict(Afns)
    level21_B = level21_to_dict(Bfns)

    ##### Run L2.2 processing to create L2.2 dictionary
    L22_dict = level21_dict_to_level22_dict(level21_A, level21_B, sph_asym_thresh=sph_asym_thresh,
                                                                  time_start=time_start, time_stop=time_stop)
    
    ##### Save L2.2 data to file
    L22_fn = save_nc_level22(L22_path, L22_dict, data_revision=data_revision)

    return L22_fn






def level21_to_level22(info_fn):
    '''
    Highest-level function to apply the Level-2.1-to-Level-2.2 algorithm. Inputs are specified via an information file.
    Files should be listed in the information file for a 24 hour period (0 UT to 0 UT), plus >15 minutes of files on either
    side. If files from green and red are both specified (as expected when run at the Science Data Center), they will be split
    and run separately by this function. The output Level 2.2 file(s) will be saved to the "Output" folder in the directory
    specified by the information file. Summary plots will also be created in that directory.
    
    INPUTS:
    
      * info_fn  -- TYPE:str.  Full path to an ASCII file in the following format:
      
                                        [PARAMETERS]
                                        Revision=001
                                        Directory=/path/to/wherever/
                                        <other parameters>

                                        [FILES]
                                        ICON_L2_MIGHTI-A_Line-of-Sight-Wind-Green_2017-05-27_000027_v01r001.NC
                                        ICON_L2_MIGHTI-A_Line-of-Sight-Wind-Green_2017-05-27_000057_v01r001.NC
                                        ICON_L2_MIGHTI-A_Line-of-Sight-Wind-Green_2017-05-27_000128_v01r001.NC
                                        etc... including files from MIGHTI-B
                                      
    OUTPUTS:
    
      *  ret     -- TYPE:str. '' if everything worked. If not, a human-readable error message for each file that failed
      
    '''
    
    info, L21_fns = read_info_file(info_fn)
    L21_fns.sort()
    
    # Parse the info
    # (0) Make sure there's a trailing "/" on the directory
    direc = info['Directory']
    if direc[-1] != '/':
        direc += '/'
    # (1) Add the directory to all the L2.1 files, which are in the Input folder
    L21_full_fns = []
    for L21_fn in L21_fns:
        L21_full_fns.append(direc + 'Input/' + L21_fn)
    # (2) Parse list of data revision numbers
    s = info['Revision'].split(',')
    data_revision = [int(x) for x in s]
    # For L2.2, we only expect a single revision number
    assert len(data_revision)==1, "Multiple revision numbers not supported for Level 2.2 processing"
    data_revision = data_revision[0]
    
    # For both red and green, call lower-level function which does all the real work
    L22_fns = []
    failure_messages = []
    for emission_color in ['red','green']:
        # Extract L2.1 files with this color
        L21_fns_color = [fn for fn in L21_full_fns if emission_color.lower() in fn.lower()]
        try:
            L22_fn = level21_to_level22_without_info_file(L21_fns_color, direc + 'Output/', data_revision=data_revision)       
            L22_fns.append(L22_fn)
        except Exception as e:
            failure_messages.append('Failed processing:\n\tColor   = %s\n%s\n'%(emission_color,traceback.format_exc()))
                
    # For both red and green, create summary plots
    for fn in L22_fns:
        try:
            # Create a 24-hour plot and 6 4-hour plots. This may change.
            plot_level22(fn, direc + 'Output/', nfigs=1)
            plot_level22(fn, direc + 'Output/', nfigs=6)
        except Exception as e:
            failure_messages.append('Failed creating summary plot for %s:\n%s\n'%(fn, traceback.format_exc()))
    
    if not failure_messages: # Everything worked
        return ''
    
    else:
        s = '\n' + '\n'.join(failure_messages)
        return s




    
def level22_to_dict(L22_fn):
    '''
    Read a Level 2.2 file and store its contents in a dictionary
    INPUTS:
    
      *  L22_fn          -- TYPE:str.   Full path to L2.2 file to be read
      
    OUTPUTS:
    
      *  L22_dict        -- TYPE:dict.  A dictionary containing the outputs of the Level 2.2 processing.
                                        See documentation for level21_dict_to_level22_dict(...) for details
                                        of the contents.
    '''
    
    L22_dict = {}
    f = netCDF4.Dataset(L22_fn)
    emission_color = L22_fn.split('/')[-1].split('_')[3].split('-')[-1].capitalize() # Red or Green

    L22_dict['lat']    = f['ICON_L2_MIGHTI_%s_Latitude'  % emission_color][:,:].T
    L22_dict['lon']    = f['ICON_L2_MIGHTI_%s_Longitude' % emission_color][:,:].T

    # Unwrap the sample longitude to avoid 0/360 jumps
    lonu = L22_dict['lon'].copy()
    lonu[0,:] = fix_longitudes(lonu[0,:], lonu[0,-1]) # Use top first longitude as target
    for i in range(np.shape(lonu)[1]):
        lonu[:,i] = fix_longitudes(lonu[:,i], lonu[0,i])
    L22_dict['lon_unwrapped'] = lonu

    L22_dict['alt']    = f['ICON_L2_MIGHTI_%s_Altitude'  % emission_color][:].T/1e3
    L22_dict['u']      = f['ICON_L2_MIGHTI_%s_Zonal_Wind' % emission_color][:,:].T
    L22_dict['v']      = f['ICON_L2_MIGHTI_%s_Meridional_Wind' % emission_color][:,:].T
    L22_dict['u_error']= f['ICON_L2_MIGHTI_%s_Zonal_Wind_Error' % emission_color][:,:].T
    L22_dict['v_error']= f['ICON_L2_MIGHTI_%s_Meridional_Wind_Error' % emission_color][:,:].T
    e                  = f['ICON_L2_MIGHTI_%s_Error_Flag' % emission_color][:,:,:]
    L22_dict['error_flags'] = np.transpose(e, (1,0,2))
    L22_dict['epoch_ms']      = f['Epoch'][:]
    L22_dict['epoch_full_ms'] = f['Epoch_Full'][:,:].T
    L22_dict['time_delta']    = f['ICON_L2_MIGHTI_%s_Time_Delta' % emission_color][:,:].T
    L22_dict['fringe_amp']    = f['ICON_L2_MIGHTI_%s_Fringe_Amplitude' % emission_color][:,:].T
    L22_dict['fringe_amp_A']  = f['ICON_L2_MIGHTI_%s_Fringe_Amplitude_A' % emission_color][:,:].T 
    L22_dict['fringe_amp_B']  = f['ICON_L2_MIGHTI_%s_Fringe_Amplitude_B' % emission_color][:,:].T
    L22_dict['fringe_amp_rel_diff'] = f['ICON_L2_MIGHTI_%s_Fringe_Amplitude_Relative_Difference' % emission_color][:,:].T
    L22_dict['wind_quality']  = f['ICON_L2_MIGHTI_%s_Wind_Quality' % emission_color][:,:].T
    L22_dict['emission_color'] = emission_color
    L22_dict['source_files'] = f.Parents
    L22_dict['acknowledgement'] = f.Acknowledgement
    L22_dict['time_start'] = datetime.strptime(f.Date_Start[-27:-4], '%Y-%m-%dT%H:%M:%S.%f')
    L22_dict['time_stop']  = datetime.strptime(f.Date_Stop [-27:-4], '%Y-%m-%dT%H:%M:%S.%f')
                                  
    
    # Convert times to datetime
    epoch = np.ma.empty(len(L22_dict['epoch_ms']), dtype=datetime)
    for i in range(len(epoch)):
        if isinstance(L22_dict['epoch_ms'],np.ma.MaskedArray) and L22_dict['epoch_ms'].mask[i]:
            epoch[i] = datetime(1989,4,18) # arbitrary fill value
            epoch[i] = np.ma.masked # mask it
        else:
            epoch[i] = datetime(1970,1,1) + timedelta(seconds = L22_dict['epoch_ms'][i]/1e3)
            
    epoch_full = np.ma.empty(np.shape(L22_dict['epoch_full_ms']), dtype=datetime)
    for i in range(np.shape(epoch_full)[0]):
        for j in range(np.shape(epoch_full)[1]):
            if isinstance(L22_dict['epoch_full_ms'],np.ma.MaskedArray) and L22_dict['epoch_full_ms'].mask[i,j]:
                epoch_full[i] = datetime(1989,4,18) # arbitrary fill value
                epoch_full[i] = np.ma.masked # mask it
            else:
                epoch_full[i] = datetime(1970,1,1) + timedelta(seconds = L22_dict['epoch_full_ms'][i,j]/1e3)
                
    L22_dict['epoch'] = epoch
    L22_dict['epoch_full'] = epoch_full
    
    f.close()
    
    return L22_dict




################################################################################################################
##########################################   TOHBAN PLOTS    ###################################################
################################################################################################################
# TODO for Tohban plots:
# - specify corners of pixels for pcolormesh instead of middle
# - break up by SLT instead of longitude or time
# - add map?
# - L2.1: identify gaps and mark as missing? Right now pcolormesh just interpolates over gaps
# - L2.2: What's the best x-axis for L2.2? Time seems good, but there are masked data points.
#         These could probably be interpolated over reasonably well. Perhaps this should be 
#         done as part of the L2.2 product?

def plot_level21(L21_fns, pngpath, timechunk=2., v_max = 200., ve_min = 1., ve_max = 100., amp_min = 1., amp_max = 1000.):
    '''
    Create Tohban plots for Level 2.1 data. Files from either MIGHTI-A or MIGHTI-B should
    be provided. The files will be broken into separate plots of length "timechunk" hours.
    
    INPUTS:
    
      *  L21_fns      --TYPE:list of str, Full paths to L2.1 files to be plotted. A or B.
      *  pngpath      --TYPE:str,         Directory to save the resulting pngs to
      
    OPTIONAL INPUTS:
    
      *  timechunk    --TYPE:float,       Time period to plot in a single file. [hr]
                                          (TODO: eventually this should be based on SLT).
      *  v_max        --TYPE:float,       Maximum wind velocity for colorbar [m/s]
      *  ve_min       --TYPE:float,       Minimum wind error for colorbar [m/s]
      *  ve_max       --TYPE:float,       Maximum wind error for colorbar [m/s]
      *  amp_min       --TYPE:float,      Minimum fringe amplitude for colorbar [m/s]
      *  amp_max       --TYPE:float,      Maximum fringe amplitude for colorbar [m/s]
                                   
    OUTPUT:
    
      *  L21_pngs     --TYPE:list of str,  Full path to the saved png files
    '''
    

    assert len(L21_fns)>0, "No files specified"

    if pngpath[-1] != '/':
        pngpath += '/'

    # Split up red/green and read
    L21_fns.sort()
    d = {} # keys are 'red' or 'green
    for color in ['red','green']:
        fns_short = [fn for fn in L21_fns if color in fn.lower()]
        if len(fns_short)>0:
            d[color] = level21_to_dict(fns_short)

    # Make nans into masked arrays so pcolormesh looks nicer.
    # (Not necessary if we upgrade to matplotlib 2.1.0)
    for color in ['red','green']:
        for v in ['los_wind','los_wind_error','amp']:
            d[color][v] = np.ma.masked_array(d[color][v], np.isnan(d[color][v]))

    # For each, add a variable giving the time as a matrix (to be used for plotting)
    for color in ['red','green']:
        if color in d.keys():
            dk = d[color]
            Nz,Nx = np.shape(dk['alt'])
            time_mat = np.empty((Nz,Nx),dtype=object)
            for j in range(Nx):
                time_mat[:,j] = dk['time'][j]
            dk['time_mat'] = time_mat    

    nx,ny = 2,3

    # Determine how many figures need to be made
    t_g = d['green']['time'] if 'green' in d.keys() else []
    t_r = d['red'  ]['time'] if 'red'   in d.keys() else []
    t = np.concatenate((t_g, t_r))
    tallmin = t.min()
    tallmax = t.max()
    if tallmax==tallmin: tallmax += timedelta(seconds=0.0001) # corner case: only one file
    nfigs = int(np.ceil((tallmax-tallmin).total_seconds()/3600./timechunk))

    # Determine if this is MIGHTI-A or MIGHTI-B
    first_sensor = L21_fns[0].split('/')[-1].split('_')[2] # MIGHTI-A or MIGHTI-B
    for fn in L21_fns:
        sensor = fn.split('/')[-1].split('_')[2]
        assert sensor == first_sensor, "Files from purely MIGHTI-A or purely MIGHTI-B must be specified."

    L21_pngs = []
    for n in range(nfigs):
        try:

            # Define time range for this plot
            tmin = tallmin + timedelta(hours=n*timechunk)
            tmax = tmin + timedelta(hours=timechunk)


            fig = plt.figure(figsize=(16,10))
            tmin_title = datetime(2100,1,1) # for specifying the title of each plot (updated below)
            tmax_title = datetime(1970,1,1)
            ################## Plot ###################
            colors = ['red', 'green']
            cmaps = ['inferno', 'viridis']
            for i,color in enumerate(colors):
                if color in d.keys():
                    dk = d[color]
                    i1 = bisect.bisect(dk['time'],tmin) - 1
                    i2 = bisect.bisect(dk['time'],tmax)
                    if dk['time'][i1] < tmin_title:
                        tmin_title = dk['time'][i1]
                    if dk['time'][i2-1] > tmax_title:
                        tmax_title = dk['time'][i2-1]

                    plt.subplot(ny,nx,i+1)
                    plt.title('%s LoS Wind'% (color.capitalize()))
                    try: # These are wrapped in try blocks so that if something unexpected happens, a plot is still created.
                        plt.pcolormesh(dk['time_mat'][:,i1:i2], dk['alt'][:,i1:i2], dk['los_wind'][:,i1:i2],
                                       vmin=-v_max, vmax=v_max, cmap='bwr')
                        cb = plt.colorbar()
                        cb.set_label('m/s',rotation=-90,labelpad=10)
                    except:
                        print 'Failed to produce wind plot:'
                        print traceback.format_exc()
                        pass

                    plt.subplot(ny,nx,i+3)
                    plt.title('%s LoS Wind Error'% (color.capitalize()))
                    try:
                        plt.pcolormesh(dk['time_mat'][:,i1:i2], dk['alt'][:,i1:i2],dk['los_wind_error'][:,i1:i2],
                                       norm=LogNorm(vmin=ve_min, vmax=ve_max), cmap=cmaps[i]+'_r')
                        cb = plt.colorbar()
                        cb.set_label('m/s',rotation=-90,labelpad=10)
                    except:
                        print 'Failed to produce wind error plot:'
                        print traceback.format_exc()
                        pass

                    plt.subplot(ny,nx,i+5)
                    plt.title('%s Fringe Amplitude'% (color.capitalize()))
                    try:
                        amp = dk['amp']
                        amp[amp<amp_min] = amp_min # so it's plotted as dark instead of missing
                        plt.pcolormesh(dk['time_mat'][:,i1:i2], dk['alt'][:,i1:i2], amp[:,i1:i2],
                                       norm=LogNorm(vmin=amp_min, vmax=amp_max), cmap=cmaps[i])
                        cb = plt.colorbar()
                        cb.set_label('arb',rotation=-90,labelpad=10)
                    except:
                        print 'Failed to produce amplitude plot:'
                        print traceback.format_exc()
                        pass

                    ###### Prettify the axes, labels, etc #####
                    for sp in [i+1,i+3,i+5]:
                        plt.subplot(ny,nx,sp)
                        plt.gca().patch.set_facecolor('gray')
                        plt.gca().xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
                        plt.ylabel('Altitude [km]')
                        plt.gca().yaxis.set_ticks([100,200,300])
                        plt.ylim((85,320))
                        plt.yscale('log')
                        plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
                        plt.gca().yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
                        plt.xlim((tmin,tmax))

            # Determine min and max time for the title
            if tmin_title.date() == tmax_title.date():
                datestr = '%s - %s' % (tmin_title.strftime('%Y-%m-%d %H:%M:%S'), tmax_title.strftime('%H:%M:%S'))
            else:
                datestr = '%s - %s' % (tmin_title.strftime('%Y-%m-%d %H:%M:%S'), tmax_title.strftime('%Y-%m-%d %H:%M:%S'))
            fig.suptitle('%s  %s' % (sensor, datestr), fontsize=28)

            fig.autofmt_xdate()
            fig.tight_layout(rect=[0,0.03,1,0.95]) # leave room for suptitle

            #### Save
            pngfn = pngpath + 'ICON_L2_%s_Line-of-Sight-Wind_%s_to_%s.png' % (sensor, tmin_title.strftime('%Y-%m-%d_%H%M%S'),
                                                                                      tmax_title.strftime('%Y-%m-%d_%H%M%S'))
            plt.savefig(pngfn, dpi=70)
            plt.close()
            L21_pngs.append(pngfn)
        except:
            # Should we catch this or let it fly? If a plot was not created, that's a failure.
            #print 'Failed to produce plot %i (out of %i):' % (n+1,nfigs)
            #print traceback.format_exc()
            raise # This should never happen, so let it crash

    return L21_pngs



def plot_level22(L22_fn, pngpath, nfigs=6, v_max = 200., ve_min = 1., ve_max = 100.):
    '''
    Create Tohban plots for Level 2.2 data.
    
    INPUTS:
    
      *  L22_fn       --TYPE:str,   Full path to a L2.2 file
      *  pngpath      --TYPE:str,   Directory to save the resulting png to
   
   OPTIONAL INPUTS:
   
      *  nfigs        --TYPE:int,   How many figures to break the data into
                                    (TODO: eventually this should be based on SLT).
      *  v_max        --TYPE:float, Maximum wind velocity for colorbar [m/s]
      *  ve_min       --TYPE:float, Minimum wind error for colorbar [m/s]
      *  ve_max       --TYPE:float, Maximum wind error for colorbar [m/s]
                                   
    OUTPUT:
    
      *  L22_pngs     --TYPE:list of str,  Full path to the saved png files
      
    '''
    
    if pngpath[-1] != '/':
        pngpath += '/'
    
    L22_dict = level22_to_dict(L22_fn)

    # Error flag meanings:
    error_str = {
        0: 'Error: missing MIGHTI-A file',
        1: 'Error: missing MIGHTI-B file',
        2: 'Error: MIGHTI-A weak signal',
        3: 'Error: MIGHTI-B weak signal',
        4: 'Error: MIGHTI-A missing this altitude',
        5: 'Error: MIGHTI-B missing this altitude',
        6: 'Error: Spherical asymmetry detected',
        7: 'Error: Unknown' 
    }

    Nalt, Nlon, Nflags = np.shape(L22_dict['error_flags'])
    lonu = L22_dict['lon_unwrapped']
    alt = L22_dict['alt']
    alt_mat = np.tile(alt, (Nlon,1)).T
    L22_pngs = []
    
    for n in range(nfigs):

        # Define grid range for this plot
        i1 = int(np.floor(1.0*Nlon/nfigs*n))
        i2 = int(np.floor(1.0*Nlon/nfigs*(n+1)))

        fig = plt.figure(figsize=(16,(Nflags/2+2)*3))
        
        #### WINDS
        for i,W in enumerate([L22_dict['u'],L22_dict['v']]):
            plt.subplot(Nflags/2 + 2, 2, i+1)
            plt.pcolormesh(lonu   [:,i1:i2], 
                           alt_mat[:,i1:i2],
                           W      [:,i1:i2],
                           cmap='bwr')
            plt.clim((-v_max,v_max))
            if i==0:
                plt.title('Zonal Wind')
            else:
                plt.title('Meridional Wind')
            cb = plt.colorbar()
            cb.set_label('Wind [m/s]',rotation=-90,labelpad=15)

        #### WIND ERROR
        for i,W in enumerate([L22_dict['u_error'],L22_dict['v_error']]):
            plt.subplot(Nflags/2 + 2, 2, i+3)
            plt.pcolormesh(lonu   [:,i1:i2], 
                           alt_mat[:,i1:i2],
                           W      [:,i1:i2],
                           cmap='viridis_r', norm=LogNorm(vmin=ve_min, vmax=ve_max))
            if i==0:
                plt.title('Zonal Wind Error')
            else:
                plt.title('Meridional Wind Error')
            cb = plt.colorbar()
            cb.set_label('Wind Error [m/s]',rotation=-90,labelpad=15)
    
        #### ERROR FLAGS
        for i in range(Nflags):
            plt.subplot(Nflags/2 + 2,2,i+5)
            # define the bins and normalize
            cmap=cmap_custom
            bounds = [-0.5,0.5,1.5]
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)            
            e = L22_dict['error_flags'][:,i1:i2,i].astype(float)
            # For spherical symmetry flag, mask it to look like data so we can see where
            # this flag actually matters
            if i==6 and isinstance(L22_dict['u'],np.ma.MaskedArray):
                e = np.ma.masked_array(e, L22_dict['u'].mask[:,i1:i2])
                
            plt.pcolormesh(lonu   [:,i1:i2],
                           alt_mat[:,i1:i2],
                           e,
                           cmap=cmap, norm=norm)
            plt.title(error_str[i])
            plt.axis('tight')
            cb = plt.colorbar(ticks=[0,1])

        #### Make plots look nice
        for sp in range(Nflags+4):
            plt.subplot(Nflags/2 + 2,2, sp+1)
            plt.gca().patch.set_facecolor('gray')
            plt.ylabel('Altitude [km]')
            plt.xlabel('Longitude [deg]')
            plt.gca().yaxis.set_ticks([100,200,300])
            plt.ylim((85,320))
            plt.yscale('log')
            plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
            plt.gca().yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
            # Replace unwrapped longitude with mod 360
            xt = plt.gca().get_xticks()
            plt.gca().set_xticklabels(['%.0f' % np.mod(x,360) for x in xt])

        plt.tight_layout(rect=[0,0.02,1,0.97]) # leave room for suptitle
        datestr = L22_fn.split('/')[-1].split('_')[-2]
        fig.suptitle('%s %s Plot %i/%i' % (datestr, L22_dict['emission_color'], n+1,nfigs), fontsize=28)
        
        #### Save
        pngfn = pngpath + L22_fn.split('/')[-1][:-3] + '_p%02i_%02i.png'%(n+1,nfigs)
        plt.savefig(pngfn, dpi=70)
        plt.close()
        L22_pngs.append(pngfn)
        
    return L22_pngs
    
    
################################################################################################################
##############################   INTERFACE FOR SCIENCE DATA CENTER:    #########################################
################################################################################################################
    
def main(argv):
    '''
    Program that will run if this file is called as a script. Depending on the 
    argument, the Level 2.1 or Level 2.2 algorithms will be launched, assuming that the Information.TXT file
    is in the directory. Examples:
    
    $ python /path/to/MIGHTI_L2.py -L2.1
    $ python /path/to/MIGHTI_L2.py -L2.2
    
    The script will return with exit code 0 on success, exit code 1 on failure. If exit code 1, that usually
    means there's a bug in the code or something completely unexpected in the data.
    
    INPUT:
        argv --TYPE:str. The command line arguments, provided by sys.argv.
        
    '''
    info_fn = './Input/Information.TXT'
    
    usagestr = 'usage: python MIGHTI_L2.py -L2.x  (where x is 1 or 2)'
    if len(argv)!=2:
        print usagestr
        sys.exit(1)
        
    plt.switch_backend('Agg') # to generate pngs without opening windows
        
    # Call the L2.1 or L2.2 algorithms, as specified.
    # Also possibly call the unit test functions
    err = ''
    if argv[1] == '-L2.1':
        err = level1_to_level21(info_fn)
    elif argv[1] == '-L2.2':
        err = level21_to_level22(info_fn)
    else:
        print usagestr
        sys.exit(1)
    
    # If the processing had a problem, print it and exit with error
    if err:
        print err
        sys.exit(1)
    else:
        sys.exit(0)
    
    
if __name__ == "__main__":
    main(sys.argv)
    
    
    
    



    
    
    
    
    
    
    
    
    
    
