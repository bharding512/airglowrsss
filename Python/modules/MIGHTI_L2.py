# A module for the conversion of MIGHTI Level 1 files to Level 1.5, 2.1, and 2.2 files.
# Level 1 files - Calibrated MIGHTI interferograms
# Level 1.5 files - Consolidated MIGHTI phases
# Level 2.1 files - Line-of-sight wind profiles (this is where the onion-peeling inversion happens)
# Level 2.2 files - Vector wind profiles (this is where the A/B matchup happens)
# Altitudes and distances are expressed in km everywhere in the code


####################################### VERSION CONTROL ############################################
# These need to be manually changed, when necessary.
# NOTE: When the major version is updated, you should change the History global attribute
# in both the L2.1 and L2.2 netcdf files, to describe the change (if that's still the convention)
software_version_major = 6 # Should only be incremented on major changes
software_version_minor = 0 # [0-99], increment on ALL published changes, resetting when the major version changes
__version__ = '%i.%02i' % (software_version_major, software_version_minor) # e.g., 2.03
####################################################################################################


############################## GLOBAL PROCESSING PARAMETERS ########################################
# Unless overridden by the user, the following parameters will be used. Note that there are different
# parameters for red and green. Sometimes there are different parameters for the two MIGHTI sensors: A and B.
global_params = {}
global_params['red'] = {
    ###################### Inversion parameters #######################
    'sigma'             : 1.0/630.0304e-9, # reciprocal of center wavelength of emission [m^-1] 
                                           # (from John Harlander)  
    'bin_size'          : 4,               # The number of rows of the interferogram to bin together to 
                                           # improve statistics at the cost of altitude resolution.   
    'account_for_local_projection': True,  # Whether to account for the fact that the line of sight is not
                                           # quite horizontal everywhere along the line of sight
    'integration_order' : 0,               # 0: Use Riemann-sum rule for discretizing line-of-sight integral
                                           # 1: Use trapezoidal rule for discretizing line-of-sight integral
    'top_layer'         : 'exp',           # 'thin': assume VER goes to zero above top layer
                                           # 'exp':  assume VER falls off exponentially in altitude
    'H'                 : 40.,              # km. The VER scale height used when top_layer='exp'.
                                           # This was found by fitting many profiles for which there was significant
                                           # emission above 300 km. Profiles were generated from the Zhang/Shepherd model and
                                           # from photochemical models fed by IRI/MSIS. (See Harding et al. [2017] SSR paper 
                                           # for details on airglow models).
                                            
    
    ###################### Quality control parameters #################
    'top_layer_thresh'  : 0.4,            # Fraction of airglow above the top observed altitude. Consider the total column 
                                           # brightness (i.e, the integral of the VER profile). When a large fraction
                                           # comes from above the top observed altitude, the quality flag is raised.
                                           # This threshold is specified as a fraction of the total column brightness.
    'terminator_thresh' : 1000.,           # [km]. Consider two points along the line of sight, both X km from the tangent 
                                           # point (where X is this parameter). One is nearer to the spacecraft than the 
                                           # tangent point, and one is farther. If these two points are on opposite sides of the 
                                           # terminator, raise a quality flag. Note that this quality flag will also be raised if
                                           # any observations at higher tangent altitudes for the same observation are flagged, because
                                           # the inversion mixes information from those observations.
                                           # DEPRECATED AS OF v5.01
    'sza_term_thresh'   : 5.,              # Mask data within X degrees of 98 degrees solar zenith angle.
    'sph_asym_thresh'   : 0.4,             # The relative difference in VER estimates from A&B, beyond which the spherical asymmetry 
                                           # flag will be raised in L2.2.
    't_diff_AB'         : 20*60.,          # [sec] Maximum time difference between A and B exposures used to determine
                                           # the wind at a given location. If it's longer than this, something went wrong
                                           # with the processing, most likely.
    'post_maneuver_window': 0.5,           # [hours] How long to wait after maneuvers before calling the thermal drift good. The
                                           # thermal environment during a maneuver is different than during the calibration orbit,
                                           # so the thermal drift correction might be wrong. It is thought that this effect is 
                                           # probably small, but it's unknown, so it will be flagged.
    'jitter_thresh'       : 0.01,          # [deg]. If the jitter is larger than this, then the S/C velocity subtraction is too
                                           # uncertain, and the data are marked as "caution."
                                           # DEPRECATED
    'chi2_thresh'         : 0.20,          # [rad^2]. If the mean-square phase across a row is more than this, then this wind
                                           # measurement is definitely bad, and should be set to quality=0
    'chi2_thresh_caution' : 0.13,          # [rad^2]. If the mean-square phase across a row is more than this, then this wind
                                           # measurement should be treated with caution (quality=0.5).
    'count_thresh_day'    : 5,             # If a profile doesn't have this many samples labeled "good" or "caution", the 
                                           # entire profile is bad due to an invalid inversion. (Adjust this when binning changes).
    'count_thresh_night'  : 5,             
}

global_params['green'] = { # See above for descriptions
    'sigma'             : 1.0/557.7339e-9,
    'bin_size'          : 1,
    'account_for_local_projection': True,
    'integration_order' : 0,
    'top_layer'         : 'exp',
    'H'                 : 40.,
    'top_layer_thresh'  : 0.4,
    'terminator_thresh' : 1000.,
    'sza_term_thresh'   : 5.,
    'sph_asym_thresh'   : 0.4,
    't_diff_AB'         : 20*60.,
    'post_maneuver_window': 0.5,
    'jitter_thresh'       : 0.01,
    'chi2_thresh'         : 0.6,
    'chi2_thresh_caution' : 0.27,
    'count_thresh_day'    : 5,
    'count_thresh_night'  : 5,
}

global_params['verbose'] = False # For debugging memory issues in the SDC.
                                 # Set to True for detailed logging.


#####################################################################################################


import numpy as np
import ICON
import bisect
from scipy import integrate, optimize, interpolate
from datetime import datetime, timedelta
import netCDF4
import getpass # for determining who is running the script
import glob
import traceback # for printing detailed error traces
import sys
import copy
import xarray as xr # Added 2021 May 27 for easier handling of Level 1.5 files.
import pandas as pd # Added 2021 May 27 for easier handling of Level 1.5 files.

import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.colors import LogNorm



# Added in v1.31: pysatMagVect was renamed to OMMBV as of v0.5.1. MIGHTI code should accept versions 0.4.0 - 0.5.1
# (Used for reporting winds in magnetic coordinates)
try:
    import OMMBV as psmv
    print('Using OMMBV, version %s' % (psmv.__version__))
except:
    try:
        import pysatMagVect as psmv
        print('Using pysatMagVect, version %s' % (psmv.__version__))
    except:
        print('WARNING: Neither OMMBV nor pysatMagVect were found')


# Added in v1.20
from pyglow import pyglow # for correcting VER for temperature visibility reduction

# Some aesthetics for plots
matplotlib.rcParams['xtick.labelsize'] = 'small'
matplotlib.rcParams['ytick.labelsize'] = 'small'

# Ignore errors having to do with NaN. These clog up the log file.
np.seterr(invalid='ignore')


# TODO: REMOVE. Used for debugging
def mem_usage():
    '''
    Return a string that describes the memory usage.
    '''
    import psutil
    vm = psutil.virtual_memory()
    s = 'T-%.2fGB A-%.2fGB %.1f%%' % (vm.total/2.**30, vm.available/2.**30, vm.percent)
    return s

def timestamp():
    '''
    Return timestamp with optional mem usage string
    '''
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' ' + mem_usage()


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



def visibility_temperature(t, lat, lon, alt, f107=None, f107a=None, f107p=None, apmsis=None):
    '''
    Return the temperature to be used to compute the fringe visibility reduction. 
    Currently implemented with MSIS.
    
    INPUTS:
      *  t    -- TYPE:datetime.          Time of requested point.
      *  lat  -- TYPE:float,  UNITS:deg. Latitude of requested point.
      *  lon  -- TYPE:float,  UNITS:deg. Longitude (0-360) of requested point.
      *  alt  -- TYPE:float,  UNITS:km.  Altitude of requested point.
      *  f107 -- TYPE:float,  UNITS:sfu. F10.7 value for the date of interest (None is use values from pyglow)
      *  f107a-- TYPE:float,  UNITS:sfu. F10.7a value for the date of interest (None is use values from pyglow)
      *  f107p-- TYPE:float,  UNITS:sfu. F10.7p value for the date of interest (None is use values from pyglow)
      *  apmsis--TYPE:array(7).          Ap vector that MSIS expects (None is use values from pyglow)

    OUTPUTS:
      *  T    -- TYPE:float,  UNITS:K.   Temperature at requested point.
    '''
    
    
    if apmsis is None:
        pt = pyglow.Point(t, lat, lon, alt)
    else:
        pt = pyglow.Point(t, lat, lon, alt, user_ind=True)
        pt.f107 = f107
        pt.f107a = f107a
        pt.f107p = f107p
        pt.apmsis = apmsis
    pt.run_msis()
    T = pt.Tn_msis
    return T
    
    

def amp_to_R_factor(sigma_opd,  T):
    '''
    Return the value g which satisfies I = g*A, where A is the fringe amplitude and I is the 
    brightness in Rayleigh. The fringe amplitude should be considered on a per-pixel basis, or
    equivalently, be considered as a mean over a row of the CCD. This was derived from first 
    principles and guided by the equations in Englert et al. (2007) doi:10.1364/AO.46.007297.
    The first principles calculation is adjusted by the multiplicative CAL_FACTOR to adjust for
    information gained after on-orbit calibration efforts.
    
    UPDATE: As of May 2019, many of these corrections have been moved to the L1 process. The
    only one that remains at L2 is the temperature correction (which must be done after the
    inversion).
    
    INPUTS:
    
      *  sigma_opd       -- TYPE:float.                  sigma times opd: Optical path difference measured 
                                                         in wavelengths. If analyzing an entire row at once, 
                                                         the mean OPD of that row should be used.
      *  T               -- TYPE:float, UNITS:K.         The temperature for visibility correction.
      
      
    OUTPUTS:
      *  g               -- TYPE:float, UNITS:R/counts.  Amplitude-to-Rayleigh factor, described above.
    '''
    
    CAL_FACTOR = 1.00 # this will be informed by on-orbit calibrations

    #### Temperature visibility correction term
    # Physical constants
    c = 299792458.            # The speed of light [m/s]
    k = 1.3806503e-23         # Boltzmann's constant [J/K]
    mO = 16 * 1.66053886e-27  # mass of oxygen [kg]
    T_factor = np.exp(-2*np.pi**2 * sigma_opd**2 * k*T/c**2/mO)

    #### Put it all together
    # The 1e10 factor is from the definition of Rayleigh (1e10 photons) and so is the 1/(4pi) factor (per column, 
    # while etendue is given in units of steradians)
    R_to_counts = CAL_FACTOR * T_factor 
    return 1.0/R_to_counts




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




def bin_array(b, y, lon = False, method='mean'):
    '''
    Downsample y by binning it, improving statistics. Every b
    elements of y will be averaged together to create a new array, y_b, 
    of length ny_b = ceil(len(y)/b). Binning starts at the end of the array, 
    so the first element of y_b may not represent exactly b samples of y.
    
    INPUTS:
    
      *  b       -- TYPE:int,       The number of rows to bin together
      *  y       -- TYPE:array(ny), The array to be binned
      
    OPTIONAL INPUTS:
    
      *  lon     -- TYPE:bool,      If True, 360-deg discontinuities will
                                    be removed before averaging (e.g., for
                                    longitude binning). (Default False)
      *  method  -- TYPE:str,       'mean': average over elements in a bin (default)
                                    'max': take maximum over elements (e.g., for quality flags)
                                    'min': take minimum over elements (e.g., for quality factor)
                                    
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
        else:
            if method == 'mean':
                y_b[i_new] = np.nanmean(y_samps)
            elif method == 'max':
                y_b[i_new] = np.nanmax(y_samps)
            elif method == 'min':
                y_b[i_new] = np.nanmin(y_samps)
            else:
                raise Exception('Input method="%s" not recognized. Try "mean", "max", or "min"' % method)
        
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

        n = sum(~np.isnan(ye_samps))
        if n==0:
            ye_b[i_new] = np.nan
        else:
            ye_b[i_new] = 1.0/sum(~np.isnan(ye_samps)) * np.sqrt(np.nansum(ye_samps**2))
        
    return ye_b
    
    
    
    
def bin_image(b, I, method='mean'):
    '''
    Downsample the interferogram in altitude to improve statistics while
    degrading vertical resolution. Every b rows will be averaged together. 
    Binning starts at high altitudes, so the lower rows of I_b may not represent 
    exactly b rows of I.
    
    INPUTS:
    
      *  b           -- TYPE:int,                        The number of rows to bin together
      *  I           -- TYPE:array(ny,nx),   UNITS:arb.  The MIGHTI interferogram
      
    OPTIONAL INPUTS:
    
      *  method  -- TYPE:str,       'mean': average over elements in a bin (default)
                                    'max': take maximum over elements (e.g., for quality flags)
                                    'min': take minimum over elements (e.g., for quality factor)

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
        I_b[:,i] = bin_array(b,I[:,i], method=method)
    return I_b





def create_observation_matrix(tang_alt, icon_alt, top_layer='exp', H=26., integration_order=0,):
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
    
      *  top_layer         -- TYPE:str,   'thin': assume VER goes to zero above top layer
                                          'exp':  assume VER falls off exponentially in altitude (default)
      *  H                 -- TYPE:float, UNITS:km. The VER scale height to use if top_layer='exp' (default 26)
      *  integration_order -- TYPE:int,
                                       * 0: Use Riemann-sum rule for discretizing line-of-sight integral (default).
                                       * 1: Use trapezoidal rule for discretizing line-of-sight integral.
                                          
    OUTPUTS:
    
      *  D          -- TYPE:array(ny,ny), UNITS:km/10.   Observation matrix. Also called the "path matrix"
                                                         or "distance matrix". The factor of 10 is to account
                                                         for the conversion from VER [ph/cm^3/s] to Rayleighs.
                                                      
    '''
    
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
    
    # Scale factor, so that when you invert a vector in Rayleigh you get a result in ph/cm^3/s
    D = D * 1e-1 # convert from km to cm (*1e5) then *1e-6 accounts for "10^10 ph/m^2" in R definition.
    
    return D






def create_local_projection_matrix(tang_alt, icon_alt):
    '''
    Define the matrix B whose entries give the factor by which a horizontal wind
    would be projected onto the line of sight. This has the same shape as the
    observation matrix (i.e., distance matrix). At the tangent point, this factor
    is 1.0. Far from the tangent point, this factor is smaller. If this effect is 
    accounted for, it makes a small change in the winds (less than 5 m/s), mostly on
    the bottomside. 
    
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
     

    
def unwrap(x, start=0):
    '''
    Unwrap a monotonically increasing phase signal to remove -2*pi jumps.
    This is very similar to np.unwrap, but only unwraps negative jumps,
    and it works with nans.
    
    INPUTS:
    
      *  x     -- TYPE:array, UNITS:rad. Signal that has -2*pi jumps to remove
      
    OPTIONAL INPUTS:
    
      * start  -- TYPE:int.              The pixel at which to begin unwrapping
      
    OUTPUTS:
    
      *  xnew  -- TYPE:array, UNITS:rad. Copy of x with -2*pi jumps removed
      
    '''

    assert(isinstance(start,int)), "Input 'start' must be an integer"

    xnew = np.zeros(len(x))
    xnew[start] = x[start]

    # Go from start forwards
    dx = np.diff(x[start:])
    #idx = dx < -np.pi # This throws a warning for nans, so use a more complicated expression:
    idx = np.array([dxi < -np.pi if not np.isnan(dxi) else False for dxi in dx], dtype=bool)
    dx[idx] = dx[idx] + 2.*np.pi
    xnew[start+1:] = xnew[start] + np.cumsum(dx)

    # Go from start backwards
    dx = np.diff(x[start::-1])
    #idx = dx > np.pi # This throws a warning for nans, so use a more complicated expression:
    idx = np.array([dxi > np.pi if not np.isnan(dxi) else False for dxi in dx], dtype=bool)
    dx[idx] = dx[idx] - 2.*np.pi
    xnew[:start] = xnew[start] + np.cumsum(dx)[::-1]

    return xnew
    
    
    
def analyze_row(row, unwrap_phase=False):
    '''
    Given a 1-D interference pattern (i.e., a row of the complex intererogram), 
    analyze it to get a scalar phase value, which represents the wind, and a
    single amplitude value, which is roughly proportional to the emission rate.
    
    INPUTS:
    
      *  row               -- TYPE:array(nx), UNITS:arb.   A row of the complex-valued, MIGHTI interferogram.
      
    OPTIONAL INPUTS:
    
      *  unwrap_phase      -- TYPE:bool.                   If True, remove 2pi discontinuities before analyzing
                                                           (default False)
    
    OUTPUTS:
    
      *  phase             -- TYPE:float,     UNITS:rad.   A scalar phase which represents the observed wind.
      *  amp               -- TYPE:float,     UNITS:arb.   A scalar amplitude which represents the brightness.
      *  chi2              -- TYPE:float,     UNITS:rad^2. Variance of the row of phase.
       
    '''
    if all(np.isnan(row)):
        return np.nan, np.nan, np.nan
    
    row_phase = np.angle(row)
    
    if unwrap_phase: # BJH added 3 Feb 2020. This should (probably) never be used on finalized public data.
        nx = len(row)
        row_phase =  unwrap(row_phase, start=int(nx/2)) # Remove negative jumps
        row_phase = -unwrap(-row_phase, start=int(nx/2)) # Remove positive jumps
    
    tot_phase = np.nanmean(row_phase)
    
    resid = row_phase - tot_phase
    chi2 = np.nanmean(resid**2)
    
    # Evaluate total amplitude
    # "Mean" is used to be resistant to rows with missing columns.
    amp = np.nanmean(abs(row))
        
    return tot_phase, amp, chi2

    
    
    
    
def perform_inversion(I, tang_alt, icon_alt, ph_err_inst, ph_err_striation, ph_err_zero_phase, ph_err_signal, I_amp_uncertainty,
                      top_layer='exp', H=26., integration_order=0, account_for_local_projection=True,
                      linear_amp = True, unwrap_phase = False, chi2_thresh = np.inf, Ncorr = 7):
    '''
    Perform the onion-peeling inversion on the interferogram to return
    a new interferogram, whose rows refer to specific altitudes. In effect,
    this function undoes the integration along the line of sight.
    
    INPUTS:
    
      *  I           -- TYPE:array(ny,nx), UNITS:arb.  The complex-valued, MIGHTI interferogram.
      *  tang_alt    -- TYPE:array(ny),    UNITS:km.   Tangent altitudes of each row of interferogram.
      *  icon_alt    -- TYPE:float,        UNITS:km.   Altitude of the satellite.
      *  ph_err_inst -- TYPE:array(ny),    UNITS:rad. Instantaneous uncertainty in the unwrapped, mean phase of each row of I.
                                                           This is provided in L1 file.
      *  ph_err_striation  -- TYPE:array(ny),    UNITS:rad. Striation uncertainty
      *  ph_err_zero_phase -- TYPE:array(ny),    UNITS:rad. Zero-wind phase uncertainty
      *  ph_err_signal     -- TYPE:array(ny),    UNITS:rad. Phase uncertainty caused by signal-dependent phase shift
      *  I_amp_uncertainty -- TYPE:array(ny),    UNITS:arb. Uncertainty in the summed amplitude of each row of I.
                                                            This is provided in L1 file.
      
    OPTIONAL INPUTS:
    
      *  top_layer   -- TYPE:str,          'thin': assume VER goes to zero above top layer
                                           'exp':  assume VER falls off exponentially in altitude (default)
      *  H                 -- TYPE:float, UNITS:km. The VER scale height to use if top_layer='exp' (default 26)
      *  integration_order -- TYPE:int,
      
                      * 0: Use Riemann-sum rule for discretizing line-of-sight integral (default).
                      * 1: Use trapezoidal rule for discretizing line-of-sight integral.
                      
      *  account_for_local_projection   -- TYPE:bool.   If False, a simple inversion is used.
                                           If True, the inversion accounts for the fact that the ray is not 
                                           perfectly tangent to each shell at each point along the ray. 
                                           (default True)
      *  linear_amp        -- TYPE:bool.   If True, a linear inversion is used on fringe amplitude. If False, the 
                                           fringe amplitude is estimated during the onion-peeling. These are the
                                           same in the absence of noise, but may make a difference at low SNR,
                                           especially for computing uncertainties. (default True)
      *  H                 -- TYPE:float, UNITS:km. The VER scale height to use if top_layer='exp' (default 26)
      *  unwrap_phase      -- TYPE:bool,   If True, unwrap the onion-peeled phase to remove 2pi discontinuities (default False)
      *  chi2_thresh       -- TYPE:float, UNITS:rad^2. If the mean-square phase residual is larger than this, treat the row as 
                                           pure noise (i.e., ignore it and set it to 0)
      *  Ncorr             -- TYPE:int,   UNITS:pixels. The row lag at which the residual striation correlation drops to 0. 
                                          (default: 7, which was determined by Monte Carlo simulations for the un-binned case.)
                                           
                                               
    OUTPUTS:
    
      *  Ip                -- TYPE:array(ny,nx), UNITS:arb. The complex-valued, onion-peeled interferogram.
      *  phase             -- TYPE:array(ny),    UNITS:rad. The mean phase of each row of Ip.
      *  amp               -- TYPE:array(ny),    UNITS:arb. The amplitude of each row of Ip.
      *  phase_err_inst    -- TYPE:array(ny),    UNITS:rad. The instantaneous error of phase
      *  phase_err_striation  -- TYPE:array(ny), UNITS:rad. The striation error of phase (24h time scales)
      *  phase_err_zero_phase -- TYPE:array(ny), UNITS:rad. The zero-phase error of phase (long time scales)
      *  phase_err_signal     -- TYPE:array(ny), UNITS:rad. Phase uncertainty caused by signal-dependent phase shift
      *  amp_uncertainty   -- TYPE:array(ny),    UNITS:arb. The uncertainty of amp.
      
    '''
    
    if top_layer not in ['exp','thin']:
        raise ValueError('Argument top_layer=\'%s\' not recognized. Use \'exp\' or \'thin\'.' % top_layer)
    if integration_order not in [0,1]:
        raise ValueError('Argument integration_order=\'%s\' not recognized. Use 0 or 1')
    
    ny,nx = np.shape(I)
    
    # Create the path matrix
    if global_params['verbose']:
        print('%s:\t\t\t\t Creating observation matrix '% (timestamp()))
        import sys
        sys.stdout.flush()
    D = create_observation_matrix(tang_alt, icon_alt, top_layer=top_layer, integration_order=integration_order, H=H)
    
    # Create local horizontal projection matrix (and set it to unity if we are to ignore this effect)
    if global_params['verbose']:
        print('%s:\t\t\t\t Creating horz proj matrix '% (timestamp()))
        import sys
        sys.stdout.flush()
    B = np.ones((ny,ny))
    if account_for_local_projection:
        B = create_local_projection_matrix(tang_alt, icon_alt)
            
    ######### Onion-peeling inversion and phase extraction #########
    phase = np.zeros(ny) # phases at each altitude
    amp   = np.zeros(ny) # fringe amplitude at each altitude
    chi2  = np.zeros(ny) # variance of each row's phase
    Ip = np.zeros((ny,nx), dtype=complex) # onion-peeled interferogram

    # This code implements Eq (9) in the MIGHTI L2 Space Science Reviews
    # paper (Harding et al. 2017).

    if global_params['verbose']:
        print('%s:\t\t\t\t Beginning onion-peeling '% (timestamp()))
        import sys
        sys.stdout.flush()
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
        p,a,c = analyze_row(Ip[i,:], unwrap_phase=unwrap_phase)
        phase[i] = p
        amp[i] = a
        chi2[i] = c
        
        # Initial quality control for chi^2 (noise), to make sure we're
        # not propagating pure noise
        if c > chi2_thresh:
            Ip[i,:] = 0.0 # So as to not contaminate other rows
            amp[i] = 0.0 # This is overly conservative. We probably 
                         # trust amplitudes more than phase. Note that
                         # this will be overwritten if linear_amp = True
        
    if linear_amp: # Replace the onion-peeled amplitude with the linear inversion
        amp_L1 = np.zeros(ny) # fringe amplitude at each row of L1 interferogram
        for i in range(ny):
            _, a, _ = analyze_row(I[i,:], unwrap_phase=unwrap_phase)
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
    # Each time scale propagated separately. Instantaneous and striation uncertainty propagate
    # in the normal way assuming rows are uncorrelated.
    # Since the zero wind phase is correlated strongly across rows, the uncertainty propagation
    # is unaffected by the inversion and is thus much more straighforward (just copy the result, basically).
    # This is true to a very good approximation since (1) the inversion is nearly linear and (2) the amplitude
    # uncertainty is much smaller than the phase uncertainty.
    
    if global_params['verbose']:
        print('%s:\t\t\t\t Beginning uncertainty propagation '% (timestamp()))
    ### Step 0: Characterize L1 and L2.1 interferograms with a single amp/phase per row
    ph_L1 = np.zeros(ny)
    A_L1 = np.zeros(ny)
    for i in range(ny):
        p,a,_ = analyze_row(I[i,:], unwrap_phase=unwrap_phase)
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
    cov_real_imag_L1_inst = np.zeros((ny,2,2))
    cov_real_imag_L1_striation = np.zeros((ny,2,2))
    for m in range(ny):
        # Jacobian of transformation from ampl/phase to real/imag.
        J = np.array([[np.cos(ph_L1[m]), -A_L1[m]*np.sin(ph_L1[m])],
                      [np.sin(ph_L1[m]),  A_L1[m]*np.cos(ph_L1[m])]])
        
        # Instantaneous
        cov_amp_phase = np.diag([I_amp_uncertainty[m], ph_err_inst[m]])**2 # assuming uncorrelated
        cov_amp_phase[np.isnan(cov_amp_phase)] = 0.0 # so the few nan *columns* don't corrupt the entire image
        cov_real_imag_L1_inst[m,:,:] = J.dot(cov_amp_phase).dot(J.T) # Error propagation
        
        # Striations
        cov_amp_phase = np.diag([I_amp_uncertainty[m], ph_err_striation[m]])**2 # assuming uncorrelated
        cov_amp_phase[np.isnan(cov_amp_phase)] = 0.0 # so the few nan *columns* don't corrupt the entire image
        cov_real_imag_L1_striation[m,:,:] = J.dot(cov_amp_phase).dot(J.T) # Error propagation
        

    ### Step 2: Propagate uncertainties through the path length inversion.
    # Note that instantaneous and striation errors are enhanced by the inversion, but
    # zero wind phase is not. This is because the latter is strongly correlated across rows.
    # The signal-dependent phase shift uncertainty is also strongly correlated across rows, and since it's mostly
    # important for the dimmest rows, the inversion does not add much. This has not been tested, but since the 
    # uncertainty is so uncertain anyway, it makes little difference.
    # Treat real and imaginary parts separately. Take special care for striation errors which have some 
    # correlation across rows.
    # Build covariance matrix of vector of real parts and of imag parts.
    
    # For instantaneous error, assume rows are uncorrelated
    cov_real_L1_inst = np.diag(cov_real_imag_L1_inst[:,0,0])
    cov_imag_L1_inst = np.diag(cov_real_imag_L1_inst[:,1,1])
    
    # For striations, rows are slightly correlated. Assume a triangular correlation as a function of row lag. The width
    # of this triangular function is an input parameter.
    corr_mat = np.diag(np.ones(ny), k=0)
    for n in range(1,Ncorr):
        corr_mat += 1.0 * (Ncorr-n)/Ncorr * (np.diag(np.ones(ny-n), k=n) + np.diag(np.ones(ny-n), k=-n))
    # Assume real and imaginary parts have the same correlation matrix (technically it's just the phase, but this simplification works).
    sig_real_L1 = np.diag(np.sqrt(cov_real_imag_L1_striation[:,0,0])) # sigma (not sigma^2) on the diagonal
    sig_imag_L1 = np.diag(np.sqrt(cov_real_imag_L1_striation[:,1,1]))
    cov_real_L1_striation = np.matmul(np.matmul(sig_real_L1, corr_mat), sig_real_L1) # covariance matrix 
    cov_imag_L1_striation = np.matmul(np.matmul(sig_imag_L1, corr_mat), sig_imag_L1) # covariance matrix 
    
    # Standard uncertainty propagation for matrix multiplication. Note zero_phase is assuming rows are correlated.
    Dinv = np.linalg.inv(D)
    cov_real_L2_inst = Dinv.dot(cov_real_L1_inst).dot(Dinv.T)
    cov_imag_L2_inst = Dinv.dot(cov_imag_L1_inst).dot(Dinv.T)    
    cov_real_L2_striation = Dinv.dot(cov_real_L1_striation).dot(Dinv.T)
    cov_imag_L2_striation = Dinv.dot(cov_imag_L1_striation).dot(Dinv.T)
    sigma_real_L2_inst = np.sqrt(np.diag(cov_real_L2_inst))
    sigma_imag_L2_inst = np.sqrt(np.diag(cov_imag_L2_inst))
    sigma_real_L2_striation = np.sqrt(np.diag(cov_real_L2_striation))
    sigma_imag_L2_striation = np.sqrt(np.diag(cov_imag_L2_striation))

    ### Step 3: Transform back to amp/phase #########
    # Each row will have a 2x2 covariance matrix describing the amplitude and phase
    cov_amp_phase_L2_inst = np.zeros((ny,2,2))
    cov_amp_phase_L2_striation = np.zeros((ny,2,2))
    for m in range(ny):
        # Jacobian of transformation from ampl/phase to real/imag
        J = np.array([[np.cos(ph_L2[m]), -A_L2[m]*np.sin(ph_L2[m])],
                      [np.sin(ph_L2[m]),  A_L2[m]*np.cos(ph_L2[m])]])
        # Jacobian of transformation from real/imag to ampl/phase
        Jinv = np.linalg.inv(J)
        cov_real_imag_inst       = np.diag([sigma_real_L2_inst[m],       sigma_imag_L2_inst[m]])**2 # assume uncorrelated
        cov_real_imag_striation  = np.diag([sigma_real_L2_striation[m],  sigma_imag_L2_striation[m]])**2 # assume uncorrelated
        cov_amp_phase_L2_inst[m,:,:]       = Jinv.dot(cov_real_imag_inst).dot(Jinv.T)
        cov_amp_phase_L2_striation[m,:,:]  = Jinv.dot(cov_real_imag_striation).dot(Jinv.T)
    # Extract amplitude and phase uncertainties
    amp_uncertainty = np.sqrt(cov_amp_phase_L2_inst[:,0,0]) # amp uncertainty only includes only shot noise as an appx. Overwritten below anyway.
    phase_uncertainty_inst       = np.sqrt(cov_amp_phase_L2_inst[:,1,1])
    phase_uncertainty_striation  = np.sqrt(cov_amp_phase_L2_striation[:,1,1])
    phase_uncertainty_zero_phase = ph_err_zero_phase # Note that this is just a copy of the input and bypasses the propagation code above.
    phase_uncertainty_signal     = ph_err_signal # Note that this is just a copy of the input and bypasses the propagation code above.
    
    if linear_amp: # Use linear propagation formula to overwrite
        DD = Dinv.dot(Dinv.T)
        amp_uncertainty = I_amp_uncertainty * np.sqrt(np.diag(DD))
        
    # Places where the values are nan should have nan uncertainty
    # (It was temporarily set to zero above for the purposes of 
    # propagating uncertainty, but should now be corrected)
    amp_uncertainty[np.isnan(amp)] = np.nan
    phase_uncertainty_inst[np.isnan(phase)] = np.nan
    phase_uncertainty_striation[np.isnan(phase)] = np.nan
    phase_uncertainty_zero_phase[np.isnan(phase)] = np.nan
    phase_uncertainty_signal[np.isnan(phase)] = np.nan
    
    if global_params['verbose']:
        print('%s:\t\t\t\t Returning from perform_inversion()'% (timestamp()))
        import sys
        sys.stdout.flush()
            
    return Ip, phase, amp, phase_uncertainty_inst, phase_uncertainty_striation, phase_uncertainty_zero_phase, phase_uncertainty_signal, amp_uncertainty, chi2





def fix_longitudes(lons, lon_target=None):
    '''
    Unwrap the list of longitudes to avoid 360-deg jumps. The list will
    be fixed so that it contains a value within 180 deg of lon_target and
    is otherwise continuous.
    
    INPUTS:
    
      *  lons       -- TYPE:array, UNITS:deg. An ordered list of longitudes to be unwrapped.
      *  lon_target -- TYPE:float, UNITS:deg. See above. (defaults to 0th element)
      
    OUTPUTS:
    
      *  lons_new   -- TYPE:array, UNITS:deg. An ordered list of longitudes with jumps removed.
      
    '''
    if lon_target is None:
        lon_target = lons[0]
    
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




def fix_longitudes_mat(lons, lon_target=None):
    '''
    Unwrap the matrix of longitudes to avoid 360-deg jumps along rows or columns. This is a
    matrix version of the function "fix_longitudes(...)" which is for arrays. The matrix will
    be fixed so that it contains a value within 180 deg of lon_target and
    is otherwise as continuous as possible.
    
    INPUTS:
    
      *  lons       -- TYPE:array(m,n), UNITS:deg. An ordered matrix of longitudes to be unwrapped.
      *  lon_target -- TYPE:float,      UNITS:deg. See above. (defaults to [0,0] element)
      
    OUTPUTS:
    
      *  lons_new   -- TYPE:array(m,n), UNITS:deg. An ordered matrix of longitudes with jumps removed.
    '''
    if lon_target is None:
        lon_target = lons[0,0]
    
    x = lons.copy()
    x[0,:] = fix_longitudes(x[0,:], lon_target)
    for i in range(np.shape(x)[1]):
        x[:,i] = fix_longitudes(x[:,i], x[0,i])
    return x
    

    

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





def interpolate_linear(x, y, x0, extrapolation='hold', prop_err = False, yerr = None, prop_nan=False):
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
      *  prop_nan      -- TYPE:bool, 
                                      * True: if one of the values is nan, return nan. 
                                      * False: if one of the values is nan, return 
                                               the other value. If both are nan, return nan.
                                              (Default)

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
        # if one value is nan, and prop_nan = False, just return the non-nan value
        if not prop_nan and np.isnan(y[j0]):
            if prop_err:
                return y[j1], yerr[j1]
            else:
                return y[j1]
        if not prop_nan and np.isnan(y[j1]):
            if prop_err:
                return y[j0], yerr[j0]
            else:
                return y[j0]
            
        # Otherwise, do it the normal way
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
      *  startstop      -- TYPE:bool, If True, use the correct start and stop variables. 
                                      If False, read the midpoint value and write it to
                                      both the "start" and "stop" variables. False is useful
                                      for debugging. (default True)
        
    OUTPUTS:
    
      *  L1_dict -- TYPE:dict. A dictionary containing information needed for
                               the level 2.1 processing. See documentation for 
                               level1_dict_to_level21_dict(...) for required keys.
                               

    '''
    
    if global_params['verbose']:
        print('%s:\t\t Opening %s' % (timestamp(), L1_fn.split('/')[-1]))
        import sys
        sys.stdout.flush()
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
    L1_dict['exp_time']                    = f['ICON_L0_MIGHTI_%s_Time_Integration' % sensor][0]*1e-3
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
    # Read attitude status register and save relevant bits
    att_reg = f['ICON_L1_MIGHTI_%s_SC_Attitude_Control_Register'%sensor][:][0]
    nb = 11
    attitude = np.zeros(nb,dtype=np.int8)
    s = '{0:b}'.format(att_reg).zfill(16)
    for b in range(nb):
        attitude[b] = int(s[-b-1])
    L1_dict['att_lvlh_normal'] = attitude[0]
    L1_dict['att_lvlh_reverse'] = attitude[1]
    L1_dict['att_limb_pointing'] = attitude[2]
    L1_dict['att_conjugate'] = attitude[6]
    L1_dict['att_zero_wind'] = attitude[10]
    # Read aperture position
    ap = f['ICON_L0_MIGHTI_%s_MT%s_Aperture1_Position' % (sensor, sensor)][0]
    assert ap in [0,2], "Aperture position not understood (%s but expected 0 or 2)"%ap
    if ap == 0:
        L1_dict['aperture'] = 'night'
    else: # ap == 2
        L1_dict['aperture'] = 'day'
    # New variables in v1.19: mag lat, mag lon, sza, slt.
    # For these new variables take middle of exposure
    L1_dict['mag_lat'] = f['ICON_L1_MIGHTI_%s_%s_Tangent_Magnetic_Latitude' % (sensor, emission_color.capitalize())][0,1,:]
    L1_dict['mag_lon'] = f['ICON_L1_MIGHTI_%s_%s_Tangent_Magnetic_Longitude' % (sensor, emission_color.capitalize())][0,1,:]
    L1_dict['slt'] = f['ICON_L1_MIGHTI_%s_%s_Tangent_Local_Solar_Time' % (sensor, emission_color.capitalize())][0,1,:]
    L1_dict['sza'] = f['ICON_L1_MIGHTI_%s_%s_Tangent_Solar_Zenith_Angle' % (sensor, emission_color.capitalize())][0,1,:]
    L1_dict['I_dc'] = f['ICON_L1_MIGHTI_%s_%s_Relative_Brightness' % (sensor, emission_color.capitalize())][0,:]
    
    # Orbit number (ensure backwards compatibility for running old L1 files without this variable)
    try:
        L1_dict['orbit_number'] = f.Orbit_Number
    except Exception as e:
        print('Orbit_Number attribute not found. Error = "%s"' % e)
        L1_dict['orbit_number'] = np.nan
        
    # Jitter (ensure backwards compatibility for running old L1 files without this variable)
    try:
        L1_dict['jitter'] = f['ICON_L1_MIGHTI_%s_SC_Pointing_Jitter' % (sensor)][0]
    except Exception as e:
        print('SC_Pointing_Jitter variable not found. Setting to 0.0. Error = "%s"' % e)
        L1_dict['jitter'] = 0.0    
        
    # Signal-dependent uncertainty (ensure backwards compatibility for running old L1 files without this variable)
    try:
        L1_dict['I_phase_signal_correction'] = f['ICON_L1_MIGHTI_%s_%s_VER_Correction_Phase' % (sensor, emission_color.capitalize())][0,:]
    except Exception as e:
#         print('VER_Correction_Phase variable not found. Setting to 0.0. Error = "%s"' % e)
        L1_dict['I_phase_signal_correction'] = np.zeros_like(L1_dict['I_phase_uncertainty'])
        
    # Signal-dependent uncertainty (ensure backwards compatibility for running old L1 files without this variable)
    try:
        L1_dict['I_phase_signal_uncertainty'] = f['ICON_L1_MIGHTI_%s_%s_VER_Correction_Phase_Uncertainties' % (sensor, emission_color.capitalize())][0,:]
    except Exception as e:
        print('VER_Correction_Phase_Uncertainties variable not found. Setting to 0.0. Error = "%s"' % e)
        L1_dict['I_phase_signal_uncertainty'] = np.zeros_like(L1_dict['I_phase_uncertainty'])
        
    # Patch for bug in L1 code regarding the uncertainty. Peg it at 40% for all of the v05 run
    L1_dict['I_phase_signal_uncertainty'] = 0.4 * abs(L1_dict['I_phase_signal_correction'])
    

    # Quality factors and flags
    L1_dict['quality'] =  f['ICON_L1_MIGHTI_%s_%s_Quality_Factor' % (sensor, emission_color.capitalize())][0]
    L1_dict['quality_flags'] = np.zeros((ny,5)) # Create 5 spaces (used to be 6 in v04)
    # Code that will let me know when Ken changes the format of the error flags    
    assert f['ICON_L1_MIGHTI_%s_Quality_Flag_Low_Signal_To_Noise_%s'%(sensor, emission_color.capitalize())].shape == (1,ny)
    assert f['ICON_L1_MIGHTI_%s_Quality_Flag_SAA'%sensor].shape == (1,)
    assert f['ICON_L1_MIGHTI_%s_Quality_Flag_Bad_Calibration'%sensor].shape == (1,)
    assert f['ICON_L1_MIGHTI_%s_%s_Quality_Factor' % (sensor, emission_color.capitalize())].shape == (1,ny)
    # Read quality flags into my array
    L1_dict['quality_flags'][:,0] = f['ICON_L1_MIGHTI_%s_Quality_Flag_Low_Signal_To_Noise_%s'%(sensor, emission_color.capitalize())][0,:]
    L1_dict['quality_flags'][:,1] = f['ICON_L1_MIGHTI_%s_Quality_Flag_SAA'%sensor][0]
    L1_dict['quality_flags'][:,2] = f['ICON_L1_MIGHTI_%s_Quality_Flag_Bad_Calibration'%sensor][0]
    # Not technically a L1 flag, but it fits here:
    L1_dict['quality_flags'][:,3] = bool(f[u'ICON_L0_MIGHTI_%s_Calibration_Lamp_1' % sensor][...].item()) | \
                                    bool(f[u'ICON_L0_MIGHTI_%s_Calibration_Lamp_2' % sensor][...].item())
    
    # Sun/moon flag (ensure backwards compatibility for running old L1 files without this variable)
    try:
        L1_dict['quality_flags'][:,4] = bool(f['ICON_L1_MIGHTI_%s_Quality_Flag_Sun_Moon_in_FoV'%sensor][0])
        assert f['ICON_L1_MIGHTI_%s_Quality_Flag_Sun_Moon_in_FoV'%sensor].shape == (1,), "Sun/Moon flag not the expected shape"
    except Exception as e:
        print('Quality_Flag_Sun_Moon_in_FoV variable not found. Setting to 0. Error = "%s"' % e)
        L1_dict['quality_flags'][:,4] = 0 

    
    # Dummy placeholder code for reading global attributes, if that matters
    nc_attrs = f.ncattrs()
    
    # Make compatible with both netCDF4-python v1.3.0 and v1.4.0:
    # Convert masked arrays to np.arrays, filling with nans.
    for v in L1_dict.keys():
        if isinstance(L1_dict[v], np.ma.masked_array):
            L1_dict[v] = L1_dict[v].filled(np.nan)
            
    if global_params['verbose']:
        print('%s:\t\t Closing %s' % (timestamp(), L1_fn.split('/')[-1]))
        import sys
        sys.stdout.flush()
    f.close()
    
    return L1_dict




def level1_dict_to_level15_dataset(L1_dict, emission_color, unwrap_phase=False):
    '''
    Starting with the level 1 data, remove S/C velocity and return 1D phase (vs altitude), with other info. This
    copies a lot from the first part of the Level 1 to Level 2.1 processing. No zero wind correction is made, and no
    smoothing is done, because the purpose of the L1.5 file is to determine what those corrections should be. This is run
    on a single file, and it is expected they will be combined by a separate function.
    
    INPUTS:
    
      *  L1_dict       -- TYPE:dict.  A dictionary containing variables needed for
                                      the level 2.1 processing. See level1_dict_to_level21_dict(...)
                                      
    OUTPUTS:
      *  L15_ds        -- TYPE:xarray.Dataset. Contains relevant variables for the Level 1.5 file.
    
    '''

    Ny, Nx = np.shape(L1_dict['I_phase'])
    sigma = global_params[emission_color]['sigma']

    ####  Load parameters from input dictionary
    Iraw = L1_dict['I_amp']*np.exp(1j*L1_dict['I_phase'])
    I_amp_uncertainty = L1_dict['I_amp_uncertainty']
    I_phase_uncertainty = L1_dict['I_phase_uncertainty']
    exp_time = L1_dict['exp_time']
    L1_fn = L1_dict['L1_fn']
    opd = L1_dict['optical_path_difference']
    sigma_opd = sigma * opd # Optical path difference, in units of wavelengths
    mighti_ecef_vectors = L1_dict['mighti_ecef_vectors']
    icon_velocity_vector = L1_dict['icon_velocity_vector']
    aperture = L1_dict['aperture']
    L1_quality_flags = L1_dict['quality_flags']
    L1_quality       = L1_dict['quality']
    I_dc = L1_dict['I_dc']
    sensor = L1_dict['sensor']

    # Load parameters which are averaged from start to stop of exposure.
    icon_alt = (L1_dict['icon_alt_start'] + L1_dict['icon_alt_stop'])/2
    icon_lat = (L1_dict['icon_lat_start'] + L1_dict['icon_lat_stop'])/2
    icon_lon = circular_mean(L1_dict['icon_lon_start'], L1_dict['icon_lon_stop'])
    tang_alt = (L1_dict['tang_alt_start'] + L1_dict['tang_alt_stop'])/2
    tang_lat = (L1_dict['tang_lat_start'] + L1_dict['tang_lat_stop'])/2
    tang_lon = circular_mean(L1_dict['tang_lon_start'], L1_dict['tang_lon_stop'])
    tmid     = L1_dict['time_start'] + (L1_dict['time_stop'] - L1_dict['time_start'])/2

    #### Remove spacecraft velocity
    icon_latlonalt = np.array([icon_lat, icon_lon, icon_alt])
    I = remove_satellite_velocity(Iraw, icon_latlonalt, icon_velocity_vector, mighti_ecef_vectors, sigma_opd)

    #### Determine az angle
    az = los_az_angle(icon_latlonalt, tang_lat, tang_lon, tang_alt)

    ### Extract 1D phase
    phL1 = np.zeros(Ny)
    ampL1 = np.zeros(Ny)
    chi2 = np.zeros(Ny)
    for i in range(Ny):
        p,a,c = analyze_row(I[i,:], unwrap_phase=unwrap_phase)
        phL1[i] = p
        ampL1[i] = a
        chi2[i] = c

    # phase to wind factor and mode
    f = phase_to_wind_factor(sigma_opd.mean())
    mode = 'Day'
    if exp_time > 45.:
        mode = 'Night'

        
    d = {
        'phase':(['time','row'], [phL1]),
        'Iamp': (['time','row'], [ampL1]),
        'Idc' : (['time','row'], [I_dc]),
        'alt' : (['time','row'], [tang_alt]),
        'lat' : (['time','row'], [tang_lat]),
        'lon' : (['time','row'], [tang_lon]),
        'az' : (['time','row'], [az]),
        'slt': (['time','row'], [L1_dict['slt']]),
        'chi2':(['time','row'], [chi2]),
        'qual':(['time','row'], [L1_quality]),
        'phase_to_wind_factor': f,
        'color':emission_color,
        'mode':(['time'], [mode]),
        'sensor':sensor,
        'cal_lamp':L1_quality_flags[0,3],
        'low_signal_corr':(['time','row'], [L1_dict['I_phase_signal_correction']]),
        }

    coords={
            'time':('time',[tmid]),
            'row' :('row',np.arange(Ny)),
        }
    
    
    ds = xr.Dataset(d, coords)
    
    # Adjust time base to match ICON convention (this doesn't affect anything in xarray but does affect the NC file)
    ds.time.encoding['units'] = 'milliseconds since 1970-01-01'
    
    return ds
    




def level21_quality(L1_quality_flags, L21_dict, L1_quality, top_layer_thresh=1.0, terminator_thresh = 0.0, 
                    chi2_thresh = np.inf, chi2_thresh_caution = np.inf, count_thresh = 0):
    '''
    Assess the quality of the L2.1 data product. This function generates quality flags and overall 
    quality factor (0-1) for wind and VER.
    
    INPUTS:
    
      * L1_quality_flags  -- TYPE:array(ny,ne). The quality flags from the L1 file, but binned in altitude to match the L2.1 data.
                                                This should have the same format as described in the documentation for
                                                level1_dict_to_level21_dict(...).
      * L21_dict          -- TYPE:dict.         A dictionary containing output variables of the level 2.1 processing.
                                                See documentation for level1_dict_to_level21_dict(...) for required keys.
      * L1_quality        -- TYPE:array(ny).    The quality factosr from the L1 processing. The L2.1 wind quality flag should
                                                be <= this.
                                                
    OPTIONAL INPUTS:
    
      * top_layer_thresh.   -- TYPE:float.      Fraction of airglow above the top observed altitude. Consider the total column 
                                                brightness (i.e, the integral of the VER profile). When a large fraction
                                                comes from above the top observed altitude, the quality flag is raised. (default 1.0)
      * terminator_thresh   -- TYPE:float, UNITS:km. Consider two points along the line of sight, both X km from the tangent 
                                                point (where X is this parameter). One is nearer to the spacecraft than the 
                                                tangent point, and one is farther. If these two points are on opposite sides of the 
                                                terminator, raise a quality flag. Note that this quality flag will also be raised if
                                                any observations at higher tangent altitudes for the same observation are flagged, because
                                                the inversion mixes information from those observations. (default 0.0)
      *  chi2_thresh         -- TYPE:float, UNITS:rad^2.  If the mean-square residual of phase within a row is larger than this,
                                                          the sample is discarded, and isn't used in the inversion. (default np.inf)
      *  chi2_thresh_caution -- TYPE:float, UNITS:rad^2.  If the mean-square residual of phase within a row is larger than this,
                                                          the sample should be treated with caution but is used in the inversion. This
                                                          affects the quality control only. (default np.inf)
      *  count_thresh        -- TYPE:int,                 A profile needs at least this many good/caution samples to pass
                                   
    OUTPUTS:
    
      * wind_quality      --TYPE:array(ny).     A number from 0 (Bad) to 1 (Good) for each altitude, quantifying 
                                                the quality of the wind data.
      * ver_quality       --TYPE:array(ny).     A number from 0 (Bad) to 1 (Good) for each altitude, quantifying
                                                the quality of the volume emission rate data.
      * quality_flags     --TYPE:array(ny,ne).  For each altitude, multiple flags exist, each of which is either 0 (False) or 1 (True). 
                                                Some flags are propaged directly from the L1 file. See documentation for 
                                                level1_dict_to_level21_dict(...) for flag definitions.
    '''
    
    if global_params['verbose']:
        print('%s:\t\t\t\t Entered level21_quality(...) with inputs:'% (timestamp()))
        print('L1_quality_flags=\n%s\n'% repr(L1_quality_flags))
        print('L21_dict = (see below)\n')
        print('L1_quality=\n%s\n'% repr(L1_quality))
        print('top_layer_thresh=%s\n'% repr(top_layer_thresh))
        print('terminator_thresh=%s\n'% repr(terminator_thresh))
        print('chi2_thresh=%s\n'% repr(chi2_thresh))
        print('chi2_thresh_caution=%s\n'% repr(chi2_thresh_caution))
        for var in ['los_wind', 'fringe_amplitude', 'H', 'alt', 'mighti_ecef_vectors', 'lat', 'lon', 'alt', 'time', 'chi2', 'jitter', 'ver']:
            print('L21_dict[\'%s\']=\n%s\n'% (var, repr(L21_dict[var])))
        import sys
        sys.stdout.flush()
        
    # Note that the "thermal drift correction is uncertain" flag is calculated after a full day of data is run, in
    # level1_to_level21_without_info_file(...). This is necessary so gaps can be detected (i.e., it can't be calculated 
    # for an individual profile.)
    
    ny = len(L21_dict['los_wind'])
    
    ############################# Create error flags (including those copied from L1) #########################
    quality_flags = np.zeros((ny, 12))
    quality_flags[:,:5] = L1_quality_flags[:,:] # copied from L1
    
    #### Significant airglow above top altitude
    if global_params['verbose']:
        print('%s:\t\t\t\t Checking airglow above top altitude'% (timestamp()))
        import sys
        sys.stdout.flush()
    a = L21_dict['fringe_amplitude'] # VER profile
    H = L21_dict['H'] #  scale height above the top altitude that was assumed in the inversion
    dz = np.diff(L21_dict['alt'])
    f = a[-1]*H / (a[-1]*H + np.sum(a[:-1]*dz)) # fraction of column brightness above top altitude
    if f > top_layer_thresh:
        quality_flags[:,7] = 1
    
    #### Line of sight crosses the terminator
    if global_params['verbose']:
        print('%s:\t\t\t\t Checking LoS crossing terminator'% (timestamp()))
        import sys
        sys.stdout.flush()
        
    ## BJH commented 2022 Mar 16. This old code explicitly computes the LoS and sees if it crosses the terminator.
    ## It may be what's causing the hangs due to the call to ephem.
#     RE = 6371. # km, earth radius
#     h_sc = 80. # km, EUV screening height
#     for i in range(ny):
#         look_xyz = L21_dict['mighti_ecef_vectors'][i,:]
#         tang_xyz = ICON.wgs84_to_ecef([ L21_dict['lat'][i], L21_dict['lon'][i], L21_dict['alt'][i] ])
#         xyz0 = tang_xyz - terminator_thresh*look_xyz # ECEF coords test point 0
#         xyz1 = tang_xyz + terminator_thresh*look_xyz # ECEF coords test point 1
#         sza0 = ICON.get_solar_zenith_angle(L21_dict['time'], *ICON.ecef_to_wgs84(xyz0)) # solar zenith angle 0
#         sza1 = ICON.get_solar_zenith_angle(L21_dict['time'], *ICON.ecef_to_wgs84(xyz1)) # solar zenith angle 1
#         sza_term = 180./np.pi * (np.pi - np.arcsin((RE + h_sc)/(RE + L21_dict['alt'][i]))) # sza of termiantor
#         if np.sign(sza0-sza_term) != np.sign(sza1-sza_term): # then the points straddle the terminator
#             quality_flags[i,8] = 1
    ## New code 2022 Mar 16, to just use SZA
    sza_bad = abs(L21_dict['sza'] - 98.) < global_params[L21_dict['emission_color']]['sza_term_thresh'] # This is a messy handling of global params
    quality_flags[sza_bad,8] = 1.
    
    # If any points above this point are flagged, flag this point too.
    if global_params['verbose']:
        print('%s:\t\t\t\t Masking points above the LoS-terminator-flagged point'% (timestamp()))
        import sys
        sys.stdout.flush()
    if any(quality_flags[:,8]):
        i = np.where(quality_flags[:,8] > 0)[0].max()
        quality_flags[:i,8] = 1
        
    #### Low SNR after inversion (this captures high chi^2 and VER < 0)
    if global_params['verbose']:
        print('%s:\t\t\t\t Checking for low SNR after inversion'% (timestamp()))
        import sys
        sys.stdout.flush()
    # Use two thresholds of chi^2.
    # 2022 May 31 BJH: Changed the ver<0 flag to 6 instead of 11. (If VER<0 then there is no way the winds or VER are actually good).
    quality_flags[:,6]  = (L21_dict['chi2'] > chi2_thresh) | (L21_dict['ver'] < 0)
    quality_flags[:,11] = (L21_dict['chi2'] > chi2_thresh_caution) & (L21_dict['chi2'] <= chi2_thresh)
    
    
    #### S/C pointing is not stable -- Jitter is too large
    if global_params['verbose']:
        print('%s:\t\t\t\t Checking for jitter'% (timestamp()))
        import sys
        sys.stdout.flush()
    if L21_dict['jitter'] > global_params[L21_dict['emission_color']]['jitter_thresh']:
        quality_flags[:,10] = 1
        
#     #### Bottom row of green is untrustworthy. Set this using the L1 quality flag
#     if global_params['verbose']:
#         print('%s:\t\t\t\t Masking bottom two rows of green'% (timestamp()))
#         import sys
#         sys.stdout.flush()
#     # TODO: remove this once we figure it out.
#     if L21_dict['emission_color'] == 'green':
#         quality_flags[:2,2] = 1 # "Bad calibration"
    # 2021 Feb 8: Keeping two bottom rows for v05 tests. Let's see how the v05 run turns out.
    
    ########################### Calculate overall quality (floor set by L1 quality) ############################
    if global_params['verbose']:
        print('%s:\t\t\t\t Beginning overall quality calculation'% (timestamp()))
        import sys
        sys.stdout.flush()
    # There are a lot of opinions built in to this section of the code.
    wind_quality = np.ones(ny)
    ver_quality = np.ones(ny)
    for i in range(ny):
        # For each possible contributing factor, collect a set of ratings or "dings" against this point.
        # The final quality factor should be the minimum of all the ratings. 
        wind_ratings = [1.0] # start with 1.0 -- if there are no dings, the quality factor is 1.0
        ver_ratings = [1.0]
#         if quality_flags[i,1]: # SAA  # Commented 2020-10-15 BJH: chi^2 filtering is replacing this.
#             wind_ratings.append(0.5)
#             ver_ratings.append(0.5)
        if quality_flags[i,3]: # Cal lamps on (TODO: Remove this once we trust it better)
            wind_ratings.append(0.5)
            ver_ratings.append(0.5)
        if quality_flags[i,6]: # SNR too low
            wind_ratings.append(0.0) # phase is definitely bad (?)
            ver_ratings.append(0.0) # in practice, it seems the VERs are bad here too, especially in the SAA
        if quality_flags[i,11]: # SNR possibly too low ("Caution")
            wind_ratings.append(0.5) # phase is maybe bad
            ver_ratings.append(0.5) # VER is probably fine, but still questionable, especially in SAA
        if quality_flags[i,7]: # airglow above 300 km
            wind_ratings.append(0.5)
            ver_ratings.append(0.5)
        if quality_flags[i,8]: # proximity to terminator
            wind_ratings.append(0.5) 
            ver_ratings.append(0.5)
        if quality_flags[i,10]: # S/C pointing is not stable:
            wind_ratings.append(0.5) # phase is not trustworthy because S/C velocity subtraction is uncertain
            ver_ratings.append(1.0) # but there's no reason to distrust VER
#         if quality_flags[i,2]: # Bad calibration, including the problem at bottom two rows of green
#             wind_ratings.append(0.0) # don't use wind but VER is fine. 2021 Feb 9: Removed this for the v05 run. Trust the L1 code.
        if np.isnan(L21_dict['ver'][i]): # If VER is nan, the quality should be set to 0 
                                         # (I'm not sure why this is happening but suspect it's from L1)
            ver_ratings.append(0.0)        
        if np.isnan(L21_dict['los_wind'][i]): # If wind is nan, the quality should be set to 0 
                                         # (I'm not sure why this is happening but suspect it's from L1)
            wind_ratings.append(0.0)
        # Lastly, append the L1 quality factor because that should be the maximum allowed rating
        wind_ratings.append(L1_quality[i])
        ver_ratings.append(L1_quality[i]) # Should I do this?
        # Compile final quality factor
        wind_quality[i] = min(wind_ratings)
        ver_quality[i]  = min(ver_ratings)
        
    # Final check to make sure there are enough continguous samples for the inversion to be trustworthy
#     count = (np.array(wind_quality) > 0).sum() # This computes total # of samples
    countvec = np.zeros(len(wind_quality))
    countvec[0] = wind_quality[0] > 0
    for i in range(1,len(wind_quality)):
        if wind_quality[i] > 0:
            countvec[i] = countvec[i-1] + 1
    count = countvec.max()
    if count < count_thresh:
        wind_quality[:] = 0
        ver_quality[:] = 0
        quality_flags[:,5] = 1
        
    if global_params['verbose']:
        print('%s:\t\t\t\t Ending overall quality calculation'% (timestamp()))
        import sys
        sys.stdout.flush()
    
    return wind_quality, ver_quality, quality_flags


    


def level1_dict_to_level21_dict(L1_dict, linear_amp = True, sigma = None, top_layer = None, H = None,
                                integration_order = None, account_for_local_projection = None, 
                                bin_size = None,
                                top_layer_thresh = None, terminator_thresh = None, zero_wind_ref = None,
                                corr_notch_drift = None, chi2_thresh = None, chi2_thresh_caution = None,
                                count_thresh = None,
                                f107 = None, f107a = None, f107p = None, apmsis = None, Tmult=1.,
                                unwrap_phase = False, zero_wind_phase = None, 
                                err_striation = None, err_zero_phase = None, Ncorr = 7):
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
                                                                      Uncertainty (random/statistical) in the mean phase of each row of interferogram
                                      * I_phase_signal_uncertainty -- TYPE:array(ny),    UNITS:rad. 
                                                                      Another component of phase uncertainty, coming from the signal-dependent phase 
                                      * I_phase_signal_correction  -- TYPE:array(ny),    UNITS:rad. 
                                                                      The actual correction used for the signal-dependent phase shift
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
                                      * acknowledgement            -- TYPE:str.
                                                                      The value in the global attribute "Acknowledgement" which will be propagated
                                      * att_lvlh_normal            -- TYPE:int. 0 or 1
                                                                      Attitude register bit 0: LVLH Normal
                                      * att_lvlh_reverse           -- TYPE:int. 0 or 1
                                                                      Attitude register bit 1: LVLH Reverse
                                      * att_limb_pointing          -- TYPE:int. 0 or 1
                                                                      Attitude register bit 2: Earth Limb Pointing
                                      * att_conjugate              -- TYPE:int. 0 or 1
                                                                      Attitude register bit 6: Conjugate Maneuver
                                      * att_zero_wind              -- TYPE:int. 0 or 1
                                                                      Attitude register bit 10: Zero Wind Maneuver
                                      * aperture                   -- TYPE:str. 'day' or 'night'. Day has a smaller aperture than night.
                                      * mag_lat                    -- TYPE:array(ny)     UNITS:deg.
                                                                      Magnetic quasi-dipole latitude
                                      * mag_lon                    -- TYPE:array(ny)     UNITS:deg.
                                                                      Magnetic quasi-dipole longitude
                                      * sza                        -- TYPE:array(ny)     UNITS:deg.
                                                                      Solar zenith angle
                                      * slt                        -- TYPE:array(ny)     UNITS:hour.
                                                                      Solar local time
                                      * orbit_number               -- TYPE:int.
                                                                      ICON orbit number
                                      * jitter                     -- TYPE:float.        UNITS:deg
                                                                      How much the S/C pointing is changing.
                                      * I_dc                       -- TYPE:array(ny)     UNITS:R.
                                                                      Brightness observed at each row. Like I_amp, but generated from the 
                                                                      DC value of the measured interferogram, not the fringe amplitude.
                                      * quality                    -- TYPE:array(ny).
                                                                      The overall quality, between 0 (bad) and 1 (good) for each altitude
                                      * quality_flags              -- TYPE:array(ny,ne).
                                                                      A register of flags (0=False, 1=True) that indicate various
                                                                      causes for concern in the L1 file. These flags are taken directly
                                                                      from the L1 file. The definitions for each position are:
                                                                       * 0: High if SNR is low enough to cause possible systematic errors
                                                                       * 1: Proximity to South Atlantic Anomaly
                                                                       * 2: Bad calibration detected
                                                                       * 3: Calibration lamps are on
                                                                       * 4: Solar/lunar contamination
                                                                                  
    OPTIONAL INPUTS - If None, defaults from MIGHTI_L2.global_params will be used 
    
      *  linear_amp          -- TYPE:bool.  If True, a linear inversion is used on fringe amplitude. If False, the 
                                            fringe amplitude is estimated during the onion-peeling. These are the
                                            same in the absence of noise, but different with noise, especially
                                            for computing uncertainties. (default True)
      *  sigma               -- TYPE:float, UNITS:m^-1. The wavenumber of the emission (1/wavelength)
      *  top_layer           -- TYPE:str,   'thin': assume VER goes to zero above top layer
                                            'exp':  assume VER falls off exponentially in altitude
      *  H                   -- TYPE:float, UNITS:km. The VER scale height to use if top_layer='exp'
      *  integration_order   -- TYPE:int,   0: Use Riemann-sum rule for discretizing line-of-sight integral
                                            1: Use trapezoidal rule for discretizing line-of-sight integral
      *  account_for_local_projection -- TYPE:bool. If False, a simple inversion is used.
                                            If True, the inversion accounts for the fact that the ray is not 
                                            perfectly tangent to each shell at each point along the ray.
      *  bin_size            -- TYPE:int,   The number of rows of the interferogram to bin together to 
                                            improve statistics at the cost of altitude resolution.
      *  top_layer_thresh.   -- TYPE:float. Fraction of airglow above the top observed altitude. Consider the total column 
                                            brightness (i.e, the integral of the VER profile). When a large fraction
                                            comes from above the top observed altitude, the quality flag is raised.
      *  terminator_thresh   -- TYPE:float, UNITS:km. Consider two points along the line of sight, both X km from the tangent 
                                            point (where X is this parameter). One is nearer to the spacecraft than the 
                                            tangent point, and one is farther. If these two points are on opposite sides of the 
                                            terminator, raise a quality flag. Note that this quality flag will also be raised if
                                            any observations at lower zenith angles for the same observation are flagged, because
                                            the inversion mixes information from those observations.
      *  zero_wind_ref       DEPRECATED
      *  corr_notch_drift    DEPRECATED   
      *  zero_wind_phase     -- TYPE:array(ny), UNITS:rad. The zero wind phase to remove from the L1 data before processing. In operational use,
                                                           this will be informed by the Zero-Phase file. If None, then zero will be used.
      *  err_striation       -- TYPE:array(ny), UNITS:rad. An estimate of the error introduced by daily calibrations, which are correlated for an 
                                                           entire 24-hour period (00:00 - 23:59 UT).
      *  err_zero_phase      -- TYPE:array(ny), UNITS:rad. An estimate of the error introduced by the zero-wind phase correction, which is correlated
                                                           over long time periods (several weeks or more).
      *  chi2_thresh         -- TYPE:float, UNITS:rad^2.  If the mean-square residual of phase within a row is larger than this,
                                                          the sample is discarded, and isn't used in the inversion.
      *  chi2_thresh_caution -- TYPE:float, UNITS:rad^2.  If the mean-square residual of phase within a row is larger than this,
                                                          the sample should be treated with caution but is used in the inversion. This
                                                          affects the quality control only.
      *  count_thresh        -- TYPE:int,                 A profile needs at least this many good/caution samples to pass                        
      *  f107                -- TYPE:float, UNITS:sfu. F10.7 value for the date of interest (None is use values from pyglow)
      *  f107a               -- TYPE:float, UNITS:sfu. F10.7a value for the date of interest (None is use values from pyglow)
      *  f107p               -- TYPE:float, UNITS:sfu. F10.7p value for the date of interest (None is use values from pyglow)
      *  apmsis              -- TYPE:array(7), Ap vector that MSIS expects (None is use values from pyglow)      
      *  Tmult               -- TYPE:float, Artificial multiplier to apply to MSIS, for sensitivity studies (default 1)
      *  unwrap_phase        -- TYPE:bool.  If True, unwrap the onion-peeled phase to remove 2pi discontinuities (default False)
      *  Ncorr               -- TYPE:int,   UNITS:pixels. The row lag at which the residual striation correlation drops to 0. 
                                          (default: 7, which was determined by Monte Carlo simulations for the un-binned case.)
                                  
    OUTPUTS:
    
      *  L21_dict            -- TYPE:dict. A dictionary containing output variables of the Level 2.1 processing:

                    * los_wind                  -- TYPE:array(ny),   UNITS:m/s.   Line-of-sight wind profile (+ towards MIGHTI)
                    * los_wind_error_inst       -- TYPE:array(ny),   UNITS:m/s.   "Instantaneous" uncertainty of los_wind
                    * los_wind_error_striation  -- TYPE:array(ny),   UNITS:m/s.   "Striation" uncertainty of los_wind
                    * los_wind_error_zero_phase -- TYPE:array(ny),   UNITS:m/s.   Accuracy of los_wind (determined by zero wind phase estimate)
                    * los_wind_error_signal     -- TYPE:array(ny),   UNITS:m/s.   Uncertainty due to signal-dependent phase shift effect
                    * los_wind_error            -- TYPE:array(ny),   UNITS:m/s.   The RSS of inst and striation errors. (backwards compatibility)
                    * lat                       -- TYPE:array(ny),   UNITS:deg.   Latitude of each point in profile
                    * lon                       -- TYPE:array(ny),   UNITS:deg.   Longitude of each point in profile
                    * alt                       -- TYPE:array(ny),   UNITS:alt.   Altitude of each point in profile
                    * time_start                -- TYPE:datetime (timezone naive) Time at start of exposure in UTC
                    * time_stop                 -- TYPE:datetime (timezone naive) Time at end of exposure in UTC
                    * time                      -- TYPE:datetime (timezone naive) Time at midpoint of exposure in UTC
                    * exp_time                  -- TYPE:float,       UNITS:s.     Exposure time
                    * az                        -- TYPE:array(ny),   UNITS:deg.   The azimuth angle of the line of sight
                                                                                  at the tangent point (deg East of North)
                    * emission_color            -- TYPE:str.                      'red' or 'green'
                    * sensor                    -- TYPE:str.                      'A' or 'B'
                    * icon_alt                  -- TYPE:float,       UNITS:km.    Spacecraft altitude
                    * icon_lat                  -- TYPE:float,       UNITS:deg.   Spacecraft latitude
                    * icon_lon                  -- TYPE:float,       UNITS:deg.   Spacecraft longitude [0,360]
                    * ver                       -- TYPE:array(ny),   UNITS:ph/cm^3/s. Volume emission rate.
                    * ver_error                 -- TYPE:array(ny),   UNITS:ph/cm^3/s. Statistical uncertainty in ver.
                    * fringe_amplitude          -- TYPE:array(ny),   UNITS:arb.   The fringe contrast, a proxy for volume emission rate
                    * fringe_amplitude_error    -- TYPE:array(ny),   UNITS:arb.   Uncertainty in fringe_amplitude (1-sigma)
                    * mighti_ecef_vectors       -- TYPE:array(ny,3).              ECEF unit vector for each line of sight
                    * icon_velocity_ecef_vector -- TYPE:array(3).    UNITS:m/s.   ECEF vector of spacecraft velocity
                    * file_creation_time        -- TYPE:datetime (timezone naive) Time this processing was run in UTC
                    * source_files              -- TYPE:list of str.              All science files that went into creating this file
                    * bin_size                  -- TYPE:int.                      Bin size used in the processing
                    * top_layer                 -- TYPE:str.                      How the top layer was handled: 'thin' or 'exp'
                    * zero_wind_ref             -- DEPRECATED
                    * corr_notch_drift          -- DEPRECATED                  
                    * H                         -- TYPE:float.       UNITS:km.    The VER scale height used if top_layer='exp' 
                    * integration_order         -- TYPE:int.                      Order of integration used in inversion: 0 or 1
                    * I                         -- TYPE:array(ny,nx) UNITS:arb.   The complex-valued, onion-peeled interferogram
                    * chi2                      -- TYPE:array(ny)    UNITS:rad^2. The normalized chi^2 from the phase
                                                                                  extraction from each row of the interferogram.
                    * wind_quality              -- TYPE:array(ny)    UNITS:arb.   A number from 0 (Bad) to 1 (Good) for each altitude.
                    * ver_quality               -- TYPE:array(ny)    UNITS:arb.   A number from 0 (Bad) to 1 (Good) for each altitude.
                    * quality_flags             -- TYPE:array(ny,ne).             For each altitude, multiple flags exist, 
                                                                                  each of which is either 0 (False) or 1 (True). 
                                                                                  Some flags are propaged directly from the L1 file.
                                                                                  The definition for each bit is:
                                                                                   * 0 : (From L1) SNR too low to reliably perform L1 processing
                                                                                   * 1 : (From L1) Proximity to South Atlantic Anomaly
                                                                                   * 2 : (From L1) Bad calibration 
                                                                                   * 3 : (From L1) Calibration lamps are on
                                                                                   * 4 : (From L1) Solar/lunar contamination
                                                                                   * 5 : Not enough valid points in profile
                                                                                   * 6 : SNR too low after inversion
                                                                                   * 7 : Significant airglow above 300 km
                                                                                   * 8 : Line of sight crosses the terminator
                                                                                   * 9 : Thermal drift correction is uncertain
                                                                                   * 10: S/C pointing is not stable
                                                                                   * 11: SNR possibly too low
                    * acknowledgement           -- TYPE:str.                      A copy of the Acknowledgement attribute in the L1 file.
                    * att_lvlh_normal           -- TYPE:int. 0 or 1               Attitude register bit 0: LVLH Normal
                    * att_lvlh_reverse          -- TYPE:int. 0 or 1               Attitude register bit 1: LVLH Reverse
                    * att_limb_pointing         -- TYPE:int. 0 or 1               Attitude register bit 2: Earth Limb Pointing
                    * att_conjugate             -- TYPE:int. 0 or 1               Attitude register bit 6: Conjugate Maneuver
                    * att_zero_wind             -- TYPE:int. 0 or 1               Attitude register bit 10: Zero Wind Maneuver
                    * mag_lat                   -- TYPE:array(ny),    UNITS:deg.  Magnetic quasi-dipole latitude
                    * mag_lon                   -- TYPE:array(ny),    UNITS:deg.  Magnetic quasi-dipole longitude
                    * sza                       -- TYPE:array(ny),    UNITS:deg.  Solar zenith angle
                    * slt                       -- TYPE:array(ny),    UNITS:hour. Solar local time
                    * ver_dc                    -- TYPE:array(ny),    UNITS:ph/cm^3/s. Like ver, but generated from the DC value of the
                                                                                  interferogram, not the fringe amplitude.
                    * orbit_node                -- TYPE:int. 0 or 1               0 = Latitude of ICON is increasing
                                                                                  1 = Latitude of ICON is decreasing
                    * orbit_number              -- TYPE:int.                      ICON orbit number
                    * jitter                    -- TYPE:float.        UNITS:deg.  How much the S/C pointing is changing. This is the stddev 
                                                                                  of the residual of a linear fit to the Euler angles in a 
                                                                                  60-sec window.
                    * zero_wind_phase           -- TYPE:array(ny),   UNITS:rad.   The zero wind phase used to process this exposure.   
                    * low_signal_corr           -- TYPE:array(ny),   UNITS:m/s.   The correction used for the low signal effect (not fed through inversion)
    '''
    #### Parse input parameters and load defaults
    if global_params['verbose']:
        print('%s:\t\t\t Parsing input params'% (timestamp()))
        import sys
        sys.stdout.flush()
    emission_color = L1_dict['emission_color']
    sensor = L1_dict['sensor']
    params = global_params[emission_color]
    if sigma is None:
        sigma = params['sigma']
    if top_layer is None:
        top_layer = params['top_layer']
    if H is None:
        H = params['H']
    if integration_order is None:
        integration_order = params['integration_order']
    if account_for_local_projection is None:
        account_for_local_projection = params['account_for_local_projection']
    if bin_size is None:
        bin_size = params['bin_size']
    bin_size = int(bin_size)
    if top_layer_thresh is None:
        top_layer_thresh = params['top_layer_thresh']
    if terminator_thresh is None:
        terminator_thresh = params['terminator_thresh']
    if chi2_thresh is None:
        chi2_thresh = params['chi2_thresh']
    if chi2_thresh_caution is None:
        chi2_thresh_caution = params['chi2_thresh_caution']
    if count_thresh is None:
        count_thresh = params['count_thresh_night']
        if L1_dict['exp_time'] < 45.: # Day mode
            count_thresh = params['count_thresh_day']
    if zero_wind_phase is None:
        zero_wind_phase = np.zeros_like(L1_dict['tang_alt_start'])
    if err_striation is None:
        err_striation = np.zeros_like(L1_dict['tang_alt_start'])
    if err_zero_phase is None:
        err_zero_phase = np.zeros_like(L1_dict['tang_alt_start'])
        
    # v5.00+: check for deprecated inputs related to zero wind and notches
    if zero_wind_ref is not None:
        print('WARNING: The "zero_wind_ref" input is deprecated and is being ignored.')
    if corr_notch_drift is not None:
        print('WARNING: The "corr_notch_drift" input is deprecated and is being ignored.')
        
            
    ####  Load parameters from input dictionary
    if global_params['verbose']:
        print('%s:\t\t\t Loading input dict'% (timestamp()))
        import sys
        sys.stdout.flush()
    Iraw = L1_dict['I_amp']*np.exp(1j*L1_dict['I_phase'])
    I_amp_uncertainty = L1_dict['I_amp_uncertainty']
    err_inst = L1_dict['I_phase_uncertainty']
    err_signal = L1_dict['I_phase_signal_uncertainty']
    source_files = L1_dict['source_files']
    exp_time = L1_dict['exp_time']
    L1_fn = L1_dict['L1_fn']
    opd = L1_dict['optical_path_difference']
    sigma_opd = sigma * opd # Optical path difference, in units of wavelengths
    mighti_ecef_vectors = L1_dict['mighti_ecef_vectors']
    icon_velocity_vector = L1_dict['icon_velocity_vector']
    aperture = L1_dict['aperture']
    L1_quality_flags = L1_dict['quality_flags']
    L1_quality       = L1_dict['quality']
    I_dc = L1_dict['I_dc']
    mag_lat = L1_dict['mag_lat']
    mag_lon = L1_dict['mag_lon']
    sza     = L1_dict['sza']
    slt     = L1_dict['slt']
    low_signal_phase_corr = L1_dict['I_phase_signal_correction']
    # Load parameters which are averaged from start to stop of exposure.
    icon_alt = (L1_dict['icon_alt_start'] + L1_dict['icon_alt_stop'])/2
    icon_lat = (L1_dict['icon_lat_start'] + L1_dict['icon_lat_stop'])/2
    icon_lon = circular_mean(L1_dict['icon_lon_start'], L1_dict['icon_lon_stop'])
    tang_alt = (L1_dict['tang_alt_start'] + L1_dict['tang_alt_stop'])/2
    tang_lat = (L1_dict['tang_lat_start'] + L1_dict['tang_lat_stop'])/2
    tang_lon = circular_mean(L1_dict['tang_lon_start'], L1_dict['tang_lon_stop'])
    tmid     = L1_dict['time_start'] + (L1_dict['time_stop'] - L1_dict['time_start'])/2
        
    # Zero wind phase adjustment
    if global_params['verbose']:
        print('%s:\t\t\t Adjusting zero wind (%s)'% (timestamp(), zero_wind_ref))
        import sys
        sys.stdout.flush()
    L1_dict['quality_flags'][np.isnan(zero_wind_phase), 2] = 1 # mark nans as "Bad Calibration" in L1
    ny0, nx0 = Iraw.shape
    for j in range(nx0):
        Iraw[:,j] = Iraw[:,j]*np.exp(-1j*zero_wind_phase)
    
    
    #### Remove Satellite Velocity
    if global_params['verbose']:
        print('%s:\t\t\t Removing S/C velocity '% (timestamp()))
        import sys
        sys.stdout.flush()
    icon_latlonalt = np.array([icon_lat, icon_lon, icon_alt])
    I = remove_satellite_velocity(Iraw, icon_latlonalt, icon_velocity_vector, mighti_ecef_vectors, sigma_opd)
                         
    #### Bin data: average nearby rows together
    if global_params['verbose']:
        print('%s:\t\t\t Binning inputs '% (timestamp()))
        import sys
        sys.stdout.flush()
    I        = bin_image(bin_size, I)
    I_dc     = bin_array(bin_size, I_dc)
    tang_lat = bin_array(bin_size, tang_lat)
    tang_lon = bin_array(bin_size, tang_lon, lon=True)
    tang_alt_lo = bin_array(bin_size, tang_alt, method='min')
    tang_alt_hi = bin_array(bin_size, tang_alt, method='max')
    tang_alt    = bin_array(bin_size, tang_alt)
    ny, nx = np.shape(I)
    mighti_ecef_vectors_center = mighti_ecef_vectors[:,nx/2,:] # For reporting in output file, determine ecef vector at center of row
    mighti_ecef_vectors_center = bin_image(bin_size, mighti_ecef_vectors_center) # bin each component separately
    L1_quality_flags = bin_image(bin_size, L1_quality_flags, method='max') # bin each flag separately, taking *max* over the bin
    L1_quality       = bin_array(bin_size, L1_quality,       method='min') # bin quality factor, taking *min* over the bin
    I_amp_uncertainty   = bin_uncertainty(bin_size, I_amp_uncertainty)
    err_inst       = bin_uncertainty(bin_size, err_inst)
    err_striation  = bin_array(bin_size, err_striation)  # Note error is not reduced by averaging (this is not exactly correct but is conservative).
    err_zero_phase = bin_array(bin_size, err_zero_phase) # Note error is not reduced by averaging
    err_signal     = bin_array(bin_size, err_signal)     # Note error is not reduced by averaging
    low_signal_phase_corr = bin_array(bin_size, low_signal_phase_corr)
    mag_lat = bin_array(bin_size, mag_lat)
    mag_lon = bin_array(bin_size, mag_lon, lon=True)
    sza     = bin_array(bin_size, sza)
    slt     = bin_array(bin_size, slt*15., lon=True)/15. # Avoid discontinuity problem by multiplying by 15 and treating like lon
    Ncorr   = int(np.ceil(1.0*Ncorr/bin_size)) # Correlation length drops (ceil is used for conservatism because an int is needed)
    
    
    #### Determine geographical locations of inverted wind
    lat, lon, alt = attribute_measurement_location(tang_lat, tang_lon, tang_alt,
                                                   integration_order=integration_order)
    
    #### Onion-peel interferogram
    if global_params['verbose']:
        print('%s:\t\t\t Entering perform_inversion '% (timestamp()))
        import sys
        sys.stdout.flush()
    Ip, phase, amp, phase_err_inst, phase_err_striation, phase_err_zero_phase, phase_err_signal, amp_uncertainty, chi2 = perform_inversion(I, 
                           tang_alt, icon_alt, 
                           err_inst, err_striation, err_zero_phase, err_signal, I_amp_uncertainty,
                           top_layer=top_layer, integration_order=integration_order,
                           account_for_local_projection=account_for_local_projection, linear_amp=linear_amp, H=H,
                           unwrap_phase = unwrap_phase, chi2_thresh = chi2_thresh, Ncorr = Ncorr)
    
    #### Transform from phase to wind
    f = phase_to_wind_factor(np.mean(sigma_opd)) # Use average OPD to analyze entire row
    v             = f * phase
    v_err_inst       = f * phase_err_inst
    v_err_striation  = f * phase_err_striation
    v_err_zero_phase = f * phase_err_zero_phase
    v_err_signal     = f * phase_err_signal
    low_signal_wind_corr = f * low_signal_phase_corr
    
    #### Handle uncertainty reporting
    v_err = np.sqrt(v_err_inst**2 + v_err_striation**2) # Included for backwards compatibility and for recommended error for tidal analysis
    
    #### Transform from fringe amplitude to VER
    r = np.zeros(ny) # amp-to-R factor, altitude dependent
    for i in range(ny):
        T = visibility_temperature(tmid, lat[i], lon[i], alt[i], f107=f107, f107a=f107a, f107p=f107p, apmsis=apmsis)
        T = T*Tmult
        r[i] = amp_to_R_factor(np.mean(sigma_opd), T)
    ver = r * amp
    ver_uncertainty = r * amp_uncertainty    
    
    #### Invert DC value of interferogram to get alternative VER product
    D = create_observation_matrix(tang_alt, icon_alt, top_layer=top_layer, integration_order=integration_order, H=H)
    ver_dc = np.linalg.solve(D, I_dc)
    
    #### Calculate azimuth angles at measurement locations
    az = los_az_angle(icon_latlonalt, lat, lon, alt)
        
    
    # Make a L2.1 dictionary
    L21_dict = {
             'los_wind'                     : v,
             'los_wind_error_inst'          : v_err_inst,
             'los_wind_error_striation'     : v_err_striation,
             'los_wind_error_zero_phase'    : v_err_zero_phase,
             'los_wind_error_signal'        : v_err_signal,
             'los_wind_error'               : v_err,
             'lat'                          : lat,
             'lon'                          : lon,
             'alt'                          : alt,
             'time_start'                   : L1_dict['time_start'],
             'time_stop'                    : L1_dict['time_stop'],
             'time'                         : tmid,
             'exp_time'                     : exp_time,
             'az'                           : az,
             'emission_color'               : emission_color,
             'sensor'                       : sensor,
             'icon_alt'                     : icon_alt,
             'icon_lat'                     : icon_lat,
             'icon_lon'                     : icon_lon,
             'ver'                          : ver,
             'ver_error'                    : ver_uncertainty,
             'fringe_amplitude'             : amp,
             'fringe_amplitude_error'       : amp_uncertainty,
             'mighti_ecef_vectors'          : mighti_ecef_vectors_center,
             'icon_velocity_ecef_vector'    : icon_velocity_vector,
             'file_creation_time'           : datetime.now(),
             'source_files'                 : [L1_fn],
             'bin_size'                     : bin_size,
             'top_layer'                    : top_layer,
             'H'                            : H,
             'integration_order'            : integration_order,
             'I'                            : Ip,
             'chi2'                         : chi2,
             'wind_quality'                 : np.nan, # A placeholder; this will be filled in below
             'ver_quality'                  : np.nan, # A placeholder; this will be filled in below
             'quality_flags'                : None,   # A placeholder; this will be filled in below
             'acknowledgement'              : L1_dict['acknowledgement'],
             'att_lvlh_normal'              : L1_dict['att_lvlh_normal'],
             'att_lvlh_reverse'             : L1_dict['att_lvlh_reverse'],
             'att_limb_pointing'            : L1_dict['att_limb_pointing'],
             'att_conjugate'                : L1_dict['att_conjugate'],
             'att_zero_wind'                : L1_dict['att_zero_wind'],
             'mag_lat'                      : mag_lat,
             'mag_lon'                      : mag_lon,
             'sza'                          : sza,
             'slt'                          : slt,
             'ver_dc'                       : ver_dc,
             'orbit_node'                   : int(icon_velocity_vector[2] < 0), # ECEF_Z < 0 --> descending
             'orbit_number'                 : L1_dict['orbit_number'],
             'jitter'                       : L1_dict['jitter'],
             'zero_wind_phase'              : zero_wind_phase,
             'low_signal_corr'              : low_signal_wind_corr,
    }
    
    #### Quality control and flagging
    if global_params['verbose']:
        print('%s:\t\t\t Computing quality flags'% (timestamp()))
        import sys
        sys.stdout.flush()
    if global_params['verbose']:
        print('%s:\t\t\t\t Entering level21_quality(...) '% (timestamp()))
        import sys
        sys.stdout.flush()
    wind_quality, ver_quality, quality_flags = level21_quality(L1_quality_flags, L21_dict, L1_quality=L1_quality,
                                                               top_layer_thresh=top_layer_thresh, 
                                                               terminator_thresh = terminator_thresh, 
                                                               chi2_thresh = chi2_thresh, chi2_thresh_caution=chi2_thresh_caution,
                                                               count_thresh = count_thresh)
    if global_params['verbose']:
        print('%s:\t\t\t\t Exited level21_quality(...)'% (timestamp()))
        import sys
        sys.stdout.flush()
    L21_dict['wind_quality'] = wind_quality
    L21_dict['ver_quality'] = ver_quality
    L21_dict['quality_flags'] = quality_flags
    
    ### Mask out points with very bad quality
    if global_params['verbose']:
        print('%s:\t\t\t\t Finding indices of low quality wind data'% (timestamp()))
        import sys
        sys.stdout.flush()
    idx = L21_dict['wind_quality'] == 0.0
    if global_params['verbose']:
        print('%s:\t\t\t\t Masking those indices'% (timestamp()))
        import sys
        sys.stdout.flush()
    L21_dict['los_wind'][idx] = np.nan
    L21_dict['los_wind_error_inst'][idx] = np.nan
    L21_dict['los_wind_error_striation'][idx] = np.nan
    L21_dict['los_wind_error_zero_phase'][idx] = np.nan
    L21_dict['los_wind_error_signal'][idx] = np.nan
    L21_dict['los_wind_error'][idx] = np.nan
    L21_dict['low_signal_corr'][idx] = np.nan
    
    if global_params['verbose']:
        print('%s:\t\t\t\t Finding indices of low quality VER data'% (timestamp()))
        import sys
        sys.stdout.flush()
    idx = L21_dict['ver_quality'] == 0.0
    if global_params['verbose']:
        print('%s:\t\t\t\t Masking those indices'% (timestamp()))
        import sys
        sys.stdout.flush()
    L21_dict['fringe_amplitude'][idx] = np.nan
    L21_dict['fringe_amplitude_error'][idx] = np.nan
    L21_dict['ver'][idx] = np.nan
    L21_dict['ver_error'][idx] = np.nan
    L21_dict['ver_dc'][idx] = np.nan
    
    if global_params['verbose']:
        print('%s:\t\t\t Returning from L2.1 core'% (timestamp()))
        import sys
        sys.stdout.flush()
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
        raise Exception('field_name="%s" is too long (%i chars). Shorten to 30 characters.' % (field_name,len(field_name)))
    if len(label_axis) > 10:
        raise Exception('label_axis="%s" is too long (%i chars). Shorten to 10 characters.' % (label_axis,len(label_axis)))
    if monoton not in ['Increase', 'Decrease', None]:
        raise Exception('Input "monoton" must be either "Increase", "Decrease" or None.')
    
    # If fill value was not specified, use the default value, if it exists.
    # It will not exist for strings, for example, for which fill values
    # cannot be set.
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
            value[np.isnan(value)]  = var._FillValue
        # For non-sequences and strings:
        elif np.isnan(value):
            value = var._FillValue
    
    # Assign value
    var[...] = value
    
    return var

   
    
    
    
def save_nc_level21(path, L21_dict, data_revision=0):
    '''
    Take the output of the Level 2.1 processing for many files, and save it as a NetCDF4 file in the official format.
    NetCDF4 file conventions are taken from "Science Operations Center Data Product Conventions" Rev 0.5.
    
    INPUTS:
    
      *  path        -- TYPE:str.  The directory the file will be saved in, including trailing "/" (e.g., '/home/user/')
      *  L21_dict    -- TYPE:dict. A dictionary containing output variables of the Level 2.1 processing. Many exposures can be 
                                   saved in the same file, but all should come from the same channel (red or green) of the same 
                                   sensor (A or B), and have times during the same date.
                                   The variables in this dictionary should match those described in the documentation for 
                                   level1_dict_to_level21_dict(...) but should have an extra dimension (on axis 0) for time,
                                   and optionally includes the "zero_wind_file" key to be saved in the "Calibration_File" attribute.
                                   Variables which do not change in time (e.g., emission color) should not have this extra 
                                   dimension. In v1.31 release, 2 extra keys are possible: ph0_day and ph0_night. If these
                                   are included, they will be saved to the file as global attributes.
                                   
    OPTIONAL INPUTS:
    
      *  data_revision -- TYPE:int,  The revision number for the data. This will be put in the filename (v01r###.NC) (default 0)
                                   
    OUTPUTS:
    
      *  L21_fn      -- TYPE:str.  The full path to the saved file.
      
    TO-DO:
    
      * How can we append to global attributes History and MODS when the processing is re-run?
      
    '''
    
    assert np.ndim(L21_dict['los_wind']) > 1, "Are you accidentally saving one inversion instead of many?"
    
    L21_dict = copy.deepcopy(L21_dict) # because netCDF4 seems to change the input when nans are involved
    data_version_major = software_version_major # enforced as per Data Product Conventions Document
    
    nt, ny, ne = np.shape(L21_dict['quality_flags']) # grab all dimensions
    _,  nyu,   = np.shape(L21_dict['zero_wind_phase']) # grab row dimension
    
    #################### Compile variables to write in file ######################
    ### Sensor:
    sensor = L21_dict['sensor']
    ### Timing:
    # 2019-12-16: This section was commented out because the SDC sometimes wants to include 
    # 1 time sample from the previous day.
#     date0 = L21_dict['time'][0].date()
#     for i in range(nt):
#         date1 = L21_dict['time'][i].date()
#         assert date0 == date1, 'Files from different dates: %s %s' % (date0, date1)
    date_ref = L21_dict['time'][nt/2].date() # Assume middle time provides the definition of "today's date"
    t_start_msec = np.array([(t - datetime(1970,1,1)).total_seconds()*1e3 for t in L21_dict['time_start']]).astype(np.int64) # milliseconds since epoch
    t_stop_msec  = np.array([(t - datetime(1970,1,1)).total_seconds()*1e3 for t in L21_dict['time_stop']]).astype(np.int64)
    t_mid_msec   = np.array([(t - datetime(1970,1,1)).total_seconds()*1e3 for t in L21_dict['time']]).astype(np.int64)
    t_file  = datetime.now()   # time this file was created
    ### Who's running this process
    user_name = getpass.getuser()
    ### Parent files
    parents = [] # This will go in global attr Parents
    for source_fn in L21_dict['source_files'].flatten():
        s = source_fn.split('/')[-1].split('.')
        pre = '.'.join(s[:-1])
        post = s[-1].upper()
        parents.append('%s > %s' % (post, pre))
    ### Calibration files
    calibration_file = ''
    if 'zero_wind_file' in L21_dict:
        calibration_file = L21_dict['zero_wind_file']


    L21_fn = 'ICON_L2-1_MIGHTI-%s_LOS-Wind-%s_%s_v%02ir%03i.NC' % (sensor,L21_dict['emission_color'].capitalize(),
                                                           date_ref.strftime('%Y-%m-%d'),
                                                           data_version_major, data_revision)
    L21_full_fn = '%s%s'%(path, L21_fn)
    ncfile = netCDF4.Dataset(L21_full_fn,mode='w',format='NETCDF4') 

    try:
        ########################## Global Attributes #################################
        ncfile.setncattr_string('Acknowledgement',                L21_dict['acknowledgement'])
        ncfile.setncattr_string('ADID_Ref',                       'NASA Contract > NNG12FA45C')
        ncfile.setncattr_string('Calibration_File',               calibration_file)
        ncfile.setncattr_string('Conventions',                    'SPDF ISTP/IACG Modified for NetCDF')
        ncfile.setncattr_string('Data_Level',                     'L2.1')
        ncfile.setncattr_string('Data_Type',                      'DP21 > Data Product 2.1: Line-of-sight Wind Profiles')
        ncfile.Data_Version_Major =                               np.uint16(data_version_major)
        ncfile.Data_Revision =                                    np.uint16(data_revision)
        ncfile.Data_Version =                                     data_version_major + 0.001 * data_revision
        ncfile.setncattr_string('Date_End',                       L21_dict['time'][-1].strftime ('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Date_Start',                     L21_dict['time'][0].strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Description',                    'ICON MIGHTI Line-of-sight Winds (DP 2.1)')
        ncfile.setncattr_string('Descriptor',                     'MIGHTI-%s > Michelson Interferometer for Global High-resolution ' % sensor+\
                                                                  'Thermospheric Imaging, Sensor %s' % sensor)
        ncfile.setncattr_string('Discipline',                     'Space Physics > Ionospheric Science')
        ncfile.setncattr_string('File',                           L21_fn)
        ncfile.setncattr_string('File_Date',                      t_file.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Generated_By',                   'ICON SDC > ICON UIUC MIGHTI L2.1 Processor v%s, B. J. Harding' % __version__)
        ncfile.setncattr_string('Generation_Date',                t_file.strftime('%Y%m%d'))
        ncfile.setncattr_string('History',                        ['v1.0: First release of MIGHTI L2.1/L2.2 software, B. J. Harding, 05 Mar 2018', 
                                                                   'v2.0: First run of on-orbit data, using external zero wind reference and smooth daily-averaged profiles, B. J. Harding, 01 May 2020',
                                                                   'v3.0: Correction for long-term mechanical drift, B. J. Harding, 04 Jun 2020',
                                                                   'v4.0: Updated correction for long-term mechanical drift to handle settling after ~May 2020 and precession cycle variation. LoS winds have changed by a bulk offset of up to 30 m/s. Studies using only perturbations from the mean (e.g., non-migrating tidal retrievals) are unlikely to be affected. B. J. Harding, 21 Oct 2020',
                                                                   'v5.0: The ad-hoc, HWM-based correction for the zero wind phase has been replaced with a self-calibration based on comparing data from the ascending and descending orbits (see the notes for the wind variables below for details). Long-term trends in the zero wind phase degraded the accuracy of version 4 over time, and the accuracy of version 4 data was tied to the accuracy of HWM. In version 5, errors on these long time scales (>100-150 days) are now accounted for, improving the accuracy to 10-25 m/s (see the "Accuracy" variable for more details) and removing the dependence on external models. A long window of data is required to implement this self-calibration, so v05 data are processed at least 100 days behind real time. For errors on precession-cycle time scales (48 days), the previous ad-hoc correction using initial red-vs-green comparisons has been replaced with a more comprehensive red-green cross-calibration (165-185 km altitude during the day) that accounts for the time-dependence of mechanical drifts of the optics. This result is consistent with a first-principles analysis of the fiducial notch positions (see Marr et al., 2020 and subsequent publications). A mission-average fiducial notch analysis is also used to correct mechanical drifts on an orbital time scale (97 minutes, or 24 hours of local time), which could affect migrating tide estimates. The RMS difference due to this effect is estimated at 5-10 m/s (root mean square). Analysis of waves with periods that do not coincide with these new corrections are not likely to be different than in version 4 (e.g., nonmigrating tides, planetary waves, and gravity waves). New variables related to error (i.e., uncertainty) estimates from various sources are now included, whereas version 4 error estimates only included the effects of shot, read, and dark noise. MIGHTI-A and MIGHTI-B variables related to emission brightness are now cross-calibrated, though not absolutely calibrated. Exposures affected by solar and lunar stray light are now flagged. Some data during periods in May and July 2020 when the sun approached the MIGHTI field of view was marked as unavailable in v04, but is now available in v05 with the exception of a few days. The data from the second row in the green channel (~91 km) is now available when the signal strength permits a wind estimate. An error in the local time calculation has been corrected, which changes the local time by up to 20 min. A new algorithm to identify cosmic ray spikes has been implemented, improving precision. A preliminary algorithm has been implemented to to correct a wind bias associated with low signal levels, and associated uncertainties are estimated (see the "Precision_Low_Signal_Effect" variable for more details). Finally, various quality control parameters have been optimized. More description is provided in the notes below. A full history of software changes can be found on Github: https://github.com/bharding512/airglowrsss/commits/master/Python/modules/MIGHTI_L2.py  B. J. Harding 08 Sep 2022'
                                                                  ])
        ncfile.setncattr_string('HTTP_LINK',                      'http://icon.ssl.berkeley.edu/Instruments/MIGHTI')
        ncfile.setncattr_string('Instrument',                     'MIGHTI-%s' % sensor)
        ncfile.setncattr_string('Instrument_Type',                'Imagers (space)')
        ncfile.setncattr_string('Link_Text',                      'MIGHTI Line-of-sight Wind Profiles (DP 2.1)')
        ncfile.setncattr_string('Link_Title',                     'ICON MIGHTI')
        ncfile.setncattr_string('Logical_File_ID',                L21_fn[:-3])
        ncfile.setncattr_string('Logical_Source',                 'ICON_L2_MIGHTI-%s_' % (sensor,))
        ncfile.setncattr_string('Logical_Source_Description',     'MIGHTI Sensor %s - Line-of-sight Wind Profiles' % (sensor,))
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
        text_supp = ["This data product contains altitude profiles of the line-of-sight winds (inverted wind profiles in the direction of the sensor's line of "
                    "sight) for 24 hours of data taken by MIGHTI. In addition to the line-of-sight wind data and the corresponding ancillary data, "
                    "such as time and location, this product contains supporting data, such as fringe amplitude profiles and "
                    "relative volume emission rate profiles. Absolute calibration and MIGHTI-A/B cross calibration of these data "
                    "is not necessary to obtain the wind data, and therefore any direct analysis of these supporting data requires caution. ",
                     
                    "There is one file for each sensor (A or B), for each color (red or green) and for each day. The profile "
                    "spans from an altitude of ~90 km (for green) or ~150 km (for red) to ~300 km, though altitudes with low signal levels are masked out. This data product is "
                    "generated from the Level 1 MIGHTI product, which comprises calibrated interference fringe amplitudes and phases. The effect of "
                    "spacecraft velocity is removed from the "
                    "interferogram phase, then (optionally) the data are binned from their native altitude sampling (~2.5 km) to improve statistics. "
                    "An onion-peeling inversion is performed to remove the effect of the line-of-sight integration. After the inversion, each row (i.e., altitude) "
                    "is analyzed to extract the phase, and thus the line-of-sight wind. Level 2.1 files from MIGHTI-A and MIGHTI-B are combined during the Level 2.2 "
                    "processing (not discussed here). See Harding et al. [2017, doi:10.1007/s11214-017-0359-3] for more details of the inversion algorithm. "
                    "Further discussion of the calibration and performance of MIGHTI after launch can be found in a forthcoming paper in Space Science "
                    "Reviews [Englert et al., 2022, in preparation].",
                    
                     "Known issues with the v05 data release are listed below. Work is in progress to resolve or mitigate these issues in future data releases. ",

                     "Known issues with v05:<br/>"
                     " * When ICON is in the South Atlantic Anomaly (SAA), radiation effects on the detector cause poor data quality. The quality control algorithm "
                     "adequately flags and masks most of the affected samples, but some outliers remain, especially near the edge of the SAA. Other uncaught outliers "
                     "are rare but can occur due to cosmic rays, stars in the field of view, moonlight, etc. <br/>"
                     " * The bottom row of data (corresponding to an altitude of ~88 km) is masked out as the signal is rarely strong enough to permit a wind "
                     "estimate, and calibrations have large uncertainties. It is unlikely but possible that this altitude will be reported in future releases, "
                     "pending further investigation. <br/>"
                     " * Airglow brightness observations are not a required mission product, and no effort was yet made to absolutely calibrate "
                     "the brightness observations for MIGHTI-A and MIGHTI-B, and thus the Relative_VER "
                     "variable should be treated with caution. In v05, MIGHTI-A and B are cross-calibrated using a conversion factor derived from on-orbit data. "
                     "However, there are some indications that this cross-calibration may be changing with time, which is not accounted for in v05. <br/>"
                     
                     " * As discussed in the variable notes below, a new zero wind phase determination has been implemented in v05. However, during the period "
                     "from 2021 Apr 26 to Aug 14, data gaps and one period of southward (\"Reverse LVLH\") pointing cause errors in this determination. The "
                     "accuracy is estimated to be degraded by a factor of two. See the *_Accuracy variable. <br/>"
                     " * During the one orbit per day when the calibration lamp is on, the wind data can be noisier and have a slight bias. Although this issue is "
                     "much improved since v04, for the sake of conservatism, these orbits are still labeled with quality=0.5 (i.e., caution). <br/>"
                     " * Some data gaps appear on days when the sun passes near MIGHTI's field of view. Most of these gaps are located near the terminator, "
                     "but some are longer lasting.<br/>"
                     " * In some cases, there are indications that the *_Precision_1_Sample variables are underestimating the true sample-to-sample noise, "
                     "suggesting that, in addition to shot noise, there is a second source of noise. It is recommend that any quantitative use of the reported "
                     "error estimates (i.e., precision and accuracy) should treat those estimates as uncertain. It is believed that most error estimates are "
                     "correct to within a factor of 2. The largest problems with error reporting occur where the airglow signal is weakest. <br/>"
                     
                     " * Imperfect daily calibration data lead to small discontinuities in the zero wind phase at the boundaries between days (i.e., "
                     " between 23:59:59 and 00:00:00 UT), which are not accounted for by the reported error variables. This was estimated to be a 2-5 "
                     "m/s (root-mean-square) error early in the mission, but is growing over time, possibly reaching 5-10 m/s by mid-2022.<br/>"
                     " * A signal-dependent phase shift is seen in atmospheric and calibration lamp fringes, possibly caused by "
                     "a charge trapping effect in the CCD. This is the subject of ongoing investigation, but a first-order correction "
                     "is implemented in the v05 dataset. The correction increases linearly with time to match the effect seen in on-orbit calibration data. "
                     "The variable *_Precision_Low_Signal_Effect is an estimate of the remaining uncertainty due to an imperfect "
                     "correction. Where this uncertainty is large, caution is recommended. "
                     "For example, for winds in the core science region (90-105 km altitude, away from the terminator), the magnitude of the correction "
                     "is small or zero, but data in the red "
                     "channel during the night and twilight are subject to a large correction (many tens of m/s) and the uncertainty is correspondingly large. "
                     "A goal for future releases is to characterize and correct this effect more accurately. <br/>"
                     " * Data near the solar terminators are subject to a variety of errors, including those described above and others related to the rapidly "
                     "varying illumination. Not all errors near the terminator are accounted for by the reported error. Users are encouraged to "
                     "use extra caution with these data.<br/>"
                     "See the documentation below for more information.<br/>"
                    ]

        ncfile.setncattr_string('Text_Supplement',                text_supp)
        ncfile.setncattr_string('Time_Resolution',                '30 - 60 seconds')
        ncfile.setncattr_string('Title',                          'ICON MIGHTI Line-of-sight Wind Profiles (DP 2.1)')
        ncfile.setncattr_string('Zero_Wind_Ref',                  'DEPRECATED - see Zero_Wind_Phase variable') # deprecated as of v05, may be removed in future
        # Added in v1.31, but may be removed in the future once smoothing is no longer necessary
        ncfile.setncattr_string('Smooth_Profile',                 'DEPRECATED') # deprecated as of v05, may be removed in future
        ncfile.setncattr_string('Smooth_Profile_Night',           'DEPRECATED') # deprecated as of v05, may be removed in future
        ncfile.setncattr_string('Smooth_Profile_Day',             'DEPRECATED') # deprecated as of v05, may be removed in future
        # Added in v2.03, but may be removed in future once notch drift is accounted for at Level 1
        ncfile.setncattr_string('Corr_Notch_Drift',               'DEPRECATED') # deprecated as of v05, may be removed in future
        ncfile.setncattr_string('Low_Signal_Patch',               'Automatic') # Indicating this file was created after the low-signal bug was identified/fixed, but was not manually patched.


        ################################## Dimensions ########################################
        prefix = 'ICON_L21' # prefix of each variable
        ncfile.createDimension('Epoch',None)
        ncfile.createDimension('Altitude', ny)
        ncfile.createDimension('Vector',3)
        ncfile.createDimension('Start_Mid_Stop',3)
        ncfile.createDimension('N_Flags', ne)
        ncfile.createDimension('Row', nyu)
        
        
        ################################## Variables #########################################
        
        
        ######### Timing Variables #########

        # Time midpoint (the official required "Epoch" variable)
        var = _create_variable(ncfile, 'Epoch', t_mid_msec, 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i8', format_fortran='I', desc='Sample time, midpoint of exposure. Number of msec since Jan 1, 1970.', 
                              display_type='time_series', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
                              units='ms', valid_min=0, valid_max=1000*365*86400*1000, var_type='support_data', chunk_sizes=[nt],
                              notes="This variable contains the time corresponding to the wind profiles reported in this file, taken at the "
                                    "midpoint of the exposure time. It is in UTC and has units of milliseconds since "
                                    "Jan 1, 1970. A human-readable version of the time can be found in the variable ICON_..._UTC_Time"
                              )

        # Time start/mid/stop
        var = _create_variable(ncfile, '%s_Time'%prefix, np.stack([t_start_msec, t_mid_msec, t_stop_msec], axis=1),
                              dimensions=('Epoch','Start_Mid_Stop'), depend_0 = 'Epoch',
                              format_nc='i8', format_fortran='I', desc='Sample time at start, mid, stop of exposure. Number of msec since Jan 1, 1970.', 
                              display_type='time_series', field_name='Time', fill_value=None, label_axis='Time', bin_location=[0.0,0.5,1.0],
                              units='ms', valid_min=0, valid_max=1000*365*86400*1000, var_type='support_data', chunk_sizes=[nt,3],
                              notes="This variable is the same as Epoch, except it has another dimension which holds the start time, middle time, "
                                    "and stop time of each exposure."
                              )
        
        # Human readable version of midpoint time
        # YYYY-MM-DD HH:mm:ss.xxxZ as per DP Conventions 0.7
        var = _create_variable(ncfile, '%s_UTC_Time'%prefix, np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3] + 'Z' for t in L21_dict['time']]),
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc=str, format_fortran='A', desc='Sample time, midpoint of exposure.', 
                              display_type='time_series', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='support_data', chunk_sizes=[nt],
                              notes="This variable is the same as Epoch but is formatted as a human-readable string.")
                               
            
        ######### Data Variables #########

        # Line-of-sight wind profile -- this one is special because it has var notes that vary depending on the zero_wind_ref
        varnotes = ["The wind is the primary data product in this file. This variable contains the projection of the horizontal wind "
                  "(at the tangent point) onto the line of sight direction. An entire altitude profile is observed simultaneously. An onion-peeling "
                  "inversion is used on the raw observations to remove the effects of the integration along the line of sight. The "
                  "line-of-sight wind is defined such that a positive value indicates motion towards the spacecraft. This direction is given "
                  "by the Line_of_Sight_Azimuth variable. It is assumed that the vertical wind is zero, but even large vertical winds "
                  "(~100 m/s) do not significantly affect this data product, since the line of sight is nearly horizontal everywhere. It "
                  "should be noted that while this measurement is ascribed to a particular latitude, longitude and altitude, it is actually "
                  "an average over many hundreds of kilometers horizontally, and 2.5-30 kilometers vertically (depending on the binning). "
                  "See Harding et al. [2017, doi:10.1007/s11214-017-0359-3] for a more complete discussion of the inversion algorithm."]


        varnotes.extend(["Knowledge of the \"zero wind phase\" is needed for any instrument using Doppler shifts to determine winds. "
                         "The zero wind phase is defined as the measured interference fringe phase that corresponds to the rest "
                         "wavelength of the emission. For the v05 data release, the zero wind phase has been determined by considering a window of LoS wind data "
                         "spanning two precession cycles (96 days). Assuming that on average the real zonal and meridional winds do not depend on the aspect "
                         "angle with which MIGHTI observes the atmosphere (an angle which is significantly different between the ascending and descending portions of the orbit), "
                         "a matrix equation can be constructed which combines data from both MIGHTI-A and "
                         "MIGHTI-B and both the ascending and descending orbits. This equation is solved for the average zonal and meridional wind, and the "
                         "zero wind phase for MIGHTI-A and MIGHTI-B. This window is moved in time to determine the appropriate zero wind phase for each date. "
                         "The value of the zero wind phase depends on emission color (red or green), aperture mode (day or night), calibration lamp status "
                         "(on or off) and row (i.e., altitude). "
                         "An additional zero-mean signal is added to this result to ensure that 48-day (i.e., 1 precession cycle) average winds are smooth "
                         "in altitude. Adjustments are smaller than the reported accuracy, so this adjustment is not "
                         "expected to change any scientific conclusions, although it does ensure more realistic wind profiles. "
                         "This is a less restrictive assumption than the smoothness criterion used in v04, which "
                         "relied on the Horizontal Wind Model 2014 and also enforced smoothness on 1-day averages. It is thus expected that the amplitude of tidal structures "
                         "in the lower thermosphere are subject to less suppression in v05+ than in v04. "
#                          "The results of this approach are consistent with the results of ICON's \"zero wind maneuver.\" " # 2022 Sep 9: We don't know if this is true yet.
                         "This version of the MIGHTI zero wind phase is independent of any external data or models (such as "
                         "the Horizontal Wind Model 2014, which was used in v04 and earlier versions). "
                         "The zero wind phase used for each wind sample is saved in the _Zero_Wind_Phase variable below. "
                         "The 1-sigma uncertainty in the winds incurred by the inaccuracy in the zero wind phase "
                         "is estimated and reported in the *_Wind_Accuracy variable below."])
                        
        var = _create_variable(ncfile, '%s_Line_of_Sight_Wind'%prefix, L21_dict['los_wind'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Line-of-sight horizontal wind profile. A positive wind is towards MIGHTI.', 
                              display_type='altitude_profile', field_name='Line-of-sight Wind', fill_value=None, label_axis='LoS Wind', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[nt,ny],
                              notes=varnotes)
        
        # Line-of-sight wind errors
        
        err_prefix = ["Various sources of error in MIGHTI winds are quantified with 1-sigma estimates and organized by their temporal persistence. " 
                      "These error sources are nearly uncorrelated with each other and can thus be added in quadrature. Users are encouraged to "
                      "contact the MIGHTI team for assistance with error propagation."]
        
        # Precision_1_Sample
        var = _create_variable(ncfile, '%s_Line_of_Sight_Wind_Precision_1_Sample'%prefix, L21_dict['los_wind_error_inst'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Line-of-sight wind precision (1 sample)', 
                              display_type='altitude_profile', field_name='LOS Wind Precision (1 sample)', fill_value=None, label_axis='Wind Prec', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nt,ny],
                              notes=[err_prefix[0],
                                     "The \"1 sample\" error variable quantifies errors that are uncorrelated from one exposure to the next, "
                                     "dominated by shot and dark noise in the detectors.  The correlation time of this error source is 30-60 "
                                     "seconds (i.e., the measurement cadence). The reported error is estimated from the fringe intensity and "
                                     "background. This is the recommended variable to use for analyses of wind fluctuations within a single day "
                                     "and a single altitude (e.g., gravity waves). Because the Level 2.2 data include interpolation of Level 2.1 "
                                     "data, some correlation remains between consecutive samples. Errors are slightly correlated across small "
                                     "altitude gaps as a result of the inversion."]
                              )
        
        # Precision_1_Day
        var = _create_variable(ncfile, '%s_Line_of_Sight_Wind_Precision_1_Day'%prefix, L21_dict['los_wind_error_striation'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Line-of-sight wind precision (1 day)', 
                              display_type='altitude_profile', field_name='LOS Wind Precision (1 day)', fill_value=None, label_axis='Wind Prec', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nt,ny],
                              notes=[err_prefix[0],
                                     "The \"1 Day\" error variable quantifies the error introduced by daily calibrations, which is correlated "
                                     "for an entire 24-hour period (00:00 - 23:59 UT). This is estimated from the magnitude of fluctuations "
                                     "in the daily-averaged phase, propagated through the inversion. Errors in day mode and night mode are "
                                     "nearly uncorrelated. For studies pertaining to "
                                     "atmospheric tidal modes that combine data from many days, this error can be treated as uncorrelated "
                                     "across time."]
                              )        
        
        # Precision_low_signal_effect
        var = _create_variable(ncfile, '%s_Line_of_Sight_Wind_Precision_Low_Signal_Effect'%prefix, L21_dict['los_wind_error_signal'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Line-of-sight wind precision (from low signal effect)', 
                              display_type='altitude_profile', field_name='LOS Wind Precision (Low Sig)', fill_value=None, label_axis='Wind Prec', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nt,ny],
                              notes=[err_prefix[0],
                                     "The \"Low Signal Effect\" error variable quantifies the error introduced by the imperfect correction for the "
                                     "signal-dependent phase shift, which is an effect seen in atmospheric and calibration-lamp fringes where "
                                     "the phase of the fringes is biased at very low signal levels. This is under investigation but could "
                                     "be caused by a charge trapping effect in the CCD. A "
                                     "correction has been implemented based upon the empirical relationship between measured phase "
                                     "and signal level of the calibration lamps for the first ~30 months of the mission. However, especially for "
                                     "cases with low signal levels, "
                                     "this correction is uncertain. The uncertainty in the resulting winds is estimated from the signal level "
                                     "and reported in this variable. It is likely to be correlated across samples nearby in time and space, but "
                                     "the correlation between different channels (red and green), sensors (MIGHTI-A and MIGHTI-B), and operating "
                                     "modes (Day and Night) is not known. Depending on the analysis being used, it could be treated as a systematic "
                                     "error or as a statistical error. Where this uncertainty is large, caution is recommended. "
                                     "For example, for winds in the core science region (90-105 km altitude), the magnitude of the correction is small or zero, "
                                     "but data in the red channel during the night and twilight are subject to a large correction (many tens of m/s) "
                                     "and the uncertainty is correspondingly large. "
                                     "A goal for future releases is to characterize and correct this effect more accurately."]
                              )
        
        # Accuracy
        var = _create_variable(ncfile, '%s_Line_of_Sight_Wind_Accuracy'%prefix, L21_dict['los_wind_error_zero_phase'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Line-of-sight wind accuracy', 
                              display_type='altitude_profile', field_name='LOS Wind Precision (1 day)', fill_value=None, label_axis='Wind Prec', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nt,ny],
                              notes=[err_prefix[0],
                                     "The \"Accuracy\" variable quantifies the error introduced by the zero-wind phase estimate. It is "
                                     "strongly correlated across time lags of days to weeks and becomes increasingly decorrelated for time "
                                     "lags longer than 2 precession cycles (96 days). This error is estimated from the discrepancy between "
                                     "various techniques of determining the zero-wind phase. This error source is irrelevant for most users "
                                     "studying perturbations from the mean (e.g., tides, waves), but may be important for studies of "
                                     "zonal mean winds, point comparisons with other data sets, and seasonal/long-term trends thereof. Errors are "
                                     "moderately correlated across "
                                     "small altitude gaps. Errors in day mode and night mode are nearly uncorrelated, implying there could "
                                     "be different offsets for day mode and night mode. This could be important for error propagation of "
                                     "odd-numbered migrating tides (e.g., DW1)."]
                              )
        
        # Error: Total/Summary
        var = _create_variable(ncfile, '%s_Line_of_Sight_Wind_Error'%prefix, L21_dict['los_wind_error'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Line-of-sight horizontal wind error profile', 
                              display_type='altitude_profile', field_name='Line-of-sight Wind Error', fill_value=None, label_axis='Wind Error', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nt,ny],
                              notes="For robust error propagation, users are encouraged to consider the individual error variables: "
                               "\"Precision_1_Sample,\" \"Precision_1_Day,\" and \"Accuracy.\" The \"Wind_Error\" variable is "
                               "included for backwards compatibility and is equal to the quadrature sum of the \"1 Sample\" error "
                               "and the \"1 Day\" error. This is the recommended uncertainty to use for analyses that collect data "
                               "from several weeks and compute perturbations from the mean (e.g., for estimating tides and planetary "
                               "waves). This error is uncorrelated across time lags larger than 24 hours. Errors are slightly correlated "
                               "across small altitude gaps. Errors in day mode and night mode are nearly uncorrelated."
)
        
        # Quality code
        var = _create_variable(ncfile, '%s_Wind_Quality'%prefix, L21_dict['wind_quality'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='A quantification of the wind quality, from 0 (Bad) to 1 (Good)', 
                              display_type='altitude_profile', field_name='Wind Quality', fill_value=None, label_axis='Quality', bin_location=0.5,
                              units='', valid_min=0.0, valid_max=1.0, var_type='data', chunk_sizes=[nt,ny],
                              notes=["A quantification of the overall quality of the wind data. "
                                     "While the intent is that the variable ICON_...Line_of_Sight_Wind_Error accurately characterizes the statistical "
                                     "error in the wind data, it is possible that systematic errors are present, or that the statistical error "
                                     "estimation is not accurate. If it is suspected that this is the case, the quality will be less than 1.0. If "
                                     "the data are definitely unusable, the quality will be 0.0 and the sample will be masked. Users should "
                                     "exercise caution when the quality is less than 1.0.",
                                     "This parameter can currently take 3 values: 0.0 (Bad), 0.5 (Caution), 1.0 (Good)"
                                     ])
        
        # Fringe amplitude profile 
        var = _create_variable(ncfile, '%s_Fringe_Amplitude'%prefix, L21_dict['fringe_amplitude'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Fringe amplitude profile', 
                              display_type='altitude_profile', field_name='Fringe Amplitude', fill_value=None, label_axis='Fringe Amp', bin_location=0.5,
                              units='arb', valid_min=-1e30, valid_max=1e30, var_type='data', chunk_sizes=[nt,ny],
                              notes="An approximate volume emission rate (VER) profile in arbitrary units. Technically this a profile of the "
                              "amplitude of the fringes, which has a dependence on thermospheric temperature and background emission. Thus, it does not "
                              "truly represent volume emission rate. However, it is a useful proxy. The units are arbitrary, but an attempt has "
                              "been made to cross-calibrate MIGHTI-A with MIGHTI-B. In contrast to the wind inversion, which is nonlinear due to the "
                              "phase extraction step, the amplitude inversion is purely linear. The Level 1 interferogram is analyzed to obtain a single "
                              "brightness value per zenith angle, and this is inverted with the distance matrix to obtain a value of the amplitude "
                              "per altitude."
                              )

        # Fringe amplitude error profile
        var = _create_variable(ncfile, '%s_Fringe_Amplitude_Error'%prefix, L21_dict['fringe_amplitude_error'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Fringe amplitude error profile', 
                              display_type='altitude_profile', field_name='Fringe Amplitude Error', fill_value=None, label_axis='Amp Err', bin_location=0.5,
                              units='arb', valid_min=0.0, valid_max=1e30, var_type='data', chunk_sizes=[nt,ny],
                              notes="The statistical (1-sigma) error in the fringe amplitude. As with the wind, systematic errors are not "
                              "included, but can arise from sources such as horizontal gradients and inaccurate calibration."
                              )
        
        # Relative VER 
        var = _create_variable(ncfile, '%s_Relative_VER'%prefix, L21_dict['ver'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Relative volume emission rate profile', 
                              display_type='altitude_profile', field_name='Relative VER', fill_value=None, label_axis='VER', bin_location=0.5,
                              units='ph/cm^3/s', valid_min=0, valid_max=1e30, var_type='data', chunk_sizes=[nt,ny],
                              notes= "The volume emission rate (VER) obtained by scaling the fringe amplitude by a calibration factor. "
                               "Pre-flight calibrations and on-orbit comparisons with ground-based instruments are used to determine the "
                               "best possible calibration. The fringe amplitude has a dependence on temperature, which is corrected using "
                               "the MSIS model. Because the on-orbit calibration is uncertain, and because the MSIS "
                               "temperature correction is not perfect, caution should be exercised when absolute "
                               "calibration is required, or when comparisons are being made between samples at different temperatures. "
                               "Please contact the MIGHTI team before performing any studies that require absolute calibration. "
                               "The statistical (1-sigma) error for this variable is provided in "
                               "the variable ICON_..._Relative_VER_Error, though it is expected that systematic calibration errors dominate "
                               "the total error. See the Fringe_Amplitude variable for a discussion of the inversion."
                              )
        
        # Relative VER Error 
        var = _create_variable(ncfile, '%s_Relative_VER_Error'%prefix, L21_dict['ver_error'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Relative volume emission rate error profile', 
                              display_type='altitude_profile', field_name='Relative VER', fill_value=None, label_axis='VER', bin_location=0.5,
                              units='ph/cm^3/s', valid_min=0, valid_max=1e30, var_type='data', chunk_sizes=[nt,ny],
                              notes= "The statistical (1-sigma) error in the relative VER estimate. This error arises mostly from shot "
                              "noise. Importantly, it is expected that systematic errors (e.g., calibration errors) dominate the total "
                              "error, but they are not included in this variable."
                              )
        
        # Quality code
        var = _create_variable(ncfile, '%s_VER_Quality'%prefix, L21_dict['ver_quality'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='A quantification of the VER quality, from 0 (Bad) to 1 (Good)', 
                              display_type='altitude_profile', field_name='VER Quality', fill_value=None, label_axis='Quality', bin_location=0.5,
                              units='', valid_min=0.0, valid_max=1.0, var_type='data', chunk_sizes=[nt,ny],
                              notes=["A quantification of the overall quality of the VER data. While the intent is that the variable VER_Error "
                                     "accurately characterizes the statistical "
                                     "error in the wind data, it is possible that systematic errors are present, or that the statistical error "
                                     "estimation is not accurate. If it is suspected that this is the case, the quality will be less than 1.0. If "
                                     "the data are definitely unusable, the the quality will be 0.0 and the sample will be masked. Users should "
                                     "exercise caution when the quality is less than 1.0.",
                                     "This parameter can currently take 3 values: 0.0 (Bad), 0.5 (Caution), 1.0 (Good)"
                                     ])

        ######### Data Location and Direction Variables #########
        
        # Altitude
        var = _create_variable(ncfile, '%s_Altitude'%prefix, L21_dict['alt'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='WGS84 altitude of each wind sample', 
                              display_type='altitude_profile', field_name='Altitude', fill_value=None, label_axis='Altitude', bin_location=0.5,
                              units='km', valid_min=50., valid_max=1000., var_type='support_data', chunk_sizes=[nt,ny],
                              notes="The altitudes of each point in the wind profile, evaluated using the WGS84 ellipsoid. If the variable "
                                    "Integration_Order=0 (which is the default value), then these altitudes are one half sample above "
                                    "the tangent altitudes of each pixel's line of sight (consistent with the assumption implicit "
                                    "in the inversion that the wind and emission rate are constant within the layer between tangent "
                                    "altitudes). If Integration_Order=1, this variable contains the tangent altitudes."
                              )

        # Latitude
        var = _create_variable(ncfile, '%s_Latitude'%prefix, L21_dict['lat'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='WGS84 latitude of each wind sample', 
                              display_type='altitude_profile', field_name='Latitude', fill_value=None, label_axis='Latitude', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='support_data', chunk_sizes=[nt,ny],
                              notes="The latitudes of each point in the wind profile, evaluated using the WGS84 ellipsoid. The latitude only "
                              "varies by several degrees from the bottom of the profile to the top. It should be noted that while a single latitude "
                              "value (the tangent latitude) is given for each point, the observation is inherently a horizontal average "
                              "over many hundreds of kilometers."
                              )

        # Longitude
        var = _create_variable(ncfile, '%s_Longitude'%prefix, L21_dict['lon'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='WGS84 longitude of each wind sample', 
                              display_type='altitude_profile', field_name='Longitude', fill_value=None, label_axis='Longitude', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[nt,ny],
                              notes="The longitudes (0-360) of each point in the wind profile, evaluated using the WGS84 ellipsoid. The longitude only "
                              "varies by several degrees from the bottom of the profile to the top. It should be noted that while a single longitude "
                              "value (the tangent longitude) is given for each point, the observation is inherently a horizontal average "
                              "over many hundreds of kilometers."
                              )
        
        # Magnetic Latitude
        var = _create_variable(ncfile, '%s_Magnetic_Latitude'%prefix, L21_dict['mag_lat'], 
                              dimensions=('Epoch', 'Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Magnetic quasi-dipole latitude of each wind sample', 
                              display_type='image', field_name='Mag Lat', fill_value=None, label_axis='Mag Lat', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='support_data', chunk_sizes=[nt,ny],
                              notes="A two-dimensional array defining the magnetic quasi-dipole latitude of the two-dimensional data grid. "
                              "The latitude varies only slightly (a few deg) with altitude, but this variation is included. "
                              "It should be noted that while a single latitude value is given for each point, the observation is "
                              "inherently a horizontal average over many hundreds of kilometers. "
                              "Quasi-dipole latitude and longitude are calculated using the fast implementation developed by "
                              "Emmert et al. (2010, doi:10.1029/2010JA015326) and the Python wrapper apexpy "
                              "(doi.org/10.5281/zenodo.1214207). "
                              )
        
        # Magnetic Longitude
        var = _create_variable(ncfile, '%s_Magnetic_Longitude'%prefix, L21_dict['mag_lon'], 
                              dimensions=('Epoch', 'Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Magnetic quasi-dipole longitude of each wind sample', 
                              display_type='image', field_name='Mag Lon', fill_value=None, label_axis='Mag Lon', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[nt,ny],
                              notes="A two-dimensional array defining the magnetic quasi-dipole longitude of the two-dimensional data grid. "
                              "The longitude varies only slightly (a few deg) with altitude, but this variation is included. "
                              "It should be noted that while a single longitude value is given for each point, the observation is "
                              "inherently a horizontal average over many hundreds of kilometers. "
                              "Quasi-dipole latitude and longitude are calculated using the fast implementation developed by "
                              "Emmert et al. (2010, doi:10.1029/2010JA015326) and the Python wrapper apexpy "
                              "(doi.org/10.5281/zenodo.1214207). Quasi-dipole longitude is defined such that zero "
                              "occurs where the geodetic longitude is near 285 deg east (depending on latitude)."
                              )
        

        # Azimuth angle of line of sight
        var = _create_variable(ncfile, '%s_Line_of_Sight_Azimuth'%prefix, L21_dict['az'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Azimuth angle of the line of sight at the tangent point. Deg East of North.', 
                              display_type='altitude_profile', field_name='Line-of-sight Azimuth', fill_value=None, label_axis='Azimuth', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[nt,ny],
                              notes="Consider the vector pointing from the spacecraft to the tangent point (i.e., the line of sight). At the tangent "
                              "point, this vector is parallel to the ground. This variable contains the azimuth angle of this vector, evaluated at "
                              "the tangent point. It follows the typical geophysical convention of degrees East of North (North=0, East=90, South=180, "
                              "West=270). It can vary by a few degrees from the top of the profile to the bottom, so one value is reported per altitude. "
                              "MIGHTI-A and MIGHTI-B will have values approximately 90 degrees apart."
                              )

        
        # The actual value used to correct the low signal effect
        var = _create_variable(ncfile, '%s_Low_Signal_Effect_Correction'%prefix, L21_dict['low_signal_corr'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Correction for low-signal effect', 
                              display_type='altitude_profile', field_name='Low-signal correction', fill_value=None, label_axis='Wind Corr', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='support_data', chunk_sizes=[nt,ny],
                              notes=["This is the correction used for the \"low signal effect\" in the lower-level processing. This correction has "
                                     "already been applied to the data, but is included here for reference. It was taken directly from the Level 1 file but "
                                     "has been converted from rad to m/s. For the red channel, it has also been binned to match the data. It has not been modified "
                                     "by the "
                                     "inversion. The uncertainty of this correction is captured by the Line_of_Sight_Wind_Precision_Low_Signal_Effect variable. See "
                                     "the notes for that variable for more information."]
                              )
                              
        

        ######### Other Metadata Variables #########
        
        # Solar zenith angle
        var = _create_variable(ncfile, '%s_Solar_Zenith_Angle'%prefix, L21_dict['sza'],
                              dimensions=('Epoch', 'Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Solar zenith angle of each wind sample', 
                              display_type='image', field_name='SZA', fill_value=None, label_axis='SZA', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=180., var_type='support_data', chunk_sizes=[nt,ny],
                              notes="Angle between the vectors towards the sun and towards zenith, at the location of each wind sample."
                              )
                               
        # Solar local time                      
        var = _create_variable(ncfile, '%s_Local_Solar_Time'%prefix, L21_dict['slt'], 
                              dimensions=('Epoch', 'Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Local solar time of each wind sample', 
                              display_type='image', field_name='LST', fill_value=None, label_axis='LST', bin_location=0.5,
                              units='hour', valid_min=0., valid_max=24., var_type='support_data', chunk_sizes=[nt,ny],
                              notes="Local solar time at the location and time of each wind sample, calculated using the equation of time."
                              )

        # Exposure Time
        var = _create_variable(ncfile, '%s_Exposure_Time'%prefix, L21_dict['exp_time'],
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='f8', format_fortran='F', desc='The exposure time for each profile', 
                              display_type='time_series', field_name='Exposure Time', fill_value=None, label_axis='Exp Time', bin_location=0.5,
                              units='s', valid_min=0., valid_max=120., var_type='metadata', chunk_sizes=[nt],
                              notes="The exposure time (i.e., integration time) for each sample. Nominally this is 30 seconds during the day "
                                    "and 60 seconds at night.")
                               

        # Chi^2 variability of phase in each row
        var = _create_variable(ncfile, '%s_Chi2'%prefix, L21_dict['chi2'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Variance of the phase in each (unwrapped) row: (std of phase)^2', 
                              display_type='altitude_profile', field_name='Variance of each phase row', fill_value=None, label_axis='Chi^2', bin_location=0.5,
                              units='rad^2', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=[nt, ny],
                              notes="In consolidating each row of the unwrapped interferogram into a single phase value, the variance of the "
                              "phase is saved in this variable. Ideally this should provide no new information beyond what is provided "
                              "by the wind uncertainty, but it is a useful diagnostic."
                              )
        
        # ICON velocity vector
        var = _create_variable(ncfile, '%s_Observatory_Velocity_Vector'%prefix, L21_dict['icon_velocity_ecef_vector'], 
                              dimensions=('Epoch','Vector'), depend_0 = 'Epoch',
                              format_nc='f8', format_fortran='F', desc="ICON S/C velocity vector in Earth-centered, Earth-fixed coordinates", 
                              display_type='time_series', field_name='ICON Velocity Vector', fill_value=None, label_axis='S/C Vel', bin_location=0.5,
                              units='m/s', valid_min=-100e6, valid_max=100e6, var_type='metadata', chunk_sizes=[nt,3],
                              notes="At each time, this is a length-3 vector [vx,vy,vz] of the ICON spacecraft's velocity in Earth-centered Earth-fixed (ECEF) "
                              "coordinates at the "
                              "midpoint time of the observation. The effect of spacecraft velocity has already been removed from the "
                              "ICON_..._Line_of_Sight_Wind variable.")

        # ICON latitude
        var = _create_variable(ncfile, '%s_Observatory_Latitude'%prefix, L21_dict['icon_lat'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='f8', format_fortran='F', desc='The WGS84 latitude of the ICON S/C', 
                              display_type='time_series', field_name='Spacecraft Latitude', fill_value=None, label_axis='S/C Lat', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='metadata', chunk_sizes=[nt],
                              notes='The latitude of the ICON spacecraft at the midpoint time of the observation, using the WGS84 ellipsoid.')

        # ICON longitude
        var = _create_variable(ncfile, '%s_Observatory_Longitude'%prefix, L21_dict['icon_lon'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='f8', format_fortran='F', desc='The WGS84 longitude of the ICON S/C', 
                              display_type='time_series', field_name='Spacecraft Longitude', fill_value=None, label_axis='S/C Lon', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='metadata', chunk_sizes=[nt],
                              notes='The longitude (0-360) of the ICON spacecraft at the midpoint time of the observation, using the WGS84 ellipsoid.')

        # ICON altitude
        var = _create_variable(ncfile, '%s_Observatory_Altitude'%prefix, L21_dict['icon_alt'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='f8', format_fortran='F', desc='The WGS84 altitude of the ICON S/C', 
                              display_type='time_series', field_name='Spacecraft Altitude', fill_value=None, label_axis='S/C Alt', bin_location=0.5,
                              units='km', valid_min=100., valid_max=2000., var_type='metadata', chunk_sizes=[nt],
                              notes='The altitude of the ICON spacecraft at the midpoint time of the observation, using the WGS84 ellipsoid.')

#         # Along-track resolution
#         val = L21_dict['resolution_along_track']
#         var = _create_variable(ncfile, '%s_Resolution_Along_Track'%prefix, val, 
#                               dimensions=(),
#                               format_nc='f8', format_fortran='F', desc='The horizontal resolution in the spacecraft velocity direction', 
#                               display_type='scalar', field_name='Along-Track Resolution', fill_value=None, label_axis='Hor Res AT', bin_location=0.5,
#                               units='km', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
#                               notes=["NOT YET IMPLEMENTED",]
#                                      )

#         # Cross-track resolution
#         val = L21_dict['resolution_cross_track']
#         var = _create_variable(ncfile, '%s_Resolution_Cross_Track'%prefix, val, 
#                               dimensions=(),
#                               format_nc='f8', format_fortran='F', desc='The horizontal resolution perpendicular to the spacecraft velocity direction', 
#                               display_type='scalar', field_name='Cross-Track Resolution', fill_value=None, label_axis='Hor Res CT', bin_location=0.5,
#                               units='km', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
#                               notes=["NOT YET IMPLEMENTED",]
#                                      )

#         # Altitude resolution
#         var = _create_variable(ncfile, '%s_Resolution_Altitude'%prefix, L21_dict['resolution_alt']
#                               dimensions=(),
#                               format_nc='f8', format_fortran='F', desc='The vertical resolution', 
#                               display_type='scalar', field_name='Vertical Resolution', fill_value=None, label_axis='Vert Res', bin_location=0.5,
#                               units='km', valid_min=0., valid_max=1e30, var_type='metadata', chunk_sizes=1,
#                               notes=["NOT YET IMPLEMENTED",]
#                                      )

        # MIGHTI ECEF vectors
        var = _create_variable(ncfile, '%s_Line_of_Sight_Vector'%prefix, L21_dict['mighti_ecef_vectors'], 
                              dimensions=('Epoch','Altitude','Vector'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='The look direction of each MIGHTI line of sight, as a vector in ECEF', 
                              display_type='altitude_profile', field_name='Line-of-sight Vector', fill_value=None, label_axis='LoS Vec', bin_location=0.5,
                              units='', valid_min=-1., valid_max=1., var_type='metadata', chunk_sizes=[nt,ny,3],
                              notes="The vector from the spacecraft to the tangent point (i.e., along MIGHTI's line of sight), as a unit "
                              "vector in Earth-centered Earth-fixed (ECEF) coordinates. A vector is provided for each tangent point for each time. If this "
                              "vector is transformed to an azimuth and zenith angle at the tangent point, the zenith angle will be 90 deg, and "
                              "the azimuth angle will be the same as the ICON_..._Line_of_Sight_Azimuth variable."
                              )
        
        # Orbit number
        var = _create_variable(ncfile, '%s_Orbit_Number'%prefix, L21_dict['orbit_number'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i8', format_fortran='I', desc='Orbit Number', 
                              display_type='time_series', field_name='Orb_Num', fill_value=None, label_axis='Orb Num', bin_location=0.5,
                              units='', valid_min=0, valid_max=1000000, var_type='metadata', chunk_sizes=[nt],
                              notes="Integer ICON orbit number"
                              )
        
        # Orbit node ascending/descending
        var = _create_variable(ncfile, '%s_Orbit_Node'%prefix, L21_dict['orbit_node'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i1', format_fortran='I', desc='Orbit Ascending/Descending Node Flag', 
                              display_type='time_series', field_name='Orb_Node', fill_value=None, label_axis='Orb Node', bin_location=0.5,
                              units='', valid_min=np.int8(0), valid_max=np.int8(1), var_type='metadata', chunk_sizes=[nt],
                              notes=["Orbit Ascending/Descending Node Flag.",
                                     "0 = Latitude of ICON is increasing.",
                                     "1 = Latitude of ICON is decreasing.",]
                              )
        
        # Bin Size
        var = _create_variable(ncfile, '%s_Bin_Size'%prefix, L21_dict['bin_size'], 
                              dimensions=(),
                              format_nc='i1', format_fortran='I', desc='How many raw samples were binned vertically for each reported sample', 
                              display_type='scalar', field_name='Bin Size', fill_value=None, label_axis='Bin Size', bin_location=0.5,
                              units='', valid_min=np.int8(1), valid_max=np.int8(100), var_type='metadata', chunk_sizes=1,
                              notes="To improve statistics, adjacent rows of the interferogram can be averaged together before the inversion. "
                              "This improves precision at the cost of vertical resolution. If no binning is performed, this value will be 1, "
                              "corresponding to ~2.5 km sampling. A value of 2 corresponds to ~5 km sampling, etc."
                              )

        # Integration order
        var = _create_variable(ncfile, '%s_Integration_Order'%prefix, L21_dict['integration_order'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i1', format_fortran='I', desc='Order used to discretize the integral for inversion: 0=Riemann, 1=Trapezoidal', 
                              display_type='time_series', field_name='Order', fill_value=None, label_axis='Order', bin_location=0.5,
                              units='', valid_min=np.int8(0), valid_max=np.int8(1), var_type='metadata', chunk_sizes=[nt],
                              notes="In formulating the inversion, an assumption must be made regarding the choice of basis functions, "
                              "which can be thought of as an assumption regarding the behavior of the wind and fringe amplitude (airglow volume emission rate) "
                              "within each altitude "
                              "layer. The most basic assumption is that these quantities are constant within each altitude layer, which corresponds "
                              "to Integration_Order=0. However, if it is assumed that the variation within each layer is linear, "
                              "Integration_Order=1. This sacrifices precision to improve vertical resolution."
                              )

        # How the top layer was handled in the inversion
        var = _create_variable(ncfile, '%s_Top_Layer_Model'%prefix, L21_dict['top_layer'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc=str, format_fortran='A', desc='How the top altitudinal layer is handled in the inversion: "exp" or "thin"', 
                              display_type='time_series', field_name='Top Layer', fill_value=None, label_axis='Top Layer', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='metadata', chunk_sizes=[nt],
                              notes="In formulating the inversion, an assumption must be made about the shape of the emission rate profile "
                              "above the top measured altitude, since this shape is not measured. It can be assumed to go to zero "
                              "(Top_Layer_Model=\"thin\") or assumed to fall off exponentially with a scale height of 26 km, a value extracted "
                              "from running a variety of airglow models (Top_Layer_Model=\"exp\"). Usually this choice will not affect the "
                              "inversion significantly. In cases where it does, the quality variable will be decreased."
                              )    
        
        # Attitude information bits
        var = _create_variable(ncfile, '%s_Attitude_LVLH_Normal'%prefix, L21_dict['att_lvlh_normal'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i1', format_fortran='I', desc='Attitude status bit 0: LVLH Normal', 
                              display_type='time_series', field_name='LVLH Norm', fill_value=None, label_axis='LVLH Norm', bin_location=0.5,
                              units='', valid_min=np.int8(0), valid_max=np.int8(1), var_type='metadata', chunk_sizes=[nt],
                              notes="LVLH Normal pointing. This variable is taken from bit 0 of the Level 1 variable "
                              "ICON_L1_MIGHTI_X_SC_Attitude_Control_Register. 0=False, 1=True"
                              )

        var = _create_variable(ncfile, '%s_Attitude_LVLH_Reverse'%prefix, L21_dict['att_lvlh_reverse'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i1', format_fortran='I', desc='Attitude status bit 1: LVLH Reverse', 
                              display_type='time_series', field_name='LVLH Rev', fill_value=None, label_axis='LVLH Rev', bin_location=0.5,
                              units='', valid_min=np.int8(0), valid_max=np.int8(1), var_type='metadata', chunk_sizes=[nt],
                              notes="LVLH Reverse pointing. This variable is taken from bit 1 of the Level 1 variable "
                              "ICON_L1_MIGHTI_X_SC_Attitude_Control_Register. 0=False, 1=True"
                              )

        var = _create_variable(ncfile, '%s_Attitude_Limb_Pointing'%prefix, L21_dict['att_limb_pointing'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i1', format_fortran='I', desc='Attitude status bit 2: Earth Limb Pointing', 
                              display_type='time_series', field_name='Att Limb', fill_value=None, label_axis='Att Limb', bin_location=0.5,
                              units='', valid_min=np.int8(0), valid_max=np.int8(1), var_type='metadata', chunk_sizes=[nt],
                              notes="Earth limb pointing. This variable is taken from bit 2 of the Level 1 variable "
                              "ICON_L1_MIGHTI_X_SC_Attitude_Control_Register. 0=False, 1=True"
                              )

        var = _create_variable(ncfile, '%s_Attitude_Conjugate_Maneuver'%prefix, L21_dict['att_conjugate'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i1', format_fortran='I', desc='Attitude status bit 6: Conjugate Maneuver', 
                              display_type='scalar', field_name='Conj. Man.', fill_value=None, label_axis='Conj. Man.', bin_location=0.5,
                              units='', valid_min=np.int8(0), valid_max=np.int8(1), var_type='metadata', chunk_sizes=[nt],
                              notes="Conjugate Maneuver. This variable is taken from bit 6 of the Level 1 variable "
                              "ICON_L1_MIGHTI_X_SC_Attitude_Control_Register. 0=False, 1=True. If it is 1, then the S/C "
                              "is performing a conjugate maneuver during this exposure."
                              )
        
        var = _create_variable(ncfile, '%s_Attitude_Zero_Wind_Maneuver'%prefix, L21_dict['att_zero_wind'], 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc='i1', format_fortran='I', desc='Attitude status bit 10: Zero Wind Maneuver', 
                              display_type='scalar', field_name='ZW Man.', fill_value=None, label_axis='ZW Man.', bin_location=0.5,
                              units='', valid_min=np.int8(0), valid_max=np.int8(1), var_type='metadata', chunk_sizes=[nt],
                              notes="Zero Wind Maneuver. This variable is taken from bit 10 of the Level 1 variable "
                              "ICON_L1_MIGHTI_X_SC_Attitude_Control_Register. 0=False, 1=True. If it is 1, then the S/C "
                              "is performing a zero wind maneuver during this exposure."
                              )
        
        
        # Quality flags      
        var = _create_variable(ncfile, '%s_Quality_Flags'%prefix, L21_dict['quality_flags'], 
                              dimensions=('Epoch','Altitude', 'N_Flags'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='i8', format_fortran='I', desc='Quality flags', 
                              display_type='image', field_name='Quality Flags', fill_value=-1, label_axis='Qual Flags', bin_location=0.5,
                              units='', valid_min=0, valid_max=1, var_type='metadata', chunk_sizes=[nt,ny,12],
                              notes=["This variable provides information on why the Wind_Quality and VER_Quality variables are reduced "
                              "from 1.0. Many quality flags can "
                              "exist for each grid point, each either 0 or 1. More than one flag can be raised per point. This variable "
                              "is a two-dimensional array with dimensions of altitude and number of flags.",
                                       "0 : (From L1) SNR too low to reliably perform L1 processing",
                                       "1 : (From L1) Proximity to South Atlantic Anomaly",
                                       "2 : (From L1) Bad calibration",
                                       "3 : (From L1) Calibration lamps are on",
                                       "4 : (From L1) Solar/lunar contamination",
                                       "5 : Not enough valid points in profile",
                                       "6 : SNR is very low after inversion",
                                       "7 : Significant airglow above 300 km",
                                       "8 : Line of sight crosses the terminator",
                                       "9 : Thermal drift correction is uncertain",
                                       "10: S/C pointing is not stable",
                                       "11: SNR is low after inversion, but maybe still usable",
                                    ])
        
        # Alternative VER product derived from DC value of interferogram
        var = _create_variable(ncfile, '%s_Relative_VER_DC'%prefix, L21_dict['ver_dc'], 
                              dimensions=('Epoch','Altitude'), depend_0 = 'Epoch', depend_1 = 'Altitude',
                              format_nc='f8', format_fortran='F', desc='Relative volume emission rate profile derived from DC value', 
                              display_type='altitude_profile', field_name='Relative VER', fill_value=None, label_axis='VER', bin_location=0.5,
                              units='ph/cm^3/s', valid_min=0, valid_max=1e30, var_type='metadata', chunk_sizes=[nt,ny],
                              notes= "The MIGHTI team recommends that users utilize the Relative_VER variable instead of this variable. This "
                               "is the same as Relative_VER, except it is derived from the DC value of the interferogram rather than the fringe "
                               "amplitude. The DC value is susceptible to contamination by stray light and background emission, but is "
                               "not sensitive to atmospheric temperature like the fringe amplitude."
                              )

        # Zero Wind Phase
        var = _create_variable(ncfile, '%s_Zero_Wind_Phase'%prefix, L21_dict['zero_wind_phase'], 
                              dimensions=('Epoch','Row'), depend_0 = 'Epoch', depend_1 = 'Row',
                              format_nc='f8', format_fortran='F', desc='The phase subtracted from Level 1 data', 
                              display_type='altitude_profile', field_name='Zero wind phase', fill_value=None, label_axis='Zero phase', bin_location=0.5,
                              units='rad', valid_min=-2*np.pi, valid_max=2*np.pi, var_type='metadata', chunk_sizes=[nt,ny],
                              notes="See notes for _Line_of_Sight_Wind above. "
                             "This variable is reported as a function of Row (Level 1 coordinates) instead of Altitude (Level 2 coordinates). "
                              )
        
        ncfile.close()

    except: # make sure the file is closed
        try:
            ncfile.close()
        except RuntimeError:
            pass
        raise
    
    return L21_full_fn
    
    
    
    
    
def combine_level21(L21_dicts):
    '''
    Combine many L2.1 dictionaries (from the output of the L1-to-L2.1 processing) into a single dictionary with the 
    same keys, but with the variables that change with time having an extra dimension (dimension 0) for time. This
    new dictionary is the proper input for the save_nc_level21 function.
    
    INPUTS:
    
      * L21_dicts  --TYPE:list of dict.   Each dictionary should have form described by level1_dict_to_level21_dict(...)
      
    OUTPUTS:
    
      * L21_dict   --TYPE:dict.           The variables in this dictionary match those described in the documentation for 
                                          level1_dict_to_level21_dict(...) but have an extra dimension (on axis 0) for time.
                                          Variables which do not change in time (e.g., emission color) do not have this extra 
                                          dimension. The variable "I" is also omitted, for speed, because it is not saved
                                          in the L2.1 file.
    '''
    N = len(L21_dicts)

    # The parameters to be saved fall into two categories: those that change with time and those that don't.
    # Define keys which should *not* change in time. All the rest will have an extra dimension 0.
    common_keys = ['acknowledgement',
                   'sensor',
                   'emission_color',
                   'bin_size', # if this varies within a day, we can't save all profiles in one file
                  ]
    changing_keys = L21_dicts[0].keys()
    for k in common_keys:
        changing_keys.remove(k)
    if 'I' in changing_keys:
        changing_keys.remove('I') # because it is large and not saved in the 2.1 file anyway.

    # Check if all of the common keys are in fact common
    for key in common_keys:
        for n in range(N):
            assert L21_dicts[n][key] == L21_dicts[0][key], 'The variable "%s" is changing with time (%s vs %s)' % (key, L21_dicts[n][key], L21_dicts[0][key])

    # Create one large dictionary
    L21_dict = {}
    for key in common_keys:
        L21_dict[key] = L21_dicts[0][key]
    for key in changing_keys:
        z = [L21_dicts[n][key] for n in range(N)]
        v = np.stack(z, axis=0) # Make time the first dimension
        L21_dict[key] = v    
    return L21_dict


    

def get_GPI(dn, year_day, f107, f107a, ap, ap3):
    '''
    Operational Code to pull the geophysical indices from the standard GPI
    ancillary file.
    INPUTS:
        dn - datetime to be used
        year_day - vector containing the year/day of the GPIs (yyyddd)
        f107 - vector containing the f107 index
        f107a - vector containing the proxy for f107a average
        ap - vector containing the daily ap index
        ap3 - array containing the 3-hour ap index
    OUTPUT:
        returns the f107, f107a, f107p, ap, and ap3 for the requested time (dn)
    NOTES:
        Can be used in conjunction with get_msisGPI to generate the msis ap vector.
    HISTORY:
        02-Jun-2017: Written by Jonathan Makela (jmakela@illinois.edu)
        14-Aug-2019: Copied from FUV_L2.py to MIGHTI_L2.py
    '''

    yd = dn.year*1000+dn.timetuple().tm_yday

    # Check bounds
    if (yd < year_day.min()) or (yd > year_day.max()):
        return np.nan, np.nan, np.nan, np.nan

    # Get the index into the year_day dimension
    i = np.argwhere(year_day == yd).flatten()[0]

    # Get the index into the 3 hour dimension
    j = dn.hour/3

    return f107[i], f107a[i], ap[i], ap3[j,i]



    
def get_msisGPI(dn, year_day, f107, f107a, ap, ap3):
    '''
    Operational Code to pull the geophysical indices from the standard GPI
    ancillary file and generate the ap vector required by MSIS.
    INPUTS:
        dn - datetime to be used
        year_day - vector containing the year/day of the GPIs (yyyddd)
        f107 - vector containing the f107 index
        f107a - vector containing the proxy for f107a average
        ap - vector containing the daily ap index
        ap3 - array containing the 3-hour ap index
    OUTPUT:
        returns the f107, f107a, msis ap for the requested time (dn)
    NOTES:
        As defined in the MSIS code the msis ap is a vector containing:
             (1) DAILY AP
             (2) 3 HR AP INDEX FOR CURRENT TIME
             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS
                 PRIOR   TO CURRENT TIME
             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS
                 PRIOR  TO CURRENT TIME
    HISTORY:
        02-Jun-2017: Written by Jonathan Makela (jmakela@illinois.edu)
        01-Nov-2017: Added the prior-day f107 value needed by MSIS (jmakela@illinois.edu)
        14-Aug-2019: Copied with superficial edits from FUV_L2.py to MIGHTI_L2.py
    '''

    # Where to store the
    my_apmsis = np.zeros(7)
    
    # Get f107, f107a, daily ap, and 3-hour ap
    my_f107, my_f107a, my_apmsis[0], my_apmsis[1] = get_GPI(dn+timedelta(hours=0),year_day,f107,f107a,ap,ap3)
    
    # Get the previous day f107
    my_f107p, _, _, _ = get_GPI(dn+timedelta(hours=0)+timedelta(days=-1),year_day,f107,f107a,ap,ap3)
    
    # Get 3-hour ap for -3, -6, and -9 hours
    _, _, _, my_apmsis[2] = get_GPI(dn+timedelta(hours = -3),year_day,f107,f107a,ap,ap3)
    _, _, _, my_apmsis[3] = get_GPI(dn+timedelta(hours = -6),year_day,f107,f107a,ap,ap3)
    _, _, _, my_apmsis[4] = get_GPI(dn+timedelta(hours = -9),year_day,f107,f107a,ap,ap3)
    
    # Get average 3-hour ap for -12 to -33 hours
    temp_ap = 0
    for delta in range(-12,-33-1,-3):
        _, _, _, temp = get_GPI(dn+timedelta(hours = delta),year_day,f107,f107a,ap,ap3)
        temp_ap += temp
    my_apmsis[5] = temp_ap/8.
    
    # Get average 3-hour ap for -36 to -57 hours
    temp_ap = 0
    for delta in range(-36,-57-1,-3):
        _, _, _, temp = get_GPI(dn+timedelta(hours = delta),year_day,f107,f107a,ap,ap3)
        temp_ap += temp
    my_apmsis[6] = temp_ap/8.
    
    return my_f107, my_f107a, my_f107p, my_apmsis
    
    
    
    
def level1_to_level15_without_info_file(L1_fns, emission_color, L15_path, vers, data_revision):
    '''
    High-level function to apply the Level-1-to-Level-1.5 algorithm to a series of Level 1 files. The input files
    must all come from the same sensor (MIGHTI-A or B, not both) and from the same date. It is expected that this 
    function will be run twice: once for red and once for green. This version
    of the function requires the user to input the arguments instead of specifying them with an 
    Information.TXT file, as will be done in the Science Data Center. The version of the file is
    copied from the first L1 file in the set, in contrast with the L2 code which specifies its own.
    
    INPUTS:
    
      *  L1_fns              -- TYPE:str.  The full path to the Level 1 files to be processed. All files must be from
                                           the same sensor (A or B) on the same date.
      *  emission_color      -- TYPE:str.  'red' or 'green'
      *  L15_path            -- TYPE:str.  The directory the Level 1.5 file will be saved in, including trailing "/"
                                           (e.g., '/home/user/')
      *  vers                -- TYPE:int.  The version number for the filename (this will be checked to see it matches L1)
      *  data_revision       -- TYPE:int.  The revision number for the filename (this need not match L1s)
   
    OUTPUTS:      
    
      *  L15_fn              -- TYPE:str.   The full path to the saved L1.5 file.
    
    '''

    # Parse inputs
    if emission_color not in ['red','green']:
        raise ValueError('Argument emission_color=\'%s\' not recognized. Use \'red\' or \'green\'.' % emission_color)
    if len(L1_fns) < 2:
        raise ValueError('Need more than 1 L1 file. len(L1_fns) = %i' % (len(L1_fns)))

    vers_L1s = L1_fns[0][-10:-7] # e.g., "v02"
    vers_commanded = 'v%02i' % vers
    assert vers_L1s == vers_commanded, "Version mismatch: commanded (%s) and input L1s (%s)" % (vers_commanded, vers_L1s)

    L1_fns.sort() # To make sure it's in the right time order

    L15s = []
    for L1_fn in L1_fns:
        L1_dict = level1_to_dict(L1_fn, emission_color)
        ds = level1_dict_to_level15_dataset(L1_dict, emission_color)
        L15s.append(ds)

    L15 = xr.concat(L15s, dim='time')
    
    # Adjust time base to match ICON convention (this doesn't affect anything in xarray but does affect the NC file)
    L15.time.encoding['units'] = 'milliseconds since 1970-01-01'

    assert (L15['sensor'] == L15['sensor'][0]).all(), "All files need to come from the same sensor (A or B)"
    dates = L15['time'].dt.floor(freq='1d')
#     assert (dates == dates[0]).all(), "All files need to come from the same date" # Commented 16 Aug 2022 because sometimes 1-2 files from prev date
    assert (L15['color'] == L15['color'][0]).all(), "Samples from different colors are mixed (this should never happen)"
    
    # Add parent files
    # Don't list every file, but rather, use wildcards
    L1_fns_short = [s.split('/')[-1] for s in L1_fns]
    assert np.std([len(fn) for fn in L1_fns_short]) == 0., "L1 fns have different lengths. Code needs to be made smarter."
    parent = ''
    for i in range(len(L1_fns_short[0])): # Loop over positions in string
        s = [fn[i] for fn in L1_fns_short]
        su = pd.unique(s)
        if len(su) == 1:
            parent += su[0]
        else:
            parent += '?'
    L15 = L15.assign_attrs(Parents=parent)
    
    # Grab date string from middle file 
    N = len(dates)
    tstr = pd.to_datetime(dates[N//2].item()).strftime('%Y-%m-%d') # Changed 16 Aug 2022 to use middle time
    sens = L15['sensor'][0].item()
    L15_fn = '%s/ICON_L1-5_MIGHTI-%s_Consolidated-Phase-%s_%s_v%02ir%03i.NC' % (L15_path, sens, emission_color.capitalize(), tstr, vers, data_revision)
    
    L15.to_netcdf(L15_fn)
    return L15_fn
    

    
def zero_wind_phase_from_L1(dsz, L1_dict):
    '''
    Given an already-loaded zero-wind-phase file and a L1 dictionary, determine the row-dependent zero wind phase
    to be applied. This function combines the "zero_phase" and "striation" variables to produce one profile. It also
    implements the error estimates, and reports them separately as striation error and zero-wind phase error.
    
    INPUTS:
          *  dsz     -- TYPE:xarray.Dataset. From a loaded zero wind file, i.e.,:
                                   dsz = xr.open_dataset('/path/to/zero_wind_file.NC')
          *  L1_dict -- TYPE:dict. L1 dictionary in the typical format. See documentation for 
                                   level1_dict_to_level21_dict(...).
                               
    OUTPUT:
          *  z0       -- TYPE:array(ny), UNITS:rad. The zero wind phase for this exposure. This is an estimate
                                                    of the phase that corresponds to the rest wavelength of the 
                                                    emission. This includes striations and zero wind phase.
          *  z0_err_striation   -- TYPE:array(ny), UNITS:rad. An estimate of the striation uncertainty of z0. This is based on the
                                                              magnitude of measured striations.
          *  z0_err_zero_phase  -- TYPE:array(ny), UNITS:rad. An estimate of the zero-wind phase uncertainty of z0. This is based on 
                                                              simulations using different techniques of estimating the zero wind phase.
          *  Ncorr_striation    -- TYPE:int.       UNITS:pixels. The row lag at which the residual striation correlation drops to 0.
                                                    
    '''
    
    # Because of the off-by-1 errors, and the start/stop time issue, the relevant time here will
    # be extracted from the filename instead of the timestamp
    tstr = L1_dict['L1_fn'].split('/')[-1].split('_')[-3]
    
    # Raise an error if the time requested is not available.
    try:
        dszt = dsz.sel(time = tstr)
    except KeyError:
        raise KeyError('"%s" is not even a timestamp in the zero wind phase file' % tstr)
        
    
    
    # Load the rest of the dimensions
    sensor = L1_dict['sensor']
    color = L1_dict['emission_color']
    mode = 'Day'
    if L1_dict['exp_time'] > 45.:
        mode = 'Night'
    cal_lamp = L1_dict['quality_flags'][0,3] # 1 = on, 0 = off
    rows = np.arange(len(L1_dict['tang_alt_start'])) # This is to shorten the array to the red rows
    
    #### Zero wind 
    z = dszt.sel(sensor=sensor, color=color, mode=mode, cal_lamp=cal_lamp, row=rows)
    
    if z.zero_phase_with_notch.isnull().all():
        raise ValueError('Zero wind phase is NaN for all rows for this case: (time=%s, sensor=%s, color=%s, mode=%s, cal_lamp=%i)' % \
                         (tstr, sensor, color, mode, cal_lamp))
    z0 = z.zero_phase_with_notch.values
    
    #### Zero wind and striation uncertainty estimate
        
    # Error in the zero wind. This was computed from the results of a (zero-accurate) HWM-based simulation, as 
    # well as estimated errors from the choice of window size and of gap matching.
    # See "zero_wind_quicklook.ipynb" for the full computation.
    zero_err_all = {
    ('A','Day','green',0):  + (5.93000830e-08 * rows**4) + (-5.72614662e-06 * rows**3) + (2.12845532e-04 * rows**2) + (-2.98831057e-03 * rows**1) + (3.27268781e-02 * rows**0) , 
    ('B','Day','green',0):  + (2.57960085e-07 * rows**4) + (-2.10005975e-05 * rows**3) + (5.15652876e-04 * rows**2) + (-3.38634469e-03 * rows**1) + (2.93805092e-02 * rows**0) , 
    ('A','Night','green',0):  + (-7.21343700e-05 * rows**4) + (1.91364270e-03 * rows**3) + (-1.41907358e-02 * rows**2) + (3.39027228e-02 * rows**1) + (5.04453365e-03 * rows**0) , 
    ('B','Night','green',0):  + (3.78402405e-05 * rows**4) + (-4.18706725e-04 * rows**3) + (1.43559276e-03 * rows**2) + (-5.41962320e-04 * rows**1) + (1.56490337e-02 * rows**0) , 
    ('A','Day','green',1):  + (2.67738566e-08 * rows**4) + (-1.79203323e-06 * rows**3) + (4.44976330e-05 * rows**2) + (-3.31112254e-04 * rows**1) + (2.89893278e-02 * rows**0) , 
    ('B','Day','green',1):  + (2.95518174e-07 * rows**4) + (-2.44738229e-05 * rows**3) + (6.11868032e-04 * rows**2) + (-4.19663948e-03 * rows**1) + (3.07019680e-02 * rows**0) , 
    ('A','Night','green',1):  + (1.46847963e-04 * rows**4) + (-2.17862878e-03 * rows**3) + (1.07264115e-02 * rows**2) + (-2.18593127e-02 * rows**1) + (4.63102603e-02 * rows**0) , 
    ('B','Night','green',1):  + (2.56661735e-05 * rows**4) + (-1.67666172e-05 * rows**3) + (-2.20354449e-03 * rows**2) + (8.48224029e-03 * rows**1) + (1.89545016e-02 * rows**0) , 
    ('A','Day','red',0):  + (-9.76674545e-09 * rows**4) + (1.30500921e-06 * rows**3) + (-4.75048443e-05 * rows**2) + (6.62966021e-04 * rows**1) + (1.95195345e-02 * rows**0) , 
    ('B','Day','red',0):  + (2.97802180e-08 * rows**4) + (-3.21995706e-06 * rows**3) + (1.28014677e-04 * rows**2) + (-1.94112095e-03 * rows**1) + (3.16355624e-02 * rows**0) , 
    ('A','Night','red',0):  + (-5.98222589e-09 * rows**4) + (3.57720281e-07 * rows**3) + (2.04770593e-05 * rows**2) + (-1.52638212e-03 * rows**1) + (4.72305507e-02 * rows**0) , 
    ('B','Night','red',0):  + (-1.26352938e-07 * rows**4) + (1.85948309e-05 * rows**3) + (-9.80868502e-04 * rows**2) + (2.22889169e-02 * rows**1) + (-1.47603371e-01 * rows**0) , 
    ('A','Day','red',1):  + (-2.18675583e-09 * rows**4) + (3.69517192e-07 * rows**3) + (-1.60967941e-05 * rows**2) + (3.07153953e-04 * rows**1) + (1.66545383e-02 * rows**0) , 
    ('B','Day','red',1):  + (2.70159995e-08 * rows**4) + (-3.09243693e-06 * rows**3) + (1.24646005e-04 * rows**2) + (-1.87330570e-03 * rows**1) + (3.14155971e-02 * rows**0) , 
    ('A','Night','red',1):  + (-3.91129461e-08 * rows**4) + (4.58251674e-06 * rows**3) + (-1.57623342e-04 * rows**2) + (1.35875707e-03 * rows**1) + (3.50490539e-02 * rows**0) , 
    ('B','Night','red',1):  + (-2.00259817e-09 * rows**4) + (2.65047803e-07 * rows**3) + (-5.26592996e-06 * rows**2) + (-3.22778413e-06 * rows**1) + (2.65827182e-02 * rows**0) , 
    }
    z0_err_zero_phase = zero_err_all[(sensor, mode, color, cal_lamp)]
    # Manually degrade the accuracy estimate for Green Night rows 8-9 (2022 Sep 8, based on preliminary run of v05 compared to HWM)
    if (mode=='Night') and (color=='green'):
        z0_err_zero_phase[8:10] *= 2
    
    # If we are in a "bad zero phase" period (as defined in zero_wind.py) then degrade the accuracy.
    bad_zero_phase_periods = [ # (start, stop)
        (pd.to_datetime('2021-04-26'), pd.to_datetime('2021-08-14')), # Windows that touch the June 2021 outage/rLVLH
    ]
    bad = False
    t = pd.to_datetime(tstr)
    for tstart, tstop in bad_zero_phase_periods:
        if (t>=tstart) & (t<=tstop):
            bad=True
    if bad:
        z0_err_zero_phase *= 2
    
    # Estimate striation error from fluctuations in daily-averaged phase
    p = z.phase # Doesn't matter if we use phase or phase_with_notch here
    phpf = p - p.rolling(row=5, center=True, min_periods=1).median() # High-pass median filtered version
    q = phpf.to_series().rolling(3, min_periods=2, center=True).std().fillna(method='ffill').fillna(method='bfill') # rolling stddev with extrapolation
    q = q.to_xarray().where(~p.isnull()) # Only evaluate where the data actually exist
    z0_err_striation = q.values
    Ncorr_striation = 1 # Unlike in preliminary v05 runs, the striation-correction bandaid is no longer used, and thus they are uncorrelated with other rows.
    
    return z0, z0_err_striation, z0_err_zero_phase, Ncorr_striation
    
        
        
        
def level1_to_level21_without_info_file(L1_fns, emission_color, L21_path, data_revision=0,
                                        gpi_yearday=None, gpi_f107=None, gpi_f107a=None, gpi_ap=None, gpi_ap3=None,
                                        sigma=None, top_layer=None, 
                                        H = None, integration_order=None, account_for_local_projection=None, 
                                        bin_size=None, top_layer_thresh=None, terminator_thresh=None, zero_wind_ref=None,
                                        corr_notch_drift=None, chi2_thresh=None, chi2_thresh_caution=None, zero_wind_file=None):
    '''
    High-level function to apply the Level-1-to-Level-2.1 algorithm to a series of Level 1 files. The input files
    must all come from the same sensor (MIGHTI-A or B, not both) and from the same date. It is expected that this 
    function will be run twice: once for red and once for green. This version
    of the function requires the user to input the arguments instead of specifying them with an 
    Information.TXT file, as will be done in the Science Data Center.
    
    INPUTS:
    
      *  L1_fns              -- TYPE:str.  The full path to the Level 1 files to be processed. All files must be from
                                           the same sensor (A or B) on the same date.
      *  emission_color      -- TYPE:str.  'red' or 'green'
      *  L21_path            -- TYPE:str.  The directory the Level 2.1 file will be saved in, including trailing "/"
                                           (e.g., '/home/user/')
      
    OPTIONAL INPUT:
    
      *  data_revision       -- TYPE:int,  The minor version of the data [0-999]. The major version is set
                                           by this software's major version. (default 0)
      *  gpi_yearday         -- TYPE:array(N), YYYYDDD ints, taken from GPI file (None is use values from pyglow)
      *  gpi_f107            -- TYPE:array(N), UNITS:sfu. F10.7 value for each day (None is use values from pyglow)
      *  gpi_f107a           -- TYPE:array(N), UNITS:sfu. F10.7a value for each day (None is use values from pyglow)
      *  gpi_ap              -- TYPE:array(N), Daily Ap for each day (None is use values from pyglow)
      *  gpi_ap3             -- TYPE:array(8,N), 3-hourly Ap value for each day (None is use values from pyglow)

    MORE OPTIONAL INPUTS - If None, defaults from MIGHTI_L2.global_params will be used 
    
      *  sigma               -- TYPE:float, UNITS:m^-1. The wavenumber of the emission (1/wavelength)
      *  top_layer           -- TYPE:str,  'thin': assume VER goes to zero above top layer
                                           'exp':  assume VER falls off exponentially in altitude
      *  H                   -- TYPE:float, UNITS:km. The VER scale height to use if top_layer='exp'
      *  integration_order   -- TYPE:int,   0: Use Riemann-sum rule for discretizing line-of-sight integral
                                            1: Use trapezoidal rule for discretizing line-of-sight integral
      *  account_for_local_projection -- TYPE:bool. If False, a simple inversion is used.
                                            If True, the inversion accounts for the fact that the ray is not 
                                            perfectly tangent to each shell at each point along the ray.
      *  bin_size            -- TYPE:int,   The number of rows of the interferogram to bin together to 
                                            improve statistics at the cost of altitude resolution.
      *  top_layer_thresh.   -- TYPE:float. Fraction of airglow above the top observed altitude. Consider the total column 
                                            brightness (i.e, the integral of the VER profile). When a large fraction
                                            comes from above the top observed altitude (~300 km) the quality flag is raised.
      *  terminator_thresh   -- TYPE:float, UNITS:km. Consider two points along the line of sight, both X km from the tangent 
                                            point (where X is this parameter). One is nearer to the spacecraft than the 
                                            tangent point, and one is farther. If these two points are on opposite sides of the 
                                            terminator, raise a quality flag. Note that this quality flag will also be raised if
                                            any observations at higher tangent altitudes for the same observation are flagged, because
                                            the inversion mixes information from those observations.
      *  zero_wind_ref       -- DEPRECATED
      *  corr_notch_drift    -- DEPRECATED
      *  chi2_thresh         -- TYPE:float, UNITS:rad^2.  If the mean-square residual of phase within a row is larger than this,
                                                          the sample is discarded, and isn't used in the inversion.
      *  chi2_thresh_caution -- TYPE:float, UNITS:rad^2.  If the mean-square residual of phase within a row is larger than this,
                                                          the sample should be treated with caution but is used in the inversion. This
                                                          affects the quality control only.
      *  zero_wind_file      -- TYPE:str,                 Full path to a "Zero-Phase-Notch" file, which contains the variables
                                                          "zero_phase" and "zero_phase_with_notch" which will be used to apply the zero wind
                                                          phase correction.
                                                          
                                           
    OUTPUTS:      
    
      *  L21_fn              -- TYPE:str.   The full path to the saved L2.1 file.
      *  msg                 -- TYPE:list.  A list of strings, one for each failure during processing of a L1 file. If no 
                                            failures occurred, this is an empty list.

    '''
    
    if global_params['verbose']:
        print('\n%s: Starting Overall Run \n' % (timestamp()))
        import sys
        sys.stdout.flush()
        
    # v5.00+: check for deprecated inputs related to zero wind and notches
    if zero_wind_ref is not None:
        print('WARNING: The "zero_wind_ref" input is deprecated and is being ignored.')
    if corr_notch_drift is not None:
        print('WARNING: The "corr_notch_drift" input is deprecated and is being ignored.')

    assert len(L1_fns)>0, "No files specified."
    x = [gpi_yearday is None, gpi_f107 is None, gpi_f107a is None, gpi_ap is None, gpi_ap3 is None]
    assert (all(x)) or (not any(x)), "All GPI inputs must be specified, or all must be None."

    L1_fns.sort() # To make sure it's in the right time order
        
    # Parse inputs
    if emission_color not in ['red','green']:
        raise ValueError('Argument emission_color=\'%s\' not recognized. Use \'red\' or \'green\'.' % emission_color)
    # For other inputs, just pass them to the lower-level function. It will replace the Nones.
    
    # Load zero wind file if it exists
    if zero_wind_file is not None:
        dsz = xr.open_dataset(zero_wind_file)
        

    #################### Run L2.1 processing for all files ###########################
    L21_dicts = []
    failure_msg = []
    for L1_fn in L1_fns:
        
        if global_params['verbose']:
            print('%s: Running L2.1 %s %s' % (timestamp(), emission_color, L1_fn.split('/')[-1]))
            import sys
            sys.stdout.flush()
        
        try:
            # Read L1 file into a dictionary
            L1_dict = level1_to_dict(L1_fn, emission_color)

            # If GPI inputs are specified, calculate the relevant inputs for MSIS:
            f107, f107a, f107p, apmsis = None, None, None, None
            if gpi_yearday is not None:
                f107, f107a, f107p, apmsis = get_msisGPI(L1_dict['time_start'], gpi_yearday, gpi_f107, gpi_f107a, gpi_ap, gpi_ap3)

            # Zero wind phase (including striations)
            z, z_err_striation, z_err_zero_phase, Ncorr = None, None, None, 7 # Arbitrary defaults
            if zero_wind_file is not None:
                z, z_err_striation, z_err_zero_phase, Ncorr = zero_wind_phase_from_L1(dsz, L1_dict)
            
            # Perform L1 to L2.1 processing
            if global_params['verbose']:
                print('%s: \tCalling level1_dict_to_level21_dict %s %s' % (timestamp(), emission_color, L1_fn.split('/')[-1]))
                import sys
                sys.stdout.flush()
            L21_dict = level1_dict_to_level21_dict(L1_dict, sigma = sigma, top_layer = top_layer, H = H,
                                                   integration_order = integration_order, 
                                                   account_for_local_projection = account_for_local_projection, 
                                                   bin_size = bin_size, 
                                                   top_layer_thresh = top_layer_thresh, terminator_thresh = terminator_thresh,
                                                   chi2_thresh = chi2_thresh, chi2_thresh_caution = chi2_thresh_caution,
                                                   f107=f107, f107a=f107a, f107p=f107p, apmsis=apmsis,
                                                   zero_wind_phase = z, err_striation = z_err_striation, err_zero_phase = z_err_zero_phase,
                                                   Ncorr = Ncorr)
            L21_dicts.append(L21_dict)
            
        except Exception as e:
            failure_msg.append('Failed processing:\t%s\n%s\n' % (L1_fn, traceback.format_exc()))
        
    # If all files failed processing, crash.
    if len(L21_dicts)==0:
        raise Exception('All files failed. Tracebacks:\n%s'%('\n'.join(failure_msg)))
    
    
    ########## Concatenate all of the L2.1 inversions into one dictionary and save ##############
    if global_params['verbose']:
        print('%s: \tCombining L2.1 dictionaries' % (timestamp()))
        import sys
        sys.stdout.flush()
    L21_dict = combine_level21(L21_dicts)
    # Add the zero wind file if it exists
    L21_dict['zero_wind_file'] = ''
    if zero_wind_file is not None:
        L21_dict['zero_wind_file'] = zero_wind_file.split('/')[-1]
    
    
    ########## Compute quality flag corrections that require a time sequence
    # Thermal drift correction is uncertain too soon after maneuvers.
    # (1) Identify gaps in the data longer than 8 minutes as maneuvers
    # (2) Flag data for a specified amount of time after the maneuver
    dt = np.diff(L21_dict['time'])
    dt = np.array([x.total_seconds() for x in dt])
    ivec = np.where(dt > 8*60.)[0] + 1 # index of first time sample after gap
    for i in ivec:
        dti = np.array([x.total_seconds() for x in L21_dict['time'] - L21_dict['time'][i]]) # seconds since maneuver ended
        # Flag all samples that are too soon after maneuver ended
        j = (dti >=0) & (dti <= global_params[emission_color]['post_maneuver_window']*3600.)
        L21_dict['quality_flags'][j,:,9] = 1
        
    # If all winds are NaN, then crash so we can figure out what happened. There's no need to save an empty file.
    if np.isnan(L21_dict['los_wind']).all():
        raise Exception('Processing completed, but all winds are NaN')
        
    
    
    ########## Save L2.1 file
    if global_params['verbose']:
        print('%s: \tSaving L2.1 file' % (timestamp()))
        import sys
        sys.stdout.flush()
    L21_fn = save_nc_level21(L21_path, L21_dict, data_revision)
    if global_params['verbose']:
        print('%s: \tCompleted level1_to_level21_without_info_file()' % (timestamp()))
        import sys
        sys.stdout.flush()

    return L21_fn, failure_msg
    
    
    
    
    
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



    
    
def level1_to_level15(info_fn):
    '''
    Highest-level function to apply the Level-1-to-Level-1.5 algorithm to a series of Level 1 files, with 
    input arguments specified via an information file. For each file, the processing is run twice:
    once for 'red' and once for 'green'. The output L1.5 files will be saved to the Output folder in
    the directory specified in the information file.
    
    INPUTS:
    
      * info_fn  -- TYPE:str.  Full path to an ASCII file in the following format:
      
                                        [PARAMETERS]
                                        Version=03
                                        Revision=001
                                        Directory=/path/to/wherever/
                                        GPI=/path/to/gpi_file.NC   (This line is optional)
                                        <other parameters>

                                        [FILES]
                                        ICON_L1_MIGHTI-A_Science_2017-03-03_191803_v04r006.NC
                                        
    OUTPUTS:   
    
      *  ret     -- TYPE:str. '' if everything worked. If not, a human-readable error message for each file that failed
    
    '''    
    # Read information.txt file
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
    # After change to multiple profiles in one file, only a single rev number is supported.
    # In the future we could split this up into A and B rev numbers, but I'm not sure how that would
    # be specified.
    assert len(data_revision)==1, "Multiple revision numbers not supported for Level 1.5 processing"
    data_revision = data_revision[0]
    # (3) Parse the version number
    vers = int(info['Version'])
    
    # Loop and call the lower-level function which does all the real work.
    # Split files into A and B sets, and run red and green for each set.
    L15_fns = []
    failure_messages = []
    for sensor in ['A', 'B']:
        L1AB = [fn for fn in L1_full_fns if 'MIGHTI-%s'%sensor in fn]
        for emission_color in ['red','green']:
            try:
                L15_fn = level1_to_level15_without_info_file(L1AB, emission_color, direc + 'Output/', vers, data_revision)
                L15_fns.append(L15_fn)
            except Exception as e:
                failure_messages.append('Failed processing:\n\tsensor  = %s\n\tcolor   = %s\n%s\n'%(sensor, emission_color, traceback.format_exc()))
        
    if not failure_messages: # Everything worked
        return ''
    
    else:
        s = '\n' + '\n'.join(failure_messages)
        return s
    
    

    
    
def level1_to_level21(info_fn):
    '''
    Highest-level function to apply the Level-1-to-Level-2.1 algorithm to a series of Level 1 files, with 
    input arguments specified via an information file. For each file, the processing is run twice:
    once for 'red' and once for 'green'. The output L2.1 files will be saved to the Output folder in
    the directory specified in the information file. Summary plot(s) will also be created in the Output
    folder.
    
    INPUTS:
    
      * info_fn  -- TYPE:str.  Full path to an ASCII file in the following format:
      
                                        [PARAMETERS]
                                        Revision=001
                                        Directory=/path/to/wherever/
                                        Zero_Phase=/path/to/zero_phase_file.NC
                                        GPI=/path/to/gpi_file.NC   (This line is optional)
                                        <other parameters>

                                        [FILES]
                                        ICON_L1_MIGHTI-A_Science_2017-03-03_191803_v04r006.NC
                                        
    OUTPUTS:   
    
      *  ret     -- TYPE:str. '' if everything worked. If not, a human-readable error message for each file that failed
    
    '''    
    # Read information.txt file
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
    # After change to multiple profiles in one file, only a single rev number is supported.
    # In the future we could split this up into A and B rev numbers, but I'm not sure how that would
    # be specified.
    assert len(data_revision)==1, "Multiple revision numbers not supported for Level 2.1 processing"
    data_revision = data_revision[0]
    # (3) Get GPI filename
    gpi_fn = None
    if 'GPI' in info.keys():
        gpi_fn = './Input/'+info['GPI']
    
    # Read GPI file, if it exists
    year_day = None
    f107 = None
    f107a = None
    ap = None
    ap3 = None
    if gpi_fn:
        gpi = netCDF4.Dataset(gpi_fn,mode='r')
        # Read the geophysical indeces
        ap3 = gpi['ap3'][:]
        ap = gpi['ap'][:]
        year_day = gpi['year_day'][:]
        f107 = gpi['f107d'][:]
        # Make sure this GPI has the average f107 in it
        if 'f107a' in gpi.variables.keys():
            f107a = gpi['f107a'][:]
        else:
            print('Cannot find f107a in provided GPI file. Using daily f107 instead')
            f107a = gpi['f107d'][:]
        gpi.close()
        
    # Read zero wind phase file from Calibrations/ directory
    fnzs = glob.glob('Calibrations/*')
    assert len(fnzs)<2, "More than one file found in Calibrations/ directory. Not sure how to proceed"
    assert len(fnzs)>0, "No files found in Calibrations/ directory."
    zero_wind_file = fnzs[0]
    assert "Zero-Phase-Notch" in zero_wind_file, "Calibration file not recognized (%s). A 'Zero-Phase-Notch' file is required" % zero_wind_file
    
    # Loop and call the lower-level function which does all the real work.
    # Split files into A and B sets, and run red and green for each set.
    L21_fns = []
    failure_messages = []
    for sensor in ['A', 'B']:
        L1AB = [fn for fn in L1_full_fns if 'MIGHTI-%s'%sensor in fn]
        for emission_color in ['red','green']:
            try:
                L21_fn, msg = level1_to_level21_without_info_file(L1AB, emission_color, direc + 'Output/', 
                                              gpi_yearday = year_day, gpi_f107=f107, gpi_f107a=f107a, gpi_ap=ap, gpi_ap3=ap3,
                                              data_revision = data_revision, zero_wind_file = zero_wind_file)
                L21_fns.append(L21_fn)
                failure_messages.extend(msg)
            except Exception as e:
                failure_messages.append('Failed processing:\n\tsensor  = %s\n\tcolor   = %s\n%s\n'%(sensor, emission_color, traceback.format_exc()))
    
    # Summary plots: One plot for each file
    if global_params['verbose']:
        print('%s: \tCreating summary plots...' % (timestamp()))
        import sys
        sys.stdout.flush()
    for fn in L21_fns:
        try:
            plot_level21(fn, direc + 'Output/')
        except Exception as e:
            failure_messages.append('Failed creating summary plot for file:\n\t%s\n%s\n'%(fn, traceback.format_exc()))

    if global_params['verbose']:
        print('%s: \tDone.' % (timestamp()))
        import sys
        sys.stdout.flush()
        
    if not failure_messages: # Everything worked
        return ''
    
    else:
        s = '\n' + '\n'.join(failure_messages)
        return s
    
    
    
    

    
    
def level21_to_dict(L21_fn, skip_att=[], keep_att=[], tstartstop = None, skip_bad_alts=False):
    ''' 
    Load a Level 2.1 file, which contains wind profiles:
        * from a single sensor (A or B)
        * from a single channel (red or green)
        * all on the same date
    Return relevant variables in a dictionary that resembles the dictionary created by the L2.1 processing, but
    with an extra dimension for time. The entire file will be read, unless one of the optional inputs is specified.
    At most one of these can be specified.
    
    INPUTS:
    
      *  L21_fn   -- TYPE:str.  The path to the Level 2.1 file to be loaded.
      
    OPTIONAL INPUTS:
      
      *  skip_att -- TYPE: list of str.  A list of attitude status variables. If any evaluate
                                         to True for a given exposure, this exposure will be 
                                         skipped and not included in the output. Specify zero or one of
                                         the inputs: skip_att, keep_att, or tstartstop.
                                         Possible values: ['att_conjugate', 'att_lvlh_reverse', 'att_zero_wind']
                                         Default: []
      *  keep_att -- TYPE: list of str.  A list of attitude status variables. Only keep data points for
                                         which these evaluate to True.  Specify zero or one of
                                         the inputs: skip_att, keep_att, or tstartstop.
                                         Possible values: ['att_conjugate', 'att_lvlh_reverse', 'att_zero_wind']
                                         Default: []
      *  tstartstop -- TYPE:2-tuple of datetimes: (tstart, tstop). Only keep data points for which the midpoint
                                         of the exposure time is between tstart and tstop, inclusive.
                                         Specify zero or one of the inputs: skip_att, keep_att, or tstartstop.
                                         Default: None
      * skip_bad_alts -- TYPE: bool.     If True, skip samples where the altitudes are outliers, which can occur
                                         for cases where there are slews labeled as science data. Thresholds are
                                         defined in the code below.
                                
                                         
      
    OUTPUTS:
    
      *  L21_dict -- TYPE:dict. A dictionary containing the following variables. Most 
                                are provided as arrays of shape (ny,nt), where ny is the number
                                of altitude samples and nt is the number of time samples.
                                
                  * lat             -- TYPE:array(ny,nt), UNITS:deg. Sample latitudes.
                  * lon             -- TYPE:array(ny,nt), UNITS:deg. Sample longitudes.
                  * alt             -- TYPE:array(ny,nt), UNITS:km.  Sample altitudes.
                  * los_wind        -- TYPE:array(ny,nt), UNITS:m/s. Line-of-sight wind component towards MIGHTI.
                  * los_wind_precision_1_sample -- TYPE:array(ny,nt), UNITS:m/s. 1-sample precision. See Var_notes
                  * los_wind_precision_1_day    -- TYPE:array(ny,nt), UNITS:m/s. 1-day precision. See Var_notes
                  * los_wind_precision_low_signal--TYPE:array(ny,nt), UNITS:m/s. Low signal precision. See Var_notes
                  * los_wind_accuracy           -- TYPE:array(ny,nt), UNITS:m/s. Accuracy. See Var_notes
                  * los_wind_error              -- TYPE:array(ny,nt), UNITS:m/s. For backwards compatbility. See Var_Notes.
                  * local_az        -- TYPE:array(ny,nt), UNITS:deg. Azimuth angle of vector pointing from 
                                       MIGHTI towards the sample location, at the sample location (deg E of N).
                  * amp             -- TYPE:array(ny,nt), UNITS:arb. Fringe amplitude at sample locations.
                  * amp_error       -- TYPE:array(ny,nt), UNITS:arb. Error in amp variable.
                  * ver             -- TYPE:array(ny,nt), UNITS:ph/cm^3/s. Volume emission rate at sample locations.
                  * ver_error       -- TYPE:array(ny,nt), UNITS:ph/cm^3/s. Error in ver variable. 
                  * time            -- TYPE:array(nt).               Array of datetime objects, one per file.
                  * icon_lat        -- TYPE:array(nt),    UNITS:deg. Spacecraft latitude.
                  * icon_lon        -- TYPE:array(nt),    UNITS:deg. Spacecraft longitude.
                  * icon_alt        -- TYPE:array(nt),    UNITS:km.  Spacecraft altitude
                  * exp_time        -- TYPE:array(nt),    UNITS:sec. Exposure time of each sample.
                  * sensor          -- TYPE:str.                     'A' or 'B'
                  * source_files    -- TYPE:list of str,             A length-1 list containing the name of the input file.
                  * version         -- TYPE:int.                     Version number of input file.
                  * revision        -- TYPE:int.                     Revision number of input file.
                  * wind_quality    -- TYPE:array(ny,nt),            Value between 0 and 1 for each sample
                  * ver_quality     -- TYPE:array(ny,nt),            Value between 0 and 1 for each sample
                  * quality_flags   -- TYPE:array(ny,nt,ne).         Quality flags for each point. See 
                                                                     documentation for level1_dict_to_level21_dict(...)
                                                                     for meaning of each flag.
                  * emission_color  -- TYPE:str,                     'red' or 'green'.
                  * acknowledgement -- TYPE:str.                     The Acknowledgment attribute in the file
                  * zero_wind_ref   -- TYPE:str.                     The Zero_Wind_ref attribute in the file
                  * att_lvlh_normal   -- TYPE:array(nt)              Attitude register bit 0: LVLH Normal
                  * att_lvlh_reverse  -- TYPE:array(nt)              Attitude register bit 1: LVLH Reverse
                  * att_limb_pointing -- TYPE:array(nt)              Attitude register bit 2: Earth Limb Pointing
                  * att_conjugate     -- TYPE:array(nt)              Attitude register bit 6: Conjugate Maneuver
                  * att_zero_wind     -- TYPE:array(nt)              Attitude register bit 10: Zero Wind Maneuver
                  * mag_lat           -- TYPE:array(ny,nt), UNITS:deg.     Magnetic quasi-dipole latitude
                  * mag_lon           -- TYPE:array(ny,nt), UNITS:deg.     Magnetic quasi-dipole longitude
                  * sza               -- TYPE:array(ny,nt), UNITS:deg.     Solar zenith angle
                  * slt               -- TYPE:array(ny,nt), UNITS:hour.    Solar local time
                  * ver_dc            -- TYPE:array(ny,nt), UNITS:ph/cm^3/s. Volume emission rate derived from 
                                                                           interferogram DC value.
                  * bin_size          -- TYPE:int.                    Number of pixels binned together in altitude dimension.
                  * orbit_node        -- TYPE:array(nt).              Descending/Ascending orbit flag
                                                                      0 = ICON latitude is increasing (ascending)
                                                                      1 = ICON latitude is decreasing (descending)
                  * orbit_number      -- TYPE:array(nt).              ICON orbit number
                  * chi2              -- TYPE:array(ny,nt). UNITS:rad^2. (std of unwrapped phase)^2 for each row
                  * zero_wind_phase   -- TYPE:array(ny,nt). UNITS:rad. The zero wind phase used for each exposure.
    '''
    
    assert bool(keep_att) + bool(skip_att) + bool(tstartstop) <= 1, "More than one optional input specified"
    
    d = netCDF4.Dataset(L21_fn,'r')

    sens  = d.Instrument[-1] # 'A' or 'B'
    color = L21_fn.split('/')[-1].split('_')[3].split('-')[2] # Red or Green
    versrev = L21_fn.split('/')[-1].split('_')[-1].split('.')[0] # e.g., v01r001
    vers = int(versrev[1:3])
    rev = int(versrev[4:])
    prefix = 'ICON_L21'
    N = len(d.variables['Epoch'])
    assert N>0, "Epoch has length 0: No samples in file %s"%L21_fn

    # Helper function to load data, replace masked array with regular array (nan-filled), etc.
    def read(v):
        ''' Read the variable "v" from the netCDF file and return it in a friendly format '''
        return np.array(np.ma.filled(d.variables[v][...], np.nan))

    # Load up a few variables and find which time indices should be loaded
    att = {}
    att['att_lvlh_normal']   = read('%s_Attitude_LVLH_Normal' % prefix)
    att['att_lvlh_reverse']  = read('%s_Attitude_LVLH_Reverse' % prefix)
    att['att_limb_pointing'] = read('%s_Attitude_Limb_Pointing' % prefix)
    att['att_conjugate']     = read('%s_Attitude_Conjugate_Maneuver' % prefix)
    try: # This was added in v05
        att['att_zero_wind']     = read('%s_Attitude_Zero_Wind_Maneuver' % prefix)
    except KeyError:
        att['att_zero_wind'] = np.zeros_like(att['att_conjugate'])
    time_msec                = read('Epoch')
    time                     = np.array([datetime(1970,1,1) + timedelta(seconds = 1e-3*s) for s in time_msec])
    idx_good = np.ones(N, dtype=bool) # boolean index
    if skip_att:
        for k in skip_att:
            assert k in att.keys(), '"%s" is not a recognized attitude status name' % k
            for i in range(N):
                if att[k][i]:
                    idx_good[i] = False
    if keep_att:
        idx_good = np.zeros(N, dtype=bool) # boolean index
        for k in keep_att:
            assert k in att.keys(), '"%s" is not a recognized attitude status name' % k
            for i in range(N):
                if att[k][i]:
                    idx_good[i] = True
    if tstartstop:
        idx_good = np.zeros(N, dtype=bool) # boolean index
        tstart, tstop = tstartstop
        for i in range(N):
            if (tstart <= time[i]) and (time[i] <= tstop):
                idx_good[i] = True
    if skip_bad_alts:
        alt = read('%s_Altitude' % prefix)[:,0]
        # The following thresholds were based off of an analysis of 2 years of data and exclude <0.01% of the data. 
        # In spot checks of excluded points, they all occurred during recoveries from slews, when data are labeled
        # as science data too early (BJH 28 Oct 2022).
        if color == 'Red':
            altmin, altmax = 157, 163 
        elif color == 'Green':
            altmin, altmax = 87, 89.5
        alt_good = (alt > altmin) & (alt < altmax)
        idx_good = idx_good & alt_good
        
    
    assert idx_good.any(), "All samples have attitude status bits in the disallowed list: %s" % skip_att   
    
    # Load variables
    z = {}
    z['lat']              = read('%s_Latitude' % prefix)[idx_good,:]
    z['lon']              = read('%s_Longitude' % prefix)[idx_good,:]
    z['alt']              = read('%s_Altitude' % prefix)[idx_good,:]
    z['los_wind']         = read('%s_Line_of_Sight_Wind' % prefix)[idx_good,:]
    z['los_wind_error']   = read('%s_Line_of_Sight_Wind_Error' % prefix)[idx_good,:]
    z['local_az']         = read('%s_Line_of_Sight_Azimuth' % prefix)[idx_good,:]
    z['amp']              = read('%s_Fringe_Amplitude' % prefix)[idx_good,:]
    z['amp_error']        = read('%s_Fringe_Amplitude_Error' % prefix)[idx_good,:]
    z['ver']              = read('%s_Relative_VER' % prefix)[idx_good,:]
    z['ver_error']        = read('%s_Relative_VER_Error' % prefix)[idx_good,:]
    z['wind_quality']     = read('%s_Wind_Quality' % prefix)[idx_good,:]
    z['ver_quality']      = read('%s_VER_Quality' % prefix)[idx_good,:]
    z['time']             = time[idx_good]
    z['icon_lat']         = read('%s_Observatory_Latitude' % prefix)[idx_good]
    z['icon_lon']         = read('%s_Observatory_Longitude' % prefix)[idx_good]
    z['icon_alt']         = read('%s_Observatory_Altitude' % prefix)[idx_good]
    for k in att.keys():
        z[k] = att[k][idx_good]
    z['exp_time']         = read('%s_Exposure_Time' % prefix)[idx_good]
    z['sensor']           = sens
    z['version']          = vers
    z['revision']         = rev
    z['source_files']     = [L21_fn]
    z['quality_flags']    = read('%s_Quality_Flags' % prefix)[idx_good,:]
    z['acknowledgement']  = d.Acknowledgement
    z['emission_color']   = color.lower()
    z['sza']              = read('%s_Solar_Zenith_Angle' % prefix)[idx_good]
    try:
        z['slt']              = read('%s_Solar_Local_Time' % prefix)[idx_good]
    except:
        z['slt']              = read('%s_Local_Solar_Time' % prefix)[idx_good]
    z['mag_lat']          = read('%s_Magnetic_Latitude' % prefix)[idx_good]
    z['mag_lon']          = read('%s_Magnetic_Longitude' % prefix)[idx_good]
    z['ver_dc']           = read('%s_Relative_VER_DC' % prefix)[idx_good]
    z['bin_size']         = read('%s_Bin_Size' % prefix)[...].item()
    z['chi2']             = read('%s_Chi2' % prefix)[idx_good,:]
    # Load new variables in a backwards compatible way
    try:
        z['orbit_node'] = read('%s_Orbit_Node' % prefix)[idx_good]
    except KeyError:
        print('WARNING: Orbit_Node variable not found.')
        z['orbit_node'] = np.nan * np.ones_like(z['icon_lat']) # fill in with nan
    try:
        z['orbit_number'] = read('%s_Orbit_Number' % prefix)[idx_good]
    except KeyError:
        print('WARNING: Orbit_Number variable not found.')
        z['orbit_number'] = np.nan * np.ones_like(z['icon_lat']) # fill in with nan
    try:
        z['zero_wind_ref'] = d.Zero_Wind_Ref
    except Exception as e:
        print('WARNING: Zero_Wind_Ref attribute not found.')
        z['zero_wind_ref'] = ''
    try:
        z['corr_notch_drift'] = bool(d.Corr_Notch_Drift)
    except Exception as e:
        print('WARNING: Corr_Notch_Drift attribute not found.')
        z['corr_notch_drift'] = False
    try:
        z['zero_wind_phase'] = read('%s_Zero_Wind_Phase' % prefix)[idx_good]
    except KeyError:
        print('WARNING: Zero_Wind_Phase variable not found.')
        z['zero_wind_phase'] = np.nan * np.ones_like(z['los_wind']) # fill in with nan
    try:
        z['los_wind_precision_1_sample'] = read('%s_Line_of_Sight_Wind_Precision_1_Sample' % prefix)[idx_good]
    except KeyError:
        print('WARNING: Wind_Precision_1_Sample variable not found.')
        z['los_wind_precision_1_sample'] = np.nan * np.ones_like(z['los_wind']) # fill in with nan
    try:
        z['los_wind_precision_1_day'] = read('%s_Line_of_Sight_Wind_Precision_1_Day' % prefix)[idx_good]
    except KeyError:
        print('WARNING: Wind_Precision_1_Day variable not found.')
        z['los_wind_precision_1_day'] = np.nan * np.ones_like(z['los_wind']) # fill in with nan
    try:
        z['los_wind_precision_low_signal'] = read('%s_Line_of_Sight_Wind_Precision_Low_Signal_Effect' % prefix)[idx_good]
    except KeyError:
        print('WARNING: Wind_Precision_Low_Signal variable not found.')
        z['los_wind_precision_low_signal'] = np.nan * np.ones_like(z['los_wind']) # fill in with nan
    try:
        z['los_wind_accuracy'] = read('%s_Line_of_Sight_Wind_Accuracy' % prefix)[idx_good]
    except KeyError:
        print('WARNING: Wind_Accuracy variable not found.')
        z['los_wind_accuracy'] = np.nan * np.ones_like(z['los_wind']) # fill in with nan


    for v in z.keys():
        if np.ndim(z[v])>1:
            # Swap so time is second dimension, not first
            z[v] = np.swapaxes(z[v], 0, 1) # This is the same as z[v].T for 2D arrays

    d.close()
    
    return z

    
    
    
    
############################################################################################################
##########################################       Level 2.2       ###########################################
############################################################################################################
   

def geog_to_mag_coords(u, v, w, lat, lon, alt, t, fast=False):
    '''
    Convert a wind (u,v,w) in geographic (east, north, up) coordinates to geomagnetic (field-aligned, meridional, 
    zonal) coordinates, using the coordinates defined by pysatMagVect (https://github.com/rstoneback/pysatMagVect)
    
    INPUTS:
    
     *  u              -- TYPE:array(n)               Geographic zonal (eastward) component of wind
     *  v              -- TYPE:array(n)               Geographic northward component of wind
     *  w              -- TYPE:array(n)               Geographic upward component of wind
     *  lat            -- TYPE:array(n) -- UNITS:deg. Geographic latitude
     *  lon            -- TYPE:array(n) -- UNITS:deg. Geographic longitude
     *  alt            -- TYPE:array(n) -- UNITS:km.  Altitude
     *  t              -- TYPE:array(n) of datetimes. Time
     
    OPTIONAL INPUTS:
     
     *  fast           -- TYPE:bool                   If True (default) then compute the rotation matrix for only
                                                      the first and last elements of the array and linearly interpolate 
                                                      for the rest.
     
    OUTPUTS:
     
     *  wind_mag_fa    -- TYPE:array(n)               Magnetic field-aligned component of wind
     *  wind_mag_mer   -- TYPE:array(n)               Magnetic meridional component of wind
     *  wind_mag_zon   -- TYPE:array(n)               Magnetic zonal component of wind
     
    '''
    # Ensure array input, not list
    u = np.array(u)
    v = np.array(v)
    w = np.array(w)
    lat = np.array(lat)
    lon = np.mod(np.array(lon), 360) # Russell's code needs -180 to 360
    alt = np.array(alt)
    t = np.array(t)
    
    n = len(u)
    ze = np.nan * np.zeros(n)
    zn = np.nan * np.zeros(n)
    zu = np.nan * np.zeros(n)
    be = np.nan * np.zeros(n)
    bn = np.nan * np.zeros(n)
    bu = np.nan * np.zeros(n)
    me = np.nan * np.zeros(n)
    mn = np.nan * np.zeros(n)
    mu = np.nan * np.zeros(n)
    
    if np.isnan(u).all():
        return np.nan*np.zeros(n), np.nan*np.zeros(n), np.nan*np.zeros(n)
    
    # If only zero or one element is non-nan, then use fast=False to avoid corner cases with fast=True
    if np.isfinite(u).sum() <= 1:
        fast = False
    
    if not fast: # Compute the entire array
        
        # Skip elements that are nan
        i = np.isfinite(u)
        try: # pySatMagVect v0.4.0 uses the ref_height parameter.
             # The new (as yet unversioned) development code does not support this parameter.
            zx,zy,zz,bx,by,bz,mx,my,mz = psmv.calculate_mag_drift_unit_vectors_ecef(lat[i], lon[i], alt[i], t[i], ref_height=80.)
        except:
            zx,zy,zz,bx,by,bz,mx,my,mz = psmv.calculate_mag_drift_unit_vectors_ecef(lat[i], lon[i], alt[i], t[i])
            
        
        # Convert unit vectors from ECEF to ENU using pysatMagVect calculation
        ze0, zn0, zu0 = psmv.ecef_to_enu_vector(zx, zy, zz, lat[i], lon[i])
        be0, bn0, bu0 = psmv.ecef_to_enu_vector(bx, by, bz, lat[i], lon[i])
        me0, mn0, mu0 = psmv.ecef_to_enu_vector(mx, my, mz, lat[i], lon[i])   
        
        ze[i] = ze0
        zn[i] = zn0
        zu[i] = zu0
        be[i] = be0
        bn[i] = bn0
        bu[i] = bu0
        me[i] = me0
        mn[i] = mn0
        mu[i] = mu0
        
        
    else: # Fast: compute the first and last elements and linearly interpolate the rest
        # Be careful not to use the first or last element if they are nan.
        i = np.where(np.isfinite(u))[0]
        i0 = i[0]
        i1 = i[-1]
        
        try: # pySatMagVect v0.4.0 uses the ref_height parameter.
             # The new (as yet unversioned) development code does not support this parameter.
            zx,zy,zz,bx,by,bz,mx,my,mz = psmv.calculate_mag_drift_unit_vectors_ecef([lat[i0], lat[i1]], 
                                                                                        [lon[i0], lon[i1]],
                                                                                        [alt[i0], alt[i1]],
                                                                                        [t[i0], t[i1]],
                                                                                        ref_height=80.)
        except:
            zx,zy,zz,bx,by,bz,mx,my,mz = psmv.calculate_mag_drift_unit_vectors_ecef([lat[i0], lat[i1]], 
                                                                                        [lon[i0], lon[i1]],
                                                                                        [alt[i0], alt[i1]],
                                                                                        [t[i0], t[i1]])
        
        # Convert unit vectors from ECEF to ENU using pysatMagVect calculation
        ze0, zn0, zu0 = psmv.ecef_to_enu_vector(zx, zy, zz, [lat[i0], lat[i1]], [lon[i0], lon[i1]])
        be0, bn0, bu0 = psmv.ecef_to_enu_vector(bx, by, bz, [lat[i0], lat[i1]], [lon[i0], lon[i1]])
        me0, mn0, mu0 = psmv.ecef_to_enu_vector(mx, my, mz, [lat[i0], lat[i1]], [lon[i0], lon[i1]])
        
        ze[i0:i1+1] = np.linspace(ze0[0], ze0[-1], i1-i0+1)
        zn[i0:i1+1] = np.linspace(zn0[0], zn0[-1], i1-i0+1)
        zu[i0:i1+1] = np.linspace(zu0[0], zu0[-1], i1-i0+1)
        be[i0:i1+1] = np.linspace(be0[0], be0[-1], i1-i0+1)
        bn[i0:i1+1] = np.linspace(bn0[0], bn0[-1], i1-i0+1)
        bu[i0:i1+1] = np.linspace(bu0[0], bu0[-1], i1-i0+1)
        me[i0:i1+1] = np.linspace(me0[0], me0[-1], i1-i0+1)
        mn[i0:i1+1] = np.linspace(mn0[0], mn0[-1], i1-i0+1)
        mu[i0:i1+1] = np.linspace(mu0[0], mu0[-1], i1-i0+1)
        
    
    # Perform coordinate rotation
    wind_mag_fa  = u*be + v*bn + w*bu
    wind_mag_zon = u*ze + v*zn + w*zu
    wind_mag_mer = u*me + v*mn + w*mu
    
    return wind_mag_fa, wind_mag_mer, wind_mag_zon



def level22_alt_grid(emission_color, bin_size=1):
    '''
    Return the L2.2 grid, onto which the L2.1 data will be interpolated and combined. This grid was determined separately using
    2 months of ICON data. The grid was defined such that the bottom grid point matches the median sampled bottom altitude.
    This altitude is fairly steady (+/- about 1 km) since the ICON pointing is commanded to keep this value constant. The 
    grid spacing was chosen to faithfully represent the grid spacing in the data (~2.9 km at 90 km altitude, and 
    ~2.2 km at 300 km altitude).
    
    INPUT:
    
      * emission_color                               'red' or 'green'
      * bin_size       -- TYPE:int.                  The desired binning (presumably matching the L2.1 data).
                                                     (Default 1)
      
    OUTPUT:
    
      * z              -- TYPE:array(ny), UNITS:km.  The 2.2 altitude grid.
    '''
    
    assert emission_color in ['red','green'], "Use emission_color = 'red' or 'green'"

    # See Jupyter notebook L22_alt_grid.ipynb
    zg = np.array([ 88.38205433,  91.3341588 ,  94.27706017,  97.21078712,
                   100.13536825, 103.05083208, 105.95720702, 108.85452142,
                   111.74280351, 114.62208146, 117.49238334, 120.35373713,
                   123.20617072, 126.04971193, 128.88438847, 131.71022798,
                   134.52725801, 137.33550602, 140.1349994 , 142.92576542,
                   145.70783131, 148.48122418, 151.24597107, 154.00209893,
                   156.74963464, 159.48860498, 162.21903665, 164.94095627,
                   167.65439038, 170.35936542, 173.05590778, 175.74404374,
                   178.42379951, 181.09520121, 183.75827488, 186.41304649,
                   189.05954192, 191.69778697, 194.32780736, 196.94962873,
                   199.56327664, 202.16877657, 204.76615393, 207.35543403,
                   209.93664212, 212.50980335, 215.07494283, 217.63208555,
                   220.18125644, 222.72248037, 225.25578209, 227.78118631,
                   230.29871765, 232.80840065, 235.31025979, 237.80431944,
                   240.29060392, 242.76913748, 245.23994428, 247.7030484 ,
                   250.15847386, 252.60624459, 255.04638446, 257.47891726,
                   259.90386671, 262.32125643, 264.73111001, 267.13345094,
                   269.52830262, 271.91568843, 274.29563161, 276.66815539,
                   279.03328289, 281.39103716, 283.7414412 , 286.08451792,
                   288.42029015, 290.74878068, 293.07001219, 295.38400733,
                   297.69078865, 299.99037865, 302.28279972, 304.56807424])
    
    zr = np.array([151.99429417, 154.75833368, 157.51385622, 160.26088806,
                   162.99945534, 165.72958414, 168.45130048, 171.16463027,
                   173.86959935, 176.56623349, 179.25455837, 181.93459959,
                   184.60638268, 187.26993308, 189.92527616, 192.57243722,
                   195.21144145, 197.84231401, 200.46507993, 203.07976421,
                   205.68639175, 208.28498736, 210.8755758 , 213.45818175,
                   216.03282979, 218.59954446, 221.15835019, 223.70927136,
                   226.25233225, 228.7875571 , 231.31497005, 233.83459516,
                   236.34645644, 238.85057781, 241.34698311, 243.83569613,
                   246.31674056, 248.79014004, 251.25591812, 253.71409828,
                   256.16470394, 258.60775845, 261.04328505, 263.47130696,
                   265.8918473 , 268.30492911, 270.71057538, 273.10880903,
                   275.49965289, 277.88312974, 280.25926227, 282.62807312,
                   284.98958484, 287.34381993, 289.6908008 , 292.03054982,
                   294.36308927, 296.68844135, 299.00662822, 301.31767196,
                   303.62159457, 305.91841801, 308.20816413, 310.49085476])

    if emission_color == 'red':
        z = zr
    else: # emission_color == 'green':
        z = zg

    
    # Handle binning. 
    # Although not immediately obvious, the correct thing to do here is to take the max
    # of the values in each bin. This is because the inversion attributes the measurement
    # location to a half-grid higher altitude, and binning takes the mean altitude.
    # L2_unbinned --> L1_unbinned --> L1_binned --> L2_binned
    # NOTE: integration_order = 0 is assumed
    zb = bin_array(bin_size, z, method='max')
    return zb


    
def level21_dict_to_level22_dict(L21_A_dict, L21_B_dict, sph_asym_thresh = None, time_start = None, time_stop = None,
                                 t_diff_AB = None):
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
                                           the start time defaults to 0 UT on the date of the midpoint between the
                                           first and last input times.
      *  time_stop       -- TYPE:datetime. A timezone-naive datetime in UT, specifying the end of the interval
                                           which the Level 2.2 data product should span. (Some L2.1 files from after
                                           the start time are needed to handle the crossover). If None (default), 
                                           the stop time defaults to 24 hours after the start time.
      *  t_diff_AB       -- TYPE:float.    Maximum time difference between A and B exposures used to determine
                                           the wind at a single point in space. If it's longer than this, there
                                           is likely a problem with the processing. If None (default), the default 
                                           value from MIGHTI_L2.global_params will be used.
                                        
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
                    * u_precision_1_sample   -- TYPE:array(ny,nx),    UNITS:m/s.  1-sample precision in u. See Var_Notes.
                    * v_precision_1_sample   -- TYPE:array(ny,nx),    UNITS:m/s.  1-sample precision in v. See Var_Notes.
                    * u_precision_1_day      -- TYPE:array(ny,nx),    UNITS:m/s.  1-day precision in u. See Var_Notes.
                    * v_precision_1_day      -- TYPE:array(ny,nx),    UNITS:m/s.  1-day precision in v. See Var_Notes.
                    * u_precision_low_signal -- TYPE:array(ny,nx),    UNITS:m/s.  Low-signal precision in u. See Var_Notes.
                    * v_precision_low_signal -- TYPE:array(ny,nx),    UNITS:m/s.  Low-signal precision in v. See Var_Notes.
                    * u_accuracy      -- TYPE:array(ny,nx),    UNITS:m/s.  Accuracy in u. See Var_Notes.
                    * v_accuracy      -- TYPE:array(ny,nx),    UNITS:m/s.  Accuracy in v. See Var_Notes.
                    * u_error         -- TYPE:array(ny,nx),    UNITS:m/s.  "Total" uncertainty in u. See Var_Notes.
                    * v_error         -- TYPE:array(ny,nx),    UNITS:m/s.  "Total" uncertainty in v. See Var_Notes.
                    * quality_flags   -- TYPE:array(ny,nx,ne), UNITS:none. The quality flags (either 0 or 1) for each point
                                                                           in the grid. Each point has a number of flags, 
                                                                           which are set to 1 under the following 
                                                                           circumstances
                                                                           
                                                                           * 0 : (From L1 A) SNR too low to reliably perform L1 processing
                                                                           * 1 : (From L1 A) Proximity to South Atlantic Anomaly
                                                                           * 2 : (From L1 A) Bad calibration 
                                                                           * 3 : (From L1 A) Calibration lamps are on
                                                                           * 4 : (From L1 A) Solar/lunar contamination
                                                                           * 5 : (From L2.1 A) Not enough valid points in profile
                                                                           * 6 : (From L2.1 A) SNR too low after inversion
                                                                           * 7 : (From L2.1 A) Significant airglow above 300 km
                                                                           * 8 : (From L2.1 A) Line of sight crosses the terminator
                                                                           * 9 : (From L2.1 A) Thermal drift correction is uncertain
                                                                           * 10: (From L2.1 A) S/C pointing is not stable
                                                                           * 11: (From L2.1 A) SNR is low after inversion, but maybe usable
                                                                           * 12: (From L1 B) SNR too low to reliably perform L1 processing
                                                                           * 13: (From L1 B) Proximity to South Atlantic Anomaly
                                                                           * 14: (From L1 B) Bad calibration 
                                                                           * 15: (From L1 B) Calibration lamps are on
                                                                           * 16: (From L1 B) Solar/lunar contamination
                                                                           * 17: (From L2.1 B) Not enough valid points in profile
                                                                           * 18: (From L2.1 B) SNR too low after inversion
                                                                           * 19: (From L2.1 B) Significant airglow above 300 km
                                                                           * 20: (From L2.1 B) Line of sight crosses the terminator
                                                                           * 21: (From L2.1 B) Thermal drift correction is uncertain
                                                                           * 22: (From L2.1 B) S/C pointing is not stable
                                                                           * 23: (From L2.1 B) SNR is low after inversion, but maybe usable
                                                                           * 24: Missing MIGHTI-A file
                                                                           * 25: Missing MIGHTI-B file
                                                                           * 26: MIGHTI-A did not sample this altitude
                                                                           * 27: MIGHTI-B did not sample this altitude
                                                                           * 28: Spherical asymmetry: A&B VER estimates disagree
                                                                           * 29: Mixing Normal and Reverse LVLH
                                                                           * 30: Unused
                                                                           * 31: Unused
                                                                           * 32: Unused
                                                                           * 33: Unknown error (notify MIGHTI team)
                                                                           
                                                                            
                   * epoch_full      -- TYPE:array(ny,nx),    UNITS:none. The average between the time of the MIGHTI A 
                                                                          and B measurements that contribute to this 
                                                                          grid point, given as a datetime object. 
                   * epoch           -- TYPE:array(nx),       UNITS:none. mean(epoch_full,axis=0), accounting for missing values, and
                                                                          given as a datetime object. This is useful because 
                                                                          while the time varies along a column, it doesn't vary 
                                                                          much (+/- 30 seconds or so).
                   * time_start      -- TYPE:datetime                     The start time for defining the reconstruction grid (inclusive)
                   * time_stop       -- TYPE:datetime                     The stop time for defining the reconstruction grid (exclusive)
                   * time_delta      -- TYPE:array(ny,nx),    UNITS:s.    The difference between the time of the MIGHTI
                                                                          A and B measurements that contribute to this 
                                                                          grid point. (MIGHTI-B time minus MIGHTI-A time)
                   * fringe_amp      -- TYPE:array(ny,nx),    UNITS:none. The fringe amplitude, estimated by combining (averaging)
                                                                          estimates from A and B.
                   * fringe_amp_error -- TYPE:array(ny,nx),   UNITS:none. Statistical (1-sigma) error in fringe_amp
                   * fringe_amp_A    -- TYPE:array(ny,nx),    UNITS:none. The fringe amplitude measured by MIGHTI A (roughly
                                                                          proportional to VER) interpolated to the L2.2 grid
                   * fringe_amp_B    -- TYPE:array(ny,nx),    UNITS:none. The fringe amplitude measured by MIGHTI B (roughly
                                                                          proportional to VER) interpolated to the L2.2 grid
                   * fringe_amp_rel_diff -- TYPE:array(ny,nx),UNITS:none. The difference between the fringe amplitude 
                                                                          of the MIGHTI A and B measurements that 
                                                                          contribute to this grid point, divided by the
                                                                          mean. When this is high, it indicates that 
                                                                          spherical asymmetry may be a problem.
                   * ver           -- TYPE:array(ny,nx), UNITS:ph/cm^3/s. Volume emission rate, estimated by combining (averaging) 
                                                                          VER estimates from MIGHTI A and B.
                   * ver_error     -- TYPE:array(ny,nx), UNITS:ph/cm^3/s. Statistical (1-sigma) error in ver
                   * ver_A         -- TYPE:array(ny,nx), UNITS:ph/cm^3/s. VER measured by MIGHTI-A, interpolated to the L2.2 grid
                   * ver_B         -- TYPE:array(ny,nx), UNITS:ph/cm^3/s. VER measured by MIGHTI-B, interpolated to the L2.2 grid
                   * ver_rel_diff  -- TYPE:array(ny,nx), UNITS:none.      The difference between the VER estimates 
                                                                          from MIGHTI A and B, divided by the
                                                                          mean. When this is high, it indicates that 
                                                                          spherical asymmetry may be a problem.
                   * wind_quality  -- TYPE:array(ny,nx),                  A number between 0 (bad) and 1 (good) for each grid point
                   * ver_quality   -- TYPE:array(ny,nx),                  A number between 0 (bad) and 1 (good) for each grid point.
                   * emission_color  -- TYPE:str,                         'red' or 'green'        
                   * source_files    -- TYPE:list of str,                 All the files used to create the data product,
                                                                          including the full paths
                   * acknowledgement -- TYPE:str,                         The Acknowledgement copied from the L2.1 files
                   * N_used_A        -- TYPE:array(nyA, n_files_A)        How many times each input MIGHTI-A file was used
                   * N_used_B        -- TYPE:array(nyB, n_files_B)        How many times each input MIGHTI-B file was used
                   * wind_mag_fa     -- TYPE:array(ny,nx), UNITS:m/s.     The component of the wind in the direction of the
                                                                          magnetic field.
                   * wind_mag_mer    -- TYPE:array(ny,nx), UNITS:m/s.     The component of the wind in the magnetic meridional 
                                                                          direction (as defined by pysatMagVect)
                   * wind_mag_zon    -- TYPE:array(ny,nx), UNITS:m/s.     The component of the wind in the magnetic zonal 
                                                                          direction (as defined by pysatMagVect)
                   * mag_lat         -- TYPE:array(ny,nx), UNITS:deg.     Magnetic quasi-dipole latitude
                   * mag_lon         -- TYPE:array(ny,nx), UNITS:deg.     Magnetic quasi-dipole longitude
                   * sza             -- TYPE:array(ny,nx), UNITS:deg.     Solar zenith angle
                   * slt             -- TYPE:array(ny,nx), UNITS:hour.    Solar local time
                   * orbit_node      -- TYPE:array(ny,nx).                Descending/Ascending orbit flag. This is usually either of:
                                                                          0 = Ascending: ICON latitude is increasing
                                                                          1 = Descending: ICON latitude is decreasing
                                                                          But when samples from both ascending and descending
                                                                          parts of the orbit are used, it is an interpolated value.
                   * orbit_number    -- TYPE:array(ny,nx).                Orbit number. This is usually an integer but is fractional
                                                                          when samples from different orbits are used.
                   * zero_wind_ref   -- TYPE:str.                         The zero wind reference used to create the Level 2.1 files.
                                                                          'external': uses HWM14 or some outside reference
                                                                          'internal': uses on-orbit value with on correction at Level 2
                   * corr_notch_drift-- TYPE:bool.                        Whether a correction for notch drift was used to create the Level 2.1 files.
    '''
            
    N_flags = 34 # Update this if the number of quality flags changes.

    ################################## Parse Inputs ####################################
    emission_color = L21_A_dict['emission_color']
    params = global_params[emission_color]
    if sph_asym_thresh is None:
        sph_asym_thresh = params['sph_asym_thresh']
    if t_diff_AB is None:
        t_diff_AB = params['t_diff_AB']

    assert L21_A_dict['emission_color'] == L21_B_dict['emission_color'], "Files for A and B are for different emissions"
#     assert L21_A_dict['zero_wind_ref'] == L21_B_dict['zero_wind_ref'], "Files for A and B use different zero wind references"
#     assert L21_A_dict['corr_notch_drift'] == L21_B_dict['corr_notch_drift'], "Files for A and B use different corr_notch_drift"
    # Make sure a LVLH Normal/Reverse maneuver doesn't happen in the middle of the dataset.
    # 2021 Sep 14 BJH: This is commented out because it didn't actually cause a problem with the geometry
    #for d in [L21_A_dict, L21_B_dict]:
    #    assert (all(d['att_lvlh_normal'])  and not any(d['att_lvlh_reverse'])) or \
    #           (all(d['att_lvlh_reverse']) and not any(d['att_lvlh_normal']))

    lat_A      =       L21_A_dict['lat']
    lon_raw_A  =       L21_A_dict['lon']
    alt_A      =       L21_A_dict['alt']
    los_wind_A =       L21_A_dict['los_wind']
    los_wind_A_precision_1_sample =   L21_A_dict['los_wind_precision_1_sample']
    los_wind_A_precision_1_day    =   L21_A_dict['los_wind_precision_1_day']
    los_wind_A_precision_low_signal=  L21_A_dict['los_wind_precision_low_signal']
    los_wind_A_accuracy           =   L21_A_dict['los_wind_accuracy']
    los_wind_A_err =   L21_A_dict['los_wind_error']
    local_az_A =       L21_A_dict['local_az']
    amp_A      =       L21_A_dict['amp']
    amp_A_err  =       L21_A_dict['amp_error']
    ver_A      =       L21_A_dict['ver']
    ver_A_err  =       L21_A_dict['ver_error']
    time_A     =       L21_A_dict['time']
    exp_time_A =       L21_A_dict['exp_time']
    emission_color =   L21_A_dict['emission_color']
    icon_lon_A =       L21_A_dict['icon_lon']
    wind_quality_A =   L21_A_dict['wind_quality']
    ver_quality_A =    L21_A_dict['ver_quality']
    qflags_A   =       L21_A_dict['quality_flags']
    mlat_A     =       L21_A_dict['mag_lat']
    sza_A      =       L21_A_dict['sza']
    orb_num_A  =       L21_A_dict['orbit_number']
    orb_desc_A =       L21_A_dict['orbit_node']
    N_alts_A, N_times_A = np.shape(lat_A)

    lat_B      =       L21_B_dict['lat']
    lon_raw_B  =       L21_B_dict['lon']
    alt_B      =       L21_B_dict['alt']
    los_wind_B =       L21_B_dict['los_wind']
    los_wind_B_precision_1_sample =   L21_B_dict['los_wind_precision_1_sample']
    los_wind_B_precision_1_day    =   L21_B_dict['los_wind_precision_1_day']
    los_wind_B_precision_low_signal=  L21_B_dict['los_wind_precision_low_signal']
    los_wind_B_accuracy           =   L21_B_dict['los_wind_accuracy']
    los_wind_B_err =   L21_B_dict['los_wind_error']
    local_az_B =       L21_B_dict['local_az']
    amp_B      =       L21_B_dict['amp']
    amp_B_err  =       L21_B_dict['amp_error']
    ver_B      =       L21_B_dict['ver']
    ver_B_err  =       L21_B_dict['ver_error']
    time_B     =       L21_B_dict['time']
    exp_time_B =       L21_B_dict['exp_time']
    emission_color =   L21_B_dict['emission_color']
    icon_lon_B =       L21_B_dict['icon_lon']
    wind_quality_B =   L21_B_dict['wind_quality']
    ver_quality_B =    L21_B_dict['ver_quality']
    qflags_B   =       L21_B_dict['quality_flags']
    mlat_B     =       L21_B_dict['mag_lat']
    sza_B      =       L21_B_dict['sza']
    orb_num_B  =       L21_B_dict['orbit_number']
    orb_desc_B =       L21_B_dict['orbit_node']
    N_alts_B, N_times_B = np.shape(lat_B)
    
    # Unwrap mag lon and SLT to avoid averaging over discontinuities. They will be "unwrapped"
    # everywhere in the following code, but will be modded right before saving.
    mlon_A = fix_longitudes_mat(L21_A_dict['mag_lon'])
    mlon_B = fix_longitudes_mat(L21_B_dict['mag_lon'])
    # SLT can be unwrapped by leveraging longitude unwrapping code (x15)
    slt_A  = fix_longitudes_mat(L21_A_dict['slt']*15.)/15.
    slt_B  = fix_longitudes_mat(L21_B_dict['slt']*15.)/15.
    
    if time_start is None:
        # Default to 0 UT on date of midpoint of first and last times in the file
        t_min = min(min(time_A), min(time_B))
        t_max = max(max(time_A), max(time_B))
        dt = (t_max - t_min).total_seconds()
        assert (dt > 3600.) or (t_max.day == t_min.day), "Less than 1 hour of data, spanning 2 dates. Unsure which day to process."
        t_mid = t_min + timedelta(seconds=dt/2)
        time_start = datetime(t_mid.year, t_mid.month, t_mid.day)

    if time_stop is None:
        # Default to 24 hours after the start time
        time_stop = time_start + timedelta(hours=24)

    assert any(time_A < time_stop ), "No MIGHTI-A files found before time_stop=%s"%(time_stop)
    assert any(time_B < time_stop ), "No MIGHTI-B files found before time_stop=%s"%(time_stop)
    assert any(time_A >= time_start), "No MIGHTI-A files found after time_start=%s"%(time_start)
    assert any(time_B >= time_start), "No MIGHTI-B files found after time_start=%s"%(time_start)
    assert any((time_A >= time_start) & (time_A < time_stop)), "No MIGHTI-A files found between time_start=%s and time_stop=%s"%\
                                                            (time_start,time_stop)
    assert any((time_B >= time_start) & (time_B < time_stop)), "No MIGHTI-B files found between time_start=%s and time_stop=%s"%\
                                                            (time_start,time_stop)
    assert len(time_A) > 1, "Need at least 2 files from MIGHTI-A"
    assert len(time_B) > 1, "Need at least 2 files from MIGHTI-B"

    ####################### Define reconstruction grid: lon/alt ########################

    # Unwrap longitudes. This must be done carefully to handle cases where large 
    # gaps due to calibration maneuvers may exist. Assume the icon longitude vs
    # time plot is approximately linear, and use this to unwrap ICON longitude. 
    # Then add/subtract 360 deg from tangent longitudes so that the tangent
    # longitudes are "close" (less than 180 deg) from ICON longitude at 
    # the same time.

    # Fit for dlon/dt (basically the speed of ICON, or the frequency of the lon vs time sawtooth)
    # Use MIGHTI-A and MIGHTI-B samples combined (since the ICON position is common)
    time_sec_A = np.array([(t-time_start).total_seconds() for t in time_A])
    time_sec_B = np.array([(t-time_start).total_seconds() for t in time_B])
    time_sec = np.concatenate((time_sec_A, time_sec_B))
    icon_lon = np.concatenate((icon_lon_A, icon_lon_B))
    def err(dlon_dt, *params):
        time_sec, lon = params
        lon_fit = np.mod(lon[0] + dlon_dt*(time_sec-time_sec[0]), 360.)
        return sum((lon_fit-lon)**2)
    dlon_dt_min = 360./(3*3600.)   # orbit every 3 hours
    dlon_dt_max = 360./(0.5*3600.) # orbit every 0.5 hours
    dlon_dt = optimize.brute(err, ranges=((dlon_dt_min,dlon_dt_max),), Ns=500, 
                             args=(time_sec, icon_lon), finish = optimize.fmin)[0]

    # Approximate for ICON's unwrapped longitude as a function of time,
    # using the dlon_dt determined above.
    def icon_lon_linear(t):
        '''
        t - datetime, or array of datetimes
        '''
        # Use A as to define the start time. This shouldn't matter, since
        # the icon_lon should trace the same path for A and B, and we're only
        # using it for 0/360 adjustments anyway.
        t0 = time_A[0]
        lon0 = icon_lon_A[0]
        # Helper function to vectorize timedelta.total_seconds()
        total_seconds = np.vectorize(lambda x: x.total_seconds())
        # Linear approx
        lon = lon0 + dlon_dt * total_seconds(t-t0)
        return lon

    icon_lon_linear_A = icon_lon_linear(time_A)
    icon_lon_linear_B = icon_lon_linear(time_B)
    jumpA = (icon_lon_linear_A - icon_lon_A)/360.
    jumpB = (icon_lon_linear_B - icon_lon_B)/360.
    njumpA = np.around(jumpA)
    njumpB = np.around(jumpB)
    djumpA = jumpA - njumpA
    djumpB = jumpB - njumpB
    assert all(abs(djumpA) < 0.2), "Contact MIGHTI Team - something went wrong in ICON longitude unwrapping."
    assert all(abs(djumpB) < 0.2), "Contact MIGHTI Team - something went wrong in ICON longitude unwrapping."
    # Apply jumps to ICON longitude to unwrap it
    icon_lonu_A = icon_lon_A + njumpA*360.
    icon_lonu_B = icon_lon_B + njumpB*360.

    # Find and apply jumps to tangent longitudes
    lon_A = np.zeros(np.shape(lon_raw_A))
    for i in range(N_alts_A):
        jump = (icon_lon_linear_A - lon_raw_A[i,:])/360.
        njump = np.around(jump)
        djump = jump - njump
        assert all(abs(djump)<0.2), "Contact MIGHTI Team - something went wrong in MIGHTI-A tang lon unwrapping."
        lon_A[i,:] = lon_raw_A[i,:] + njump*360.

    lon_B = np.zeros(np.shape(lon_raw_B))
    for i in range(N_alts_B):
        jump = (icon_lon_linear_B - lon_raw_B[i,:])/360.
        njump = np.around(jump)
        djump = jump - njump
        assert all(abs(djump)<0.2), "Contact MIGHTI Team - something went wrong in MIGHTI-B tang lon unwrapping."
        lon_B[i,:] = lon_raw_B[i,:] + njump*360.

    # Do some sanity checks
    assert (abs(np.diff(lon_A,axis=0))).max() < 5., "Large jump detected in longitude profile (MIGHTI-A)"
    assert (abs(np.diff(lon_B,axis=0))).max() < 5., "Large jump detected in longitude profile (MIGHTI-B)"

    # Determine start and stop longitude to define grid.
    # Use the full range of sampled unwrapped longitudes,
    # including previous and next days, then trim later based on time.
    lon_start = min(lon_A.min(), lon_B.min())
    lon_stop  = max(lon_A.max(), lon_B.max())

    # Define longitude grid in such a way that it respects the L2.1 cadence.
    dlon_A = np.diff(lon_A[0,:])
    dlon_A = np.concatenate((dlon_A,[dlon_A[-1]])) # Add last sample again to make same length
    dlon_B = np.diff(lon_B[0,:])
    dlon_B = np.concatenate((dlon_B,[dlon_B[-1]]))

    dlon_A_interp = interpolate.interp1d(lon_A[0,:], dlon_A, kind='zero', bounds_error=False, fill_value=np.inf)
    dlon_B_interp = interpolate.interp1d(lon_B[0,:], dlon_B, kind='zero', bounds_error=False, fill_value=np.inf)

    loni = lon_start
    dlon_last = 4.0 # default value to start with
    dlon_max = 6.0 # the maximum permittable gap 
    lon_vec = []
    # Build longitude grid iteratively, using the L2.1 cadence to define the L2.2 cadence
    while loni <= lon_stop:
        # Evaluate delta lon for A and B at this longitude, and take min
        dA = dlon_A_interp(loni).item()
        dB = dlon_B_interp(loni).item()
        dlon = min(dA,dB)
        if np.isinf(dlon): # This should not happen
            print('WARNING: possibly invalid value longitude delta definition (dlon=%s)'%dlon)
            # Should this crash, or should we let it fly? UPDATE 2020/04/01: Let it fly. This sometimes occurs if the edges
            # of the time sample occur near a maneuver, or for missing data from A or B.
#             raise Exception('WARNING: invalid value longitude delta definition (dlon=%s)'%dlon) 
        if dlon <= 0.0: # This is a problem
            print('WARNING: invalid value longitude delta definition (dlon=%s)'%dlon)
            # Should this crash, or should we let it fly? I think we have to crash here since the longitude grid doesn't make sense.
            # UPDATE 2022 Nov 1: Let it fly but use the last lon gap. This only happened 3 times in 2 years: with rLVLH transitions
            # and on a mislabeled slew which is now caught.
#             raise Exception('WARNING: invalid value longitude delta definition (dlon=%s)'%dlon) 
            dlon = dlon_last
        if dlon > dlon_max: # if there are really large gaps (e.g., EUV calibration), then use a reasonable cadence
            dlon = dlon_last
        loni += dlon
        lon_vec.append(loni)
#         dlon_last = dlon # TODO: Enable this in v06 (It was a mistake to omit it, but it only slightly changes the grid.)
    lon_vec = np.array(lon_vec)
    
    assert np.std(lon_vec) > 0, "Contact MIGHTI Team - longitude grid not changing"
    assert len(lon_vec) > 0, "Contact MIGHTI Team - zero length longitude grid"

    # Call function to get altitude grid (changed in v1.31). Respect binning.
    b = np.min([L21_A_dict['bin_size'], L21_B_dict['bin_size']])
    alt = level22_alt_grid(emission_color, b) # See function in this module.

    # Create grid for longitude (since we might change how this is done in the future)
    lon,_ = np.meshgrid(lon_vec, alt)
    N_alts, N_lons = np.shape(lon)
    
    # Sanity check before possibly running for hours and hours
    assert N_alts < 100, 'ERROR? Unreasonable number of altitudes in grid (N_alts = %i)' % N_alts
    assert N_lons < 86400./30., 'ERROR? Unreasonable number of longitudes in grid (N_lons = %i)' % N_lons

    ############### Interpolate values to reconstruction grid ##################
    # Use bilinear interpolation. This is somewhat complicated because
    # the sample grid is not exactly regular, but it's close, and 
    # we are approximating it as such. We're implementing our own
    # bilinear interpolation so we can control extrapolation in 
    # longitude and altitude as desired. Bilinear interpolation is 
    # used because it is insensitive to the units used to define
    # the sample grid.
    # This proceeds in steps:
    # 1) Setup
    # 2) Quality flagging
    # 3) Interpolation
    # 4) Inverting (i.e., rotating LoS winds to cardinal)
    # 5) More quality flagging
    # 6) Quality factor calculation
    # 7) Magnetic coordinate conversion
    # 8) Trimming data points
    
    # Output variables, which will be defined on the reconstruction grid
    U = np.nan*np.zeros((N_alts, N_lons))                # zonal wind
    V = np.nan*np.zeros((N_alts, N_lons))                # meridional wind
    U_precision_1_sample   = np.nan*np.zeros((N_alts, N_lons)) # Precision 1 sample
    V_precision_1_sample   = np.nan*np.zeros((N_alts, N_lons))
    U_precision_1_day      = np.nan*np.zeros((N_alts, N_lons)) # Precision 1 day
    V_precision_1_day      = np.nan*np.zeros((N_alts, N_lons))
    U_precision_low_signal = np.nan*np.zeros((N_alts, N_lons)) # Precision low signal
    V_precision_low_signal = np.nan*np.zeros((N_alts, N_lons))
    U_accuracy             = np.nan*np.zeros((N_alts, N_lons)) # Accuracy
    V_accuracy             = np.nan*np.zeros((N_alts, N_lons))
    U_err = np.nan*np.zeros((N_alts, N_lons))            # zonal wind "total" uncertainty
    V_err = np.nan*np.zeros((N_alts, N_lons))            # meridional wind "total" uncertainty
    lat = np.nan*np.zeros((N_alts, N_lons))              # latitude
    time = np.empty((N_alts, N_lons), dtype=object)      # time ascribed to L2.2 data point (as datetime objects)
    time_delta = np.nan*np.zeros((N_alts, N_lons))       # difference between A and B times used (seconds)
    ver_intp_A = np.nan*np.zeros((N_alts, N_lons))       # VER from MIGHTI-A interpolated to grid
    ver_intp_B = np.nan*np.zeros((N_alts, N_lons))       # VER from MIGHTI-B interpolated to grid
    ver   = np.nan*np.zeros((N_alts, N_lons))            # combined VER (mean of A and B)
    ver_err = np.nan*np.zeros((N_alts, N_lons))          # VER uncertainty
    ver_rel_diff = np.nan*np.zeros((N_alts, N_lons))     # relative difference in A and B VER
    amp_intp_A = np.nan*np.zeros((N_alts, N_lons))       # fringe amplitude from A (related to VER) interpolated to grid
    amp_intp_B = np.nan*np.zeros((N_alts, N_lons))       # fringe amplitude from B (related to VER) interpolated to grid
    amp   = np.nan*np.zeros((N_alts, N_lons))            # combined fringe amplitude (mean of A and B)
    amp_err = np.nan*np.zeros((N_alts, N_lons))          # fringe amplitude uncertainty
    amp_rel_diff = np.nan*np.zeros((N_alts, N_lons))     # relative difference in A and B fringe amplitudes
    wind_quality = np.zeros((N_alts, N_lons))            # Wind quality factor, in (0, 1.0)
    ver_quality  = np.zeros((N_alts, N_lons))            # VER quality factor, in (0, 1.0)
    qflags = np.zeros((N_alts, N_lons, N_flags), dtype=bool) # Quality flags, one set per grid point. See above for definition.
    N_used_A = np.zeros((N_alts_A, N_times_A))           # Number of times each L2.1 MIGHTI-A file was used
    N_used_B = np.zeros((N_alts_B, N_times_B))           # Same for MIGHTI-B
    mlat = np.nan*np.zeros((N_alts, N_lons))                    # Mag lat
    mlon = np.nan*np.zeros((N_alts, N_lons))                    # Mag lon
    sza  = np.nan*np.zeros((N_alts, N_lons))                    # solar zenith angle
    slt  = np.nan*np.zeros((N_alts, N_lons))                    # solar local time
    orb_num  = np.nan*np.zeros((N_alts, N_lons))         # Orbit number
    orb_desc = np.nan*np.zeros((N_alts, N_lons))         # Orbit descending flag
    
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


            ##################### Initial Quality Flagging ##########################
            # Never extrapolate in longitude. If this flag is raised here, it probably indicates an 
            # extended calibration routine like the EUV lunar cal that crossed a day boundary.
            # Mark as missing, and continue.
            if kA0 < 0 or kA1 >= N_times_A or kB0 < 0 or kB1 >= N_times_B:
                if kA0 < 0 or kA1 >= N_times_A:
                    qflags[i,k,24] = 1
                if kB0 < 0 or kB1 >= N_times_B:
                    qflags[i,k,25] = 1
                continue

            # Determine if there are "missing" files by checking the time between the straddling
            # files we just found and comparing to the exposure time of the files.
            # If so, throw quality flag and continue. This means that we never interpolate over gaps, 
            # and that one missing A file will likely create 2 missing grid points. This is, however, the
            # right thing to do when bilinear interpolation is being used (as opposed to NN).
            # Note that it's the greater of the previous and next exposure times that matter.
            exp_time_A_thresh = 1.1 * max(exp_time_A[kA0], exp_time_A[kA1]) # Extra 10% included for buffer
            exp_time_B_thresh = 1.1 * max(exp_time_B[kB0], exp_time_B[kB1])
            missing_A = abs((time_A[kA1] - time_A[kA0]).total_seconds()) > exp_time_A_thresh
            missing_B = abs((time_B[kB1] - time_B[kB0]).total_seconds()) > exp_time_B_thresh
            if missing_A or missing_B:
                if missing_A:
                    qflags[i,k,24] = 1
                if missing_B:
                    qflags[i,k,25] = 1
                continue

            # If the desired altitude is outside the range of altitudes sampled by the 
            # instruments, throw quality flag and continue. For example, this flag is 
            # raised at 90 km if MIGHTI-B's lowest tangent point is 95 km.
            # For this, allow some wiggle room to handle case where MIGHTI A samples
            # at 90.01 km but we wanted 90.00 km.
            altmin_A = max(min(alt_A[:,kA0]), min(alt_A[:,kA1])) - (alt[1]  - alt[0])/2.
            altmin_B = max(min(alt_B[:,kB0]), min(alt_B[:,kB1])) - (alt[1]  - alt[0])/2.
            altmax_A = min(max(alt_A[:,kA0]), max(alt_A[:,kA1])) + (alt[-1] - alt[-2])/2.
            altmax_B = min(max(alt_B[:,kB0]), max(alt_B[:,kB1])) + (alt[-1] - alt[-2])/2.
            if alt_pt > min(altmax_A, altmax_B) or alt_pt < max(altmin_A, altmin_B):
                if alt_pt > altmax_A or alt_pt < altmin_A:
                    qflags[i,k,26] = 1
                if alt_pt > altmax_B or alt_pt < altmin_B:
                    qflags[i,k,27] = 1
                continue
                
            # Make sure not to mix Normal and Reverse LVLH
            if L21_A_dict['att_lvlh_normal'][kA0] + L21_A_dict['att_lvlh_normal'][kA1] + \
               L21_B_dict['att_lvlh_normal'][kB0] + L21_B_dict['att_lvlh_normal'][kB1] not in [0,4]:
                qflags[i,k,29] = 1
                continue


            ######################## Interpolating ############################
            # If it passed all the error checks, perform bilinear interpolation (altitude, then longitude).
            # Variables to interpolate to this point:
            #   - los_wind (A and B)
            #   - az       (A and B)
            #   - lat      (A and B, to be averaged)
            #   - time     (A and B, to be averaged and subtracted)
            #   - ver      (A and B, to be compared)
            #   - quality factors (A and B, to be used to compute L2.2 quality factors)
            #   - qflags   (A and B, to be passed through)

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
                    # Do interpolation of value to the desired altitude, for each longitude
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

            # Interpolate A&B los_wind, errors, temperature, lat, az, VER, etc. to this grid point
            los_wind_A_pt, los_wind_A_pt_err = \
                            bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], los_wind_A[:,kA0:kA1+1],
                                            prop_err = True, valerr = los_wind_A_err[:,kA0:kA1+1])
            los_wind_B_pt, los_wind_B_pt_err = \
                            bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], los_wind_B[:,kB0:kB1+1],
                                            prop_err = True, valerr = los_wind_B_err[:,kB0:kB1+1])
            _, los_wind_A_pt_precision_1_sample = \
                            bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], los_wind_A[:,kA0:kA1+1],
                                            prop_err = True, valerr = los_wind_A_precision_1_sample[:,kA0:kA1+1])
            _, los_wind_B_pt_precision_1_sample = \
                            bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], los_wind_B[:,kB0:kB1+1],
                                            prop_err = True, valerr = los_wind_B_precision_1_sample[:,kB0:kB1+1])
            _, los_wind_A_pt_precision_1_day = \
                            bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], los_wind_A[:,kA0:kA1+1],
                                            prop_err = True, valerr = los_wind_A_precision_1_day[:,kA0:kA1+1])
            _, los_wind_B_pt_precision_1_day = \
                            bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], los_wind_B[:,kB0:kB1+1],
                                            prop_err = True, valerr = los_wind_B_precision_1_day[:,kB0:kB1+1])
            _, los_wind_A_pt_precision_low_signal = \
                            bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], los_wind_A[:,kA0:kA1+1],
                                            prop_err = True, valerr = los_wind_A_precision_low_signal[:,kA0:kA1+1])
            _, los_wind_B_pt_precision_low_signal = \
                            bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], los_wind_B[:,kB0:kB1+1],
                                            prop_err = True, valerr = los_wind_B_precision_low_signal[:,kB0:kB1+1])
            _, los_wind_A_pt_accuracy = \
                            bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], los_wind_A[:,kA0:kA1+1],
                                            prop_err = True, valerr = los_wind_A_accuracy[:,kA0:kA1+1])
            _, los_wind_B_pt_accuracy = \
                            bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], los_wind_B[:,kB0:kB1+1],
                                            prop_err = True, valerr = los_wind_B_accuracy[:,kB0:kB1+1])
            ver_A_pt, ver_A_pt_err = \
                            bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], ver_A[:,kA0:kA1+1],
                                            prop_err = True, valerr = ver_A_err[:,kA0:kA1+1])
            ver_B_pt, ver_B_pt_err = \
                            bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], ver_B[:,kB0:kB1+1],
                                            prop_err = True, valerr = ver_B_err[:,kB0:kB1+1])
            amp_A_pt, amp_A_pt_err = \
                            bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], amp_A[:,kA0:kA1+1],
                                            prop_err = True, valerr = amp_A_err[:,kA0:kA1+1])
            amp_B_pt, amp_B_pt_err = \
                            bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], amp_B[:,kB0:kB1+1],
                                            prop_err = True, valerr = amp_B_err[:,kB0:kB1+1])
            local_az_A_pt = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], local_az_A[:,kA0:kA1+1])
            local_az_B_pt = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], local_az_B[:,kB0:kB1+1])
            lat_A_pt      = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], lat_A     [:,kA0:kA1+1])
            lat_B_pt      = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], lat_B     [:,kB0:kB1+1])
            mlat_A_pt     = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], mlat_A[:,kA0:kA1+1])
            mlat_B_pt     = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], mlat_B[:,kB0:kB1+1])
            mlon_A_pt     = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], mlon_A[:,kA0:kA1+1])
            mlon_B_pt     = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], mlon_B[:,kB0:kB1+1])
            slt_A_pt      = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], slt_A[:,kA0:kA1+1])
            slt_B_pt      = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], slt_B[:,kB0:kB1+1])
            sza_A_pt      = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], sza_A[:,kA0:kA1+1])
            sza_B_pt      = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], sza_B[:,kB0:kB1+1])
            # Record which files were used
            nA = np.tile(np.arange(N_alts_A),(2,1)).T.astype(float)
            nB = np.tile(np.arange(N_alts_B),(2,1)).T.astype(float)
            idx_A = bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], nA)
            idx_B = bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], nB)
            idx_A_0 = int(np.floor(idx_A))
            idx_B_0 = int(np.floor(idx_B))
            idx_A_1 = int(np.ceil(idx_A))
            idx_B_1 = int(np.ceil(idx_B))
            N_used_A[idx_A_0:idx_A_1+1, kA0:kA1+1] += 1
            N_used_B[idx_B_0:idx_B_1+1, kB0:kB1+1] += 1
            # "Interpolate" quality factor and flags from A&B to this point. This is difficult. One wants to be conservative
            # and say that "if any point used in the interpolation is bad, then the result is bad." However that means that
            # an entire row of data (alt = 95 km for instance) could be marked "bad" just because the MTB sample at 92.1 km contributed
            # a little bit. So we'll use a compromise strategy: Use the interpolated quality factor and round to nearest 0.5.
            # Use same idea for quality flags.
            wind_quality_A_pt = 0.5 * round(2.*bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], wind_quality_A[:,kA0:kA1+1]))
            wind_quality_B_pt = 0.5 * round(2.*bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], wind_quality_B[:,kB0:kB1+1]))
            ver_quality_A_pt  = 0.5 * round(2.*bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], ver_quality_A[:,kA0:kA1+1]))
            ver_quality_B_pt  = 0.5 * round(2.*bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], ver_quality_B[:,kB0:kB1+1]))
            N_flags_L21 = np.shape(qflags_A)[2] # Number of flags at L2.1 (assume same for A and B)
            qflags_A_pt = np.zeros(N_flags_L21) # A 1D array of flags from A, "interpolated" to this grid point
            qflags_B_pt = np.zeros(N_flags_L21) # A 1D array of flags from B, "interpolated" to this grid point
            for nflag in range(N_flags_L21):
                # TODO: I'm worried that there will be a disconnect between how this flag ends up, and how the final quality factor ends up.
                # Double-check that it's working ok.
                qflags_A_pt[nflag] = round(bilinear_interp(lon_A[:,kA0:kA1+1], alt_A[:,kA0:kA1+1], qflags_A[:,kA0:kA1+1,nflag]))
                qflags_B_pt[nflag] = round(bilinear_interp(lon_B[:,kB0:kB1+1], alt_B[:,kB0:kB1+1], qflags_B[:,kB0:kB1+1,nflag]))
            
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
            # Interpolate orbit_number and orbit_descending flag
            orb_num_A_pt  = interpolate_linear([lon_A_0, lon_A_1], orb_num_A[kA0:kA1+1],  lon_pt)
            orb_num_B_pt  = interpolate_linear([lon_B_0, lon_B_1], orb_num_B[kB0:kB1+1],  lon_pt)
            orb_desc_A_pt = interpolate_linear([lon_A_0, lon_A_1], orb_desc_A[kA0:kA1+1], lon_pt)
            orb_desc_B_pt = interpolate_linear([lon_B_0, lon_B_1], orb_desc_B[kB0:kB1+1], lon_pt)
            
            # Calculate time difference and verify it is small. If it's large,
            # this implies something went wrong with the A/B matchup (i.e., longitude unwrapping).
            time_delta_pt = (t_B - t_A).total_seconds()
            time_pt = t_A + timedelta(seconds=time_delta_pt/2.)
            assert abs(time_delta_pt) < t_diff_AB, \
                   "A/B match-up using files %.0f sec apart. Max allowed is %.0f sec" % (time_delta_pt, t_diff_AB)
            
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
            # propagate uncertainties, separately for each kind of error
            # Precision 1 sample
            Sig_y = np.array([[los_wind_A_pt_precision_1_sample**2, 0.0],
                              [0.0, los_wind_B_pt_precision_1_sample**2]]) # covariance matrix of y
            Sig_x = invA.dot(Sig_y.dot(invA.T)) # standard linear error propagation
            u_precision_1_sample = np.sqrt(Sig_x[0,0])
            v_precision_1_sample = np.sqrt(Sig_x[1,1])
            # Precision 1 day
            Sig_y = np.array([[los_wind_A_pt_precision_1_day**2, 0.0],
                              [0.0, los_wind_B_pt_precision_1_day**2]]) # covariance matrix of y
            Sig_x = invA.dot(Sig_y.dot(invA.T)) # standard linear error propagation
            u_precision_1_day = np.sqrt(Sig_x[0,0])
            v_precision_1_day = np.sqrt(Sig_x[1,1])
            # Precision low signal
            Sig_y = np.array([[los_wind_A_pt_precision_low_signal**2, 0.0],
                              [0.0, los_wind_B_pt_precision_low_signal**2]]) # covariance matrix of y
            Sig_x = invA.dot(Sig_y.dot(invA.T)) # standard linear error propagation
            u_precision_low_signal = np.sqrt(Sig_x[0,0])
            v_precision_low_signal = np.sqrt(Sig_x[1,1])
            # Accuracy
            Sig_y = np.array([[los_wind_A_pt_accuracy**2, 0.0],
                              [0.0, los_wind_B_pt_accuracy**2]]) # covariance matrix of y
            Sig_x = invA.dot(Sig_y.dot(invA.T)) # standard linear error propagation
            u_accuracy = np.sqrt(Sig_x[0,0])
            v_accuracy = np.sqrt(Sig_x[1,1])
            # Total
            Sig_y = np.array([[los_wind_A_pt_err**2, 0.0],
                              [0.0, los_wind_B_pt_err**2]]) # covariance matrix of y
            Sig_x = invA.dot(Sig_y.dot(invA.T)) # standard linear error propagation
            u_err = np.sqrt(Sig_x[0,0])
            v_err = np.sqrt(Sig_x[1,1])
            # Now "invert" VER (simply average the two, and propagate uncertainty)
            ver_pt = (ver_A_pt + ver_B_pt)/2.
            ver_err_pt = np.sqrt(ver_A_pt_err**2 + ver_B_pt_err**2)/2.
            ver_rel_diff_pt = abs(ver_A_pt - ver_B_pt)/ver_pt # this is used for quality flagging
            # And same for fringe amplitude
            amp_pt = (amp_A_pt + amp_B_pt)/2.
            amp_err_pt = np.sqrt(amp_A_pt_err**2 + amp_B_pt_err**2)/2.
            amp_rel_diff_pt = abs(amp_A_pt - amp_B_pt)/amp_pt # this is used for quality flagging

            # Fill in all the relevant variables at this grid point
            U[i,k] = u
            V[i,k] = v
            U_precision_1_sample[i,k] = u_precision_1_sample
            V_precision_1_sample[i,k] = v_precision_1_sample
            U_precision_1_day[i,k] = u_precision_1_day
            V_precision_1_day[i,k] = v_precision_1_day
            U_precision_low_signal[i,k] = u_precision_low_signal
            V_precision_low_signal[i,k] = v_precision_low_signal
            U_accuracy[i,k] = u_accuracy
            V_accuracy[i,k] = v_accuracy
            U_err[i,k] = u_err
            V_err[i,k] = v_err
            lat[i,k] = (lat_A_pt + lat_B_pt)/2.
            time[i,k] = time_pt
            time_delta[i,k] = time_delta_pt
            ver_intp_A[i,k] = ver_A_pt
            ver_intp_B[i,k] = ver_B_pt
            ver[i,k]   = ver_pt
            ver_err[i,k] = ver_err_pt
            ver_rel_diff[i,k] = ver_rel_diff_pt    
            amp_intp_A[i,k] = amp_A_pt
            amp_intp_B[i,k] = amp_B_pt
            amp[i,k]   = amp_pt
            amp_err[i,k] = amp_err_pt
            amp_rel_diff[i,k] = amp_rel_diff_pt
            mlat[i,k] = (mlat_A_pt + mlat_B_pt)/2.
            mlon[i,k] = circular_mean(mlon_A_pt, mlon_B_pt) # even though A and B are individually unwrapped, circular mean is needed 
                                                            # since A and B might be one or more cycles apart.
            slt[i,k]  = circular_mean(15.*slt_A_pt, 15.*slt_B_pt)/15.   # same as above
            sza[i,k]  = (sza_A_pt + sza_B_pt)/2.
            orb_num[i,k] = (orb_num_A_pt + orb_num_B_pt)/2. # is this the best way to do it?
            orb_desc[i,k] = (orb_desc_A_pt + orb_desc_B_pt)/2. # is this the best way to do it?
            
            ###################### More quality flagging #######################            
            # Pass through all L2.1 flags
            qflags[i,k,:N_flags_L21] = qflags_A_pt
            qflags[i,k,N_flags_L21:2*N_flags_L21] = qflags_B_pt

            # Check spherical symmetry
            if ver_rel_diff_pt > sph_asym_thresh: # Note this is comparing VER not fringe amplitude
                qflags[i,k,28] = 1

            # Unknown error - this should never happen
            if (np.isnan(u) or np.isnan(v)) and all(qflags[i,k,:] == 0): # Unknown error
                qflags[i,k,-1] = 1
                
            #################### Quality factor calculation ###################
            # Construct a list of ratings, and take the min
            wind_ratings = [1.0] # start with 1.0. If there are no dings, the quality factor is 1.0
            ver_ratings = [1.0]
            # Note that this L2.2 quality code only considers the newly-introduced L2.2 quality flags. It is up
            # to the L1 and L2.1 code to adjust the quality factor based on its flags. Of course, this L2.2 quality 
            # code considers the L2.1 and L1 (via the L2.1) quality factors as "ceilings" for its factor.
            # Missing MIGHTI-A or MIGHTI-B files: automatic 0
            if any(qflags[i,k,[24,25]]):
                wind_ratings.append(0.0)
                ver_ratings.append(0.0)
            # Didn't sample this altitude: automatic 0
            if any(qflags[i,k,[26,27]]):
                wind_ratings.append(0.0)
                ver_ratings.append(0.0)
            # Spherical asymmetry:
            if qflags[i,k,28]:
                # BJH 2020 Feb 3: temporarily disabled this flag until we fix VER estimates
                # BJH 2022 Jun 7: re-enabled after A/B cross calibration in v05
                wind_ratings.append(0.5)
                ver_ratings.append(0.5)
            # Quality flag automatically set to 0 if it's nan
            if np.isnan(u) or np.isnan(v):
                wind_ratings.append(0.0)
            if np.isnan(ver_pt):
                ver_ratings.append(0.0)
            
            # Quality flag for time_start/time_stop will be handled below
            # Add L2.1 quality factors and determine final rating.
            wind_ratings.append(wind_quality_A_pt)
            wind_ratings.append(wind_quality_B_pt)
            ver_ratings.append(ver_quality_A_pt)
            ver_ratings.append(ver_quality_B_pt)
            wind_quality[i,k] = min(wind_ratings)
            ver_quality[i,k] = min(ver_ratings)
            
    ####################### Magnetic coordinates ######################
    # This section comprises about 25% of the runtime of the L2.2 processing
    # (which is not a problem - it's already pretty fast)
    wind_mag_fa  = np.nan*np.zeros((N_alts, N_lons))
    wind_mag_mer = np.nan*np.zeros((N_alts, N_lons))
    wind_mag_zon = np.nan*np.zeros((N_alts, N_lons))
    for k in range(N_lons):
        # Assume vertical wind is zero.
        # Use time_start as IGRF time reference.
        # Use fast formulation to save time (accurate to within 0.08% and much faster)
        wfa, wmer, wzon = geog_to_mag_coords(U[:,k], V[:,k], np.zeros(N_alts), lat[:,k], lon[:,k], 
                                             alt, [time_start]*N_alts, fast=True)
        wind_mag_fa[:,k] = wfa
        wind_mag_mer[:,k] = wmer
        wind_mag_zon[:,k] = wzon

    ############### Trimming data points at beginning and end ##################### 
    # Determine the "average" time per column to be reported in file.
    epoch_sec = np.nan * np.zeros(N_lons) # referenced to start time
    for k in range(N_lons):
        dt_k = [(ti - time_start).total_seconds() for ti in time[:,k] if ti is not None]
        if len(dt_k) > 0:
            epoch_sec[k] = np.mean(dt_k)

    # Identify start and stop columns. Don't keep times outside time start/stop, and 
    # trim empty samples at the beginning and end.
    # Find the first index after time_start. This is done manually because other built-in
    # options (bisect.bisect, np.searchsorted, etc.) do not play nicely with nans.
    k0 = 0
    for k in range(len(epoch_sec)):
        if np.isfinite(epoch_sec[k]) and epoch_sec[k] >= 0.0:
            k0 = k
            break    # Find the last index before time_stop (if it exists)
    sec_stop = (time_stop - time_start).total_seconds()
    k_save = -1
    for k in range(k0,len(epoch_sec)):
        if np.isfinite(epoch_sec[k]) and epoch_sec[k] < sec_stop:
            # The last time this test passes is the k1 we want (minus 1)
            k_save = k
    k1 = k_save + 1

    # Fill in gaps. This is not technically necessary but allows us to adhere to 
    # a netCDF "best practice" of coordinate variables being non-masked.
    x = np.arange(N_lons)
    good = np.isfinite(epoch_sec)
    ep_int = interpolate.InterpolatedUnivariateSpline(x[good],epoch_sec[good], ext=0, k=1)
    epoch_sec[~good] = ep_int(x[~good])

    # Convert to datetime
    epoch = np.array([time_start + timedelta(seconds=s) for s in epoch_sec])
    
    # Make sure that all data with quality 0 are masked.
    # This might happen on the edges of good regions, where the algorithm is
    # interpolating between good and bad points.
    i = wind_quality == 0.0
    for v in [U, V, wind_mag_fa, wind_mag_mer, wind_mag_zon, 
              U_err, V_err, 
              U_precision_1_sample, U_precision_1_day, U_precision_low_signal, U_accuracy,
              V_precision_1_sample, V_precision_1_day, V_precision_low_signal, V_accuracy,]:
        v[i] = np.nan
    i = ver_quality == 0.0
    for v in [ver, amp, ver_err, amp_err]:
        v[i] = np.nan    

    ############## Output section ##################
    # Note that the before/after trimming is done here.
    L22_dict = {}
    L22_dict['lat']                 = lat  [:,k0:k1]
    L22_dict['lon']                 = np.mod(lon[:,k0:k1], 360.)
    L22_dict['lon_unwrapped']       = lon  [:,k0:k1]
    L22_dict['alt']                 = alt
    L22_dict['u']                   = U    [:,k0:k1]
    L22_dict['v']                   = V    [:,k0:k1]
    L22_dict['u_precision_1_sample']   = U_precision_1_sample[:,k0:k1]
    L22_dict['v_precision_1_sample']   = V_precision_1_sample[:,k0:k1]
    L22_dict['u_precision_1_day']      = U_precision_1_day[:,k0:k1]
    L22_dict['v_precision_1_day']      = V_precision_1_day[:,k0:k1]
    L22_dict['u_precision_low_signal'] = U_precision_low_signal[:,k0:k1]
    L22_dict['v_precision_low_signal'] = V_precision_low_signal[:,k0:k1]
    L22_dict['u_accuracy']             = U_accuracy[:,k0:k1]
    L22_dict['v_accuracy']             = V_accuracy[:,k0:k1]
    L22_dict['u_error']             = U_err[:,k0:k1]
    L22_dict['v_error']             = V_err[:,k0:k1]
    L22_dict['quality_flags']       = qflags[:,k0:k1,:]
    L22_dict['epoch_full']          = time [:,k0:k1] # 2D
    L22_dict['epoch']               = epoch[k0:k1] # 1D
    L22_dict['time_delta']          = time_delta[:,k0:k1]
    L22_dict['time_start']          = time_start
    L22_dict['time_stop']           = time_stop
    L22_dict['fringe_amp']          = amp[:,k0:k1]
    L22_dict['fringe_amp_error']    = amp_err[:,k0:k1]
    L22_dict['fringe_amp_A']        = amp_intp_A[:,k0:k1]
    L22_dict['fringe_amp_B']        = amp_intp_B[:,k0:k1]
    L22_dict['fringe_amp_rel_diff'] = amp_rel_diff[:,k0:k1]
    L22_dict['ver']                 = ver[:,k0:k1]
    L22_dict['ver_error']           = ver_err[:,k0:k1]
    L22_dict['ver_A']               = ver_intp_A[:,k0:k1]
    L22_dict['ver_B']               = ver_intp_B[:,k0:k1]
    L22_dict['ver_rel_diff']        = ver_rel_diff[:,k0:k1]
    L22_dict['emission_color']      = emission_color
    L22_dict['source_files']        = np.concatenate((L21_A_dict['source_files'], L21_B_dict['source_files']))
    L22_dict['acknowledgement']     = L21_A_dict['acknowledgement']
    L22_dict['wind_quality']        = wind_quality[:,k0:k1]
    L22_dict['ver_quality']         = ver_quality [:,k0:k1]
    L22_dict['N_used_A']            = N_used_A
    L22_dict['N_used_B']            = N_used_B
    L22_dict['wind_mag_fa']         = wind_mag_fa[:,k0:k1]
    L22_dict['wind_mag_mer']        = wind_mag_mer[:,k0:k1]
    L22_dict['wind_mag_zon']        = wind_mag_zon[:,k0:k1]
    L22_dict['mag_lat']             = mlat[:,k0:k1]
    L22_dict['mag_lon']             = mlon[:,k0:k1]
    L22_dict['sza']                 = sza[:,k0:k1]
    L22_dict['slt']                 = slt[:,k0:k1]
    L22_dict['orbit_number']        = orb_num[:,k0:k1]
    L22_dict['orbit_node']          = orb_desc[:,k0:k1]
    L22_dict['zero_wind_ref']       = L21_A_dict['zero_wind_ref']
    L22_dict['corr_notch_drift']    = L21_A_dict['corr_notch_drift']
    
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
    L22_fn = 'ICON_L2-2_MIGHTI_Vector-Wind-%s_%s_v%02ir%03i.NC' % (L22_dict['emission_color'].capitalize(),
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
        ncfile.setncattr_string('Date_End',                       t_stop.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Date_Start',                     t_start.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Description',                    'ICON MIGHTI Cardinal Vector Winds (DP 2.2)')
        ncfile.setncattr_string('Descriptor',                     'MIGHTI > Michelson Interferometer for Global High-resolution Thermospheric Imaging')
        ncfile.setncattr_string('Discipline',                     'Space Physics > Ionospheric Science')
        ncfile.setncattr_string('File',                           L22_fn)
        ncfile.setncattr_string('File_Date',                      t_file.strftime('%a, %d %b %Y, %Y-%m-%dT%H:%M:%S.%f')[:-3] + ' UTC')
        ncfile.setncattr_string('Generated_By',                   'ICON SDC > ICON UIUC MIGHTI L2.2 Processor v%s, B. J. Harding' % __version__)
        ncfile.setncattr_string('Generation_Date',                t_file.strftime('%Y%m%d'))
        ncfile.setncattr_string('History',                        ['v1.0: First release of MIGHTI L2.1/L2.2 software, B. J. Harding, 05 Mar 2018', 
                                                                   'v2.0: First run of on-orbit data, using external zero wind reference and smooth daily-averaged profiles, B. J. Harding, 01 May 2020',
                                                                   'v3.0: Correction for long-term mechanical drift, B. J. Harding, 04 Jun 2020',
                                                                   'v4.0: Updated correction for long-term mechanical drift to handle settling after ~May 2020 and precession cycle variation. LoS winds have changed by a bulk offset of up to 30 m/s. Studies using only perturbations from the mean (e.g., non-migrating tidal retrievals) are unlikely to be affected. B. J. Harding 21 Oct 2020',
                                                                   'v5.0: The ad-hoc, HWM-based correction for the zero wind phase has been replaced with a self-calibration based on comparing data from the ascending and descending orbits (see the notes for the wind variables below for details). Long-term trends in the zero wind phase degraded the accuracy of version 4 over time, and the accuracy of version 4 data was tied to the accuracy of HWM. In version 5, errors on these long time scales (>100-150 days) are now accounted for, improving the accuracy to 10-25 m/s (see the "Accuracy" variable for more details) and removing the dependence on external models. A long window of data is required to implement this self-calibration, so v05 data are processed at least 100 days behind real time. For errors on precession-cycle time scales (48 days), the previous ad-hoc correction using initial red-vs-green comparisons has been replaced with a more comprehensive red-green cross-calibration (165-185 km altitude during the day) that accounts for the time-dependence of mechanical drifts of the optics. This result is consistent with a first-principles analysis of the fiducial notch positions (see Marr et al., 2020 and subsequent publications). A mission-average fiducial notch analysis is also used to correct mechanical drifts on an orbital time scale (97 minutes, or 24 hours of local time), which could affect migrating tide estimates. The RMS difference due to this effect is estimated at 5-10 m/s (root mean square). Analysis of waves with periods that do not coincide with these new corrections are not likely to be different than in version 4 (e.g., nonmigrating tides, planetary waves, and gravity waves). New variables related to error (i.e., uncertainty) estimates from various sources are now included, whereas version 4 error estimates only included the effects of shot, read, and dark noise. MIGHTI-A and MIGHTI-B variables related to emission brightness are now cross-calibrated, though not absolutely calibrated. Exposures affected by solar and lunar stray light are now flagged. Some data during periods in May and July 2020 when the sun approached the MIGHTI field of view was marked as unavailable in v04, but is now available in v05 with the exception of a few days. The data from the second row in the green channel (~91 km) is now available when the signal strength permits a wind estimate. An error in the local time calculation has been corrected, which changes the local time by up to 20 min. A new algorithm to identify cosmic ray spikes has been implemented, improving precision. A preliminary algorithm has been implemented to to correct a wind bias associated with low signal levels, and associated uncertainties are estimated (see the "Precision_Low_Signal_Effect" variable for more details). Finally, various quality control parameters have been optimized. More description is provided in the notes below. A full history of software changes can be found on Github: https://github.com/bharding512/airglowrsss/commits/master/Python/modules/MIGHTI_L2.py  B. J. Harding 08 Sep 2022'
                                                                  ])
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
        ncfile.setncattr_string('Text_Supplement',                [
        'This data product contains cardinal (i.e., zonal and meridional) thermospheric winds ' +
        'obtained by combining Level 2.1 (line-of-sight winds) from MIGHTI A and MIGHTI B. The cardinal winds are given as a function of time ' +
        '(spanning 24 hours) and altitude (spanning nominally 90-300 km). In addition to the cardinal vector wind data and the corresponding ancillary data, ' +
        'such as time and location, this product contains supporting data, such as fringe amplitude profiles and relative volume emission rate profiles. ' +
        'Absolute calibration and MIGHTI-A/B cross calibration of these data is not necessary to obtain the wind data, and therefore any direct analysis ' +
        'of these supporting data requires caution. There is one file per emission color (red or green).',
         
        'Cardinal wind observations are enabled by the ~90-degree offset between the two MIGHTI sensors. First, MIGHTI A measures a wind component along '+
        'its line of sight. Five to eight minutes later, depending on tangent point altitude, the spacecraft has moved to a position such that MIGHTI B '+
        'measures a nearly orthogonal wind component at approximately the same location. A coordinate rotation is performed on the two line-of-sight '+
        'components to obtain the northward and eastward components reported in this file. The assumption is that the thermospheric wind has not changed '+
        'during this time interval. Because the Level 2.1 data are naturally on an irregular grid, '+
        'they are first interpolated to a regular, pre-defined grid of longitude and altitude before the coordinate rotation is performed. See Harding et al. [2017, '+
        'doi:10.1007/s11214-017-0359-3] for more details of the Level 2.2 algorithm. ' +
        'Further discussion of the calibration and performance of MIGHTI after launch can be found in a forthcoming paper in Space Science ' +
        'Reviews [Englert et al., 2022, in preparation].',

             "Known issues with the v05 data release are listed below. Work is in progress to resolve or mitigate these issues in future data releases. ",

                     "Known issues with v05:<br/>"
                     " * When ICON is in the South Atlantic Anomaly (SAA), radiation effects on the detector cause poor data quality. The quality control algorithm "
                     "adequately flags and masks most of the affected samples, but some outliers remain, especially near the edge of the SAA. Other uncaught outliers "
                     "are rare but can occur due to cosmic rays, stars in the field of view, moonlight, etc. <br/>"
                     " * The bottom row of data (corresponding to an altitude of ~88 km) is masked out as the signal is rarely strong enough to permit a wind "
                     "estimate, and calibrations have large uncertainties. It is unlikely but possible that this altitude will be reported in future releases, "
                     "pending further investigation. <br/>"
                     " * Airglow brightness observations are not a required mission product, and no effort was yet made to absolutely calibrate "
                     "the brightness observations for MIGHTI-A and MIGHTI-B, and thus the Relative_VER "
                     "variable should be treated with caution. In v05, MIGHTI-A and B are cross-calibrated using a conversion factor derived from on-orbit data. "
                     "However, there are some indications that this cross-calibration may be changing with time, which is not accounted for in v05. <br/>"
                     
                     " * As discussed in the variable notes below, a new zero wind phase determination has been implemented in v05. However, during the period "
                     "from 2021 Apr 26 to Aug 14, data gaps and one period of southward (\"Reverse LVLH\") pointing cause errors in this determination. The "
                     "accuracy is estimated to be degraded by a factor of two. See the *_Accuracy variable. <br/>"
                     " * During the one orbit per day when the calibration lamp is on, the wind data can be noisier and have a slight bias. Although this issue is "
                     "much improved since v04, for the sake of conservatism, these orbits are still labeled with quality=0.5 (i.e., caution). <br/>"
                     " * Some data gaps appear on days when the sun passes near MIGHTI's field of view. Most of these gaps are located near the terminator, "
                     "but some are longer lasting.<br/>"
                     " * In some cases, there are indications that the *_Precision_1_Sample variables are underestimating the true sample-to-sample noise, "
                     "suggesting that, in addition to shot noise, there is a second source of noise. It is recommend that any quantitative use of the reported "
                     "error estimates (i.e., precision and accuracy) should treat those estimates as uncertain. It is believed that most error estimates are "
                     "correct to within a factor of 2. The largest problems with error reporting occur where the airglow signal is weakest. <br/>"
                     
                     " * Imperfect daily calibration data lead to small discontinuities in the zero wind phase at the boundaries between days (i.e., "
                     " between 23:59:59 and 00:00:00 UT), which are not accounted for by the reported error variables. This was estimated to be a 2-5 "
                     "m/s (root-mean-square) error early in the mission, but is growing over time, possibly reaching 5-10 m/s by mid-2022.<br/>"
                     " * A signal-dependent phase shift is seen in atmospheric and calibration lamp fringes, possibly caused by "
                     "a charge trapping effect in the CCD. This is the subject of ongoing investigation, but a first-order correction "
                     "is implemented in the v05 dataset. The correction increases linearly with time to match the effect seen in on-orbit calibration data. "
                     "The variable *_Precision_Low_Signal_Effect is an estimate of the remaining uncertainty due to an imperfect "
                     "correction. Where this uncertainty is large, caution is recommended. "
                     "For example, for winds in the core science region (90-105 km altitude, away from the terminator), the magnitude of the correction "
                     "is small or zero, but data in the red "
                     "channel during the night and twilight are subject to a large correction (many tens of m/s) and the uncertainty is correspondingly large. "
                     "A goal for future releases is to characterize and correct this effect more accurately. <br/>"
                     " * Data near the solar terminators are subject to a variety of errors, including those described above and others related to the rapidly "
                     "varying illumination. Not all errors near the terminator are accounted for by the reported error. Users are encouraged to "
                     "use extra caution with these data.<br/>"
                     "See the documentation below for more information.<br/>"
                    
        ])
        ncfile.setncattr_string('Time_Resolution',                '30 - 60 seconds')
        ncfile.setncattr_string('Title',                          'ICON MIGHTI Cardinal Vector Winds (DP 2.2)')
        ncfile.setncattr_string('Zero_Wind_Ref',                  'DEPRECATED')
        ncfile.setncattr_string('Corr_Notch_Drift',               'DEPRECATED')
        ncfile.setncattr_string('Low_Signal_Patch',               'Automatic') # Indicating this file was created after the low-signal bug was identified/fixed, but was not manually patched.



        ################################## Dimensions ########################################
        prefix = 'ICON_L22' # variable prefix
        var_alt_name = '%s_Altitude'%prefix
        
        ny,nx,nflags = np.shape(L22_dict['quality_flags'])
        ncfile.createDimension('Epoch',None)
        ncfile.createDimension(var_alt_name, ny)
        ncfile.createDimension('N_Flags', nflags)
        

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
                              "MIGHTI-B sample times, which differ by 5-8 minutes. Where MIGHTI-A or MIGHTI-B samples are missing, data are "
                              "reported as missing, but gaps in Epoch are interpolated over to adhere to the netCDF4 standard that coordinate "
                              "variables should have no missing values. "
                              "The matchup between MIGHTI-A and B happens at slightly different "
                              "times at different altitudes, a complication which is ignored by this variable. The effect is small (plus or minus 30-60 "
                              "seconds), but in cases where it is important, it is recommended to use the alternative time variable "
                              "Epoch_Full, which is two dimensional and captures the variation with altitude. "
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
                              "the difference between the MIGHTI-A and MIGHTI-B times that contributed to each point. Epoch_Full "
                              "contains the average time."
                              )
        
        # Human readable versions of those two variables
        utctime = []
        for t in L22_dict['epoch']:
            if t is None:
                utctime.append('')
            else:
                utctime.append(t.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3] + 'Z') # YYYY-MM-DD HH:mm:ss.xxxZ as per DP Conventions 0.7
        utctime = np.array(utctime)
        var = _create_variable(ncfile, '%s_UTC_Time'%prefix, utctime, 
                              dimensions=('Epoch'), depend_0 = 'Epoch',
                              format_nc=str, format_fortran='A', desc='Sample time, average of A and B measurements.', 
                              display_type='time_series', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
                              units='', valid_min=None, valid_max=None, var_type='support_data', chunk_sizes=[nx],
                              notes="This variable is the same as Epoch but is formatted as a human-readable string. Missing "
                              "grid points are labeled with empty strings.")       
        
#         utctime_full = np.empty((ny,nx), dtype='S23')
#         for i in range(ny):
#             for j in range(nx):
#                 t = L22_dict['epoch_full'][i,j]
#                 if t is None:
#                     utctime_full[i,j] = ''
#                 else:
#                     utctime_full[i,j] = t.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
#         #utctime_full = np.array(utctime_full)
#         var = _create_variable(ncfile, '%s_UTC_Time_Full'%prefix, utctime_full.T, 
#                               dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
#                               format_nc=str, format_fortran='A', desc='Sample time, average of A and B measurements.', 
#                               display_type='image', field_name='Time', fill_value=None, label_axis='Time', bin_location=0.5,
#                               units='', valid_min=None, valid_max=None, var_type='support_data', chunk_sizes=[nx,ny],
#                               notes="This variable is the same as Epoch_Full but is formatted as a human-readable string. Missing "
#                               "grid points are labeled with empty strings.")       

        ######### Data Variables #########
        # var notes for zonal and meridional wind are identical, and have a special case if zero_wind_ref == 'external'
        varnotes = ["The zonal (positive eastward) and meridional (positive northward) winds are the primary "
                  "data product in this file. They are defined on a grid with dimensions of time and altitude, "
                  "spanning 24 hours and nominally 90-300 km (150-300 km for the red channel). The altitude, time, "
                  "latitude and longitude corresponding to each point in the grid are given by other variables in "
                  "this file. It should be noted that while each measurement is ascribed to a particular latitude, "
                  "longitude, altitude, and time, it is actually an average over many hundreds of kilometers "
                  "horizontally and 2.5-30 kilometers vertically (depending on the binning). It also assumes stationarity "
                  "over the 5-8 minutes between the MIGHTI-A and B measurements used for each point. See Harding et "
                  "al. [2017, doi:10.1007/s11214-017-0359-3] for a more complete discussion of the inversion algorithm."]
        varnotes.extend(["Knowledge of the \"zero wind phase\" is needed for any instrument using Doppler shifts to determine winds. "
                         "The zero wind phase is defined as the measured interference fringe phase that corresponds to the rest "
                         "wavelength of the emission. For the v05 data release, the zero wind phase has been determined by considering a window of LoS wind data "
                         "spanning two precession cycles (96 days). Assuming that on average the real zonal and meridional winds do not depend on the aspect "
                         "angle with which MIGHTI observes the atmosphere (an angle which is significantly different between the ascending and descending portions of the orbit), "
                         "a matrix equation can be constructed which combines data from both MIGHTI-A and "
                         "MIGHTI-B and both the ascending and descending orbits. This equation is solved for the average zonal and meridional wind, and the "
                         "zero wind phase for MIGHTI-A and MIGHTI-B. This window is moved in time to determine the appropriate zero wind phase for each date. "
                         "The value of the zero wind phase depends on emission color (red or green), aperture mode (day or night), calibration lamp status "
                         "(on or off) and row (i.e., altitude). "
                         "An additional zero-mean signal is added to this result to ensure that 48-day (i.e., 1 precession cycle) average winds are smooth "
                         "in altitude. Adjustments are smaller than the reported accuracy, so this adjustment is not "
                         "expected to change any scientific conclusions, although it does ensure more realistic wind profiles. "
                         "This is a less restrictive assumption than the smoothness criterion used in v04, which "
                         "relied on the Horizontal Wind Model 2014 and also enforced smoothness on 1-day averages. It is thus expected that the amplitude of tidal structures "
                         "in the lower thermosphere are subject to less suppression in v05+ than in v04. "
#                          "The results of this approach are consistent with the results of ICON's \"zero wind maneuver.\" " # 2022 Sep 9: We don't know if this is true yet.
                         "This version of the MIGHTI zero wind phase is independent of any external data or models (such as "
                         "the Horizontal Wind Model 2014, which was used in v04 and earlier versions). "
                         "The zero wind phase used for each wind sample is saved in the _Zero_Wind_Phase variable below. "
                         "The 1-sigma uncertainty in the winds incurred by the inaccuracy in the zero wind phase "
                         "is estimated and reported in the *_Wind_Accuracy variable below."])
        
        # Zonal Wind
        var = _create_variable(ncfile, '%s_Zonal_Wind'%prefix, L22_dict['u'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Zonal component of the horizontal wind. Positive Eastward.', 
                              display_type='image', field_name='Zonal Wind', fill_value=None, label_axis='Zonal Wind', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=varnotes)

        # Meridional Wind
        var = _create_variable(ncfile, '%s_Meridional_Wind'%prefix, L22_dict['v'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Meridional component of the horizontal wind. Positive Northward.', 
                              display_type='image', field_name='Meridional Wind', fill_value=None, label_axis='Merid Wind', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=varnotes)    
        
        # Wind Errors
        
        err_prefix = ["Various sources of error in MIGHTI winds are quantified with 1-sigma estimates and organized by their temporal persistence. " 
                      "These error sources are nearly uncorrelated with each other and can thus be added in quadrature. Users are encouraged to "
                      "contact the MIGHTI team for assistance with error propagation."]
        
        # Zonal Wind
        var = _create_variable(ncfile, '%s_Zonal_Wind_Precision_1_Sample'%prefix, L22_dict['u_precision_1_sample'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='1-sample precision in the zonal wind estimate.', 
                              display_type='image', field_name='Zonal Wind Prec 1-Sample', fill_value=None, label_axis='Z Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=[err_prefix[0],
                                     "The \"1 sample\" error variable quantifies errors that are uncorrelated from one exposure to the next, "
                                     "dominated by shot and dark noise in the detectors.  The correlation time of this error source is 30-60 "
                                     "seconds (i.e., the measurement cadence). The reported error is estimated from the fringe intensity and "
                                     "background. This is the recommended variable to use for analyses of wind fluctuations within a single day "
                                     "and a single altitude (e.g., gravity waves). Because the Level 2.2 data include interpolation of Level 2.1 "
                                     "data, some correlation remains between consecutive samples. Errors are slightly correlated across small "
                                     "altitude gaps as a result of the inversion."]
                              )
        
        var = _create_variable(ncfile, '%s_Zonal_Wind_Precision_1_Day'%prefix, L22_dict['u_precision_1_day'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='1-day precision in the zonal wind estimate.', 
                              display_type='image', field_name='Zonal Wind Prec 1-Day', fill_value=None, label_axis='Z Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=[err_prefix[0],
                                     "The \"1 Day\" error variable quantifies the error introduced by daily calibrations, which is correlated "
                                     "for an entire 24-hour period (00:00 - 23:59 UT). This is estimated from the magnitude of fluctuations "
                                     "in the daily-averaged phase, propagated through the inversion. Errors in day mode and night mode are "
                                     "nearly uncorrelated. For studies pertaining to "
                                     "atmospheric tidal modes that combine data from many days, this error can be treated as uncorrelated "
                                     "across time."]
                              )
        
        var = _create_variable(ncfile, '%s_Zonal_Wind_Precision_Low_Signal_Effect'%prefix, L22_dict['u_precision_low_signal'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Low-signal precision in the zonal wind estimate.', 
                              display_type='image', field_name='Zonal Wind Prec Low Sig', fill_value=None, label_axis='Z Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=[err_prefix[0],
                                     "The \"Low Signal Effect\" error variable quantifies the error introduced by the imperfect correction for the "
                                     "signal-dependent phase shift, which is an effect seen in atmospheric and calibration-lamp fringes where "
                                     "the phase of the fringes is biased at very low signal levels. This is under investigation but could "
                                     "be caused by a charge trapping effect in the CCD. A "
                                     "correction has been implemented based upon the empirical relationship between measured phase "
                                     "and signal level of the calibration lamps for the first ~30 months of the mission. However, especially for "
                                     "cases with low signal levels, "
                                     "this correction is uncertain. The uncertainty in the resulting winds is estimated from the signal level "
                                     "and reported in this variable. It is likely to be correlated across samples nearby in time and space, but "
                                     "the correlation between different channels (red and green), sensors (MIGHTI-A and MIGHTI-B), and operating "
                                     "modes (Day and Night) is not known. Depending on the analysis being used, it could be treated as a systematic "
                                     "error or as a statistical error. Where this uncertainty is large, caution is recommended. "
                                     "For example, for winds in the core science region (90-105 km altitude), the magnitude of the correction is small or zero, "
                                     "but data in the red channel during the night and twilight are subject to a large correction (many tens of m/s) "
                                     "and the uncertainty is correspondingly large. "
                                     "A goal for future releases is to characterize and correct this effect more accurately."]
                              )
        
        var = _create_variable(ncfile, '%s_Zonal_Wind_Accuracy'%prefix, L22_dict['u_accuracy'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Accuracy of the zonal wind estimate.', 
                              display_type='image', field_name='Zonal Wind Prec 1-Day', fill_value=None, label_axis='Z Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=[err_prefix[0],
                                     "The \"Accuracy\" variable quantifies the error introduced by the zero-wind phase estimate. It is "
                                     "strongly correlated across time lags of days to weeks and becomes increasingly decorrelated for time "
                                     "lags longer than 2 precession cycles (96 days). This error is estimated from the discrepancy between "
                                     "various techniques of determining the zero-wind phase. This error source is irrelevant for most users "
                                     "studying perturbations from the mean (e.g., tides, waves), but may be important for studies of "
                                     "zonal mean winds, point comparisons with other data sets, and seasonal/long-term trends thereof. Errors are "
                                     "moderately correlated across "
                                     "small altitude gaps. Errors in day mode and night mode are nearly uncorrelated, implying there could "
                                     "be different offsets for day mode and night mode. This could be important for error propagation of "
                                     "odd-numbered migrating tides (e.g., DW1)."]
                              )
        
        var = _create_variable(ncfile, '%s_Zonal_Wind_Error'%prefix, L22_dict['u_error'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Error in the zonal wind estimate.', 
                              display_type='image', field_name='Zonal Wind Error', fill_value=None, label_axis='Z Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="For robust error propagation, users are encouraged to consider the individual error variables: "
                               "\"Precision_1_Sample,\" \"Precision_1_Day,\" and \"Accuracy.\" The \"Wind_Error\" variable is "
                               "included for backwards compatibility and is equal to the quadrature sum of the \"1 Sample\" error "
                               "and the \"1 Day\" error. This is the recommended uncertainty to use for analyses that collect data "
                               "from several weeks and compute perturbations from the mean (e.g., for estimating tides and planetary "
                               "waves). This error is uncorrelated across time lags larger than 24 hours. Errors are slightly correlated "
                               "across small altitude gaps. Errors in day mode and night mode are nearly uncorrelated."
                              )

        # Meridional Winds
        var = _create_variable(ncfile, '%s_Meridional_Wind_Precision_1_Sample'%prefix, L22_dict['v_precision_1_sample'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='1-sample precision in the meridional wind estimate.', 
                              display_type='image', field_name='Merid Wind Prec 1-Sample', fill_value=None, label_axis='M Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=[err_prefix[0],
                                     "The \"1 sample\" error variable quantifies errors that are uncorrelated from one exposure to the next, "
                                     "dominated by shot and dark noise in the detectors.  The correlation time of this error source is 30-60 "
                                     "seconds (i.e., the measurement cadence). The reported error is estimated from the fringe intensity and "
                                     "background. This is the recommended variable to use for analyses of wind fluctuations within a single day "
                                     "and a single altitude (e.g., gravity waves). Because the Level 2.2 data include interpolation of Level 2.1 "
                                     "data, some correlation remains between consecutive samples. Errors are slightly correlated across small "
                                     "altitude gaps as a result of the inversion."]
                              )
        
        var = _create_variable(ncfile, '%s_Meridional_Wind_Precision_1_Day'%prefix, L22_dict['v_precision_1_day'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='1-day precision in the meridional wind estimate.', 
                              display_type='image', field_name='Meridional Wind Prec 1-Day', fill_value=None, label_axis='M Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=[err_prefix[0],
                                     "The \"1 Day\" error variable quantifies the error introduced by daily calibrations, which is correlated "
                                     "for an entire 24-hour period (00:00 - 23:59 UT). This is estimated from the magnitude of fluctuations "
                                     "in the daily-averaged phase, propagated through the inversion. Errors in day mode and night mode are "
                                     "nearly uncorrelated. For studies pertaining to "
                                     "atmospheric tidal modes that combine data from many days, this error can be treated as uncorrelated "
                                     "across time."]
                              )        
        
        var = _create_variable(ncfile, '%s_Meridional_Wind_Precision_Low_Signal_Effect'%prefix, L22_dict['v_precision_low_signal'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Low-signal precision in the meridional wind estimate.', 
                              display_type='image', field_name='Merid Wind Prec Low Sig', fill_value=None, label_axis='M Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=[err_prefix[0],
                                     "The \"Low Signal Effect\" error variable quantifies the error introduced by the imperfect correction for the "
                                     "signal-dependent phase shift, which is an effect seen in atmospheric and calibration-lamp fringes where "
                                     "the phase of the fringes is biased at very low signal levels. This is under investigation but could "
                                     "be caused by a charge trapping effect in the CCD. A "
                                     "correction has been implemented based upon the empirical relationship between measured phase "
                                     "and signal level of the calibration lamps for the first ~30 months of the mission. However, especially for "
                                     "cases with low signal levels, "
                                     "this correction is uncertain. The uncertainty in the resulting winds is estimated from the signal level "
                                     "and reported in this variable. It is likely to be correlated across samples nearby in time and space, but "
                                     "the correlation between different channels (red and green), sensors (MIGHTI-A and MIGHTI-B), and operating "
                                     "modes (Day and Night) is not known. Depending on the analysis being used, it could be treated as a systematic "
                                     "error or as a statistical error. Where this uncertainty is large, caution is recommended. "
                                     "For example, for winds in the core science region (90-105 km altitude), the magnitude of the correction is small or zero, "
                                     "but data in the red channel during the night and twilight are subject to a large correction (many tens of m/s) "
                                     "and the uncertainty is correspondingly large. "
                                     "A goal for future releases is to characterize and correct this effect more accurately."]
                              )
        
        var = _create_variable(ncfile, '%s_Meridional_Wind_Accuracy'%prefix, L22_dict['v_accuracy'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Accuracy of the meridional wind estimate.', 
                              display_type='image', field_name='Meridional Wind Prec 1-Day', fill_value=None, label_axis='M Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes=[err_prefix[0],
                                     "The \"Accuracy\" variable quantifies the error introduced by the zero-wind phase estimate. It is "
                                     "strongly correlated across time lags of days to weeks and becomes increasingly decorrelated for time "
                                     "lags longer than 2 precession cycles (96 days). This error is estimated from the discrepancy between "
                                     "various techniques of determining the zero-wind phase. This error source is irrelevant for most users "
                                     "studying perturbations from the mean (e.g., tides, waves), but may be important for studies of "
                                     "zonal mean winds, point comparisons with other data sets, and seasonal/long-term trends thereof. Errors are "
                                     "moderately correlated across "
                                     "small altitude gaps. Errors in day mode and night mode are nearly uncorrelated, implying there could "
                                     "be different offsets for day mode and night mode. This could be important for error propagation of "
                                     "odd-numbered migrating tides (e.g., DW1)."]
                              )
        
        var = _create_variable(ncfile, '%s_Meridional_Wind_Error'%prefix, L22_dict['v_error'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Error in the meridional wind estimate.', 
                              display_type='image', field_name='Meridional Wind Error', fill_value=None, label_axis='M Wind Err', bin_location=0.5,
                              units='m/s', valid_min=0.0, valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="For robust error propagation, users are encouraged to consider the individual error variables: "
                               "\"Precision_1_Sample,\" \"Precision_1_Day,\" and \"Accuracy.\" The \"Wind_Error\" variable is "
                               "included for backwards compatibility and is equal to the quadrature sum of the \"1 Sample\" error "
                               "and the \"1 Day\" error. This is the recommended uncertainty to use for analyses that collect data "
                               "from several weeks and compute perturbations from the mean (e.g., for estimating tides and planetary "
                               "waves). This error is uncorrelated across time lags larger than 24 hours. Errors are slightly correlated "
                               "across small altitude gaps. Errors in day mode and night mode are nearly uncorrelated."
                              )
        

        # Quality factor
        var = _create_variable(ncfile, '%s_Wind_Quality'%prefix, L22_dict['wind_quality'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='A quantification of the quality, from 0 (Bad) to 1 (Good)', 
                              display_type='image', field_name='Wind Quality', fill_value=None, label_axis='Wind Qual', bin_location=0.5,
                              units='', valid_min=0.0, valid_max=1.0, var_type='data', chunk_sizes=[nx,ny],
                              notes=["A quantification of the overall quality of the wind data. "
                              "While the intent is that the XXX_Wind_Error variable accurately characterizes the statistical "
                              "error in the wind data, it is possible that systematic errors are present, or that the statistical error "
                              "estimation is not accurate. If this is suspected to be the case, the quality will be less than 1.0. If "
                              "the data are definitely unusable, the quality will be 0.0 and the sample will be masked. Users should "
                              "exercise caution when the quality is less than 1.0.",
                              "Currently, the quality can take values of 0 (Bad), 0.5 (Caution), or 1 (Good)."
                                     ])

        # Fringe amplitude (combining MIGHTI-A and MIGHTI-B)
        var = _create_variable(ncfile, '%s_Fringe_Amplitude'%prefix, L22_dict['fringe_amp'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Fringe Amplitude', 
                              display_type='image', field_name='Fringe Amplitude', fill_value=None, label_axis='Fringe Amp', bin_location=0.5,
                              units='arb', valid_min=0, valid_max=1e30, var_type='data', chunk_sizes=[nx,ny],
                              notes="An approximate volume emission rate (VER) profile in arbitrary units, estimated by combining "
                              "MIGHTI-A and MIGHTI-B data. Technically this is not the VER, but rather the amplitude of the fringes, "
                              "which has a dependence on thermospheric temperature and background emission. Thus, it does not truly "
                              "represent volume emission rate. However, it is a useful proxy. The units are arbitrary, as the "
                              "fringe amplitudes are not calibrated. See also variables Fringe_Amplitude_Relative_Difference, "
                              "Fringe_Amplitude_A, and Fringe_Amplitude_B."
                              )
        
        # Fringe amplitude error
        var = _create_variable(ncfile, '%s_Fringe_Amplitude_Error'%prefix, L22_dict['fringe_amp_error'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Error in the fringe amplitude estimate', 
                              display_type='image', field_name='Fringe Amplitude Error', fill_value=None, label_axis='Amp Error', bin_location=0.5,
                              units='arb', valid_min=0, valid_max=1e30, var_type='data', chunk_sizes=[nx,ny],
                              notes="The statistical (1-sigma) error in the fringe amplitude estimate, propagated from error in the "
                              "MIGHTI-A and MIGHTI-B inversions. The units are arbitrary, as the fringe amplitudes are not absolutely "
                              "calibrated. Systematic errors, such as those arising from airglow gradients or cross-calibration, are not "
                              "included in this variable, but are probably the dominant source of total error."
                              )
        
        # VER (combining MIGHTI-A and MIGHTI-B)
        var = _create_variable(ncfile, '%s_Relative_VER'%prefix, L22_dict['ver'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Relative volume emission rate', 
                              display_type='image', field_name='Relative VER', fill_value=None, label_axis='Rel VER', bin_location=0.5,
                              units='ph/cm^3/s', valid_min=0, valid_max=1e30, var_type='data', chunk_sizes=[nx,ny],
                              notes= "The volume emission rate (VER) obtained by averaging the VER from MIGHTI-A and MIGHTI-B, which is "
                               "obtained by scaling the fringe amplitude by a calibration factor, as described in Data Product 2.1. "
                               "Pre-flight calibrations and on-orbit comparisons with ground-based instruments are used to determine the "
                               "best possible calibration. The fringe amplitude has a dependence on temperature, which is corrected using "
                               "the MSIS model. Because the on-orbit calibration is uncertain, and because the MSIS "
                               "temperature correction is not perfect, caution should be exercised when absolute "
                               "calibration is required, or when precise comparisons are being made between samples at very "
                               "different temperatures. "
                               "Please contact the MIGHTI team before performing any studies that require absolute calibration. "
                               "The statistical (1-sigma) error for this variable is provided in "
                               "the variable ICON_..._Relative_VER_Error, though it is expected that systematic calibration errors dominate "
                               "the total error."
                              )
        
        # VER Error
        var = _create_variable(ncfile, '%s_Relative_VER_Error'%prefix, L22_dict['ver_error'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Error in VER estimate (statistical)', 
                              display_type='image', field_name='Relative VER Error', fill_value=None, label_axis='VER Err', bin_location=0.5,
                              units='ph/cm^3/s', valid_min=0, valid_max=1e30, var_type='data', chunk_sizes=[nx,ny],
                              notes= "The statistical (1-sigma) error in the relative VER estimate, propagated from error in the "
                              "MIGHTI-A and MIGHTI-B inversions. This error arises mostly from shot noise. Importantly, it is "
                              "expected that systematic errors (e.g., calibration errors) dominate the total error, but they are "
                              "not included in this variable."
                              )
        
        # Quality factor for VER/Fringe amplitude
        var = _create_variable(ncfile, '%s_VER_Quality'%prefix, L22_dict['ver_quality'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='A quantification of the quality, from 0 (Bad) to 1 (Good)', 
                              display_type='image', field_name='VER Quality', fill_value=None, label_axis='VER Qual', bin_location=0.5,
                              units='', valid_min=0.0, valid_max=1.0, var_type='data', chunk_sizes=[nx,ny],
                              notes=["A quantification of the overall quality of the VER data. "
                              "While the intent is that the XXX_VER_Error variable accurately characterizes the statistical "
                              "error in the VER data, it is possible that systematic errors are present, or that the statistical error "
                              "estimation is not accurate. If it is suspected that this is the case, the quality will be less than 1.0. If "
                              "the data are definitely unusable, the quality will be 0.0 and the sample will be masked. Users should "
                              "exercise caution when the quality is less than 1.0.",
                              "Currently, the quality can take values of 0 (Bad), 0.5 (Caution), or 1 (Good)."
                               ])       
        
        # Magnetic field-aligned/meridional/zonal wind
        var = _create_variable(ncfile, '%s_Magnetic_Field_Aligned_Wind'%prefix, L22_dict['wind_mag_fa'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Magnetic field-aligned component of the wind', 
                              display_type='image', field_name='MgFA Wind', fill_value=None, label_axis='MgFA Wind', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="The component of the wind in the direction of the magnetic field line, assuming vertical winds "
                               "are negligible. This variable is calculated by taking the geographic zonal and meridional wind "
                               "(the primary data products in this file) and expressing the wind vector in a local magnetic coordinate "
                               "system defined using the Python package OMMBV (https://github.com/rstoneback/OMMBV). "
                               "The "
                               "coordinate system used here is orthogonal and is identical to the coordinate system used to express the ion "
                               "drifts in the ICON IVM data product 2.7 (i.e., the variables ICON_L27_Ion_Velocity_Meridional, ICON_L27_Ion_Velocity_Zonal, "
                               "and ICON_L27_Ion_Velocity_Field_Aligned)."
                               )
        
        var = _create_variable(ncfile, '%s_Magnetic_Meridional_Wind'%prefix, L22_dict['wind_mag_mer'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Magnetic meridional component of the wind', 
                              display_type='image', field_name='MgMerWind', fill_value=None, label_axis='MgMer Wind', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="The component of the wind in the magnetic meridional direction, assuming vertical winds are negligible. "
                               "This variable is calculated by taking the geographic zonal and meridional wind "
                               "(the primary data products in this file) and expressing the wind vector in a local magnetic coordinate "
                               "system defined using the Python package OMMBV (https://github.com/rstoneback/OMMBV). "
                               "The magnetic meridional unit vector is orthogonal to the magnetic field line but within the plane of "
                               "the magnetic meridian (defined by the apex of the field line and its footpoint ). "
                               "At the magnetic equator, the meridional direction points up, while away from the equator it has "
                               "a poleward component (north in the northern hemisphere, south in the southern hemisphere). "
                               "Note that in some ion-neutral coupling models, a definition of magnetic meridional is often used that is "
                               "horizontal (i.e., perpendicular to gravity) and generally northward. The definition used here is "
                               "perpendicular to B and thus has primarily a vertical component at ICON latitudes. "
                               "Note also that "
                               "the definition of magnetic meridional and zonal used here differs from quasi-dipole and apex coordinate bases. The "
                               "coordinate system used here is orthonormal and is identical to the coordinate system used to express the ion "
                               "drifts in the ICON IVM data product 2.7 (i.e., the variables ICON_L27_Ion_Velocity_Meridional, ICON_L27_Ion_Velocity_Zonal, "
                               "and ICON_L27_Ion_Velocity_Field_Aligned)."
                              )
        
        var = _create_variable(ncfile, '%s_Magnetic_Zonal_Wind'%prefix, L22_dict['wind_mag_zon'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Magnetic zonal component of the wind', 
                              display_type='image', field_name='MgZon Wind', fill_value=None, label_axis='MgZon Wind', bin_location=0.5,
                              units='m/s', valid_min=-4000., valid_max=4000., var_type='data', chunk_sizes=[nx,ny],
                              notes="The component of the wind in the magnetic zonal direction, assuming vertical winds are negligible. "
                               "This variable is calculated by taking the geographic zonal and meridional wind "
                               "(the primary data products in this file) and expressing the wind vector in a local magnetic coordinate "
                               "system defined using the Python package OMMBV (https://github.com/rstoneback/OMMBV). "
                               "At the magnetic equator, the zonal direction points horizontally, while away from the equator "
                               "it can have a slightly vertical component. Note that "
                               "the definition of magnetic meridional and zonal used here differs from quasi-dipole and apex coordinate bases. The "
                               "coordinate system used here is orthonormal and is identical to the coordinate system used to express the ion "
                               "drifts in the ICON IVM data product 2.7 (i.e., the variables ICON_L27_Ion_Velocity_Meridional, ICON_L27_Ion_Velocity_Zonal, "
                               "and ICON_L27_Ion_Velocity_Field_Aligned)."
                              )

        ######### Data Location Variables #########

        # Altitude
        val = L22_dict['alt']
        var = _create_variable(ncfile, var_alt_name, val, 
                              dimensions=(var_alt_name), depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='WGS84 altitude of each wind sample', 
                              display_type='image', field_name='Altitude', fill_value=None, label_axis='Altitude', bin_location=0.5,
                              units='km', valid_min=50., valid_max=1000., var_type='support_data', chunk_sizes=[ny], monoton='Increase',
                              notes="A one-dimensional array defining the altitude dimension of the data grid (the other dimension "
                              "being time). Altitude is defined using the WGS84 ellipsoid."
                              )


        # Longitude
        var = _create_variable(ncfile, '%s_Longitude'%prefix, L22_dict['lon'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='WGS84 longitude of each wind sample', 
                              display_type='image', field_name='Longitude', fill_value=None, label_axis='Longitude', bin_location=0.5,
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
                              display_type='image', field_name='Latitude', fill_value=None, label_axis='Latitude', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='support_data', chunk_sizes=[nx,ny],
                              notes="A two-dimensional array defining the latitude of the two-dimensional data grid. The latitude varies "
                              "only slightly (a few deg) with altitude, but this variation is included. Latitude is defined using the WGS84 "
                              "ellipsoid. It should be noted that while a single latitude value is given for each point, the observation is "
                              "inherently a horizontal average over many hundreds of kilometers."
                              )
        
        # Magnetic Latitude
        var = _create_variable(ncfile, '%s_Magnetic_Latitude'%prefix, L22_dict['mag_lat'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Magnetic quasi-dipole latitude of each wind sample', 
                              display_type='image', field_name='Mag Lat', fill_value=None, label_axis='Mag Lat', bin_location=0.5,
                              units='deg', valid_min=-90., valid_max=90., var_type='support_data', chunk_sizes=[nx,ny],
                              notes="A two-dimensional array defining the magnetic quasi-dipole latitude of the two-dimensional data grid. "
                              "The latitude varies only slightly (a few deg) with altitude, but this variation is included. "
                              "It should be noted that while a single latitude value is given for each point, the observation is "
                              "inherently a horizontal average over many hundreds of kilometers. "
                              "Quasi-dipole latitude and longitude are calculated using the fast implementation developed by "
                              "Emmert et al. (2010, doi:10.1029/2010JA015326) and the Python wrapper apexpy "
                              "(doi.org/10.5281/zenodo.1214207). "
                              )
        
        # Magnetic Longitude
        var = _create_variable(ncfile, '%s_Magnetic_Longitude'%prefix, L22_dict['mag_lon'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Magnetic quasi-dipole longitude of each wind sample', 
                              display_type='image', field_name='Mag Lon', fill_value=None, label_axis='Mag Lon', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=360., var_type='support_data', chunk_sizes=[nx,ny],
                              notes="A two-dimensional array defining the magnetic quasi-dipole longitude of the two-dimensional data grid. "
                              "The longitude varies only slightly (a few deg) with altitude, but this variation is included. "
                              "It should be noted that while a single longitude value is given for each point, the observation is "
                              "inherently a horizontal average over many hundreds of kilometers. "
                              "Quasi-dipole latitude and longitude are calculated using the fast implementation developed by "
                              "Emmert et al. (2010, doi:10.1029/2010JA015326) and the Python wrapper apexpy "
                              "(doi.org/10.5281/zenodo.1214207). Quasi-dipole longitude is defined such that zero "
                              "occurs where the geodetic longitude is near 285 deg east (depending on latitude)."
                              )
        
        ######### Other Metadata Variables #########
        
        # Solar zenith angle
        var = _create_variable(ncfile, '%s_Solar_Zenith_Angle'%prefix, L22_dict['sza'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Solar zenith angle of each wind sample', 
                              display_type='image', field_name='SZA', fill_value=None, label_axis='SZA', bin_location=0.5,
                              units='deg', valid_min=0., valid_max=180., var_type='support_data', chunk_sizes=[nx,ny],
                              notes="Angle between the vectors towards the sun and towards zenith, for each point in the grid."
                              )
                               
        # Solar local time                      
        var = _create_variable(ncfile, '%s_Local_Solar_Time'%prefix, L22_dict['slt'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Local solar time of each wind sample', 
                              display_type='image', field_name='LST', fill_value=None, label_axis='LST', bin_location=0.5,
                              units='hour', valid_min=0., valid_max=24., var_type='support_data', chunk_sizes=[nx,ny],
                              notes="Local solar time at each point in the grid, calculating using the equation of time."
                              )

        # Difference between the MIGHTI-A and MIGHTI-B times contributing to each point
        var = _create_variable(ncfile, '%s_Time_Delta'%prefix, L22_dict['time_delta'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Difference between MIGHTI-A and B times contributing to each point', 
                              display_type='image', field_name='Time Delta', fill_value=None, label_axis='Time Delta', bin_location=0.5,
                              units='s', valid_min=-30.*60, valid_max=30.*60, var_type='support_data', chunk_sizes=[nx,ny],
                              notes="To determine the cardinal wind at each point, a MIGHTI-A line-of-sight wind is combined "
                              "with a MIGHTI-B line-of-sight wind from several minutes later. This variable contains this time "
                              "difference for every point. During standard operations (LVLH Normal), this variable should be positive, but can "
                              "potentially become negative during conjugate operations or when ICON is observing to the south (LVLH Reverse)."
                              )

        # Fringe amplitude profile from MIGHTI-A
        var = _create_variable(ncfile, '%s_Fringe_Amplitude_A'%prefix, L22_dict['fringe_amp_A'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Fringe Amplitude from MIGHTI-A', 
                              display_type='image', field_name='Fringe Amplitude A', fill_value=None, label_axis='Amp A', bin_location=0.5,
                              units='arb', valid_min=-1e30, valid_max=1e30, var_type='metadata', chunk_sizes=[nx,ny],
                              notes="See Fringe_Amplitude. This variable contains the fringe amplitude measured "
                              "by MIGHTI-A, interpolated to the reconstruction grid. This is one of two variables "
                              "used to create Fringe_Amplitude."
                              )

        # Fringe amplitude profile from MIGHTI-B
        var = _create_variable(ncfile, '%s_Fringe_Amplitude_B'%prefix, L22_dict['fringe_amp_B'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Fringe Amplitude from MIGHTI-B', 
                              display_type='image', field_name='Fringe Amplitude B', fill_value=None, label_axis='Amp B', bin_location=0.5,
                              units='arb', valid_min=-1e30, valid_max=1e30, var_type='metadata', chunk_sizes=[nx,ny],
                              notes="See Fringe_Amplitude. This variable contains the fringe amplitude measured "
                              "by MIGHTI-B, interpolated to the reconstruction grid. This is one of two variables "
                              "used to create Fringe_Amplitude."
                              )

        # VER profile from MIGHTI-A
        var = _create_variable(ncfile, '%s_Relative_VER_A'%prefix, L22_dict['ver_A'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Relative VER from MIGHTI-A', 
                              display_type='image', field_name='VER A', fill_value=None, label_axis='VER A', bin_location=0.5,
                              units='ph/cm^3/s', valid_min=-1e30, valid_max=1e30, var_type='metadata', chunk_sizes=[nx,ny],
                              notes="See Relative_VER. This variable contains the VER measured "
                              "by MIGHTI-A, interpolated to the reconstruction grid. This is one of two variables "
                              "used to create Relative_VER. When A and B are significantly different, large horizontal "
                              "gradients are suspected, and the quality is reduced."
                              )

        # VER profile from MIGHTI-B
        var = _create_variable(ncfile, '%s_Relative_VER_B'%prefix, L22_dict['ver_B'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='Relative VER from MIGHTI-B', 
                              display_type='image', field_name='VER B', fill_value=None, label_axis='VER B', bin_location=0.5,
                              units='ph/cm^3/s', valid_min=-1e30, valid_max=1e30, var_type='metadata', chunk_sizes=[nx,ny],
                              notes="See Relative_VER. This variable contains the VER measured "
                              "by MIGHTI-B, interpolated to the reconstruction grid. This is one of two variables "
                              "used to create Relative_VER. When A and B are significantly different, large horizontal "
                              "gradients are suspected, and the quality is reduced."
                              )

        # VER relative difference (fringe amp rel diff is not saved because it is not used, and also it can be calculated from
        # Fringe_Amplitude_A and Fringe_Amplitude_B)
        var = _create_variable(ncfile, '%s_VER_Relative_Difference'%prefix, L22_dict['ver_rel_diff'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc="Difference in MIGHTI A and B's VER estimates, divided by the mean", 
                              display_type='image', field_name='VER Difference', fill_value=None, label_axis='VER Diff', bin_location=0.5,
                              units='', valid_min=0.0, valid_max=1e10, var_type='metadata', chunk_sizes=[nx,ny],
                              notes="The absolute value of the difference between Relative_VER_A and Relative_VER_B, divided "
                              "by the average. Ideally, MIGHTI A and B should measure the same VER. When they do not, this is an "
                              "indication of potential violations of the spherical symmetry assumption inherent to the inversion. This "
                              "is the parameter used to determine if the spherical asymmetry flag is raised."
                              )
        
        # Orbit number
        var = _create_variable(ncfile, '%s_Orbit_Number'%prefix, L22_dict['orbit_number'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='ICON orbit number', 
                              display_type='image', field_name='Orbit Number', fill_value=None, label_axis='Orb Num', bin_location=0.5,
                              units='', valid_min=0.0, valid_max=1e6, var_type='data', chunk_sizes=[nx,ny],
                              notes=["The ICON orbit number corresponding to each grid point. This is usually an integer, "
                              "but when samples from two different orbits are used, an interpolated (fractional) value is used."
                               ])
        
        # Orbit ascending/descending node flag
        var = _create_variable(ncfile, '%s_Orbit_Node'%prefix, L22_dict['orbit_node'].T, 
                              dimensions=('Epoch', var_alt_name), depend_0 = 'Epoch', depend_1 = var_alt_name,
                              format_nc='f8', format_fortran='F', desc='ICON orbit ascending/descending flag', 
                              display_type='image', field_name='Orbit Node', fill_value=None, label_axis='Orb Node', bin_location=0.5,
                              units='', valid_min=0.0, valid_max=1.0, var_type='data', chunk_sizes=[nx,ny],
                              notes=["A flag indicating whether ICON is in the ascending (0) or descending (1) part of the orbit. "
                              "For some grid points, samples from MIGHTI-A are on the descending part of the orbit, while samples from "
                              "MIGHTI-B are ascending. In these cases an interpolated value is used (between 0 and 1)."
                               ])

        # Quality flags      
        var = _create_variable(ncfile, '%s_Quality_Flags'%prefix, np.transpose(L22_dict['quality_flags'], (1,0,2)), # Transpose so time is first dim
                              dimensions=('Epoch', var_alt_name, 'N_Flags'), depend_0 = 'Epoch', depend_1 = var_alt_name, depend_2 = 'N_Flags',
                              format_nc='i8', format_fortran='I', desc='Quality flags', 
                              display_type='image', field_name='Quality Flags', fill_value=-1, label_axis='Qual Flags', bin_location=0.5,
                              units='', valid_min=0, valid_max=1, var_type='metadata', chunk_sizes=[nx,ny,nflags],
                              notes=["This variable provides information on why the Quality variable is reduced from 1.0. Many quality flags can "
                              "be raised for each grid point, and each flag takes values 0 or 1. More than one flag can be raised per point. "
                              "This variable is a three-dimensional array with dimensions time, altitude, and number of flags. Each entry is 0 or 1. "
                              "Most quality flags are passed through from the L1 and L2.1 algorithms (after interpolation to the L2.2 grid). "
                              "Some additional quality flags are created in L2.2. The N_Flags dimension is defined below: ",

                                       "* 0 : (From L1 A) SNR too low to reliably perform L1 processing",
                                       "* 1 : (From L1 A) Proximity to South Atlantic Anomaly",
                                       "* 2 : (From L1 A) Bad calibration",
                                       "* 3 : (From L1 A) Calibration lamps are on",
                                       "* 4 : (From L1 A) Solar/lunar contamination",
                                       "* 5 : (From L2.1 A) Not enough valid points in profile",
                                       "* 6 : (From L2.1 A) SNR too low after inversion",
                                       "* 7 : (From L2.1 A) Significant airglow above 300 km",
                                       "* 8 : (From L2.1 A) Line of sight crosses the terminator",
                                       "* 9 : (From L2.1 A) Thermal drift correction is uncertain",
                                       "* 10: (From L2.1 A) S/C pointing is not stable",
                                       "* 11: (From L2.1 A) SNR is low after inversion, but maybe usable",
                                       "* 12: (From L1 B) SNR too low to reliably perform L1 processing",
                                       "* 13: (From L1 B) Proximity to South Atlantic Anomaly",
                                       "* 14: (From L1 B) Bad calibration",
                                       "* 15: (From L1 B) Calibration lamps are on",
                                       "* 16: (From L1 B) Solar/lunar contamination",
                                       "* 17: (From L2.1 B) Not enough valid points in profile",
                                       "* 18: (From L2.1 B) SNR too low after inversion",
                                       "* 19: (From L2.1 B) Significant airglow above 300 km",
                                       "* 20: (From L2.1 B) Line of sight crosses the terminator",
                                       "* 21: (From L2.1 B) Thermal drift correction is uncertain",
                                       "* 22: (From L2.1 B) S/C pointing is not stable",
                                       "* 23: (From L2.1 B) SNR is low after inversion, but maybe usable",
                                       "* 24: Missing MIGHTI-A file",
                                       "* 25: Missing MIGHTI-B file",
                                       "* 26: MIGHTI-A did not sample this altitude",
                                       "* 27: MIGHTI-B did not sample this altitude",
                                       "* 28: Spherical asymmetry: A&B VER estimates disagree",
                                       "* 29: Mixing Normal and Reverse LVLH",
                                       "* 30: Unused",
                                       "* 31: Unused",
                                       "* 32: Unused",
                                       "* 33: Unknown error (notify MIGHTI team)",
                                    ])

        ncfile.close()
        
    except: # Make sure the file is closed
        ncfile.close()
        raise
            
    return L22_full_fn




def level21_to_level22_without_info_file(A_curr_fn, B_curr_fn, A_prev_fn, B_prev_fn, A_next_fn, B_next_fn, 
                                         L22_path, data_revision=0, sph_asym_thresh=None, time_start=None, time_stop=None,
                                         skip_att = ['att_conjugate', 'att_zero_wind']):
    '''
    High-level function to apply the Level-2.1-to-Level-2.2 algorithm to Level 2.1 files from MIGHTI-A and B from the same
    color, and create a single Level 2.2 file (in the L22_path directory). This version
    of the function requires the user to input parameters manually, instead of specifying an 
    Information.TXT file, like is done at the Science Data Center. Files taken during conjugate maneuvers are skipped.
    
    INPUTS:
      *  A_curr_fn       -- TYPE:str       The filename of the L2.1 MIGHTI-A file on the day to be processed
      *  B_curr_fn       -- TYPE:str       The filename of the L2.1 MIGHTI-B file on the day to be processed
      *  A_prev_fn       -- TYPE:str       The filename of the L2.1 MIGHTI-A file on the previous day.
                                           If this == '', it will be skipped.
      *  B_prev_fn       -- TYPE:str       The filename of the L2.1 MIGHTI-B file on the previous day.
                                           If this == '', it will be skipped.
      *  A_next_fn       -- TYPE:str       The filename of the L2.1 MIGHTI-A file on the next day.
                                           If this == '', it will be skipped.
      *  B_next_fn       -- TYPE:str       The filename of the L2.1 MIGHTI-B file on the next day.
                                           If this == '', it will be skipped.
      *  L22_path        -- TYPE:str.      The directory the L2.2 file will be saved in, including trailing "/"
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
                                           which the Level 2.2 data product should span. (Some L2.1 data from before
                                           the start time are needed to handle the crossover). If None (default), 
                                           the start time defaults to 0 UT on the date of A_curr_fn
      *  time_stop       -- TYPE:datetime. A timezone-naive datetime in UT, specifying the end of the interval
                                           which the Level 2.2 data product should span. (Some L2.1 data from after
                                           the stop time are needed to handle the crossover). If None (default), 
                                           the stop time defaults to 24 hours after the start time. 
      *  skip_att        -- TYPE: list of str.  A list of attitude status variables. If any evaluate are active, 
                                           the associated L2.1 file(s) will be skipped and not included 
                                           in L2.2 analysis. Default value: ['att_conjugate']. A possible 
                                           alternate default value: ['att_conjugate', 'att_lvlh_reverse']. See
                                           documentation for level1_dict_to_level21_dict(...) for complete list.
    OUTPUTS:
    
      *  L22_fn          -- TYPE:str.      The full path to the saved L2.2 file.
      
    '''
    
    ################## Load data and combine into one dictionary for A and one for B ##########################
    # Load entire current day and do sanity checks
    # BJH 28 Oct 2022: Setting skip_bad_alts=True
    L21_A_curr = level21_to_dict(A_curr_fn, skip_att = skip_att, skip_bad_alts=True)
    L21_B_curr = level21_to_dict(B_curr_fn, skip_att = skip_att, skip_bad_alts=True)
    assert len(L21_A_curr['time'])>0, "No MIGHTI-A samples in %s" % A_curr_fn
    assert len(L21_B_curr['time'])>0, "No MIGHTI-B samples in %s" % B_curr_fn


    def combine(d0, d1, tstart, tstop):
        '''
        Combine the two L2.1 dictionaries (d0 and d1), only keeping samples taken between
        times tstart and tstop. It's assumed that d0['time'] < d1['time'].
        '''
        d = {}
        t0 = d0['time']
        t1 = d1['time']
        i0 = (t0 >= tstart) & (t0 <= tstop)
        i1 = (t1 >= tstart) & (t1 <= tstop)
        for v in d0.keys():

            # Special case variables.
            if v in ['source_files']: 
                d[v] = np.concatenate((d0[v], d1[v]))
                continue
                
            # These do not need to match and will not be included in output file (redundant information in
            # source_files)
            if v in ['version', 'revision']:
                continue

            # All other variables can be handled with the following general code:
            ndim = np.ndim(d0[v]) # number of dimensions of the variable
            if ndim == 0: # Don't append anything
                # Make sure they are the same though
                assert d0[v] == d1[v], "Variable '%s' does not match: v = %s vs %s" % (v, d0[v], d1[v])
                d[v] = d0[v]
            elif ndim == 1: # It's the time dimension
                d[v] = np.concatenate((d0[v][i0], d1[v][i1]))
            elif ndim == 2: # Time is the second dimension
                d[v] = np.concatenate((d0[v][:,i0], d1[v][:,i1]), axis=1)
            elif ndim == 3: # Time is the second dimension
                d[v] = np.concatenate((d0[v][:,i0,:], d1[v][:,i1,:]), axis=1)
            else:
                raise Exception('4 dimensions?!')
        return d

    # Define start and stop times (based on "A_curr" data) to include a 30 minute buffer
    # before and after the current 24-hour period.
    nt = len(L21_A_curr['time'])
    t0 = L21_A_curr['time'][nt/2] # Assume midpoint time in A_curr_fn defines "today's date"
    tstart = datetime(t0.year, t0.month, t0.day) - timedelta(hours=0.5)
    tstop  = tstart + timedelta(hours=25.)

    # Combine curr, prev, and next data, if they exist.
    # BJH 28 Oct 2022: Setting skip_bad_alts=True
    if A_prev_fn:
        L21_A_prev = level21_to_dict(A_prev_fn, skip_att = skip_att, skip_bad_alts=True)
        assert L21_A_prev['time'][-1] < L21_A_curr['time'][0], "Files in reverse time order"
        L21_A_curr = combine(L21_A_prev, L21_A_curr, tstart, tstop)
    if B_prev_fn:
        L21_B_prev = level21_to_dict(B_prev_fn, skip_att = skip_att, skip_bad_alts=True)
        assert L21_B_prev['time'][-1] < L21_B_curr['time'][0], "Files in reverse time order"
        L21_B_curr = combine(L21_B_prev, L21_B_curr, tstart, tstop)
    if A_next_fn:
        L21_A_next = level21_to_dict(A_next_fn, skip_att = skip_att, skip_bad_alts=True)
        assert L21_A_curr['time'][-1] < L21_A_next['time'][0], "Files in reverse time order"
        L21_A_curr = combine(L21_A_curr, L21_A_next, tstart, tstop)
    if B_next_fn:
        L21_B_next = level21_to_dict(B_next_fn, skip_att = skip_att, skip_bad_alts=True)
        assert L21_B_curr['time'][-1] < L21_B_next['time'][0], "Files in reverse time order"
        L21_B_curr = combine(L21_B_curr, L21_B_next, tstart, tstop)



    ########################### Run L2.2 processing and save #############################
    # determine time_start
    time_start = datetime(t0.year, t0.month, t0.day)
    
    # Run processing
    L22_dict = level21_dict_to_level22_dict(L21_A_curr, L21_B_curr, time_start=time_start)

    # Save L2.2 data to file
    L22_fn = save_nc_level22(L22_path, L22_dict, data_revision=data_revision)

    return L22_fn






def level21_to_level22(info_fn):
    '''
    Highest-level function to apply the Level-2.1-to-Level-2.2 algorithm. Inputs are specified via an information file. 
    If files from green and red are both specified (as expected when run at the Science Data Center), they will be split
    and run separately by this function. The output Level 2.2 file(s) will be saved to the "Output" folder in the directory
    specified by the information file. Summary plots will also be created in that directory.
    
    INPUTS:
    
      * info_fn  -- TYPE:str.  Full path to an ASCII file in the following format:
      
                                        [PARAMETERS]
                                        Revision=001
                                        Directory=/path/to/wherever/
                                        Day=2019-01-03
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
    L21_fns.sort() # These do not include the directory
    assert len(L21_fns)>0, "No files specified"
    
    # Parse the info file
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
    # (3) Parse day
    t = datetime.strptime(info['Day'], '%Y-%m-%d')
    
    # For both red and green, separate previous, current, and next files, and then call the
    # lower-level function which does all the real work.
    L22_fns = []
    L21_fns_used = []
    failure_messages = []
    for emission_color in ['red','green']:
        try:
            # Extract the 6 (3 per sensor) L2.1 files of this color
            # Today
            z = [fn for fn in L21_full_fns if emission_color in fn.lower() and t.strftime('%Y-%m-%d') in fn and 'MIGHTI-A' in fn]
            assert len(z) == 1, "%i files found for MIGHTI-A %s on %s in the following list:\n%s" % (len(z), emission_color, t, L21_fns)
            A_curr = z[0]

            z = [fn for fn in L21_full_fns if emission_color in fn.lower() and t.strftime('%Y-%m-%d') in fn and 'MIGHTI-B' in fn]
            assert len(z) == 1, "%i files found for MIGHTI-B %s on %s in the following list:\n%s" % (len(z), emission_color, t, L21_fns)
            B_curr = z[0]

            # Yesterday - it's ok if this file is not specified
            tprev = t - timedelta(days=1)
            z = [fn for fn in L21_full_fns if emission_color in fn.lower() and tprev.strftime('%Y-%m-%d') in fn and 'MIGHTI-A' in fn]
            assert len(z) <=1, "%i files found for MIGHTI-A %s on %s in the following list:\n%s" % (len(z), emission_color, tprev, L21_fns)
            A_prev = ''
            if len(z)==1:
                A_prev = z[0]

            z = [fn for fn in L21_full_fns if emission_color in fn.lower() and tprev.strftime('%Y-%m-%d') in fn and 'MIGHTI-B' in fn]
            assert len(z) <=1, "%i files found for MIGHTI-B %s on %s in the following list:\n%s" % (len(z), emission_color, tprev, L21_fns)
            B_prev = ''
            if len(z)==1:
                B_prev = z[0]

            # Tomorrow - it's ok if this file is not specified
            tnext = t + timedelta(days=1)
            z = [fn for fn in L21_full_fns if emission_color in fn.lower() and tnext.strftime('%Y-%m-%d') in fn and 'MIGHTI-A' in fn]
            assert len(z) <=1, "%i files found for MIGHTI-A %s on %s in the following list:\n%s" % (len(z), emission_color, tnext, L21_fns)
            A_next = ''
            if len(z)==1:
                A_next = z[0]

            z = [fn for fn in L21_full_fns if emission_color in fn.lower() and tnext.strftime('%Y-%m-%d') in fn and 'MIGHTI-B' in fn]
            assert len(z) <=1, "%i files found for MIGHTI-B %s on %s in the following list:\n%s" % (len(z), emission_color, tnext, L21_fns)
            B_next = ''
            if len(z)==1:
                B_next = z[0]

            # Call lower-level processing function
            L22_fn = level21_to_level22_without_info_file(A_curr, B_curr, A_prev, B_prev, A_next, B_next, 
                                                          direc + 'Output/', data_revision=data_revision)       
            L22_fns.append(L22_fn)
            # Register which L2.1 files were used
            L21_fns_used.extend([s.split('/')[-1] for s in [A_curr, B_curr, A_prev, B_prev, A_next, B_next]]) # Not including directory
            
        except Exception as e:
            failure_messages.append('Failed processing:\n\tColor   = %s\n%s\n'%(emission_color,traceback.format_exc()))
            
    # For both red and green, create Tohban plots
    for fn in L22_fns:
        try:
            plot_level22(fn, direc + 'Output/')
        except Exception as e:
            failure_messages.append('Failed creating Tohban plot for %s:\n%s\n'%(fn, traceback.format_exc()))
    # Create summary combined plot
    assert len(L22_fns)==2, "Detected %i L2.2 files. Need 2." % len(L22_fns)
    try:
        fng, fnr = L22_fns[0], L22_fns[1]
        if 'Green' in fnr:
            fnr, fng = L22_fns[0], L22_fns[1]
        plot_level22_summary(fng, fnr, direc + 'Output/', 'u')
        plot_level22_summary(fng, fnr, direc + 'Output/', 'v')
    except Exception as e:
        failure_messages.append('Failed creating combined summary plot:\n%s\n'%(traceback.format_exc()))
    
    # Warn if input files were not used.
    for in_fn in L21_fns:
        if in_fn not in L21_fns_used:
            # 2020 Jan 21 BJH: The SDC sometimes adds extra files, e.g., the day after tomorrow. In this
            # case the code should not report a code failure, just print to the log.
#             failure_messages.append('Input file not used: %s' % in_fn)
            print('Input file not used: %s' % in_fn)
    
    
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
                                        The variables are:
                                        
                                        * All of the variables described in documentation for level21_dict_to_level22_dict(...), plus:
                                        * epoch_ms         -- TYPE:array(nx),   UNITS:ms.  Same as epoch, except in units of 
                                                                                           "ms since Jan 1, 1970"
                                        * epoch_full_ms    -- TYPE:array(ny,nx) UNITS:ms.  Same as epoch_full, except in units of 
                                                                                           "ms since Jan 1, 1970"                    
                                        * EXCEPT the following variables (which are not saved in the L2.2 NC file):
                                        * fringe_amp_rel_diff
                                        * N_used_A
                                        * N_used_B
    '''
    
    L22_dict = {}
    f = netCDF4.Dataset(L22_fn)
    emission_color = L22_fn.split('/')[-1].split('_')[3].split('-')[-1].capitalize() # Red or Green
    pre = 'ICON_L22_' # prefix for all variables
                               
    L22_dict['lat']    = f[pre + 'Latitude'][:,:].T
    L22_dict['lon']    = f[pre + 'Longitude'][:,:].T

    # Unwrap the sample longitude to avoid 0/360 jumps
    L22_dict['lon_unwrapped'] = fix_longitudes_mat(L22_dict['lon'])

    L22_dict['alt']    = f[pre + 'Altitude'][:].T
    L22_dict['u']      = f[pre + 'Zonal_Wind'][:,:].T
    L22_dict['v']      = f[pre + 'Meridional_Wind'][:,:].T
    L22_dict['u_error']= f[pre + 'Zonal_Wind_Error'][:,:].T
    L22_dict['v_error']= f[pre + 'Meridional_Wind_Error'][:,:].T
    e                  = f[pre + 'Quality_Flags'][:,:,:]
    L22_dict['quality_flags'] = np.transpose(e, (1,0,2))
    L22_dict['epoch_ms']      = f['Epoch'][:]
    L22_dict['epoch_full_ms'] = f['Epoch_Full'][:,:].T
    L22_dict['epoch']      = None # To be filled in below
    L22_dict['epoch_full'] = None # To be filled in below
    L22_dict['time_start'] = datetime.strptime(f.Date_Start[-27:-4], '%Y-%m-%dT%H:%M:%S.%f')
    try: # In v05 we changed Date_Stop to Date_End
        L22_dict['time_stop']  = datetime.strptime(f.Date_End [-27:-4], '%Y-%m-%dT%H:%M:%S.%f')
    except:
        L22_dict['time_stop']  = datetime.strptime(f.Date_Stop [-27:-4], '%Y-%m-%dT%H:%M:%S.%f')
    L22_dict['time_delta']    = f[pre + 'Time_Delta'][:,:].T
    L22_dict['fringe_amp']    = f[pre + 'Fringe_Amplitude'][:,:].T
    L22_dict['fringe_amp_error'] = f[pre + 'Fringe_Amplitude_Error'][:,:].T
    L22_dict['fringe_amp_A']  = f[pre + 'Fringe_Amplitude_A'][:,:].T 
    L22_dict['fringe_amp_B']  = f[pre + 'Fringe_Amplitude_B'][:,:].T
    L22_dict['ver'] = f[pre + 'Relative_VER'][:,:].T
    L22_dict['ver_error'] = f[pre + 'Relative_VER_Error'][:,:].T
    L22_dict['ver_A'] = f[pre + 'Relative_VER_A'][:,:].T
    L22_dict['ver_B'] = f[pre + 'Relative_VER_B'][:,:].T
    L22_dict['ver_rel_diff']  = f[pre + 'VER_Relative_Difference'][:,:].T
    L22_dict['wind_quality']  = f[pre + 'Wind_Quality'][:,:].T
    L22_dict['ver_quality']   = f[pre + 'VER_Quality'][:,:].T
    L22_dict['emission_color'] = emission_color
    L22_dict['source_files'] = f.Parents
    L22_dict['acknowledgement'] = f.Acknowledgement
    L22_dict['wind_mag_fa']   = f[pre + 'Magnetic_Field_Aligned_Wind'][:,:].T
    L22_dict['wind_mag_mer']  = f[pre + 'Magnetic_Meridional_Wind'][:,:].T
    L22_dict['wind_mag_zon']  = f[pre + 'Magnetic_Zonal_Wind'][:,:].T
    L22_dict['mag_lat']       = f[pre + 'Magnetic_Latitude'][:,:].T
    L22_dict['mag_lon']       = f[pre + 'Magnetic_Longitude'][:,:].T
    try:
        L22_dict['slt']           = f[pre + 'Solar_Local_Time'][:,:].T
    except:
        L22_dict['slt']           = f[pre + 'Local_Solar_Time'][:,:].T
    L22_dict['sza']           = f[pre + 'Solar_Zenith_Angle'][:,:].T
    try:
        L22_dict['orbit_number']     = f[pre + 'Orbit_Number'][:,:].T
    except IndexError:
        print('WARNING: Orbit_Number not found')
    try:
        L22_dict['orbit_node'] = f[pre + 'Orbit_Node'][:,:].T
    except IndexError:
        print('WARNING: Orbit_Node not found')
    
    # Convert times to datetime
    epoch = np.empty(len(L22_dict['epoch_ms']), dtype=datetime)
    for i in range(len(epoch)):
        if np.ma.is_masked(L22_dict['epoch_ms']) and L22_dict['epoch_ms'].mask[i]:
            epoch[i] = datetime(1989,4,18) # arbitrary fill value
        else:
            epoch[i] = datetime(1970,1,1) + timedelta(seconds = L22_dict['epoch_ms'][i]/1e3)
    epoch = np.ma.masked_where(L22_dict['epoch_ms'].mask, epoch)
    np.ma.set_fill_value(epoch, datetime(1989,4,18))
            
    epoch_full = np.empty(np.shape(L22_dict['epoch_full_ms']), dtype=datetime)
    for i in range(np.shape(epoch_full)[0]):
        for j in range(np.shape(epoch_full)[1]):
            if np.ma.is_masked(L22_dict['epoch_full_ms']) and L22_dict['epoch_full_ms'].mask[i,j]:
                epoch_full[i,j] = datetime(1989,4,18) # arbitrary fill value
            else:
                epoch_full[i,j] = datetime(1970,1,1) + timedelta(seconds = L22_dict['epoch_full_ms'][i,j]/1e3)
    epoch_full = np.ma.masked_where(L22_dict['epoch_full_ms'].mask, epoch_full)
        
    L22_dict['epoch'] = epoch
    L22_dict['epoch_full'] = epoch_full
    
    f.close()
    
    return L22_dict




################################################################################################################
##########################################   TOHBAN PLOTS    ###################################################
################################################################################################################

def unwrap_slt(t, slt):
    '''
    Unwrap solar local time. slt can be 1D, or 2D (in which case its rows are unwrapped)
    This code is more complicated than simply removing negative jumps, because it needs to handle
    large data gaps.
    
    INPUTS:
    
      *  t        --TYPE:array(nt) of datetime                 Time
      *  slt      --TYPE:array(nt) or array(nt,nx). UNITS:hr.  Solar Local Time to be unwrapped
      
    OUTPUTS:
      * slt_u     --TYPE:same as slt.               UNITS:hr.  Unwrapped version of input slt.
    
    '''
    
    # In case input is a masked array, fill with nans
    slt = np.ma.filled(slt, np.nan)
    
    def unwrap_slt_1D(t, slt_1D, t0, slt0):
        '''Helper function to solve the 1D unwrapping problem
        t and slt_1D are as described above
        t0 and slt0 are scalars. Array will be unwrapped such that at time t0, the unwrapped slt
        will be within 12h of slt0'''
        # Find the average slope of the SLT vs t line
        
        time_sec = np.array([(ti - t0).total_seconds() for ti in t])
        def err(dslt_dt):
            slt_fit = np.mod(slt0 + dslt_dt*time_sec, 24.)
            return np.nanmean((slt_fit-slt_1D)**2)
        dslt_dt_min = 24./(1*3600)
        dslt_dt_max = 24./(2*3600)
        dslt_dt = optimize.brute(err, ranges=((dslt_dt_min,dslt_dt_max),), Ns=500, finish = optimize.fmin)[0]

        # Add or subtract 24 from the input array to get close to the unwrapped line
        slt_fit = slt0 + dslt_dt*time_sec
        d = slt_fit - slt_1D
        n = np.round(d/24.)
        slt_u = slt_1D + n*24.
        return slt_u
    
    # If input is already 1D, then simply call lower-level function
    if np.ndim(slt) == 1:
        k = np.where(np.isfinite(slt))[0][0] # find reference that's not nan
        return unwrap_slt_1D(t, slt, t[k], slt[k])
    
    # If input is 2D, then unwrap each row individually
    # Find a reference value to use 
    ij = np.where(np.isfinite(slt))
    i = ij[0][0]
    j = ij[1][0]
    t0 = t[j]
    slt0 = slt[i,j]
    slt_u = np.zeros_like(slt)
    for i in range(np.shape(slt)[0]):
        slt_u[i,:] = unwrap_slt_1D(t, slt[i,:], t0, slt0)
    return slt_u



def fill(x):
    '''
    If there are nans in the middle of the array x, fill them in with linear interpolation.
    
    INPUT:
     * x             --TYPE:array(nx).    Array with possible nans in the middle. Can also be of 
                                          type np.ma.MaskedArray. This array must not have nans at
                                          the beginning or end.
    
    OUTPUT:
     * x_filled      --TYPE:array(nx).    Copy of x with nans replaced with interpolated values
    '''
    
    y = np.ma.filled(x, np.nan).copy()
    n = np.arange(len(y))
    i = np.isfinite(y) # which indices are good
    # Interpolate between missing samples
    f = interpolate.interp1d(n[i], y[i])
    y[~i] = f(n[~i])
    return y



def plot_level21(L21_fn, pngpath, v_max = 200., ve_min = 1., ve_max = 100., a_min = 1., a_max = 1000.,
                 ae_min = 0.1, ae_max = 100., gap_slt = 0.4, close = True, first_only = False):
    '''
    Create Tohban plots for Level 2.2 data. The nominal L2.2 file contains 24 hours of data. The 
    Tohban plots are split up by solar local time. There is one plot per noon-to-noon LT period.
    
    INPUTS:
    
      *  L21_fn       --TYPE:str,   Full path to a L2.2 file
      *  pngpath      --TYPE:str,   Directory to save the resulting png(s) to
      
   OPTIONAL INPUTS:
   
      *  v_max        --TYPE:float, Maximum wind velocity for colorbar [m/s]
      *  ve_min       --TYPE:float, Minimum wind error for colorbar [m/s]
      *  ve_max       --TYPE:float, Maximum wind error for colorbar [m/s]
      *  a_min        --TYPE:float, Minimum VER for colorbar [ph/cm^3/s]
      *  a_max        --TYPE:float, Maximum VER for colorbar [ph/cm^3/s]
      *  ae_min       --TYPE:float, Minimum VER error for colorbar [ph/cm^3/s]
      *  ae_max       --TYPE:float, Maximum VER error for colorbar [ph/cm^3/s]
      *  gap_slt      --TYPE:float, Data gaps longer than this will be shown explicitly as gaps [hours].
                                    A 60-second exposure for ICON covers ~0.25 hours of SLT.
      *  close        --TYPE:bool,  If True, close the figure after saving it.
      *  first_only   --TYPE:bool,  If True, only make the first plot (useful for debugging).
                                   
    OUTPUT:
    
      *  L21_pngs     --TYPE:list of str,  Full path to the saved png files
    '''
    from mpl_toolkits.basemap import Basemap # For putting map on Tohban plots
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.colors import LinearSegmentedColormap

    # Custom colormap for quality factor/flags
    low_rgb = [0.7, 0.7, 1.0]
    #mid_rgb = np.array([252,238,170])/256.
    mid_rgb = np.array([253,185,40])/256.
    hi_rgb  = [0.7, 0.0, 0.0]
    cdict = {'red':   ((0.0, low_rgb[0], low_rgb[0]),
                       (0.5, mid_rgb[0], mid_rgb[0]),
                       (1.0, hi_rgb [0], hi_rgb [0])),
             'green': ((0.0, low_rgb[1], low_rgb[1]),
                       (0.5, mid_rgb[1], mid_rgb[1]),
                       (1.0, hi_rgb [1], hi_rgb [1])),
             'blue':  ((0.0, low_rgb[2], low_rgb[2]),
                       (0.5, mid_rgb[2], mid_rgb[2]),
                       (1.0, hi_rgb [2], hi_rgb [2]))}
    cmap_byr = LinearSegmentedColormap('Custom', cdict)

    if pngpath[-1] != '/':
        pngpath += '/'

    L21_dict = level21_to_dict(L21_fn)

    Nalt, Nt = np.shape(L21_dict['los_wind'])
    alt = np.ma.filled(L21_dict['alt'], np.nan) # convert from masked array to array
    t = L21_dict['time']
    slt = L21_dict['slt']
    sensor = L21_dict['sensor'] # A or B

    L21_pngs = []
    nrows, ncols = 4, 4
    csize = '3%' # For adding colorbar
    cpad = 0.08  # For adding colorbar
    clabelpad = 15 # For colorbar label

    # Load magnetic coordinates
    lonlat_m10, lonlat_0, lonlat_p10 = mag_lines()

    assert any(t), "All data in the L2.2 are masked. No plot can be created."

    # Organize figures by SLT. Create one plot for each (SLT) noon-to-noon period.
    sltu_2D = unwrap_slt(t,slt)
    sltu = np.nanmedian(sltu_2D, axis=0) # Reference SLT for each column, median to avoid discontinuities
    n = np.around(sltu/24.) # for each element, the figure it should go in

    # Create istart, istop
    istart = []
    istop = []
    for ii in np.unique(n[~np.isnan(n)]):
        idx = np.where(n == ii)[0]
        istart.append(idx[0])
        istop.append(idx[-1]+1)
    nfigs = len(istart)

    if first_only:
        istop = istop[0:1]
        istart = istart[0:1]
        nfigs = 1


    # Loop and create figures
    for n, (i1, i2) in enumerate(zip(istart, istop)):

        if i2-i1 <= 1: # corner case: skip this plot
            print('Skipping plot %i for %s\n\t(Can\'t make plot with only one time sample)' % (n, L21_fn))
            continue
        
        # pcolormesh wants the coordinates of the corners (not the middle) of each pixel
        # middle values:
        xm = np.mod(sltu[i1:i2] + 12., 24) - 12.
        dx = np.diff(xm)
        Ym = alt[:,i1:i2] # this is a matrix, unlike the L2.2 case
        dY = np.diff(Ym, axis=0)
        ny,nx = np.shape(Ym)
        X = np.zeros((ny+1, nx+1))
        Y = np.zeros((ny+1, nx+1))
        # Find SLT edges
        X[:,0] = xm[0] - dx[0]
        X[:,1:-1] = xm[:-1] + dx
        X[:,-1] = xm[-1] + dx[-1]
        # Find altitude edges (column 0 is special case)
        Y[0,1:] = Ym[0,:] - dY[0,:]/2
        Y[1:-1,1:] = Ym[:-1,:] + dY/2
        Y[-1,1:] = Ym[-1,:] + dY[-1,:]/2
        Y[:,0] = Y[:,1] # Copy second column into first
        # Mask gaps
        igap = np.where(dx > gap_slt)[0]

        fig, axarr = plt.subplots(nrows, ncols, figsize=(7.3*ncols,3*nrows))
        ax_format = [] # Will hold axes containing the normal Altitude vs. Longitude vs. Color plots, for later common formatting

        ########### Wind and wind error ###########
        ax = axarr[0,0]
        C = L21_dict['los_wind'][:,i1:i2]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap='bwr', vmin=-v_max, vmax=v_max)
        ax.set_title('LoS Wind')        
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('m/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ax = axarr[1,0]
        C = L21_dict['los_wind_error'][:,i1:i2]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap='viridis_r', norm=LogNorm(vmin=ve_min, vmax=ve_max))
        ax.set_title('LoS Wind Error')        
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('m/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ############ VER and VER Error ##########
        ax = axarr[0,1]
        C = L21_dict['ver'][:,i1:i2]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap='viridis', norm=LogNorm(vmin=a_min, vmax=a_max))
        ax.set_title('Relative VER')        
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('ph/cm$^3$/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ax = axarr[1,1]
        C = L21_dict['ver_error'][:,i1:i2]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap='viridis', norm=LogNorm(vmin=ae_min, vmax=ae_max))
        ax.set_title('Relative VER Error')      
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('ph/cm$^3$/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ############ Quality Factors #############
        cmap = cmap_byr
        bounds = [-0.33,0.33,0.67,1.33]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)  

        ax = axarr[2,0]
        C = 1 - L21_dict['wind_quality'][:,i1:i2]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Wind Quality')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.5,1])
        cb.ax.set_yticklabels(['Good','Caution','Bad'], rotation=-90, va='center')
        ax_format.append(ax)

        ax = axarr[2,1]
        C = 1 - L21_dict['ver_quality'][:,i1:i2]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('VER Quality')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.5,1])
        cb.ax.set_yticklabels(['Good','Caution','Bad'], rotation=-90, va='center')
        ax_format.append(ax)

        ############ Location #############
        try: # because I'm worried something funny may happen with this hacky code
            ax = axarr[3,0]
            j = int(ny/2) # Use middle altitude for reference lat/lon
            # Annoying code to work longitude ranges for Basemap [-360, 720]
            lon = L21_dict['lon'][j,i1:i2]
            lonu = fix_longitudes(lon)
            lonoff = lonu[0] - (np.mod(lonu[0] + 180., 360.)-180.) # offset to put min lon in [-180, 180]
            lonplot = lonu - lonoff
            latplot = L21_dict['lat'][j,i1:i2]
            # Compute lonmin and lonmax for plot. This is an unfortunately hacky way to ensure
            # that the x axis on the map matches the x axis on the other plots, even for 
            # partial orbits.
            c = np.polyfit(xm, lonplot, 1)
            dlon_dslt = c[0]
            lon_slt0 = c[1]
            lonmin = dlon_dslt*-12 + lon_slt0
            lonmax = dlon_dslt*12  + lon_slt0
            m = Basemap(ax=ax, llcrnrlon=lonmin, urcrnrlon=lonmax, llcrnrlat=-50., urcrnrlat=50.,\
                        resolution='l',projection='merc',lat_ts=0.,\
                        fix_aspect=False)
            m.plot(lonplot, latplot, 'C1-', lw=3, latlon=True, label='Observations')
            xm,ym = m(lonlat_0[:,0], lonlat_0[:,1])
            ax.plot(xm,ym, 'k-', lw=1, label='MLAT=0')
            xm,ym = m(lonlat_m10[:,0], lonlat_m10[:,1])
            ax.plot(xm,ym, 'k--', lw=1, label='MLAT=+/-10 deg')
            xm,ym = m(lonlat_p10[:,0], lonlat_p10[:,1])
            ax.plot(xm,ym, 'k--', lw=1)
            cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad) # Helps alignment
            cax.set_visible(False) # Helps alignment
            m.drawparallels(np.arange(-90.,91.,30.), labels=[0,1,0,1]); # LRTB
            m.drawmeridians(np.arange(-180.,361.,90.), labels=[0,1,0,1]);
            ax.legend(loc='lower left', prop={'size':8})
            m.drawcoastlines(linewidth=0.5);
            m.fillcontinents()
        except Exception as e:
            print('Error creating map: %s' % e)

        ############ Quality flags ############
        # There isn't enough room to plot all the flags, so we'll plot the most important ones (subject to change)
        cmap = cmap_byr
        bounds = [-0.5,0.5,1.5]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)  

        ax = axarr[0,2]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,0]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Very low SNR (before inversion)')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)
        
        ax = axarr[1,2]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,6]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Very low SNR (after inversion)')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)
        
        ax = axarr[2,2]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,11]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Cautiously low SNR (after inversion)')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)

        # 2022-07-05 Replaced SAA plot with "not enough valid samples in profile" plot
        ax = axarr[3,2]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,5]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Not enough valid points in profile')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)

        ax = axarr[3,3]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,7]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Airglow above 300 km')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)

        ax = axarr[3,1]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,8]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('LoS crosses terminator')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)
        
        ax = axarr[0,3]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,9]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Thermal drift correction uncertain')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)
        
        # 2021 Feb 18: replaced S/C jitter with solar/lunar contam flag
        ax = axarr[1,3]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,4]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Solar/lunar contamination')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)
        
        ax = axarr[2,3]
        C = 1.0 * L21_dict['quality_flags'][:,i1:i2,3]
        C[:,igap] = np.nan
        h = ax.pcolormesh(X, Y, C, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Calibration lamps are on')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)

        #### Turn off unused axes
#         axarr[3,3].axis('off')

        #### Make the 2D plots look nice
        for ax in ax_format:
            plt.sca(ax)
            ax.patch.set_facecolor('gray')
            ax.set_ylabel('Altitude [km]')
            ax.yaxis.set_ticks([100,200,300])
            ax.set_ylim((85,320))
            ax.set_yscale('log')
            ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
            ax.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
            ax.set_xlim((-12,12))
            ax.set_xticks(np.arange(-12,13,3))
            ax.set_xticklabels(['%.0f' % np.mod(x,24.) for x in ax.get_xticks()])
            ax.set_xlabel('Solar LT [hr]')


        plt.tight_layout(rect=[0,0.02,1,0.96]) # leave room for suptitle
        # Create date and time string for this plot
        t0 = L21_dict['time'][i1]
        t1 = L21_dict['time'][i2-1]
        t_mid = t0 + timedelta(seconds=(t1-t0).total_seconds()/2.)
        datetime_str = t_mid.strftime('%Y-%m-%d-%H%M%S')

        fig.suptitle('%s MIGHTI-%s %s' % (datetime_str, sensor, L21_dict['emission_color'].capitalize()), fontsize=28)

        #### Save
        versrev = L21_fn.split('/')[-1].split('_')[-1].split('.')[0] # e.g., v01r001
        desc = L21_fn.split('/')[-1].split('_')[-3] #e.g., Vector-Wind-Green
        pngfn = pngpath + 'ICON_L2-1_MIGHTI-%s_Plot-%s_%s_%s.PNG'%(sensor, desc, datetime_str, versrev)
        plt.savefig(pngfn, dpi=120)
        L21_pngs.append(pngfn)
        if close: 
            # There seems to be a memory leak here, despite this:
            fig.clf()
            plt.close(fig)
            plt.close('all')
            
    # Be extra careful for memory leaks
    del L21_dict

    return L21_pngs





def plot_level22(L22_fn, pngpath, v_max = 200., ve_min = 1., ve_max = 100., 
                 a_min = 1., a_max = 1000., ae_min = 0.1, ae_max = 100., close=True,
                 first_only=False):
    '''
    Create Tohban plots for Level 2.2 data. The nominal L2.2 file contains 24 hours of data. The 
    Tohban plots are split up by solar local time. There is one plot per noon-to-noon LT period.
    
    INPUTS:
    
      *  L22_fn       --TYPE:str,   Full path to a L2.2 file
      *  pngpath      --TYPE:str,   Directory to save the resulting png(s) to.
   
    OPTIONAL INPUTS:
   
      *  v_max        --TYPE:float, Maximum wind velocity for colorbar [m/s]
      *  ve_min       --TYPE:float, Minimum wind error for colorbar [m/s]
      *  ve_max       --TYPE:float, Maximum wind error for colorbar [m/s]
      *  a_min        --TYPE:float, Minimum VER for colorbar [ph/cm^3/s]
      *  a_max        --TYPE:float, Maximum VER for colorbar [ph/cm^3/s]
      *  ae_min       --TYPE:float, Minimum VER error for colorbar [ph/cm^3/s]
      *  ae_max       --TYPE:float, Maximum VER error for colorbar [ph/cm^3/s]
      *  close        --TYPE:bool,  If True, close the figure after saving it.
      *  first_only   --TYPE:bool,  If True, only make the first plot (useful for debugging).
                                   
    OUTPUT:
    
      *  L22_pngs     --TYPE:list of str,  Full path to the saved png files
      
    '''
    from mpl_toolkits.basemap import Basemap # For putting map on Tohban plots
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.colors import LinearSegmentedColormap

    # Custom colormap for quality factor/flags
    low_rgb = [0.7, 0.7, 1.0]
    #mid_rgb = np.array([252,238,170])/256.
    mid_rgb = np.array([253,185,40])/256.
    hi_rgb  = [0.7, 0.0, 0.0]
    cdict = {'red':   ((0.0, low_rgb[0], low_rgb[0]),
                       (0.5, mid_rgb[0], mid_rgb[0]),
                       (1.0, hi_rgb [0], hi_rgb [0])),
             'green': ((0.0, low_rgb[1], low_rgb[1]),
                       (0.5, mid_rgb[1], mid_rgb[1]),
                       (1.0, hi_rgb [1], hi_rgb [1])),
             'blue':  ((0.0, low_rgb[2], low_rgb[2]),
                       (0.5, mid_rgb[2], mid_rgb[2]),
                       (1.0, hi_rgb [2], hi_rgb [2]))}
    cmap_byr = LinearSegmentedColormap('Custom', cdict)

    if pngpath[-1] != '/':
        pngpath += '/'

    L22_dict = level22_to_dict(L22_fn)

    Nalt, Nlon = np.shape(L22_dict['u'])
    lonu = L22_dict['lon_unwrapped']
    alt = np.ma.filled(L22_dict['alt'], np.nan) # convert from masked array to array
    t = L22_dict['epoch']
    
    L22_pngs = []
    nrows, ncols = 6, 4
    csize = '3%' # For adding colorbar
    cpad = 0.08  # For adding colorbar
    clabelpad = 15 # For colorbar label

    # Load magnetic coordinates
    lonlat_m10, lonlat_0, lonlat_p10 = mag_lines()

    assert any(t), "All data in the L2.2 are masked. No plot can be created."

    # Organize figures by SLT. Create one plot for each (SLT) noon-to-noon period.
    sltu_2D = unwrap_slt(t,L22_dict['slt'])
    sltu = np.nanmedian(sltu_2D, axis=0) # Reference SLT for each column, median to avoid discontinuities
    n = np.around(sltu/24.) # for each element, the figure it should go in

    # Create istart, istop (i.e., deciding which columns go in which figures)
    istart = []
    istop = []
    for ii in np.unique(n[~np.isnan(n)]):
        idx = np.where(n == ii)[0]
        istart.append(idx[0])
        istop.append(idx[-1]+1)
    nfigs = len(istart)

    if first_only:
        istop = istop[0:1]
        istart = istart[0:1]
        nfigs = 1

    # Loop and create figures
    for n, (i1, i2) in enumerate(zip(istart, istop)):

        if i2-i1 <= 1: # corner case: skip this plot
            continue

        # pcolormesh wants the coordinates of the corners (not the middle) of each pixel
        x = np.mod(sltu[i1:i2] + 12., 24) - 12.
        x = fill(x) # (SLT) Fill in nans in the middle with linear interpolation
        y = alt
        nx = len(x)
        ny = len(y)
        dx = np.diff(x)
        dy = np.diff(y)
        xo = np.concatenate([[x[0]-dx[0]/2], x[:-1] + dx/2, [x[-1]+dx[-1]/2]])
        yo = np.concatenate([[y[0]-dy[0]/2], y[:-1] + dy/2, [y[-1]+dy[-1]/2]])
        X,Y = np.meshgrid(x,y)

        fig, axarr = plt.subplots(nrows, ncols, figsize=(7.3*ncols,3*nrows))
        ax_format = [] # Will hold axes containing the normal Altitude vs. Longitude vs. Color plots, for later common formatting

        ########### Wind and wind error ###########
        ax = axarr[0,0]
        h = ax.pcolormesh(X, Y, L22_dict['u'][:,i1:i2], cmap='bwr', vmin=-v_max, vmax=v_max)
        ax.set_title('Zonal Wind')        
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('m/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ax = axarr[1,0]
        h = ax.pcolormesh(X, Y, L22_dict['v'][:,i1:i2], cmap='bwr', vmin=-v_max, vmax=v_max)
        ax.set_title('Meridional Wind')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('m/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ax = axarr[0,1]
        h = ax.pcolormesh(X, Y, L22_dict['u_error'][:,i1:i2], cmap='viridis_r', norm=LogNorm(vmin=ve_min, vmax=ve_max))
        ax.set_title('Zonal Wind Error')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('m/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ax = axarr[1,1]
        h = ax.pcolormesh(X, Y, L22_dict['v_error'][:,i1:i2], cmap='viridis_r', norm=LogNorm(vmin=ve_min, vmax=ve_max))
        ax.set_title('Meridional Wind Error')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('m/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ############ VER/Amplitude and error ############
        ax = axarr[2,0]
        h = ax.pcolormesh(X, Y, L22_dict['ver'][:,i1:i2], cmap='viridis', norm=LogNorm(vmin=a_min, vmax=a_max))
        ax.set_title('Relative VER')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('ph/cm$^3$/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)

        ax = axarr[2,1]
        h = ax.pcolormesh(X, Y, L22_dict['ver_error'][:,i1:i2], cmap='viridis', norm=LogNorm(vmin=ae_min, vmax=ae_max))
        ax.set_title('Relative VER Error')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        fig.colorbar(h, cax=cax).set_label('ph/cm$^3$/s',rotation=-90,labelpad=clabelpad)
        ax_format.append(ax)


        ############ Quality Factors #############
        cmap = cmap_byr
        bounds = [-0.33,0.33,0.67,1.33]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)  

        ax = axarr[3,0]
        h = ax.pcolormesh(X, Y, 1 - L22_dict['wind_quality'][:,i1:i2], cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Wind Quality')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.5,1])
        cb.ax.set_yticklabels(['Good','Caution','Bad'], rotation=-90, va='center')
        ax_format.append(ax)

        ax = axarr[3,1]
        h = ax.pcolormesh(X, Y, 1 - L22_dict['ver_quality'][:,i1:i2], cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('VER Quality')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.5,1])
        cb.ax.set_yticklabels(['Good','Caution','Bad'], rotation=-90, va='center')
        ax_format.append(ax)

        ############ Location #############
        try: # because I'm worried something funny may happen with this hacky code
            ax = axarr[5,0]
            j = int(Nalt/2) # Use middle altitude for reference lat/lon
            # Annoying code to work longitude ranges for Basemap [-360, 720]
            lonoff = lonu[j,i1] - (np.mod(lonu[j,i1] + 180., 360.)-180.) # offset to put min lon in [-180, 180]
            lonplot = lonu[j,i1:i2] - lonoff
            latplot = L22_dict['lat'][j,i1:i2]
            # Compute lonmin and lonmax for plot. This is an unfortunately hacky way to ensure
            # that the x axis on the map matches the x axis on the other plots, even for 
            # partial orbits.
            c = np.polyfit(x, lonplot, 1)
            dlon_dslt = c[0]
            lon_slt0 = c[1]
            lonmin = dlon_dslt*-12 + lon_slt0
            lonmax = dlon_dslt*12  + lon_slt0
            m = Basemap(ax=ax, llcrnrlon=lonmin, urcrnrlon=lonmax, llcrnrlat=-50., urcrnrlat=50.,\
                        resolution='l',projection='merc',lat_ts=0.,\
                        fix_aspect=False)
            m.plot(lonplot, latplot, 'C1-', lw=3, latlon=True, label='Observations')
            xm,ym = m(lonlat_0[:,0], lonlat_0[:,1])
            ax.plot(xm,ym, 'k-', lw=1, label='MLAT=0')
            xm,ym = m(lonlat_m10[:,0], lonlat_m10[:,1])
            ax.plot(xm,ym, 'k--', lw=1, label='MLAT=+/-10 deg')
            xm,ym = m(lonlat_p10[:,0], lonlat_p10[:,1])
            ax.plot(xm,ym, 'k--', lw=1)
            cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad) # Helps alignment
            cax.set_visible(False) # Helps alignment
            m.drawparallels(np.arange(-90.,91.,30.), labels=[0,1,0,1]); # LRTB
            m.drawmeridians(np.arange(-180.,361.,90.), labels=[0,1,0,1]);
            ax.legend(loc='lower left', prop={'size':8})
            m.drawcoastlines(linewidth=0.5);
            m.fillcontinents()
        except Exception as e:
            print('Error creating map: %s' % e)

        ############ Quality flags derived at L2.2 ############
        # There isn't enough room to plot all the flags, so we'll plot the most important ones (subject to change)
        cmap = cmap_byr
        bounds = [-0.5,0.5,1.5]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)  

        ax = axarr[0,2]
        h = ax.pcolormesh(X, Y, L22_dict['quality_flags'][:,i1:i2,28], cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Spherical asymmetry detected')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)

        ax = axarr[4,0]
        h = ax.pcolormesh(X, Y, L22_dict['quality_flags'][:,i1:i2,33], cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('Unknown Error')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,1])
        ax_format.append(ax)

        ############ Quality Flags derived at L1 and L2.1 (separate MIGHTI A and B) #############
        cmap = cmap_byr
        bounds = [-0.167, 0.167, 0.5, 0.833, 1.167]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N) 

        ax = axarr[1,2]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,24] + 0.667*L22_dict['quality_flags'][:,i1:i2,25]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Missing file')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)

        ax = axarr[2,2]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,0] + 0.667*L22_dict['quality_flags'][:,i1:i2,12]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Very low SNR before inversion')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)

        ax = axarr[3,2]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,6] + 0.667*L22_dict['quality_flags'][:,i1:i2,18]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Very low SNR after inversion')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)
        
        ax = axarr[4,2]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,11] + 0.667*L22_dict['quality_flags'][:,i1:i2,23]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Cautiously low SNR after inversion')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)

        ax = axarr[4,1]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,7] + 0.667*L22_dict['quality_flags'][:,i1:i2,19]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Airglow above 300 km')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)

        ax = axarr[0,3]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,8] + 0.667*L22_dict['quality_flags'][:,i1:i2,20]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: LoS crosses terminator')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)

        ax = axarr[1,3]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,2] + 0.667*L22_dict['quality_flags'][:,i1:i2,14]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Bad calibration')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)

        ax = axarr[2,3]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,26] + 0.667*L22_dict['quality_flags'][:,i1:i2,27]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Unobserved altitudes')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)

        # 2022 Jul 25: Replaced plot of SAA with plot of "not enough valid points in profile"
        ax = axarr[5,2]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,5] + 0.667*L22_dict['quality_flags'][:,i1:i2,17]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Not enough valid points in profile')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)
        
        ax = axarr[3,3]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,9] + 0.667*L22_dict['quality_flags'][:,i1:i2,21]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Thermal drift correction uncertain')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)
        
        # 2021 Feb 8: replaced plot of S/C jitter flag with plot of solar/lunar contamination
        ax = axarr[4,3]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,4] + 0.667*L22_dict['quality_flags'][:,i1:i2,16]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Solar/lunar contamination')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)
        
        ax = axarr[5,3]
        z = 0.333*L22_dict['quality_flags'][:,i1:i2,3] + 0.667*L22_dict['quality_flags'][:,i1:i2,15]
        h = ax.pcolormesh(X, Y, z, cmap=cmap, norm=norm, vmin=0.0, vmax=1.0)
        ax.set_title('MIGHTI A or B: Calibration lamps are on')
        cax = make_axes_locatable(ax).append_axes('right', size=csize, pad=cpad)
        cb = fig.colorbar(h, cax=cax, ticks=[0,0.333,0.667,1])
        cb.ax.set_yticklabels(['Neither','A','B','Both'], rotation=-90, va='center')
        ax_format.append(ax)

        #### Turn off unused axes
        axarr[5,1].axis('off')

        #### Make the 2D plots look nice
        for ax in ax_format:
            plt.sca(ax)
            ax.patch.set_facecolor('gray')
            ax.set_ylabel('Altitude [km]')
            ax.yaxis.set_ticks([100,200,300])
            ax.set_ylim((85,320))
            ax.set_yscale('log')
            ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
            ax.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
            ax.set_xlim((-12,12))
            ax.set_xticks(np.arange(-12,13,3))
            ax.set_xticklabels(['%.0f' % np.mod(x,24.) for x in ax.get_xticks()])
            ax.set_xlabel('Solar LT [hr]')

        plt.tight_layout(rect=[0,0.02,1,0.96]) # leave room for suptitle
        # Create date and time string for this plot
        t0 = L22_dict['epoch'][i1]
        t1 = L22_dict['epoch'][i2-1]
        t_mid = t0 + timedelta(seconds=(t1-t0).total_seconds()/2.)
        datetime_str = t_mid.strftime('%Y-%m-%d-%H%M%S')

        fig.suptitle('%s %s' % (datetime_str, L22_dict['emission_color']), fontsize=28)

        #### Save
        versrev = L22_fn.split('/')[-1].split('_')[-1].split('.')[0] # e.g., v01r001
        desc = L22_fn.split('/')[-1].split('_')[-3] #e.g., Vector-Wind-Green
        pngfn = pngpath + 'ICON_L2-2_MIGHTI_Plot-%s_%s_%s.PNG'%(desc,datetime_str,versrev)
        plt.savefig(pngfn, dpi=120)
        L22_pngs.append(pngfn)
        if close: 
            # There seems to be a memory leak here, despite this:
            fig.clf()
            plt.close(fig)
            plt.close('all')
    # Be extra careful for memory leaks
    del L22_dict

    return L22_pngs



def plot_level22_summary(L22_fng, L22_fnr, pngpath, val, vmax = 150., close=False):
    '''
    Create an overall summary plot for Level 2.2 data, combining red and green channels.
    
    INPUTS:
    
      *  L22_fng      --TYPE:str,   Full path to a L2.2 file (green channel)
      *  L22_fnr      --TYPE:str,   Full path to a L2.2 file (red channel)
      *  pngpath      --TYPE:str,   Directory to save the resulting png to.  If None, figure won't be saved.
      *  val          --TYPE:str,   'u' or 'v' for zonal wind or meridional wind
   
    OPTIONAL INPUTS:
    
      *  v_max        --TYPE:float, Maximum wind velocity for colorbar [m/s]
      *  close        --TYPE:bool,  If True, close the figure after saving it.    

                                   
    OUTPUT:
    
      *  L22_png      --TYPE:str,  Full path to the saved png file
      
    '''
    try:
        assert val in ['u','v'], 'Input "val" is "%s" but must be "u" or "v"' % val
        zon_mer = 'Zonal'
        if val=='v':
            zon_mer = 'Meridional'

        dg = level22_to_dict(L22_fng)
        dr = level22_to_dict(L22_fnr)

        # Quality control and other things
        for d in [dr, dg]:
            idx = d['wind_quality'] < 0.5
            if not np.isscalar(d['u'].mask): # If mask is not correct then skip this step. (Useful for TIEGCM or WACCM-X sampled files)
                d['u'].mask[idx] = True
                d['v'].mask[idx] = True

        # 2022 Oct 31: Check if red and green orbits are out of sync. This is very rare and only happens if the UT switch happens 
        # when ICON is crossing lon = 0.
        if dr['lon_unwrapped'][0,0] - dg['lon_unwrapped'][0,0] > 180.:
            dr['lon_unwrapped'] -= 360.
        elif dr['lon_unwrapped'][0,0] - dg['lon_unwrapped'][0,0] < -180.:
            dr['lon_unwrapped'] += 360.

        ### Process dr and dg to add new SLT-based orbit number
        for d in [dr, dg]:

            # Which indices to plot -- recalculate approximate SLT here so there are no gaps
            Nt = len(d['epoch'])
            lonu = d['lon_unwrapped'][0,:]
            t0 = pd.to_datetime(d['epoch'][-1]).date()
            thr = (pd.to_datetime(d['epoch']) - pd.to_datetime(t0)).total_seconds()/3600.
            sltu = thr + 24./360.*lonu
            orb_num_new = np.int32(np.floor(sltu/24.))
            slt = np.mod(sltu, 24.)
            d['Norbs'] = max(orb_num_new)
            d['slt_new'] = slt
            d['orbit_number_new'] = orb_num_new
        # 2022 Oct 31: Handle the very rare case where one sample is different in the red and green
        # resulting in different number of orbits. Use the maximum to be able to plot all the data.
        # This only happened a couple times in 3 years.
        Norbs = max(dg['Norbs'],dr['Norbs'])

        # Put in jmin and jmax here
        dr['jmin'] = 1
        dr['jmax'] = 16
        dg['jmin'] = 0
        dg['jmax'] = 39


        from matplotlib.ticker import NullFormatter # https://github.com/matplotlib/matplotlib/issues/8027/
        from matplotlib.ticker import FixedLocator, FixedFormatter # https://brohrer.github.io/matplotlib_ticks.html

        nx = 25 # Number of grids wide

        fig = plt.figure(figsize=(6,10), dpi = 120)

        for n in range(Norbs):

            plt.subplot2grid((Norbs+1,nx), (n,0), rowspan=1, colspan=nx-1)
            ax = plt.gca()

            for d in [dg, dr]:
                jmin = d['jmin']
                jmax = d['jmax']
                slt = d['slt_new']
                im = d['orbit_number_new'] == n
                if im.any(): # Don't try to plot something if there is no data (It's very rare to have this check fail)
                    m = ax.pcolormesh(slt[im], d['alt'][jmin:jmax], d[val][jmin:jmax,im], cmap='RdBu_r', vmin=-vmax, vmax=vmax)


            plt.subplots_adjust(hspace = .000)
            ax.set_ylim((90,350))
            ax.set_yscale('log')
            y_formatter = FixedFormatter(["90", "150", "300"])
            y_locator = FixedLocator([90, 150, 300])
            ax.yaxis.set_major_formatter(y_formatter)
            ax.yaxis.set_major_locator(y_locator)
            ax.yaxis.set_minor_formatter(NullFormatter())
            if n==0:
                plt.title('%s\n%s Wind'%(t0.strftime('%Y-%m-%d'), zon_mer))
                ax.set_ylabel('Height\n[km]')  
            else:
                ax.set_yticks([])

            ax.set_xlim((0, 24))
            ax.set_facecolor('gray')
            ax.set_xticks([])

        ##### Colorbar
        plt.subplot2grid((Norbs+1,nx), (0,nx-1), rowspan=Norbs, colspan=1)
        cax = plt.gca()
        cb = plt.colorbar(mappable=m, cax=cax, ax=None)
        cb.set_label('m/s', labelpad=-8, rotation=-90)

        ##### Latitude
        plt.subplot2grid((Norbs+1,nx), (Norbs,0), rowspan=1, colspan=nx-1)
        ax = plt.gca()
        plt.plot(d['slt_new'], d['lat'].T, 'k.', ms=1)
        plt.ylim((-18,48))
        plt.yticks([-15,0,15,30,45])
        ax.set_ylabel('Latitude\n[deg]')
        ax.set_xticks([0,6,12,18,24])
        ax.set_xlim((0, 24))
        ax.set_xlabel('Local Solar Time [hr]')    
        
        # Be extra careful for memory leaks
        del d, dr, dg
        
    except Exception as e:
        fig = plt.figure(figsize=(6,10), dpi = 120)
        plt.text(0.5, 0.5, 'Summary plots failed, but the\nunderlying .NC files\nmight be okay.\nContact MIGHTI Team\nbefore approving', va='center',ha='center', fontsize=16)
        plt.ylim((0,1))
        plt.xlim((0,1))
        raise
        
        
    #### Save
    pngfn = None
    if pngpath is not None:
        versrev = L22_fng.split('/')[-1].split('_')[-1].split('.')[0] # e.g., v01r001
        desc = L22_fng.split('/')[-1].split('_')[-3] #e.g., Vector-Wind-Green
        pngfn = pngpath + 'ICON_L2-2_MIGHTI_Plot-Vector-Wind-Combined-%s_%s_%s.PNG'%(zon_mer, t0.strftime('%Y-%m-%d'), versrev)
        plt.savefig(pngfn, dpi=120)
    
    if close: 
        # There seems to be a memory leak here, despite this:
        fig.clf()
        plt.close(fig)
        plt.close('all')

    
    return pngfn


def mag_lines():
    '''
    Helper function to define magnetic latitude -10, 0, and 10 degrees.
    
    OUTPUTS:
    
    lonlat_m10 -- TYPE:array(n,2), UNITS:deg. Geographic Longitude and Latitude of MLAT=-10
    lonlat_0   -- TYPE:array(n,2), UNITS:deg. Geographic Longitude and Latitude of MLAT=0
    lonlat_p10 -- TYPE:array(n,2), UNITS:deg. Geographic Longitude and Latitude of MLAT=+10
    '''
    
    lonlat_m10 = np.array([[-1.80000000e+02, -6.89663865e+00],
       [-1.73934866e+02, -7.89473684e+00],
       [-1.72727273e+02, -8.09094760e+00],
       [-1.65454545e+02, -9.18946721e+00],
       [-1.58181818e+02, -1.02280458e+01],
       [-1.50909091e+02, -1.12260749e+01],
       [-1.43636364e+02, -1.21902194e+01],
       [-1.36363636e+02, -1.31261068e+01],
       [-1.36114617e+02, -1.31578947e+01],
       [-1.29090909e+02, -1.40714935e+01],
       [-1.21818182e+02, -1.50878128e+01],
       [-1.14545455e+02, -1.62505655e+01],
       [-1.07272727e+02, -1.75944366e+01],
       [-1.03152466e+02, -1.84210526e+01],
       [-1.00000000e+02, -1.90623662e+01],
       [-9.27272727e+01, -2.04554277e+01],
       [-8.54545455e+01, -2.14698205e+01],
       [-7.81818182e+01, -2.18049765e+01],
       [-7.09090909e+01, -2.11985346e+01],
       [-6.36363636e+01, -1.94469971e+01],
       [-6.10635145e+01, -1.84210526e+01],
       [-5.63636364e+01, -1.64569984e+01],
       [-5.04255273e+01, -1.31578947e+01],
       [-4.90909091e+01, -1.23970405e+01],
       [-4.18181818e+01, -7.89982720e+00],
       [-4.18093917e+01, -7.89473684e+00],
       [-3.45454545e+01, -3.90585729e+00],
       [-3.16027050e+01, -2.63157895e+00],
       [-2.72727273e+01, -9.11362092e-01],
       [-2.00000000e+01,  1.01222503e+00],
       [-1.27272727e+01,  2.06884921e+00],
       [-5.45454545e+00,  2.48026461e+00],
       [ 1.81818182e+00,  2.44923610e+00],
       [ 9.09090909e+00,  2.10314163e+00],
       [ 1.63636364e+01,  1.50189804e+00],
       [ 2.36363636e+01,  6.88556237e-01],
       [ 3.09090909e+01, -2.56265826e-01],
       [ 3.81818182e+01, -1.19303805e+00],
       [ 4.54545455e+01, -1.95113074e+00],
       [ 5.27272727e+01, -2.39259438e+00],
       [ 6.00000000e+01, -2.47651452e+00],
       [ 6.72727273e+01, -2.28518665e+00],
       [ 7.45454545e+01, -1.99157011e+00],
       [ 8.18181818e+01, -1.77998248e+00],
       [ 8.90909091e+01, -1.76111782e+00],
       [ 9.63636364e+01, -1.92681549e+00],
       [ 1.03636364e+02, -2.16655645e+00],
       [ 1.10909091e+02, -2.33444107e+00],
       [ 1.18181818e+02, -2.33197268e+00],
       [ 1.25454545e+02, -2.16680375e+00],
       [ 1.32727273e+02, -1.95996319e+00],
       [ 1.40000000e+02, -1.89764262e+00],
       [ 1.47272727e+02, -2.15012851e+00],
       [ 1.52687138e+02, -2.63157895e+00],
       [ 1.54545455e+02, -2.79630350e+00],
       [ 1.61818182e+02, -3.78633647e+00],
       [ 1.69090909e+02, -4.99365044e+00],
       [ 1.76363636e+02, -6.26923698e+00],
       [ 1.83636364e+02, -7.50601752e+00],
       [ 1.86091614e+02, -7.89473684e+00],
       [ 1.90909091e+02, -8.64922866e+00],
       [ 1.98181818e+02, -9.71476394e+00],
       [ 2.05454545e+02, -1.07314437e+01],
       [ 2.12727273e+02, -1.17122790e+01],
       [ 2.20000000e+02, -1.26607199e+01],
       [ 2.23902453e+02, -1.31578947e+01],
       [ 2.27272727e+02, -1.35945186e+01],
       [ 2.34545455e+02, -1.45658774e+01],
       [ 2.41818182e+02, -1.56469785e+01],
       [ 2.49090909e+02, -1.69010974e+01],
       [ 2.56363636e+02, -1.83184386e+01],
       [ 2.56876685e+02, -1.84210526e+01],
       [ 2.63636364e+02, -1.97870406e+01],
       [ 2.70909091e+02, -2.10292572e+01],
       [ 2.78181818e+02, -2.17399031e+01],
       [ 2.85454545e+02, -2.16337604e+01],
       [ 2.92727273e+02, -2.04756936e+01],
       [ 2.99061810e+02, -1.84210526e+01],
       [ 3.00000000e+02, -1.81032581e+01],
       [ 3.07272727e+02, -1.45325622e+01],
       [ 3.09584369e+02, -1.31578947e+01],
       [ 3.14545455e+02, -1.01639914e+01],
       [ 3.18274624e+02, -7.89473684e+00],
       [ 3.21818182e+02, -5.81797757e+00],
       [ 3.28247928e+02, -2.63157895e+00],
       [ 3.29090909e+02, -2.24527946e+00],
       [ 3.36363636e+02,  1.70549785e-01],
       [ 3.43636364e+02,  1.63571825e+00],
       [ 3.50909091e+02,  2.34109777e+00],
       [ 3.58181818e+02,  2.51013362e+00],
       [ 3.65454545e+02,  2.31071579e+00],
       [ 3.72727273e+02,  1.83205057e+00],
       [ 3.80000000e+02,  1.11806273e+00],
       [ 3.87272727e+02,  2.25153663e-01],
       [ 3.94545455e+02, -7.36236307e-01],
       [ 4.01818182e+02, -1.60466114e+00],
       [ 4.09090909e+02, -2.21684555e+00],
       [ 4.16363636e+02, -2.47690612e+00],
       [ 4.23636364e+02, -2.40579198e+00],
       [ 4.30909091e+02, -2.13881868e+00],
       [ 4.38181818e+02, -1.86606516e+00],
       [ 4.45454545e+02, -1.74406975e+00],
       [ 4.52727273e+02, -1.82601168e+00],
       [ 4.60000000e+02, -2.04671554e+00],
       [ 4.67272727e+02, -2.26764156e+00],
       [ 4.74545455e+02, -2.35688517e+00],
       [ 4.81818182e+02, -2.26449603e+00],
       [ 4.89090909e+02, -2.05760473e+00],
       [ 4.96363636e+02, -1.89871579e+00],
       [ 5.03636364e+02, -1.97676564e+00],
       [ 5.10909091e+02, -2.42435614e+00],
       [ 5.12715523e+02, -2.63157895e+00],
       [ 5.18181818e+02, -3.25468679e+00],
       [ 5.25454545e+02, -4.37257051e+00],
       [ 5.32727273e+02, -5.63105858e+00],
       [ 5.40000000e+02, -6.89663865e+00]])
    
    lonlat_0 = np.array([[-1.80000000e+02,  2.90487384e+00],
       [-1.78373570e+02,  2.63157895e+00],
       [-1.72727273e+02,  1.70386852e+00],
       [-1.65454545e+02,  6.10021466e-01],
       [-1.58181818e+02, -4.05664667e-01],
       [-1.50909091e+02, -1.36953882e+00],
       [-1.43636364e+02, -2.29116986e+00],
       [-1.40831637e+02, -2.63157895e+00],
       [-1.36363636e+02, -3.17163838e+00],
       [-1.29090909e+02, -4.04274534e+00],
       [-1.21818182e+02, -4.97019120e+00],
       [-1.14545455e+02, -6.04013854e+00],
       [-1.07272727e+02, -7.30729788e+00],
       [-1.04255692e+02, -7.89473684e+00],
       [-1.00000000e+02, -8.73680021e+00],
       [-9.27272727e+01, -1.01239583e+01],
       [-8.54545455e+01, -1.11316219e+01],
       [-7.81818182e+01, -1.13846254e+01],
       [-7.09090909e+01, -1.05470471e+01],
       [-6.36363636e+01, -8.40225687e+00],
       [-6.25317751e+01, -7.89473684e+00],
       [-5.63636364e+01, -4.99024863e+00],
       [-5.23132529e+01, -2.63157895e+00],
       [-4.90909091e+01, -7.81479445e-01],
       [-4.32073017e+01,  2.63157895e+00],
       [-4.18181818e+01,  3.39111761e+00],
       [-3.45454545e+01,  6.79374775e+00],
       [-3.13756709e+01,  7.89473684e+00],
       [-2.72727273e+01,  9.20476725e+00],
       [-2.00000000e+01,  1.06952133e+01],
       [-1.27272727e+01,  1.14788784e+01],
       [-5.45454545e+00,  1.17397502e+01],
       [ 1.81818182e+00,  1.16297846e+01],
       [ 9.09090909e+00,  1.12411799e+01],
       [ 1.63636364e+01,  1.06185188e+01],
       [ 2.36363636e+01,  9.79950626e+00],
       [ 3.09090909e+01,  8.85884952e+00],
       [ 3.81818182e+01,  7.92778763e+00],
       [ 3.85017458e+01,  7.89473684e+00],
       [ 4.54545455e+01,  7.17161638e+00],
       [ 5.27272727e+01,  6.73081052e+00],
       [ 6.00000000e+01,  6.64065813e+00],
       [ 6.72727273e+01,  6.81625386e+00],
       [ 7.45454545e+01,  7.08868574e+00],
       [ 8.18181818e+01,  7.28449618e+00],
       [ 8.90909091e+01,  7.30464741e+00],
       [ 9.63636364e+01,  7.16326920e+00],
       [ 1.03636364e+02,  6.97087634e+00],
       [ 1.10909091e+02,  6.87134012e+00],
       [ 1.18181818e+02,  6.96088319e+00],
       [ 1.25454545e+02,  7.22794149e+00],
       [ 1.32727273e+02,  7.54549472e+00],
       [ 1.40000000e+02,  7.72047750e+00],
       [ 1.47272727e+02,  7.57424174e+00],
       [ 1.54545455e+02,  7.01451596e+00],
       [ 1.61818182e+02,  6.06754373e+00],
       [ 1.69090909e+02,  4.85844917e+00],
       [ 1.76363636e+02,  3.55102938e+00],
       [ 1.81627680e+02,  2.63157895e+00],
       [ 1.83636364e+02,  2.28823630e+00],
       [ 1.90909091e+02,  1.14526538e+00],
       [ 1.98181818e+02,  9.43432514e-02],
       [ 2.05454545e+02, -8.93050892e-01],
       [ 2.12727273e+02, -1.83562355e+00],
       [ 2.19152899e+02, -2.63157895e+00],
       [ 2.20000000e+02, -2.73585842e+00],
       [ 2.27272727e+02, -3.60503288e+00],
       [ 2.34545455e+02, -4.49409831e+00],
       [ 2.41818182e+02, -5.48236814e+00],
       [ 2.49090909e+02, -6.64883482e+00],
       [ 2.55768931e+02, -7.89473684e+00],
       [ 2.56363636e+02, -8.00762700e+00],
       [ 2.63636364e+02, -9.45515866e+00],
       [ 2.70909091e+02, -1.06986272e+01],
       [ 2.78181818e+02, -1.13753096e+01],
       [ 2.85454545e+02, -1.11194654e+01],
       [ 2.92727273e+02, -9.64455314e+00],
       [ 2.97308313e+02, -7.89473684e+00],
       [ 3.00000000e+02, -6.83629819e+00],
       [ 3.07272727e+02, -2.92306400e+00],
       [ 3.07754507e+02, -2.63157895e+00],
       [ 3.14545455e+02,  1.36354360e+00],
       [ 3.16840647e+02,  2.63157895e+00],
       [ 3.21818182e+02,  5.19666074e+00],
       [ 3.28428178e+02,  7.89473684e+00],
       [ 3.29090909e+02,  8.14330208e+00],
       [ 3.36363636e+02,  1.00495545e+01],
       [ 3.43636364e+02,  1.11637361e+01],
       [ 3.50909091e+02,  1.16638939e+01],
       [ 3.58181818e+02,  1.17239142e+01],
       [ 3.65454545e+02,  1.14668000e+01],
       [ 3.72727273e+02,  1.09571595e+01],
       [ 3.80000000e+02,  1.02301614e+01],
       [ 3.87272727e+02,  9.33746819e+00],
       [ 3.94545455e+02,  8.38201285e+00],
       [ 3.98657954e+02,  7.89473684e+00],
       [ 4.01818182e+02,  7.51691930e+00],
       [ 4.09090909e+02,  6.90677022e+00],
       [ 4.16363636e+02,  6.64452154e+00],
       [ 4.23636364e+02,  6.70463879e+00],
       [ 4.30909091e+02,  6.95214569e+00],
       [ 4.38181818e+02,  7.20482002e+00],
       [ 4.45454545e+02,  7.31837321e+00],
       [ 4.52727273e+02,  7.24893720e+00],
       [ 4.60000000e+02,  7.06428747e+00],
       [ 4.67272727e+02,  6.90146648e+00],
       [ 4.74545455e+02,  6.89027940e+00],
       [ 4.81818182e+02,  7.07779225e+00],
       [ 4.89090909e+02,  7.39181165e+00],
       [ 4.96363636e+02,  7.66328961e+00],
       [ 5.03636364e+02,  7.69592668e+00],
       [ 5.10909091e+02,  7.34722401e+00],
       [ 5.18181818e+02,  6.58336382e+00],
       [ 5.25454545e+02,  5.48556896e+00],
       [ 5.32727273e+02,  4.20727716e+00],
       [ 5.40000000e+02,  2.90487384e+00]])
    
    lonlat_p10 = np.array([[-1.80000000e+02,  1.31806916e+01],
       [-1.79866290e+02,  1.31578947e+01],
       [-1.72727273e+02,  1.19606722e+01],
       [-1.65454545e+02,  1.08223414e+01],
       [-1.58181818e+02,  9.75361312e+00],
       [-1.50909091e+02,  8.73139526e+00],
       [-1.44727125e+02,  7.89473684e+00],
       [-1.43636364e+02,  7.75040342e+00],
       [-1.36363636e+02,  6.81799221e+00],
       [-1.29090909e+02,  5.90682935e+00],
       [-1.21818182e+02,  4.96470576e+00],
       [-1.14545455e+02,  3.91187988e+00],
       [-1.07272727e+02,  2.69122434e+00],
       [-1.06950209e+02,  2.63157895e+00],
       [-1.00000000e+02,  1.32875181e+00],
       [-9.27272727e+01,  2.54727232e-02],
       [-8.54545455e+01, -8.77955905e-01],
       [-7.81818182e+01, -9.65469003e-01],
       [-7.09090909e+01,  1.26827273e-01],
       [-6.36363636e+01,  2.55621740e+00],
       [-6.34813079e+01,  2.63157895e+00],
       [-5.63636364e+01,  6.06466795e+00],
       [-5.30595454e+01,  7.89473684e+00],
       [-4.90909091e+01,  9.98325206e+00],
       [-4.27223458e+01,  1.31578947e+01],
       [-4.18181818e+01,  1.35742242e+01],
       [-3.45454545e+01,  1.63544287e+01],
       [-2.72727273e+01,  1.83453108e+01],
       [-2.68434008e+01,  1.84210526e+01],
       [-2.00000000e+01,  1.95634108e+01],
       [-1.27272727e+01,  2.02184891e+01],
       [-5.45454545e+00,  2.04478806e+01],
       [ 1.81818182e+00,  2.03562685e+01],
       [ 9.09090909e+00,  2.00074414e+01],
       [ 1.63636364e+01,  1.94368716e+01],
       [ 2.36363636e+01,  1.86823618e+01],
       [ 2.58551978e+01,  1.84210526e+01],
       [ 3.09090909e+01,  1.78157790e+01],
       [ 3.81818182e+01,  1.69670617e+01],
       [ 4.54545455e+01,  1.62858018e+01],
       [ 5.27272727e+01,  1.58905693e+01],
       [ 6.00000000e+01,  1.58111134e+01],
       [ 6.72727273e+01,  1.59673518e+01],
       [ 7.45454545e+01,  1.62081917e+01],
       [ 8.18181818e+01,  1.63861666e+01],
       [ 8.90909091e+01,  1.64229595e+01],
       [ 9.63636364e+01,  1.63370075e+01],
       [ 1.03636364e+02,  1.62276666e+01],
       [ 1.10909091e+02,  1.62223175e+01],
       [ 1.18181818e+02,  1.64046168e+01],
       [ 1.25454545e+02,  1.67585495e+01],
       [ 1.32727273e+02,  1.71626131e+01],
       [ 1.40000000e+02,  1.74375221e+01],
       [ 1.47272727e+02,  1.74163153e+01],
       [ 1.54545455e+02,  1.70016484e+01],
       [ 1.61818182e+02,  1.61937799e+01],
       [ 1.69090909e+02,  1.50854702e+01],
       [ 1.76363636e+02,  1.38235170e+01],
       [ 1.80174807e+02,  1.31578947e+01],
       [ 1.83636364e+02,  1.25611700e+01],
       [ 1.90909091e+02,  1.13813258e+01],
       [ 1.98181818e+02,  1.02809067e+01],
       [ 2.05454545e+02,  9.23762319e+00],
       [ 2.12727273e+02,  8.23477864e+00],
       [ 2.15281778e+02,  7.89473684e+00],
       [ 2.20000000e+02,  7.27995405e+00],
       [ 2.27272727e+02,  6.36203047e+00],
       [ 2.34545455e+02,  5.44440273e+00],
       [ 2.41818182e+02,  4.45701865e+00],
       [ 2.49090909e+02,  3.32332622e+00],
       [ 2.52968064e+02,  2.63157895e+00],
       [ 2.56363636e+02,  2.01863805e+00],
       [ 2.63636364e+02,  6.51342439e-01],
       [ 2.70909091e+02, -5.01734460e-01],
       [ 2.78181818e+02, -1.04943682e+00],
       [ 2.85454545e+02, -5.83290324e-01],
       [ 2.92727273e+02,  1.17669678e+00],
       [ 2.96225931e+02,  2.63157895e+00],
       [ 3.00000000e+02,  4.20671545e+00],
       [ 3.06980785e+02,  7.89473684e+00],
       [ 3.07272727e+02,  8.04317382e+00],
       [ 3.14545455e+02,  1.18544487e+01],
       [ 3.17357577e+02,  1.31578947e+01],
       [ 3.21818182e+02,  1.50616612e+01],
       [ 3.29090909e+02,  1.74473297e+01],
       [ 3.33444056e+02,  1.84210526e+01],
       [ 3.36363636e+02,  1.90346944e+01],
       [ 3.43636364e+02,  1.99521450e+01],
       [ 3.50909091e+02,  2.03789130e+01],
       [ 3.58181818e+02,  2.04372863e+01],
       [ 3.65454545e+02,  2.02114280e+01],
       [ 3.72727273e+02,  1.97479881e+01],
       [ 3.80000000e+02,  1.90791866e+01],
       [ 3.85837049e+02,  1.84210526e+01],
       [ 3.87272727e+02,  1.82562161e+01],
       [ 3.94545455e+02,  1.73798016e+01],
       [ 4.01818182e+02,  1.65966245e+01],
       [ 4.09090909e+02,  1.60479571e+01],
       [ 4.16363636e+02,  1.58139763e+01],
       [ 4.23636364e+02,  1.58683613e+01],
       [ 4.30909091e+02,  1.60874109e+01],
       [ 4.38181818e+02,  1.63120930e+01],
       [ 4.45454545e+02,  1.64233450e+01],
       [ 4.52727273e+02,  1.63905480e+01],
       [ 4.60000000e+02,  1.62771386e+01],
       [ 4.67272727e+02,  1.62049053e+01],
       [ 4.74545455e+02,  1.62883091e+01],
       [ 4.81818182e+02,  1.65655992e+01],
       [ 4.89090909e+02,  1.69649793e+01],
       [ 4.96363636e+02,  1.73277318e+01],
       [ 5.03636364e+02,  1.74721342e+01],
       [ 5.10909091e+02,  1.72605420e+01],
       [ 5.18181818e+02,  1.66429396e+01],
       [ 5.25454545e+02,  1.56686139e+01],
       [ 5.32727273e+02,  1.44640313e+01],
       [ 5.40000000e+02,  1.31806916e+01]])
    
    return lonlat_m10, lonlat_0, lonlat_p10
    
    
################################################################################################################
##############################   INTERFACE FOR SCIENCE DATA CENTER:    #########################################
################################################################################################################
    
def main(argv):
    '''
    Program that will run if this file is called as a script. Depending on the 
    argument, the Level 1.5, Level 2.1, or Level 2.2 algorithms will be launched, assuming that the Information.TXT file
    is in the directory. Examples:
    
    $ python /path/to/MIGHTI_L2.py -L1.5
    $ python /path/to/MIGHTI_L2.py -L2.1
    $ python /path/to/MIGHTI_L2.py -L2.2
    
    The script will return with exit code 0 on success, exit code 1 on failure. If exit code 1, that usually
    means there's a bug in the code or something completely unexpected in the data.
    
    INPUT:
        argv --TYPE:str. The command line arguments, provided by sys.argv.
        
    '''
    global_params['verbose'] = False # Temporary debugging turned off in v05
    
    info_fn = './Input/Information.TXT'
    
    usagestr = 'usage: python MIGHTI_L2.py -Lxxx  (where xxx is 1.5, 2.1, or 2.2)'
    if len(argv)!=2:
        print(usagestr)
        sys.exit(1)
        
    plt.switch_backend('Agg') # to generate pngs without opening windows
        
    # Call the L2.1 or L2.2 algorithms, as specified.
    # Also possibly call the unit test functions
    err = ''
    if argv[1] == '-L1.5':
        err = level1_to_level15(info_fn)
    elif argv[1] == '-L2.1':
        err = level1_to_level21(info_fn)
    elif argv[1] == '-L2.2':
        err = level21_to_level22(info_fn)
    else:
        print(usagestr)
        sys.exit(1)
    
    # If the processing had a problem, print it and exit with error
    if err:
        print(err)
        sys.exit(1)
    else:
        sys.exit(0)
    
    
if __name__ == "__main__":
    main(sys.argv)
    
    
    
    



    
    
    
    
    
    
    
    
    
    
