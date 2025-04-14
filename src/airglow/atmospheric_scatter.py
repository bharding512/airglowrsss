# A source of functions for radiative transfer and atmospheric scattering calculations.
# Currently used by FPI_atmospheric_scatter_inversion.ipynb

import numpy as np
from numpy import sqrt, exp, arctan2, sin, cos, linspace, zeros, pi, arange, ravel_multi_index, meshgrid
from numpy import arccos, arcsin, sum, reshape, mod, arctan, array, mean, concatenate
import sys
import time as time_mod
from IPython.display import display, clear_output
from scipy import interpolate
import scipy.sparse as sp
from scipy.sparse import linalg as splinalg
from numpy.linalg import norm
from scipy.interpolate import RectBivariateSpline, RectSphereBivariateSpline
import bisect
import warnings

RE = 6371e3
# THIS WAS COMMENTED OUT JUST TO FIX PY3 ERROR. IF USING, UNCOMMENT!
#sec = lambda(x): 1/cos(x)


# Copied from MIGHTI_L2
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
        
    # Special corner case: if x0 is x[-1], return y[-1]
    if x0==x[-1]:
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


def particle_velocity_density(v,u,T,I):
    '''
    A plot over v will show the velocity distribution of particles that are emitting
    airglow. Integrating this distribution will give you I. The particles are from 
    a gas with bulk LoS velocity u and temperature T.
    
    Basically, this function is a forward model for the spectrum in velocity space.
    '''

    c = 299792458.
    k = 1.3806503e-23
    m = 16/6.0221367e26

    # https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    g = I*sqrt(m/(2*pi*k*T))*exp(-m*(v-u)**2 / (2*k*T))
    
    return g




def spectral_intensity(dhf,x,y,v,uvw,T,h):
    '''
    The column brightness in a narrow wavelength (or equivalently, velocity) band.
    dhf gives the column brightness as a function of x,y.
    x and y are arrays which are passed to dhf(x,y).
    v is a scalar representing the velocity to calculate the spectrum at
    uvw is a 3-array of the neutral wind (assumed constant across field)
    T is a scalar of the neutral temperature (assumed constant across field)
    '''
    #th = arctan(sqrt(x**2+y**2)/h) # zenith to points
    # Convert x,y to look angles th, ph using spherical geometry
    rho = sqrt(x**2 + y**2)
    R = sqrt(RE**2 + (RE+h)**2 - 2*RE*(RE+h)*cos(rho/(RE+h)))
    th = arcsin((RE+h)/R * sin(rho/(RE+h)))
    ph = arctan2(x,y) # azimuth to points
    gamma = arcsin(RE/(RE+h) * sin(th)) # zenith angle (slant angle) at pierce point    
    loswind = uvw[0]*sin(gamma)*sin(ph) + uvw[1]*sin(gamma)*cos(ph) + uvw[2]*cos(gamma)
    #print loswind
    
    return particle_velocity_density(v, loswind, T, dhf(x,y))



def analyze_spectrum(v, spectrum):
    ''' Analyze a spectrum in velocity space, to get wind and temperature'''
    
    if any(np.isnan(spectrum)):
        return np.nan, np.nan, np.nan
    
    # Achieve this by fitting a Gaussian to the spectrum
    def my_func(p, v, spectrum):
        model = p[3] + particle_velocity_density(v,p[0],p[1],p[2])
        return model-spectrum

    from scipy.optimize import leastsq
    p0 = [0,1000,sum(spectrum)*(v[1]-v[0]),0.]
    p,flag = leastsq(my_func, p0, args=(v,spectrum))
    return p[0],p[1],p[2]



def white_light_scatter(dhf, tau0, P, h, om=1., alpha=0., M=20, N=20, R=20, N_int=30, R_int=30,
                        tol=1e-4, K=20, verbose=True, return_J=False):
    '''
    This function calculates the scattered light distribution for a given source.
    
    dhf: function that takes x,y and returns the column brightness
    tau0: optical thickness of atmosphere
    P: scattering phase function P(u0,u1,phi0,phi1)
    h: height of airglow layer [m]
    om: single-scattering albedo (default 1)
    alpha: ground albedo (0-1, default 0)
    M: number of tau grid points
    N: number of u grid points
    R: number of phi grid points
    N_int: number of u grid points to use in sky integration
    R_int: number of phi grid poitns to use in sky integration
    tol: stop iterating when the relative change in the solution is less than this
    K: maximum iterations
    verbose: whether to print progress
    return_J: optional flag (default False) which will return source function J
              instead of intensity I, for debugging purposes.
    
    Returns:
    I: (MxNxR) scattered light distribution (or, if return_J==True, returns source
               function, which is size (M-1)xNxR
    tau_bounds: (M) the bounds of the cells in tau
    u_bounds:   (N+1) the bounds of the cells in u
    phi_bounds: (R+1) the bounds of the cells in phi
    '''    
    
    # Define coordinates
    tau = linspace(0,tau0,M) # tau defined at grid walls
    tau_mid = (tau[1:] + tau[:-1])/2
    dtau = tau[1] - tau[0]

    th_bounds = linspace(pi,0,N+1)
    u_bounds = cos(th_bounds) # u defined at grid center
    u = (u_bounds[1:] + u_bounds[:-1])/2
    du = u_bounds[1:] - u_bounds[:-1]

    phi_bounds = linspace(0,2*pi,R+1) # phi defined at grid center
    phi = (phi_bounds[1:] + phi_bounds[:-1])/2
    dphi = phi[1] - phi[0]
    
    if h < 1000:
        warnings.warn('Are you sure you specified "h" in meters, not kilometers? (h=%.1f)'%h)

    ##################################################################
    # Compute matrices which will implement scattering calculations

    #### A: From J to I
    Arow = []
    Acol = []
    Aval = []
    shap = (M,N,R)
    for m in range(M):
        for r in range(R):
            # u > 0
            for n in range(N/2,N):
                mp = arange(0,m)
                row_index = ravel_multi_index((m,n,r),shap)
                col_indexes = ravel_multi_index((mp,n,r),shap)
                vals = 1/u[n] * exp(-dtau/u[n]*(m-mp-0.5)) * dtau
                Arow.extend(m*[row_index])
                Acol.extend(col_indexes)
                Aval.extend(vals)
            # u < 0
            for n in range(0,N/2):
                mp = arange(m,M-1)
                row_index = ravel_multi_index((m,n,r),shap)
                col_indexes = ravel_multi_index((mp,n,r),shap)
                vals = -1/u[n] * exp(-dtau/u[n]*(m-mp-0.5)) * dtau
                Arow.extend((M-1-m)*[row_index])
                Acol.extend(col_indexes)
                Aval.extend(vals)

        if verbose:
            clear_output(wait=True)
            time_mod.sleep(0.01)
            print('Building A: %i/%i'%(m+1,M))
            sys.stdout.flush()

    A = sp.coo_matrix((Aval, [Arow, Acol]), shape=(M*N*R,(M-1)*N*R))
    A = A.tocsr()

    #### B: From I to J
    Brow = []
    Bcol = []
    Bval = []
    shap = (M,N,R)
    for m in range(M-1):
        for n in range(N):
            for r in range(R):
                row_index = ravel_multi_index((m,n,r),shap)
                np,rp = meshgrid(range(N),range(R))
                col_indexes0 = ravel_multi_index((m,np,rp),shap)
                col_indexes1 = ravel_multi_index((m+1,np,rp),shap)
                val = om/(4*pi) * P(u[n],u[np],phi[r],phi[rp]) / 2.0 * du*dphi
                Brow.extend(2*N*R*[row_index])
                Bcol.extend(col_indexes0.flatten())
                Bcol.extend(col_indexes1.flatten())
                Bval.extend(val.flatten())
                Bval.extend(val.flatten())

        if verbose:
            clear_output(wait=True)
            time_mod.sleep(0.01)
            print('Building B: %i/%i'%(m+1,M-1))
            sys.stdout.flush()

    B = sp.coo_matrix((Bval, [Brow, Bcol]), shape=((M-1)*N*R,M*N*R))
    B = B.tocsr()

    ####################################################################
    # Compute single-scatter source function. This is the only step which
    # needs the airglow distribution

    # Initial state (iteration 0 is the initialization)
    J = zeros((M-1,N,R,K+1)) # defined on tau grid midpoints
    I = zeros((M,  N,R,K))   # defined on tau grid walls

    # Requires an integration over the airglow VER.
    # See notes page 28 for equation describing integration over sky.
    th_int_bounds = linspace(pi/2,0.0,N_int+1)
    u_int_bounds = cos(th_int_bounds)
    #u_int_bounds = linspace(0,1,N_int+1)
    u_int = (u_int_bounds[1:] + u_int_bounds[:-1])/2
    du_int = u_int_bounds[1:]-u_int_bounds[:-1]
    phi_int_bounds = linspace(0,2*pi,R_int+1)
    phi_int = (phi_int_bounds[1:] + phi_int_bounds[:-1])/2
    dphi_int = phi_int_bounds[1:]-phi_int_bounds[:-1]    

    for n in range(N):
        for r in range(R):
            # Perform calculation in a vectorized way to avoid expensive for loops.
            # Capitalized variables indicate 3-D variables (tau, u, phi).
            # In including the ground albedo effect in J0, it is assumed that
            # most scattering happens near the ground, in order to simplify
            # the trigonometry of ground reflections.
            TAU,U,PHI = meshgrid(tau_mid, u_int, phi_int, indexing='ij')
            _,DU,DPHI = meshgrid(tau_mid, du_int, dphi_int, indexing='ij')
            TH = arccos(U)
            GAMMA = arcsin(RE/(RE+h) * sin(TH))
            RHO = (TH - GAMMA) * (RE+h)
            X = -RHO * cos(PHI)
            Y = -RHO * sin(PHI)
            dJ = om/(4*pi) * (P(u[n], U, phi[r], PHI) * exp(-TAU/U) + alpha * P(u[n], -U, phi[r], PHI) * exp(-(2*tau0 - TAU)/U)) * \
                             dhf(X,Y) * sec(GAMMA) * DU*DPHI
            J[:,n,r,0] = sum(dJ, axis=(1,2))
            
        if verbose:
            clear_output(wait=True)
            time_mod.sleep(0.01)
            print('Initial source function: %i/%i'%(n+1,N))
            sys.stdout.flush()

    ######################################################################
    # Compute multiple scattering solution based on Hansen's method
    # of successive scattering.

    for k in range(K):

        # Compute I
        Iflat = A.dot(J[:,:,:,k].flatten())
        I[:,:,:,k] = reshape(Iflat, (M,N,R))
        
        # If the ground is reflective, add the extra contribution
        if alpha > 0.0:
            nlt0 = arange(N/2) # n values less than zero
            I[:,nlt0,:] += alpha*I[M-1,N-1-nlt0,:]

        # Compute J for next iteration
        Jflat = B.dot(I[:,:,:,k].flatten())
        J[:,:,:,k+1] = reshape(Jflat,(M-1,N,R))

        if verbose:
            clear_output(wait=True)
            time_mod.sleep(0.01)
            print('Scattering iteration: %i/%i'%(k+1,K))
            sys.stdout.flush()

        # See how much solution is changing
        change = norm(J[:,:,:,k+1])/norm(sum(J[:,:,:,:k+1],axis=3))
        if change < tol:
            if verbose:
                print('Halted')
            break
            
    Ifull = sum(I,axis=3) # All scattering orders
    Jfull = sum(I,axis=3) 
    
    return Ifull, tau, u_bounds, phi_bounds




def scatter_along_los(dhf, tau0, P, h, th_look, phi_look, om=1., q=None, full_fov_deg=None, M=20, N=20, R=20, N_int=30, R_int=30, Q_int=30,
                        tol=1e-4, K=20, verbose=True):
    '''
    For given lines of sight, calculate the direct and scattered brightness measurements, for both
    the in-FoV and out-of-FoV (i.e., stray light) components.
    
    dhf: function that takes x,y and returns the column brightness
    tau0: optical thickness of atmosphere
    P: scattering phase function P(u0,u1,phi0,phi1)
    h: height of airglow layer [m]
    th_look:  (array) zenith angle of observation (deg, geophysical convention)
    phi_look: (array) azimuth angle of observation (deg, geophysical convention)
    om: single-scattering albedo
    q: stray light function (None = ignore stray light)
    full_fov_deg: the field of view in degrees (only used if q is specified)
    M: number of tau grid points
    N: number of u grid points
    R: number of phi grid points
    N_int: number of u grid points to use in sky integration
    R_int: number of phi grid points to use in sky integration
    Q_int: number of u and phi grid points to use in stray light calculation
    tol: stop iterating when the relative change in the solution is less than this
    K: maximum iterations
    verbose: whether to print progress
    
    Returns:
    g_sc:  (array) measured brightness from atmospheric scattering
    g_dir: (array) measured brightness from direct ray
    s_sc:  (array) measured brightness from stray light of scattered signal
    s_dir: (array) measured brightness from stray light of direct signal
    '''


    Ifull, tau, u_bounds, phi_bounds = white_light_scatter(dhf, tau0, P, h, om=om, M=M, N=N, R=R, N_int=N_int, R_int=R_int,
                        tol=tol, K=K, verbose=verbose)

    
    # Precalculate interpolation function
    ui = (u_bounds[1:] + u_bounds[:-1])/2
    phii = (phi_bounds[1:] + phi_bounds[:-1])/2
    thi = arccos(ui)[::-1]
    
 
    
    ######## VERSION USING ROTATED COORDINATES #############
    Ii = Ifull[-1,::-1,:]
    
    # There is a discontinuity at u=0 (th=pi/2), so interpolate over either
    # the top half or bottom half plane. In practice, we only need the bottom.
    # Another detail is that we should specify a point at zenith, otherwise
    # it's extrapolating.
    zenith_val = np.mean(Ii[0,:])
    # Define interpolation for phi in 90 to 270
    Iint0 = RectSphereBivariateSpline(thi[:N/2], phii, Ii[:N/2,:], pole_values=(zenith_val, None))
    # Define interpolation for phi in 270 to 90, using a rotated coordinate system to avoid problems at phi=0
    r0 = sum(phii < pi/2) + 1 # index just after 90 deg
    r1 = sum(phii < 3*pi/2) - 1 # index just before 270 deg
    phii_rot = np.concatenate((phii[r1:]-2*pi, phii[:r0+1])) + pi
    Ii_rot = np.concatenate((Ii[:,r1:], Ii[:,:r0+1]),axis=1)
    Iint1 = RectSphereBivariateSpline(thi[:N/2], phii_rot, Ii_rot[:N/2,:], pole_values=(zenith_val, None))

    def I_sc(u,phi):
        '''
        Given a calculated, discretized, scattered intensity (Ifull), evaluate it 
        at the given u, phi. Use linear interpolation on a sphere.
        Global variables: Iint0, Iint1
        '''
        if phi > np.pi/2 and phi < 3*np.pi/2:
            return Iint0(arccos(u),phi).item()
        else:
            phi_rot = np.mod(phi+pi, 2*pi)
            return Iint1(arccos(u),phi_rot).item()
    
    
    def I_dir(dhf,tau0,u,phi):
        '''
        Evaluate the direct intensity.
        Global variables: RE,h
        '''
        if u <= 0.0:
            return 0.0 # no incident light looking downwards
        th = arccos(u)
        gamma = arcsin(RE/(RE+h) * sqrt(1-u**2))
        rho = (th-gamma) * (RE+h)
        x = -rho*cos(phi)
        y = -rho*sin(phi)
        return dhf(x,y) * sec(gamma) * exp(-tau0/u)


    # Determine correct u,phi (in math coordinates, incoming ray convention) to use
    # for direct signal
    g_sc = []
    g_dir = []
    for i in range(len(th_look)):
        u_dir = cos(th_look[i]*pi/180)
        phi_dir = mod(pi/2 - phi_look[i]*pi/180 + pi, 2*pi)
        g_sc.append(I_sc(u_dir, phi_dir))
        g_dir.append(I_dir(dhf, tau0, u_dir, phi_dir))
        
        
    g_sc  = array(g_sc)
    g_dir = array(g_dir)
        
    # Stray light calculation, if desired:
    s_sc = zeros(len(g_sc))
    s_dir = zeros(len(g_dir))
    if q:    # Grid points for stray light integral:
        if full_fov_deg is None:
            raise Exception('If "q" is specified, "full_fov_deg" needs to be specified')
        dS = 2*pi*(1-cos(full_fov_deg/2 * pi/180)) # solid angle
        
        th_bound_vec = linspace(pi/2, full_fov_deg/2*pi/180, Q_int+1) # don't include FoV
        u_bound_vec = cos(th_bound_vec)
        u_vec = (u_bound_vec[1:] + u_bound_vec[:-1])/2
        du_vec = u_bound_vec[1:] - u_bound_vec[:-1]

        phi_bound_vec = linspace(0,2*pi,Q_int+1)
        phi_vec = (phi_bound_vec[1:] + phi_bound_vec[:-1])/2
        dphi_vec = phi_bound_vec[1:] - phi_bound_vec[:-1]
        
        
        for looki in range(len(th_look)):
            # Collect integral term at each grid point
            d_gsc_stray = zeros((Q_int,Q_int))  # terms of integral at each grid point (see notes 2016-03-28)
            d_gdir_stray = zeros((Q_int,Q_int)) # terms of integral at each grid point (see notes 2016-03-28)
            for i in range(Q_int):
                for j in range(Q_int):
                    th = arccos(u_vec[i])
                    ph = phi_vec[j]
                    # th (rad) is the angle of the incoming ray relative to -z
                    # ph (rad) is the (math-convention) azimuth angle of the incoming ray
                    # th_look (deg) is the angle of the look direction relative to +z
                    # phi_look (deg) is the (geophysical-convention) azimuth angle of the look direction
                    th_look_rad = th_look[looki]*pi/180
                    phi_look_rad = phi_look[looki]*pi/180
                    # unit vector (x,y,z)=(E,N,U) of look direction
                    look = [sin(th_look_rad)*sin(phi_look_rad), sin(th_look_rad)*cos(phi_look_rad), cos(th_look_rad)]
                    # unit vector (x,y,z)=(E,N,U) of direction to source
                    src  = [-sin(th)*cos(ph), -sin(th)*sin(ph), cos(th)]
                    angle_between = arccos(np.dot(look,src))
                    d_gsc_stray[i,j] = u_vec[i] * I_sc(u_vec[i],phi_vec[j]) * q(angle_between) * \
                                       du_vec[i] * dphi_vec[j] 
                    d_gdir_stray[i,j] = u_vec[i] * I_dir(dhf,tau0, u_vec[i],phi_vec[j]) * q(angle_between) * \
                                       du_vec[i] * dphi_vec[j]     

            # Add these stray components to the previously-calculated in-FoV component, of both
            # the direct and atmosphere-scattered signal.
            # Determine correct u,phi (in math coordinates, incoming ray convention) to use
            # for direct signal.
            # Technically I should scale the in-FoV component by dS, but that's the same as
            # scaling the stray component by 1/dS. (For backwards compatibility).
            u_dir = cos(th_look[looki]*pi/180)
            phi_dir = mod(pi/2 - phi_look[looki]*pi/180 + pi, 360)
            s_sc[looki]  = 1/dS * sum(d_gsc_stray)
            s_dir[looki] = 1/dS * sum(d_gdir_stray)
    
    return g_sc, g_dir, s_sc, s_dir





def descatter_asi_data(imcut, yim, tau0, P, h,  om=1., M=10, N=20, R=10, N_int=20, R_int=20, 
                       tol=1e-4, K=20, N_iters=5, verbose=True):
    '''
    Remove the effects of atmospheric extinction and scattering from the ASI data. This assumes data
    are in the format given by Carlos (which already accounts for van Rhijn, aperture projection, 
    flat field, solid angle, vignetting, stray light, etc.)
    '''
        
    imcut_meas = imcut.copy()
    imcut_scatt = zeros(len(imcut))
    iter_change = zeros(N_iters)

    thim = abs(arctan2(yim,h)) # zenith angles for each imcut location (yim)

    for n_iter in range(N_iters):
        
        if verbose:
            clear_output(wait=True)
            time_mod.sleep(0.01)
            print('%i/%i' % (n_iter+1, N_iters))
            sys.stdout.flush()

        # Define 2D brightness distribution function from current guess (sky - scatt)
        imcut_iter = exp(tau0/cos(thim)) * (imcut_meas - cos(thim)*imcut_scatt)
        dhf_interp_iter = interpolate.interp1d(yim, imcut_iter, bounds_error = False)

        def dhf_iter(x,y):
            '''
            Airglow source VER function, integrated vertically. (i.e., thin-shell approximation).
            It is assumed that input variables are arrays of the same size.
            '''
            r_max = (RE+h) * arccos(RE/(RE+h))
            f = dhf_interp_iter(y)
            f[y < yim[0]] = imcut_iter[0] # Zero-order-hold extrapolation
            f[y > yim[-1]] = imcut_iter[-1] # Zero-order-hold extrapolation
            f[x**2 + y**2 > r_max**2] = 0.0 # Nothing beyond horizon
            return f

        Ifull, tau, u_bounds, phi_bounds = white_light_scatter(dhf_iter, tau0, P, h, om=om, 
                                    M=M, N=N, R=R, N_int=N_int, R_int=R_int, tol=tol, K=K, verbose=False)
        
        # Using interpolation
        def I_sc(u,phi):
            '''
            Given a calculated, discretized, scattered intensity (Ifull), evaluate it 
            at the given u, phi. Use linear interpolation on a sphere.
            Global variables: Ifull, u_bounds, phi_bounds
            '''
            ui = (u_bounds[1:] + u_bounds[:-1])/2
            phii = (phi_bounds[1:] + phi_bounds[:-1])/2
            thi = arccos(ui)[::-1]
            Ii = Ifull[-1,::-1,:]
            # There is a discontinuity at u=0 (th=pi/2), so interpolation over either
            # the top half or bottom half plane. In practice, we only need the bottom.
            # Another detail is that we should specify a point at zenith, otherwise
            # it's extrapolating.
            zenith_val = np.mean(Ii[0,:])
            # Because this interpolator is not ideal, create a ghost data point to
            # enforce linear interpolation across 0-2pi jump in phi
            phii_new = concatenate(([0],phii))
            Ii_new = zeros((N,R+1))
            Ii_new[:,0] = (Ii[:,-1] + Ii[:,0])/2
            Ii_new[:,1:] = Ii
            Iint = RectSphereBivariateSpline(thi[:N/2], phii_new, Ii_new[:N/2,:], pole_values=(zenith_val, None))

            return Iint(arccos(u),phi).item()
        

        # Loop over cut and remove scattered component
        imcut_scatt_new = zeros(len(imcut_meas))
        for i in range(len(yim)):
            u = h/sqrt(h**2 + yim[i]**2)
            phi = mod(arctan2(yim[i],0)+pi,2*pi) # change from ray direction to look direction
            imcut_scatt_new[i] = I_sc(u,phi)

        iter_change[n_iter] = norm(imcut_scatt - imcut_scatt_new)/norm(imcut_scatt)
        imcut_scatt = imcut_scatt_new
        
    # Create function for return
    
    # Define 2D brightness distribution function from current guess (sky - scatt)
    imcut_descatter = exp(tau0/cos(thim)) * (imcut_meas - cos(thim)*imcut_scatt)
    dhf_interp_descatter = interpolate.interp1d(yim, imcut_descatter, bounds_error = False)

    # Remove this calculated scattered component and define a new airglow source function
    def dhf_descatter(x,y):
        '''
        Airglow source VER function, integrated vertically. (i.e., thin-shell approximation).
        It is assumed that input variables are arrays of the same size.
        '''
        r_max = (RE+h) * arccos(RE/(RE+h))
        f = dhf_interp_descatter(y)
        f[y < yim[0]] = imcut_descatter[0] # Zero-order-hold extrapolation
        f[y > yim[-1]] = imcut_descatter[-1] # Zero-order-hold extrapolation
        f[x**2 + y**2 > r_max**2] = 0.0 # Nothing beyond horizon
        return f
    
    return dhf_descatter




def descatter_asi_data_2D(imd, asi_mask_d, Xd, Yd, tau0, P, h, om=1., M=10, N=20, R=10, N_int=20, R_int=20, 
                       tol=1e-4, K=20, N_iters=5, verbose=True):
    '''
    Remove the effects of atmospheric extinction and scattering from the ASI data. This assumes data (imd)
    already has all instrumental effects removed: stray light, solid angle, flat field, vignetting, aperture projection,
    and also van Rhijn.
    '''
        
    imd_meas = imd.copy()
    imd_scatt = zeros(np.shape(imd))
    iter_change = zeros(N_iters)

    thid = abs(arctan2(sqrt(Xd**2 + Yd**2),h)) # zenith angles for each imd location

    for n_iter in range(N_iters):
        
        if verbose:
            clear_output(wait=True)
            time_mod.sleep(0.01)
            print('%i/%i' % (n_iter+1, N_iters))
            sys.stdout.flush()

        # Define 2D brightness distribution function from current guess (sky - scatt)
        imd_iter = exp(tau0/cos(thid)) * (imd_meas - cos(thid)*imd_scatt)

        def dhf_iter(x,y):
            '''
            Airglow source VER function, integrated vertically. (i.e., thin-shell approximation).
            It is assumed that input variables are arrays of the same size.
            '''
            sc = 1.0 # how much more to favor y over x when interpolating
            f = interpolate.griddata((Xd[asi_mask_d].flatten(),sc*Yd[asi_mask_d].flatten()),imd_iter[asi_mask_d].flatten(),(x,sc*y),
                                     method='nearest')
            # Special case if scalar inputs.
            if not hasattr(x,'__len__'): # don't worry about horizon or whatever. It's probably fine.
                return f 

            r_max = (RE+h) * arccos(RE/(RE+h))
            f[x**2 + y**2 > r_max**2] = 0.0 # nothing beyond horizon (likely doesn't matter)
            return f
        
        
        Ifull, tau, u_bounds, phi_bounds = white_light_scatter(dhf_iter, tau0, P, h, om=om,
                                    M=M, N=N, R=R, N_int=N_int, R_int=R_int, tol=tol, K=K, verbose=False)
        
        # Using interpolation
        def I_sc(u,phi):
            '''
            Given a calculated, discretized, scattered intensity (Ifull), evaluate it 
            at the given u, phi. Use linear interpolation on a sphere.
            Global variables: Ifull, u_bounds, phi_bounds
            '''
            ui = (u_bounds[1:] + u_bounds[:-1])/2
            phii = (phi_bounds[1:] + phi_bounds[:-1])/2
            thi = arccos(ui)[::-1]
            Ii = Ifull[-1,::-1,:]
            # There is a discontinuity at u=0 (th=pi/2), so interpolation over either
            # the top half or bottom half plane. In practice, we only need the bottom.
            # Another detail is that we should specify a point at zenith, otherwise
            # it's extrapolating.
            zenith_val = np.mean(Ii[0,:])
            # Because this interpolator is not ideal, create a ghost data point to
            # enforce linear interpolation across 0-2pi jump in phi
            phii_new = concatenate(([0],phii))
            Ii_new = zeros((N,R+1))
            Ii_new[:,0] = (Ii[:,-1] + Ii[:,0])/2
            Ii_new[:,1:] = Ii
            Iint = RectSphereBivariateSpline(thi[:N/2], phii_new, Ii_new[:N/2,:], pole_values=(zenith_val, None))

            return Iint(arccos(u),phi).item()
        

        # Loop over cut and remove scattered component
        imd_scatt_new = zeros(np.shape(imd_meas))
        for i in range(np.shape(imd)[0]):
            for j in range(np.shape(imd)[1]):
                u = h/sqrt(h**2 + Xd[i,j]**2 + Yd[i,j]**2)
                phi = mod(arctan2(Yd[i,j],Xd[i,j])+pi,2*pi) # change from ray direction to look direction
                imd_scatt_new[i,j] = u * I_sc(u,phi)

        iter_change[n_iter] = norm(imd_scatt - imd_scatt_new)/norm(imd_scatt)
        imd_scatt = imd_scatt_new
        
    # Create function for return
    
    # Define 2D brightness distribution function from current guess (sky - scatt)
    imd_descatter = exp(tau0/cos(thid)) * (imd_meas - cos(thid)*imd_scatt)

    # Remove this calculated scattered component and define a new airglow source function
    def dhf_descatter(x,y):
        '''
        Array inputs are expected, but scalar inputs will work too.
        '''
        sc = 1.0 # how much more to favor y over x when interpolating
        f = interpolate.griddata((Xd[asi_mask_d].flatten(),sc*Yd[asi_mask_d].flatten()),imd_descatter[asi_mask_d].flatten(),(x,sc*y),
                                 method='nearest')
        # Special case if scalar inputs.
        if not hasattr(x,'__len__'): # don't worry about horizon or whatever. It's probably fine.
            return f 

        r_max = (RE+h) * arccos(RE/(RE+h))
        f[x**2 + y**2 > r_max**2] = 0.0 # nothing beyond horizon (likely doesn't matter)
        return f
    
    return dhf_descatter




def scattered_and_direct_spectra(dhf, tau0, P, h, th_look, phi_look, uvw, T, om=1., q=None, full_fov_deg=None, M=20, N=20, 
                                 R=20, N_int=30, R_int=30, Q_int=30, tol=1e-4, K=20, L=20, verbose=True):

    Nlook = len(th_look)
    v = linspace(-3000,3000,L)
    spectrum_direct = zeros((Nlook,L))
    spectrum_scatt = zeros((Nlook,L))
    spectrum_direct_stray = zeros((Nlook,L))
    spectrum_scatt_stray = zeros((Nlook,L))

    for l in range(L):

        # Define function that will give brightness in this wavelength band
        dhf_l = lambda x,y: spectral_intensity(dhf,x,y,v[l],uvw,T,h)

        # Calculate direct and scattered spectra
        # Use this function to evaluate white-light scatter at this wavelength
        g_sc, g_dir, s_sc, s_dir = scatter_along_los(dhf_l, tau0, P, h, th_look, phi_look, om=om, M=M, N=N, R=R, q=q,
                                          full_fov_deg=full_fov_deg,
                                          N_int=N_int, R_int=R_int, Q_int=Q_int, K=K, tol=tol, verbose=False)
        spectrum_scatt[:,l] = g_sc
        spectrum_direct[:,l] = g_dir
        spectrum_scatt_stray[:,l] = s_sc
        spectrum_direct_stray[:,l] = s_dir

        if verbose:
            clear_output(wait=True)
            time_mod.sleep(0.01)
            print('%i/%i'%(l+1,L))
            sys.stdout.flush()

    return v, spectrum_scatt, spectrum_direct, spectrum_scatt_stray, spectrum_direct_stray





