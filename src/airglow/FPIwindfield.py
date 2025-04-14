# Import a bunch of stuff into the namespace so it 
# acts like an IPython notebook. Probably not good coding 
# practice, I know...
from numpy import *
from pylab import *
#from numpy.linalg import *
import matplotlib

# Import the usual modules
import FPIprocess
import fpiinfo
from datetime import datetime, timedelta
import ASI # for cnv_azel2latlon
import pytz
import scipy.sparse as sp
from mpl_toolkits.basemap import Basemap
import os
import glob
import FPI
from mpl_toolkits.basemap import Basemap
import subprocess
from collections import deque
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


class WindField:
    '''
    A class representing a regional wind field estimate. This handles loading raw data, inverting,
    and visualizing the results.
    Visualizations:
        - Wind field with quiver: create pngs and movie with quiver and optional color for vertical wind
        - Wind field with tracers: create movie. We should probably delete the pngs.
        - Quick look summary plot
        - Wind flow movie with div and vort plots (TODO)
        - Div and Vort plots
    TODO:
        - Make eval_field do a linear interpolation instead of nearest neighbor.
          (using griddata, this might even be faster).
    '''
    
    def __init__(self, network_names, year, doy,
                 cloudthresh = -17., timestep_fit = 180., startstop_hr = None, \
                 Nx = 25, Ny = 25, ht = 250., errbarthresh = 30., oldthresh = 60., minlocs=5., \
                 chi2mult = 1.0, tol0=1e-4, tol1 = 1e-1, tol2 = 0.1, L=5, L2 = 10,\
                 printstuffinner = False, printstuffouter = False, maxiters=100, \
                 lam1small = 1e-4, lam0guess = 1e-5, lam0interior = 1e10, \
                 lam1interior = 1, distthresh = 250, Nmc = 100, estimate_vertical_wind=False, \
                 estimate_uncertainty=False, use_zenith_ref=False,):
        '''
        Initializes the Wind Field object. Loads and organizes data, but does not run inversion.
        INPUTS:
            network_names:   list of str, which networks to use in inversion. e.g., ['renoir', 'peru'] or ['nation']
            year:            int,         year of night of data to estimate wind field.
            doy:             int,         day of year of night of data to estimate wind field. This specifies a date in 
                                          local time. A night of data spans two local-time dates. 
                                          This specifies the first one.
        OPTIONAL INPUTS:
            cloudthresh:     float, deg C. above this number, it is considered cloudy and the data are ignored
                                    (default -17)
            timestep_fit:    float, seconds. time step for fitting wind fields (as small as your patience allows)
                                    (default 180)
            startstop_hr:    list of float, hours in UT. If None, run the whole night. If [h1, h2], only run samples
                             where the UT time in hours (mod 24) is > h1 and < h2 (default None)
            Nx:              int,   number of grid points in longitudinal direction (default 25)
            Ny:              int,   number of grid points in latitudinal direction (default 25)
            ht:              float, km. assumed emission altitude. (default 250)
            errbarthresh:    float, m/s. error bar above which measurement is ignored. (default 30)
            oldthresh:       float, minutes. if there is a gap longer than this, don't interpolate through it.
                                    (default 60)
            minlocs:         int.   minimum number of available locations at a given time in order to estimate a 
                                    wind field (default 5)
            chi2mult:        float. discrepancy principle factor (set errors equal to chi2mult times chi2)
                                    (default 1.0)
            tol0:            float. stop inner problem when data cost is within 100*tol1 percent of desired_cost
                                    (default 1e-4)
            tol1:            float. stop outer problem when grad cost is within 100*tol0 percent of midpoint
                                    (this is not that sensitive) (default 1e-1)
            tol2:            float. the change during outer line search to consider close enough to 0, in log10 
                                    (default 0.1)
            L:               float. multiplicative factor for inner line search (default 5)
            L2:              float. multiplicative factor for outer line search (default 10)
            printstuffinner: bool.  print output for inner problem (default False)
            printstuffouter: bool.  print output for outer problem (default False)
            maxiters:        int.   maximum number of iterations before giving up (default 100)
            lam1small:       float. starting value for search for optimal lam1 (must be small enough) (default 1e-4)
            lam0guess:       float. starting value for optimal lam0 search (this isn't important) (default 1e-5)
            lam0interior:    float. An extremely large value which, if it provides a low-data-cost solution, 
                                    indicates that a flat wind field will fit the data just fine. This should 
                                    be as large as possible without causing numerical problems in the inversion.
                                    (default 1e10)
            lam1interior:    float. The curvature weighting to accompany lam0interior. This is arbitrary 
                                    as long as it is < inf. (default 1)
            distthresh:      float, km. Points farther than this from a measurement location will not be plotted.
                                    (default 250)
            Nmc:             int.   number of Monte Carlo trials for estimating uncertainty. Accuracy goes as sqrt(1/Nmc).
                                    (default 100)
            estimate_vertical_wind: bool. If False, set the vertical wind to a constant, determined by 
                                          the zenith measurements.
            estimate_uncertainty: bool. If True, estimate the uncertainty in the resulting wind field.
                                        This is time-consuming and not recommended for Nx,Ny>20.
                                        (default False)
            use_zenith_ref:  bool.  If True, override the Doppler reference used in the npz file, and re-establish
                                    the Doppler reference using zenith reference, to force the vertical wind to be zero.
                                    If this is True, then estimate_vertical_wind will be set to False. (default False).
        '''
        
        # Record input
        if use_zenith_ref:
            estimate_vertical_wind = False
        
        self.network_names = network_names
        self.year = year
        self.doy = doy
        self.cloudthresh = cloudthresh
        self.timestep_fit = timestep_fit
        self.startstop_hr = startstop_hr
        self.Nx = Nx
        self.Ny = Ny
        self.ht = ht
        self.errbarthresh = errbarthresh
        self.oldthresh = oldthresh
        self.minlocs = minlocs
        self.chi2mult = chi2mult
        self.tol0 = tol0
        self.tol1 = tol1
        self.tol2 = tol2
        self.L = L
        self.L2 = L2
        self.printstuffinner = printstuffinner
        self.printstuffouter = printstuffouter
        self.maxiters = maxiters
        self.lam1small = lam1small
        self.lam0guess = lam0guess
        self.lam0interior = lam0interior
        self.lam1interior = lam1interior
        self.distthresh = distthresh
        self.Nmc = Nmc  
        self.estimate_vertical_wind = estimate_vertical_wind
        self.estimate_uncertainty = estimate_uncertainty
        self.use_zenith_ref = use_zenith_ref
        
                
        # Set other parameters
        self.vertweight = 1.0 # Extra weighting for vertical wind structure (1.0 = no extra weighting)
        self.message = ''
        base_dir = '/rdata/airglow/fpi/results/level3_windfield/'
        self.dir = '%s%s/' % (base_dir, self.get_stub())
        self.toffset = 20 # minutes after the first observation to start inverting
        self.success = True # True by default. Will be set to False if necessary
        
        # Write the parameters to the log, so users can see which parameters
        # were used in this inversion
        self.message += 'Initializing wind field estimation with the following parameters:\n'
        fields = self.__dict__
        keys = fields.keys()
        keys.sort()
        for param in keys:
            if param not in ['message']: # list the parameters here not to show
                self.message+='\t%28s: %s\n' % (param, fields[param])
        self.message += '\n'
        
        # For completeness, list the class variables that will be assigned later
        self.instr_names = None
        self.dat = None # cleaned level 0 data
        self.all_t = None
        self.Dc = None
        self.Dg = None
        self.windfield = None # inverted wind fields
        self.lats = None # all measurement locations used throughout the night
        self.lons = None
        self.t_divvort = None
        self.div = None
        self.vort = None

        
        # Make directory to save results in
        try:
            os.mkdir(self.dir)
        except Exception as e:
            self.message += 'Previous results will be overwritten in\n\t%s\n' % self.dir

            
        # Load the data
        year = self.year
        doy = self.doy
        dn = datetime(year,1,1) + timedelta(days=doy-1)

        instr_names = []
        for network_name in network_names:
            sites = fpiinfo.get_network_info(network_name).keys()
            for site in sites:
                instr_names.extend(fpiinfo.get_instr_at(site, dn))
        instr_names.sort()

        dat = {}
        for instr_name in instr_names:

            self.message += 'Loading %s\n' % instr_name
            try:
                r = FPIprocess.load_level0(instr_name, year, doy)
            except:
                self.message += '\t%s: No Data. Ignoring this instrument.\n' % instr_name
                continue # skip this instrument

            instrdat = {}
            fpir = r['FPI_Results']
            directions = fpir['direction']
            uniqdirecs = list(set(directions))
            losu = fpir['LOSwind']
            losue = fpir['sigma_LOSwind'] # sigma_LOSwind or sigma_fit_LOSwind?
            times = fpir['sky_times']
            cloud = fpir['Clouds']
            site_name = fpiinfo.get_site_of(instr_name, dn)
            site = fpiinfo.get_site_info(site_name)
                        

            if not cloud:
                self.message += '\t%s: No cloud data. Hoping it was clear.\n' % instr_name

            try:
                # Calculate doppler reference
                ref = fpir['reference']
                if use_zenith_ref:
                    ref = 'zenith'
                    
                if ref=='laser': # use the doppler reference function above
                    drefvec, drefevec = FPI.DopplerReference(fpir, reference='laser')
                    losu = losu - drefvec
                    losue = sqrt(losue**2 + drefevec**2)
                elif ref=='zenith': # use zenith reference and remove the zenith measurement
                    drefvec, drefevec = FPI.DopplerReference(fpir, reference='zenith')
                    losu = losu - drefvec
                    losue = sqrt(losue**2 + drefevec**2)
                    idx = array([d == 'Zenith' for d in directions]) # select zenith measurements
                    losu[idx] = nan # nan them out
                    self.message += '\t%s: Using zenith reference. '% instr_name + \
                                    'Errors are larger and zenith measurement is ignored.\n' 
            except Exception as e:
                self.message += '\t%s: Doppler reference failed: %s. Ignoring this instrument.\n' % (instr_name, e)
                continue # skip this instrument


            # Collect data from each direction
            for direc in uniqdirecs:
                if direc not in ['Unknown','Laser']:
                    du  = array([u for (u,d) in zip(losu,directions)  if d==direc])
                    due = array([u for (u,d) in zip(losue,directions) if d==direc])
                    dt  = array([t.astimezone(pytz.utc) for (t,d) in zip(times,directions) if d==direc])
                    thisdat = {}
                    # manually correct anything that needs to be manually corrected (see above)
                    du = self.quality_hack(year, doy, instr_name, dt, du, direc) 

                    az = site['Directions'][direc]['az']
                    ze = site['Directions'][direc]['ze']
                    if ze < 0:
                        az = az + 180
                        ze = -ze
                    latx, lonx = ASI.ConvertAzEl2LatLon(az, 90-ze, ht, site['Location'][0], site['Location'][1])
                    thisdat['losu'] = du
                    thisdat['losue'] = due
                    thisdat['t'] = dt
                    thisdat['az'] = az
                    thisdat['ze'] = ze
                    thisdat['lat'] = latx
                    thisdat['lon'] = lonx
                    if cloud:
                        c = array([t for (t,d) in zip(cloud['mean'],directions) if d==direc])
                        thisdat['c'] = c
                    else:
                        thisdat['c'] = nan*ones_like(dt)
                    instrdat[direc] = thisdat
                else:
                    self.message += '\t%s: Unknown direction ignored.\n' % instr_name 
            dat[instr_name] = instrdat

        # Required to pass:
        # 1) Error bar less than errbarthresh
        # 2) Either:
        #     a) no cloud data (nans)
        #     b) cloud < cloudthresh
        self.message += 'Performing quality control based on cloud data and error bars\n'
        instr_names = dat.keys()
        for iname in instr_names:
            instrdat = dat[iname]
            for direc in instrdat.keys():
                u = instrdat[direc]['losu']
                ue = instrdat[direc]['losue']
                c = instrdat[direc]['c']
                t = instrdat[direc]['t']
                missing = [isnan(x) for x in c]
                cloudidx = logical_or(missing, c < cloudthresh)
                erroridx = logical_and(ue < errbarthresh, ue > 0)
                nanidx = ~isnan(u)
                idx = logical_and(cloudidx, logical_and(erroridx, nanidx))
                u = u[idx]
                ue = ue[idx]
                c = c[idx]
                t = t[idx]
                instrdat[direc]['losu'] = u
                instrdat[direc]['losue'] = ue
                instrdat[direc]['c'] = c
                instrdat[direc]['t'] = t
                self.message += '\t%s %15s: %i/%i passed\n' % (iname, direc, sum(idx),len(idx))
                if sum(idx) == 0:
                    del instrdat[direc]
                    
        instr_names = dat.keys()

        # Gather all times
        all_t = []
        for iname in instr_names:
            for dname in dat[iname].keys():
                all_t.extend(dat[iname][dname]['t'])
        all_t.sort()
        
        # Save the important variables as class variables
        self.instr_names = instr_names
        self.dat = dat
        self.all_t = all_t

        # Write log file
        logfn = '%s%s_windfield_log.txt' % (self.dir, self.get_stub())
        f = open(logfn, 'w')
        f.write(self.message)
        f.close()

        
        
        
        
        
    def __repr__(self):
        return 'WindField Object: \n\tnetworks: %s\n\tyear: %i\n\tdoy: %i' % (self.network_names, self.year, self.doy)

        
        
        
        
        
        
    def quality_hack(self, year, doy, instr_name, dt, du, direc):
        '''
        This is a constantly-evolving function to do any ad-hoc corrections needed for certain
        data sets. Hopefully, this is not how the final version of processing will
        be implemented, but will merely act as a placeholder until more sophisticated 
        processing routines are developed, or until the database has been cleaned up.
        TODO:
            - better automatic correction for laser drift and laser outliers. Ideally this
              will be done in level 0.
        '''


        # For "wave" event on 2014 doy 127, ignore the first part of PAR data,
        # when the laser is acting up.
        # Planned ultimate resolution: Go through database and remove data
        # that coincide with large variantions in laser intensity or etalon
        # gap.
        if year == 2014 and doy == 127 and instr_name == 'minime06':
            hour = array([t.hour for t in dt])
            minute = array([t.minute for t in dt])
            delete = (hour + minute/60.) < 4.5
            du[delete] = nan
            self.message += '\tminime06: Manually deleted early-night data because of laser problems.\n'
        elif year == 2014 and doy == 142 and instr_name == 'minime06':
            hour = array([t.hour for t in dt])
            minute = array([t.minute for t in dt])
            delete = (hour + minute/60.) < 6.0
            du[delete] = nan
            self.message += '\tminime06: Manually deleted early-night data because of laser problems.\n'
        elif year == 2014 and doy == 98 and instr_name == 'minime06':
            du[:5] = nan
            self.message += '\tminime06: Manually deleted first few samples because of laser problems.\n'

        # Ignore VTI's east and west looks (should probably do this for all days).
        elif year == 2014 and doy in [58,79,84] and instr_name == 'minime09' and direc in ['East','West'] :
            du[:] = nan

        elif year == 2013 and doy == 151: # ignore impossibly large values
            idx = abs(du) > 300
            du[idx] = nan

        elif year == 2014 and doy == 84 and instr_name == 'minime05':
            du[0] = nan

        # Try to ignore large jumps in PAR wind.
        # (This is turned off for now, because it's getting rid of good data. Let it fly.
        '''
        if instr_name == 'minime06' and not (year==2013 and doy==274) and not (year==2014 and doy in [49,50] and year!=2012) \
            and not (year==2012 and doy == 318):
            ddu = diff(du)
            bad = abs(ddu) > 100
            idx = logical_or(concatenate(([False],bad)), concatenate((bad,[False])))
            du[idx] = nan
            N = sum(idx)
            if N > 0:
                self.message += '\tminime06: Ignored %i samples due to large jumps.\n' % N
        '''
        # Ignore data during weird interference on Feb 19 storm at PAR
        if instr_name == 'minime06' and year == 2014 and doy == 49:
            hour = array([t.hour for t in dt])
            minute = array([t.minute for t in dt])
            thour = hour + minute/60.
            delete = logical_and(thour > 6.9, thour < 7.3)
            du[delete] = nan
            self.message += '\tminime06: Manually deleted %i samples during interference.\n' % (sum(delete))
            
        # Bad gap variation in last few samples:
        if instr_name == 'minime08' and year == 2012 and doy == 210:
            du[-2:] = nan

        else:
            du = du

        return du
        
    
    
    
    def get_start_stop_time(self, timestep):
        '''
        Determine the starting and stopping time for the inversion and for
        plotting, based upon the sample times and the specified start stop
        times, if they exist.
        '''
        all_t = self.all_t
        t0 = all_t[0] + timedelta(minutes=self.toffset) # to make sure other sites have started
        t1 = all_t[-1] - timedelta(minutes=self.toffset)
        if self.startstop_hr is not None: # The user specified start stop times
            hourvec = array([t.hour + t.minute/60.0 for t in all_t])
            hr1 = self.startstop_hr[0]
            hr2 = self.startstop_hr[1]
            if hr2 < hr1: # the times cross the 24-hour mark
                ok = logical_or(hourvec >= hr1, hourvec <= hr2)
            else: # they don't
                ok = logical_and(hourvec >= hr1, hourvec <= hr2)
            t0 = min(array(all_t)[ok])
            t1 = max(array(all_t)[ok])
        dt = timestep/60.0
        t1 = t1 - timedelta(seconds = t1.second, minutes=mod(t1.minute, dt)) # "Floor" to whole # of minutes
        t0 = t0 + timedelta(seconds = 60-t0.second, minutes=dt-1-mod(t0.minute, dt)) # "Ceiling" to whole # of minutes
        return t0,t1
    
    
        
        
        
        
        
    def plot_los_winds(self, save=True):
        '''
        Make summary plot of the line of sight winds used in the inversion. If save=True, 
        save this plot in the appropriate folder. Return the figure.
        '''
        
        all_t = self.all_t
        instr_names = self.instr_names
        dat = self.dat
        
        Nd = max([len(dat[iname].keys()) for iname in instr_names]) # max num of looks at a site
        Ni = len(instr_names)

        fig = figure(figsize=(4*Ni,3*Nd), dpi = 150)
        instr_names.sort()
        for instri in range(Ni):
            instr_name = instr_names[instri]
            instrdat = dat[instr_name]
            direcs = instrdat.keys()
            direcs.sort()
            for direci in range(len(direcs)):
                direc = direcs[direci]
                u =  instrdat[direc]['losu']
                ue = instrdat[direc]['losue']
                t =  instrdat[direc]['t']
                c =  instrdat[direc]['c']
                subplot(Nd,Ni,1+direci*Ni + instri)
                #plot(t,u,'k.-')
                errorbar(t,u,fmt='k.-',yerr=ue)
                plot([all_t[0], all_t[-1]],[0,0],'k--', dashes=(3,3))
                plt.gca().xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H'))
                title('%s %s' % (instr_name, direc))
                ylim((-200,200))
                xlim((all_t[0], all_t[-1]))
                if direci==len(direcs)-1:
                    xlabel('UT')
                if instri==0:
                    ylabel('LoS velocity away [m/s]')

        plt.tight_layout()
        if save:
            savefig('%s%s_loswinds.png' % (self.dir, self.get_stub()))
        return fig
    
    
    
    
    
    
    
    def get_stub(self):
        '''
        Return the string used to identify this reconstruction (e.g., nation_2013_274)
        '''
        network_names_str = ''.join(self.network_names)
        return '%s_%i_%03i' % (network_names_str, self.year, self.doy)
       
        
      
      
      
      
        
    def write_ASCII(self):
        '''
        Create the ASCII text file specifying the wind field at all grid points and times.
        '''
        raise Exception('Not Implemented')
      
       
       
       
       
        
        
    def windfieldsolve(self, lam0, lam1, A, W, z):
        '''
        Find the wind field which minimizes the unconstrained form.
        '''
        Dg = self.Dg
        Dc = self.Dc
        
        Mleft = A.T*W.T*W*A
        Mright = self.Dg.T*self.Dg + lam1*self.Dc.T*self.Dc
        rhs = A.T*W.T*W*z
        M = Mleft + lam0*Mright
        uvwest = sp.linalg.spsolve(M,rhs)
        datacost = norm(W*(A*uvwest - z))**2
        gradcost = norm(Dg*uvwest)**2
        curvcost = norm(Dc*uvwest)**2

        return uvwest, datacost, gradcost, curvcost





    

    def solve_inner_problem(self, lam0guess, lam1, A, W, z):
        ''' 
        Given a curvature weighting (lam1), find the regularization weighting
        (lam0) in order to achieve the desired data cost, i.e. use the 
        discrepancy principle.
        INPUTS:
            lam0guess: starting value, to be optimized (the closer to the answer, the faster)
            lam1: constant
            A, W, z: observation matrix, weight matrix, observations
        '''
        Dg = self.Dg
        Dc = self.Dc
        L = self.L
        tol0 = self.tol0
        tol1 = self.tol1
        tol2 = self.tol2
        maxiters = self.maxiters
        desired_datacost = self.chi2mult*len(z)
        printstuffinner = self.printstuffinner
        
        # Zero-th step: Find a value of lam0 below the desired value
        numIters = 0
        datacost = inf
        lam0 = lam0guess
        while(datacost > desired_datacost):
            lam0 = lam0/L
            lastdatacost = datacost
            uvwest, datacost, gradcost, curvcost = self.windfieldsolve(lam0, lam1, A, W, z)
            regcost = norm(Dg*uvwest)**2 + lam1*norm(Dc*uvwest)**2
            if printstuffinner:
                self.message += 'DOWN:  %.3e: %.3f vs %.3e \n' % (lam0,datacost, regcost)
            numIters += 1
            if numIters > maxiters:
                uvwest[:] = nan
                raise Exception('First Line search failed')

        # First step: Find two values of c that straddle the desired cost                
        numIters = 0
        datacost = 0
        while(datacost < desired_datacost):
            lam0 = lam0*L
            lastdatacost = datacost
            uvwest, datacost, gradcost, curvcost = self.windfieldsolve(lam0, lam1, A, W, z)
            regcost = norm(Dg*uvwest)**2 + lam1*norm(Dc*uvwest)**2
            if printstuffinner:
                self.message += 'UP:    %.3e: %.3f vs %.3e \n' % (lam0,datacost, regcost)
            numIters += 1
            if numIters > maxiters:
                raise Exception('Second Line search failed')

        # Second step: Bisection Method to find optimal c
        numIters = 0
        c0 = lam0/L
        c1 = lam0
        cn = lam0 # This line is only important if no iterations are necessary
        f0 = lastdatacost - desired_datacost
        f1 = datacost - desired_datacost
        while abs((datacost-desired_datacost)/desired_datacost) > tol0: # within a tolerance
            cn = (c0 + c1)/2
            uvwest, datacost, gradcost, curvcost = self.windfieldsolve(cn, lam1, A, W, z)
            regcost = norm(Dg*uvwest)**2 + lam1*norm(Dc*uvwest)**2
            fn = datacost - desired_datacost
            if fn > 0:
                f1 = fn
                c1 = cn
            else:
                f0 = fn
                c0 = cn

            if printstuffinner:  
                self.message += 'BISEC: %.3e - %.3e: %.3f (%.3f) vs %.3e\n' % (c0,c1,datacost,desired_datacost, regcost)
            numIters +=1
            if numIters > maxiters or f0 > 0.0 or f1 < 0.0:
                raise Exception('Bisection failed')

        lam0 = cn

        return lam0, uvwest, datacost, gradcost, curvcost






    
    def solve_outer_problem(self, A, W, z):
        '''
        Find the optimal lam0 and lam1 to be used in the wind field estimate
        '''
        lam0interior = self.lam0interior
        lam1interior = self.lam1interior
        lam1small = self.lam1small
        lam0guess = self.lam0guess
        printstuffouter = self.printstuffouter
        desired_datacost = self.chi2mult*len(z)
        L2 = self.L2
        tol0 = self.tol0
        tol1 = self.tol1
        tol2 = self.tol2
        maxiters = self.maxiters

        # x == lam1
        # y == gradcost
        # w == lam0, just for the record

        # Zero-th step: Check if the solution is in the interior of the feasible set:

        uvwest, datacost, gradcost, curvcost = self.windfieldsolve(lam0interior, lam1interior, A, W, z)
        if datacost <= desired_datacost:
            # Stop. A (potentially not unique, but still smooth) solution has been found.
            if printstuffouter:
                self.message += 'Zero-structure solution fits data with datacost < desired_datacost' + \
                      ' (%f < %f). Halting.\n' % (datacost, desired_datacost)
            return lam0interior, lam1interior
        elif printstuffouter:
            self.message += 'Zero-structure solution does not fit data' + \
                      ' (%f > %f). Continuing.\n' % (datacost, desired_datacost)
        # If that passed, we know the solution is on the boundary of the feasible set.
        # Continue with the ad-hoc algorithm to find the optimal lam1.

        # First, do line search to find point where gradcost starts increasing
        xo = lam1small
        wo, _ , _ , yo, _ = self.solve_inner_problem(lam0guess, xo, A, W, z)
        lam1v = [xo]
        gradcostv = [yo]
        lam0v = [wo]
        diffgc = 0
        niters = 0
        while (diffgc < tol2):
            xn = xo * L2
            try:
                wn, _ , _ , yn, _ = self.solve_inner_problem(wo, xo, A, W, z)
            except Exception as e:
                raise
            lam1v.append(xo)
            gradcostv.append(yn)
            lam0v.append(wn)
            diffgc = log10(yn) - log10(yo)
            xo = xn
            yo = yn
            wo = wn
            if printstuffouter:
                self.message += 'LAM1UP: lam1: %e\t gradcost: %e \tgradcostdiff: %e (%e)\n' % (xn, yn, diffgc, tol2)
            if niters > maxiters:
                raise Exception('First outer line search reached maxiters')
            niters += 1
        x0 = lam1v[-2] # save the left endpoint
        y0 = gradcostv[-2]
        if printstuffouter:
            self.message += '----\n'


        # Second, continue line search to find point where gradcost stops increasing
        niters = 0
        while (diffgc >= tol2):
            xn = xo * L2
            try:
                wn, _ , _ , yn, _ = self.solve_inner_problem(wo, xo, A, W, z)
            except Exception as e:
                raise
                # TODO: fail gracefully, nan results, and print the result to the log
            lam1v.append(xo)
            gradcostv.append(yn)
            lam0v.append(wn)
            diffgc = log10(yn) - log10(yo)
            xo = xn
            yo = yn
            wo = wn # use solution for lam0 this time as initial guess for next time
            if printstuffouter:
                self.message += 'LAM1UP: lam1: %e\t gradcost: %e \tgradcostdiff: %e (%e)\n' % (xn, yn, diffgc, tol2)
            if niters > maxiters:
                raise Exception('Second outer line search reached maxiters')
            niters += 1
        x1 = xo # save right endpoint
        y1 = yo
        if printstuffouter:
            self.message += '----\n'


        # Use bisection until the desired value of lam1 (x) is found.
        # Use geometrical average in bisection.
        #y_desired = (y0+y1)/2
        y_desired = sqrt(y0*y1)

        # Set up first step of loop
        xn = sqrt(x0*x1)
        wn, _ , _ , yn , _ = self.solve_inner_problem(lam0guess, xn, A, W, z)
        lam1v.append(xn)
        lam0v.append(wn)
        gradcostv.append(yn)
        if printstuffouter:
            self.message += 'BISEC: lam1 endpoints: %e - %e. Desired gradcost: %e\n' % (x0, x1, y_desired)


        numIters = 0
        # Bisection algorithm
        while abs(yn - y_desired) > tol1*y_desired:

            if yn > y_desired:
                x1 = xn
                y1 = yn
            else:
                x0 = xn
                y0 = yn

            xn = sqrt(x0*x1)
            try:
                # Use last value of lam0 as guess to speed up inner inversion
                wn, _ , _ , yn , _ = self.solve_inner_problem(wn, xn, A, W, z)
            except Exception as e:
                raise
            lam1v.append(xn)
            lam0v.append(wn)
            gradcostv.append(yn)
            numIters += 1
            if printstuffouter:
                self.message += 'BISEC: lam1: %e\t gradcost: %e (%e) \t lam0: %e\n' % (xn, yn, y_desired, wn)

            if numIters > maxiters or y1 < y_desired or y0 > y_desired:
                raise Exception('Outer optimization (for lam1) failed') 
        lam0 = wn
        lam1 = xn 
        return lam0, lam1

   
   
   
   
    

    # Define function to estimate uncertainty
    def estimate_uncertainty(lam0, lam1, A, W, z, Dg, Dc, Nmc):
        '''
        Estimate the uncertainty in the wind field estimate. Use the
        Nicolls approach, intepreting the estimate as a MAP estimate,
        and using the a posteriori covariance as the estimate covariance.
        '''
        
        Dg = self.Dg
        Dc = self.Dc
        Nmc = self.Nmc
        Nx = self.Nx
        Ny = self.Ny
        
        Mleft = A.T*W.T*W*A
        Mright = Dg.T*Dg + lam1*Dc.T*Dc
        M = Mleft + lam0*Mright   

        Sig_uvw = linalg.inv(M.todense())
        uvwerr = diag(Sig_uvw)

        uerr = uvwerr[:Nx*Ny]
        verr = uvwerr[Nx*Ny:2*Nx*Ny]
        Uerr = reshape(uerr, (Ny,Nx))
        Verr = reshape(verr, (Ny,Nx))
        Werr = zeros((Ny,Nx))
        if estimate_vertical_wind:
            werr = uvwerr[2*Nx*Ny:]
            Werr = reshape(werr, (Ny,Nx)) 
        return Uerr, Verr, Werr






    def display_point(self, lon, lat, lonsm, latsm):
        '''
        Return whether or not the given point should be displayed.
        lon, lat: point under consideration
        lonsm, latsm: vectors of measurement locations
        distthresh: maximum distance from a measurement location
        returns True or False
        '''
        dx = (lon - lonsm)*cos(lat*pi/180)*111.
        dy = (lat - latsm)*111.
        r = sqrt(dx**2 + dy**2)
        if len(r)==0:
            return False
        return min(r) < self.distthresh






    
    def run_inversion(self):
        ''' 
        Run the inversion to estimate the regional wind field. This function doesn't plot
        anything.
        '''
        # Load necessary class variables so I don't have to type "self." for everything
        all_t = self.all_t
        Nx = self.Nx
        Ny = self.Ny
        dat = self.dat
        timestep_fit = self.timestep_fit
        instr_names = self.instr_names
        vertweight = self.vertweight
        minlocs = self.minlocs
        chi2mult = self.chi2mult
        tol0 = self.tol0
        tol1 = self.tol1
        L = self.L
        L2 = self.L2
        maxiters = self.maxiters 
        lam1small = self.lam1small
        lam0guess = self.lam0guess
        lam0interior = self.lam0interior
        lam1interior = self.lam1interior
        distthresh = self.distthresh
        oldthresh = self.oldthresh
        
        self.message += '\nStarting Inversion\n'
        
        def tand(x):
            return tan(pi/180.*x)
        def cosd(x):
            return cos(pi/180.*x)
        def sind(x):
            return sin(pi/180.*x)

        # Define grid in terms of latitude and longitude
        # TODO: Is it worth going to x,y or can we redefine gradients appropriately? Probably the latter.
        lons = array([dat[iname][dname]['lon'] for iname in instr_names for dname in dat[iname].keys()])
        lats = array([dat[iname][dname]['lat'] for iname in instr_names for dname in dat[iname].keys()])
        lonmin = lons.min() - 2.5
        lonmax = lons.max() + 2.5
        latmin = lats.min() - 2.5
        latmax = lats.max() + 2.5
        xs = linspace(lonmin, lonmax, Nx)
        ys = linspace(latmin, latmax, Ny)

        [X,Y] = meshgrid(xs,ys)
        x = reshape(X,(Nx*Ny,))
        y = reshape(Y,(Nx*Ny,))

        # Rectify to grid locations
        for iname in instr_names:
            instrdat = dat[iname]
            for direc in instrdat.keys():
                thisx = instrdat[direc]['lon']
                thisy = instrdat[direc]['lat']
                locidx = argmin((x-thisx)**2 + (y-thisy)**2)
                instrdat[direc]['rectlon'] = x[locidx]
                instrdat[direc]['rectlat'] = y[locidx]
                instrdat[direc]['grididx'] = locidx


        # Preconstruct the D matrix (regularization penalty matrix)

        # Dxx measures x curvature, Dyy measures y curvature, Dxy measures xy curvature, Dyx measures yx curvature
        # Dx measures y gradient, Dy measures y gradient
        # D matrix will comprise these in some ratio.

        # Construct Dxx
        Di = zeros(9*(Nx-2)*Ny)
        Dj = zeros(9*(Nx-2)*Ny)
        Dv = zeros(9*(Nx-2)*Ny)
        count = 0
        for idx in range(Nx*Ny):
            i,j = unravel_index(idx, (Ny,Nx))
            ileft = i
            iright = i
            jleft = j-1
            jright = j+1
            if (ileft >= 0 and ileft < Ny and iright >=0 and iright < Ny and \
                jleft >= 0 and jleft < Nx and jright >=0 and jright < Nx):
                dx = (X[i,jright] - X[i,j]) * 111. * cosd(Y[i,j])
                dy = (Y[iright,j] - Y[i,j]) * 111.        
                idxleft = ravel_multi_index((ileft,jleft),(Ny,Nx))
                idxright = ravel_multi_index((iright,jright),(Ny,Nx))
                Di[count:count+3] = [idx,idx,idx]
                Dj[count:count+3] = [idxleft,idx,idxright]
                Dv[count:count+3] = array([-1,2,-1])/dx**2
                count += 3
                Di[count:count+3] = [Nx*Ny + idx,    Nx*Ny + idx,Nx*Ny + idx]
                Dj[count:count+3] = [Nx*Ny + idxleft,Nx*Ny + idx,Nx*Ny + idxright]
                Dv[count:count+3] = array([-1,2,-1])/dx**2
                count += 3
                Di[count:count+3] = [2*Nx*Ny + idx,    2*Nx*Ny + idx,2*Nx*Ny + idx]
                Dj[count:count+3] = [2*Nx*Ny + idxleft,2*Nx*Ny + idx,2*Nx*Ny + idxright]
                Dv[count:count+3] = array([-1,2,-1])/dx**2 * vertweight
                count += 3
        Dxx = sp.coo_matrix((Dv,[Di,Dj]), (3*Nx*Ny,3*Nx*Ny))

        # Construct Dyy
        Di = zeros(9*(Nx-2)*Ny)
        Dj = zeros(9*(Nx-2)*Ny)
        Dv = zeros(9*(Nx-2)*Ny)
        count = 0
        for idx in range(Nx*Ny):
            i,j = unravel_index(idx, (Ny,Nx))
            ileft = i-1
            iright = i+1
            jleft = j
            jright = j
            if (ileft >= 0 and ileft < Ny and iright >=0 and iright < Ny and \
                jleft >= 0 and jleft < Nx and jright >=0 and jright < Nx):
                dx = (X[i,jright] - X[i,j]) * 111. * cosd(Y[i,j])
                dy = (Y[iright,j] - Y[i,j]) * 111. 
                idxleft = ravel_multi_index((ileft,jleft),(Ny,Nx))
                idxright = ravel_multi_index((iright,jright),(Ny,Nx))
                Di[count:count+3] = [idx,idx,idx]
                Dj[count:count+3] = [idxleft,idx,idxright]
                Dv[count:count+3] = array([-1,2,-1])/dy**2
                count += 3
                Di[count:count+3] = [Nx*Ny + idx,    Nx*Ny + idx,Nx*Ny + idx]
                Dj[count:count+3] = [Nx*Ny + idxleft,Nx*Ny + idx,Nx*Ny + idxright]
                Dv[count:count+3] = array([-1,2,-1])/dy**2
                count += 3
                Di[count:count+3] = [2*Nx*Ny + idx,    2*Nx*Ny + idx,2*Nx*Ny + idx]
                Dj[count:count+3] = [2*Nx*Ny + idxleft,2*Nx*Ny + idx,2*Nx*Ny + idxright]
                Dv[count:count+3] = array([-1,2,-1])/dy**2 * vertweight
                count += 3
        Dyy = sp.coo_matrix((Dv,[Di,Dj]), (3*Nx*Ny,3*Nx*Ny))

        # Construct Dxy
        Di = zeros(9*(Nx-2)*(Ny-2))
        Dj = zeros(9*(Nx-2)*(Ny-2))
        Dv = zeros(9*(Nx-2)*(Ny-2))
        count = 0
        for idx in range(Nx*Ny):
            i,j = unravel_index(idx, (Ny,Nx))
            ileft = i-1
            iright = i+1
            jleft = j-1
            jright = j+1
            if (ileft >= 0 and ileft < Ny and iright >=0 and iright < Ny and \
                jleft >= 0 and jleft < Nx and jright >=0 and jright < Nx):
                dx = (X[i,jright] - X[i,j]) * 111. * cosd(Y[i,j])
                dy = (Y[iright,j] - Y[i,j]) * 111. 
                idxleft = ravel_multi_index((ileft,jleft),(Ny,Nx))
                idxright = ravel_multi_index((iright,jright),(Ny,Nx))
                Di[count:count+3] = [idx,idx,idx]
                Dj[count:count+3] = [idxleft,idx,idxright]
                Dv[count:count+3] = array([-1,2,-1])/(dx*dy)
                count += 3
                Di[count:count+3] = [Nx*Ny + idx,    Nx*Ny + idx,Nx*Ny + idx]
                Dj[count:count+3] = [Nx*Ny + idxleft,Nx*Ny + idx,Nx*Ny + idxright]
                Dv[count:count+3] = array([-1,2,-1])/(dx*dy)
                count += 3
                Di[count:count+3] = [2*Nx*Ny + idx,    2*Nx*Ny + idx,2*Nx*Ny + idx]
                Dj[count:count+3] = [2*Nx*Ny + idxleft,2*Nx*Ny + idx,2*Nx*Ny + idxright]
                Dv[count:count+3] = array([-1,2,-1])/(dx*dy) * vertweight
                count += 3
        Dxy = sp.coo_matrix((Dv,[Di,Dj]), (3*Nx*Ny,3*Nx*Ny))

        # Construct Dyx
        Di = zeros(9*(Nx-2)*(Ny-2))
        Dj = zeros(9*(Nx-2)*(Ny-2))
        Dv = zeros(9*(Nx-2)*(Ny-2))
        count = 0
        for idx in range(Nx*Ny):
            i,j = unravel_index(idx, (Ny,Nx))
            ileft = i-1
            iright = i+1
            jleft = j+1
            jright = j-1
            if (ileft >= 0 and ileft < Ny and iright >=0 and iright < Ny and \
                jleft >= 0 and jleft < Nx and jright >=0 and jright < Nx):
                dx = (X[i,jright] - X[i,j]) * 111. * cosd(Y[i,j])
                dy = (Y[iright,j] - Y[i,j]) * 111. 
                idxleft = ravel_multi_index((ileft,jleft),(Ny,Nx))
                idxright = ravel_multi_index((iright,jright),(Ny,Nx))
                Di[count:count+3] = [idx,idx,idx]
                Dj[count:count+3] = [idxleft,idx,idxright]
                Dv[count:count+3] = array([-1,2,-1])/(dx*dy)
                count += 3
                Di[count:count+3] = [Nx*Ny + idx,    Nx*Ny + idx,Nx*Ny + idx]
                Dj[count:count+3] = [Nx*Ny + idxleft,Nx*Ny + idx,Nx*Ny + idxright]
                Dv[count:count+3] = array([-1,2,-1])/(dx*dy)
                count += 3
                Di[count:count+3] = [2*Nx*Ny + idx,    2*Nx*Ny + idx,2*Nx*Ny + idx]
                Dj[count:count+3] = [2*Nx*Ny + idxleft,2*Nx*Ny + idx,2*Nx*Ny + idxright]
                Dv[count:count+3] = array([-1,2,-1])/(dx*dy) * vertweight
                count += 3
        Dyx = sp.coo_matrix((Dv,[Di,Dj]), (3*Nx*Ny,3*Nx*Ny))


        # D1 measures x gradient, D2 measures y gradient
        # (u contains u and v)

        # Solution is xest = (AT*A + c*DT*D)^-1 * AT*z
        # Cov-transformed version:
        # xest = (AT*WT*W*A + c*DT*D)^-1 * AT*WT*z

        # Construct Dx
        Di = zeros(6*(Nx-1)*Ny)
        Dj = zeros(6*(Nx-1)*Ny)
        Dv = zeros(6*(Nx-1)*Ny)
        count = 0
        for idx in range(Nx*Ny):
            i,j = unravel_index(idx, (Ny,Nx))
            newi = i
            newj = j-1
            if (newj >= 0 and newi >= 0):
                dx = (X[i,newj] - X[i,j]) * 111. * cosd(Y[i,j])
                dy = (Y[newi,j] - Y[i,j]) * 111. 
                idxleft = ravel_multi_index((newi,newj),(Ny,Nx))
                Di[count:count+2] = [idx,idx]
                Dj[count:count+2] = [idx,idxleft]
                Dv[count:count+2] = array([1,-1])/dx
                count += 2
                Di[count:count+2] = [Nx*Ny + idx, Nx*Ny + idx]
                Dj[count:count+2] = [Nx*Ny + idx, Nx*Ny + idxleft]
                Dv[count:count+2] = array([1,-1])/dx
                count += 2
                Di[count:count+2] = [2*Nx*Ny + idx, 2*Nx*Ny + idx]
                Dj[count:count+2] = [2*Nx*Ny + idx, 2*Nx*Ny + idxleft]
                Dv[count:count+2] = array([1,-1])/dx * vertweight
                count += 2
        Dx = sp.coo_matrix((Dv,[Di,Dj]), (3*Nx*Ny,3*Nx*Ny))


        # Construct Dy
        Di = zeros(6*(Nx-1)*Ny)
        Dj = zeros(6*(Nx-1)*Ny)
        Dv = zeros(6*(Nx-1)*Ny)
        count = 0
        for idx in range(Nx*Ny):
            i,j = unravel_index(idx, (Ny,Nx))
            newi = i-1
            newj = j
            if (newj >= 0 and newi >= 0):
                dx = (X[i,newj] - X[i,j]) * 111. * cosd(Y[i,j])
                dy = (Y[newi,j] - Y[i,j]) * 111. 
                idxleft = ravel_multi_index((newi,newj),(Ny,Nx))
                Di[count:count+2] = [idx,idx]
                Dj[count:count+2] = [idx,idxleft]
                Dv[count:count+2] = array([1,-1])/dy
                count += 2
                Di[count:count+2] = [Nx*Ny + idx, Nx*Ny + idx]
                Dj[count:count+2] = [Nx*Ny + idx, Nx*Ny + idxleft]
                Dv[count:count+2] = array([1,-1])/dy
                count += 2
                Di[count:count+2] = [2*Nx*Ny + idx, 2*Nx*Ny + idx]
                Dj[count:count+2] = [2*Nx*Ny + idx, 2*Nx*Ny + idxleft]
                Dv[count:count+2] = array([1,-1])/dy * vertweight
                count += 2
        Dy = sp.coo_matrix((Dv,[Di,Dj]), (3*Nx*Ny,3*Nx*Ny))

        # Construct Iw (Only used if minimize_vertical_wind=True)
        Di = zeros(Nx*Ny)
        Dj = zeros(Nx*Ny)
        Dv = zeros(Nx*Ny)
        count = 0
        for idx in range(Nx*Ny):
            i,j = unravel_index(idx, (Ny,Nx))
            Di[count] = 2*Nx*Ny + idx
            Dj[count] = 2*Nx*Ny + idx
            Dv[count] = 1
            count += 1
        DIw = sp.coo_matrix((Dv,[Di,Dj]), (3*Nx*Ny,3*Nx*Ny))


        Dg = sp.vstack((Dx,Dy))
        Dc = sp.vstack((Dxx,Dyy,Dxy,Dyx))

        # Store the important variables
        self.Dg = Dg
        self.Dc = Dc
        
        # Now that setup is over, actually run the inversion
        windfield = {}
        windfield['t'] = []
        windfield['lat'] = []
        windfield['lon'] = []
        windfield['u'] = []
        windfield['v'] = []
        windfield['w'] = []
        windfield['lat_measured'] = []
        windfield['lon_measured'] = []
        windfield['u_sigma'] = []
        windfield['v_sigma'] = []
        windfield['w_sigma'] = []
        windfield['lam0'] = []
        windfield['lam1'] = []


        t0, t1 = self.get_start_stop_time(timestep_fit)
        timestepvec = range(0, int((t1-t0).total_seconds())+1, timestep_fit)

        numfailed = 0

        for timeindex in range(len(timestepvec)):

            # Construct the observation matrix

            ti = t0 + timedelta(seconds=timestepvec[timeindex]) # time for inversion

            Ai = []
            Aj = []
            Av = []
            z  = []
            zerr = []
            xm = []
            ym = []
            i = 0
            wvec = []
            wevec = []
            for iname in instr_names:
                instrdat = dat[iname]
                direcs = instrdat.keys()
                for direc in direcs:
                    t     = instrdat[direc]['t']
                    losu  = instrdat[direc]['losu']
                    losue = instrdat[direc]['losue']

                    # Interpolate the values to ti. Ignore this direction if it is too old.
                    # TODO: maybe split off the interpolation stuff to another function
                    si1 = sum(t < ti)
                    si0 = si1 - 1
                    if si0 >= 0 and si1 < len(t): # ti falls within start/stop times
                        tgap = (t[si1]-t[si0]).total_seconds()/60.
                        if tgap < oldthresh: # it's valid to interpolate
                            idx = instrdat[direc]['grididx']
                            az = instrdat[direc]['az']
                            ze = instrdat[direc]['ze']

                            # Interpolate
                            interplosu = losu[si0] + \
                                        (losu[si1]-losu[si0])*(ti-t[si0]).total_seconds()/(t[si1]-t[si0]).total_seconds()
                            interplosue= losue[si0] + \
                                        (losue[si1]-losue[si0])*(ti-t[si0]).total_seconds()/(t[si1]-t[si0]).total_seconds()
                            # If estimate_vertical_wind==False, ignore zenith measurements. (But save them for later)
                            ignore_this_one = not self.estimate_vertical_wind and direc=='Zenith'
                            if ignore_this_one:
                                wvec.append(interplosu)
                                wevec.append(interplosue)
                            # Ignore the measurements denoted above and also nan measurements.
                            if not isnan(interplosu) and not ignore_this_one:
                                z.append(interplosu)
                                zerr.append(interplosue)
                                xm.append(x[idx])
                                ym.append(y[idx])

                                # Add this measurement
                                # for eastward component:
                                Ai.append(i)
                                Aj.append(idx)
                                Av.append(sind(ze)*sind(az))
                                # for northward component:
                                Ai.append(i)
                                Aj.append(idx + Nx*Ny)
                                Av.append(sind(ze)*cosd(az))
                                # for zenith component:
                                Ai.append(i)
                                Aj.append(idx + 2*Nx*Ny)
                                Av.append(cosd(ze))

                                i += 1

            if len(z) > minlocs: # If there are enough measurements
                A = sp.coo_matrix((Av,[Ai,Aj]), (i, 3*Nx*Ny))
                z = array(z)
                zerr = array(zerr)
                success = True

                # Incorporate error bars
                W = sp.spdiags(1./zerr, 0, i, i) # weights, i.e. cov(z)^{-1/2}
                
                # Shave off the right third of the matrices, if we're not estimating the vertical wind.
                # Set the vertical wind to the weighted average of all vertical wind measurements
                # http://stackoverflow.com/questions/13352280/slicing-sparse-matrices-in-scipy-which-types-work-best
                wconst = zeros(Nx*Ny)
                if not self.estimate_vertical_wind:
                    # Take weighted average of vertical wind measurements
                    w = 0
                    if len(wvec) > 0:
                        weights = 1./array(wevec)**2
                        w = sum(weights*array(wvec))/sum(weights)
                    wconst = w*ones(Nx*Ny)
                    Aw = A.tocsc()[:,2*Nx*Ny:].tocoo() # Right third of A matrix, dealing with vertical wind.
                    # Move knowns to left-hand side of the equation: z = [Au, Av, Aw]*[u; v; w] ---> z - Aw*w = [Au, Av]*[u; v]
                    z = z - Aw.dot(wconst)
            
                    # Resize matrices
                    A = A.tocsc()[:,:2*Nx*Ny].tocoo()
                    Dc = Dc.tocsc()[:,:2*Nx*Ny].tocoo()
                    Dg = Dg.tocsc()[:,:2*Nx*Ny].tocoo()
                    self.Dc = Dc
                    self.Dg = Dg

                #Store matrices for debugging
                self.A = A
                self.W = W
                self.z = z
                
                 
                lam0 = nan
                lam1 = nan
                uvwest = nan*zeros(shape(A)[1])
                try:
                    # Solve the problem of finding the regularization parameters
                    lam0, lam1 = self.solve_outer_problem(A, W, z)
                    # Use the regularization parameters and solve the Tikhonov problem.
                    uvwest, datacost, gradcost, curvcost = self.windfieldsolve(lam0, lam1, A, W, z)
                    #print('%s' % ti)
                except Exception as e:
                    numfailed += 1
                    self.message += '%s ignored: Inversion failed: "%s".\n' % (ti,str(e))
                    raise

                # Estimate uncertainty
                U_sigma = nan*zeros((Ny,Nx))
                V_sigma = nan*zeros((Ny,Nx))
                W_sigma = nan*zeros((Ny,Nx))
                if self.estimate_uncertainty:
                    try:
                        U_sigma, V_sigma, W_sigma = self.estimate_uncertainty(lam0, lam1, A, W, z)
                    except Exception as e:
                        self.message += '%s uncertainty estimate failed: "%s". \n' % (ti,str(e))
                        raise

                # Rearrange
                uest = uvwest[:Nx*Ny]
                vest = uvwest[Nx*Ny:2*Nx*Ny]
                Uest = reshape(uest, (Ny,Nx))
                Vest = reshape(vest, (Ny,Nx))
                West = reshape(wconst,(Ny,Nx))
                if self.estimate_vertical_wind:
                    west = uvwest[2*Nx*Ny:]
                    West = reshape(west, (Ny,Nx)) 

                # Save
                # need: t, lat, lon, u, v, w, latmeasured, lonmeasured, 
                windfield['t'].append(ti)
                windfield['lat'].append(Y)
                windfield['lon'].append(X)
                windfield['u'].append(Uest)
                windfield['v'].append(Vest)
                windfield['w'].append(West)
                windfield['u_sigma'].append(U_sigma)
                windfield['v_sigma'].append(V_sigma)
                windfield['w_sigma'].append(W_sigma)
                windfield['lat_measured'].append(ym)
                windfield['lon_measured'].append(xm)
                windfield['lam0'].append(lam0)
                windfield['lam1'].append(lam1)


            else: # fewer than 4 measurements. append nans.
                windfield['t'].append(ti)
                windfield['lat'].append(Y)
                windfield['lon'].append(X)
                windfield['u'].append(nan*zeros((Ny,Nx)))
                windfield['v'].append(nan*zeros((Ny,Nx)))
                windfield['w'].append(nan*zeros((Ny,Nx)))
                windfield['u_sigma'].append(nan*zeros((Ny,Nx)))
                windfield['v_sigma'].append(nan*zeros((Ny,Nx)))
                windfield['w_sigma'].append(nan*zeros((Ny,Nx)))
                windfield['lat_measured'].append(ym)
                windfield['lon_measured'].append(xm)
                windfield['lam0'].append(nan)
                windfield['lam1'].append(nan)

                self.message += '%s ignored: not enough data points (%i).\n' % (ti,len(z))


        self.message += '%i/%i inversions failed.\n' % (numfailed, len(timestepvec))
        
        # Save results
        self.windfield = windfield
        self.lats = lats
        self.lons = lons
        
        if len(self.windfield['t']) == 0:
            self.success = False
            self.message += 'All inversions failed. Setting success=False.\n'

        # Write log file
        logfn = '%s%s_windfield_log.txt' % (self.dir, self.get_stub())
        f = open(logfn, 'w')
        f.write(self.message)
        f.close()





        

    def eval_field(self,lon,lat,dn, alt=250):
        '''
        Evaluate the wind model at the requested location and time
        INPUTS:
            lon - longitude of point (deg)
            lat - latitude of point (deg)
        OPTIONAL INPUTS:
            dn - datetime of interest.  Defaults to March 23, 2011 at 9:30 UT for no particular reason.
            alt - altitude of point (km).  Defaults to 250 km to simulate redline airglow layer.
            distthresh - the distance from the nearest sample point, above which a pixel should be set to nan (km)
        OUTPUTS:
            (U,V,W,UE,VE,WE,lam0,lam1,lonm,latm):
            U, V, W - the zonal, meridional, and vertical wind at the requested location and time (m/s).
            UE, VE, WE - the corresponding uncertainties (m/s)
            lam0, lam1 - the regularization parameters used for the inversion
            lonm,latm - the locations of the measurements used to generate the wind field estimate.
        HISTORY:
            Written by Jonathan J. Makela (jmakela@illinois.edu) on 17 March 2014.
            Adapted to use output from wind field estimator by Brian J. Harding on 23 May 2014.
                Using zero-th order interpolation (i.e., nearest neighbor), since the 
                calculated wind field is already dense (hopefully!)
        TODO:
            Make this faster by pre-calculating and storing a griddata interpolation object.
            This will also allow us to easily do higher-order interpolation.
        '''

        # Interpolate in time
        allt = self.windfield['t']
        # Convert to naive utc time
        allt = [t.astimezone(pytz.utc).replace(tzinfo=None) for t in allt]
        # Find nearest sample
        i = argmin([abs((t - dn).total_seconds()) for t in allt])

        # Interpolate in space
        lons = self.windfield['lon'][i]
        lats = self.windfield['lat'][i]
        idx = argmin(abs(lons-lon) + abs(lats-lat))
        j,k = unravel_index(idx, lats.shape)

        u = self.windfield['u'][i][j,k]
        v = self.windfield['v'][i][j,k]
        w = self.windfield['w'][i][j,k]
        u_sigma = self.windfield['u_sigma'][i][j,k]
        v_sigma = self.windfield['v_sigma'][i][j,k]
        w_sigma = self.windfield['w_sigma'][i][j,k]
        lam0 = self.windfield['lam0'][i]
        lam1 = self.windfield['lam1'][i]

        xm = array(self.windfield['lon_measured'][i])
        ym = array(self.windfield['lat_measured'][i])

        # Return nan if we're too far from any measurements    
        if not self.display_point(lon, lat, xm,  ym):
            u = nan
            v = nan
            w = nan
            u_sigma = nan
            v_sigma = nan
            w_sigma = nan

        return (u,v,w,u_sigma,v_sigma,w_sigma,lam0,lam1,xm,ym)
    
    
    
    
    
    
    
    def make_quicklook(self, timestep = None, save = True, show_vert_wind=False, scale=250.):
        '''
        Make a summary figure with many subplots, each showing the wind field at different 
        times throughout the night, plotted using quiver.
        OPTIONAL INPUTS:
            timestep: sec. Cadence with which to plot wind field. If None, use time step from the fit. 
                     (default None)
            save:     bool. If True, save the resulting figure to self.save_dir
            show_vert_wind: bool. If True, plot the vertical wind as a color.
            scale:    float, m/s. An arrow 1/10th of the figure width refers to this velocity
        RETURNS:
            figure, saved_filename
        '''
        if timestep is None:
            timestep = self.timestep_fit
        all_t = self.all_t
        distthresh = self.distthresh
        lons = self.lons
        lats = self.lats
        
        spnx = 4 # subplot x size (y size is automatically calculated)
        pivot = 'middle' # for quiver
        headwidth = 4 
        headlength = 4
        headaxislength= headlength-1
        minshaft = 2
        qcolor = 'k'
        Nq = 10 # approximately how many arrows across and down you want
        WN = 4 # upsample factor for vertical wind color plot
        # Figure out the color limits for vertical wind based on the data
        wmax = 30.
        for wmat in self.windfield['w']:
            mx = nanmax(abs(wmat))
            if mx > wmax:
                wmax = mx
        climits = (-wmax, wmax)
        
        lonmin = lons.min() - 2.5
        lonmax = lons.max() + 2.5
        latmin = lats.min() - 2.5
        latmax = lats.max() + 2.5
        xs = linspace(lonmin, lonmax, Nq)
        ys = linspace(latmin, latmax, Nq)
        [Xq,Yq] = meshgrid(xs,ys)

        t0, t1 = self.get_start_stop_time(timestep)
        timestepvec = range(0,int((t1-t0).total_seconds())+1, timestep)
        
        start_t = t0
        start_t = start_t.astimezone(pytz.utc).replace(tzinfo=None)
             
        spny = int(ceil(1.0*len(timestepvec)/spnx))
        figsize=(3.5*spnx+0.1, 4*spny+0.1)
        sc = figsize[0]/spnx/(10.*scale) # 1/10th of the plot width <==> 250 m/s, for quiver

        fig = figure(figsize=figsize, dpi=150)
        for timeindex in range(len(timestepvec)):

            subplot(spny,spnx,timeindex+1)
            t = start_t + timedelta(seconds = timestepvec[timeindex])

            # Load up the necessary values
            U = zeros((Nq,Nq))
            V = zeros((Nq,Nq))
            W = zeros((Nq,Nq))
            for i in range(shape(Xq)[0]):
                for j in range(shape(Xq)[1]):
                    (u,v,w,ue,ve,we,lam0,lam1,xm,ym) = self.eval_field(Xq[i,j], Yq[i,j], t)
                    U[i,j] = u
                    V[i,j] = v

            # Transform U,V to inches
            U = sc*U
            V = sc*V

            # setup Lambert Conformal basemap.
            m = Basemap(llcrnrlon=lons.min()-3,llcrnrlat=lats.min()-3,urcrnrlon=lons.max()+3,urcrnrlat=lats.max()+3,
                    projection='merc', area_thresh=1000,
                    resolution='i')#,lat_1=45.,lat_2=55,lat_0=40,lon_0=-85.)
            # draw coastlines and fill the continents (alpha value is set to partial transparancy because for some
            # reason the fill is over top the quivers used to denote the wind vectors
            m.drawcoastlines()
            #m.fillcontinents(alpha=.5)
            #m.drawstates()

            if show_vert_wind:
                # Sample the vertical wind field higher, if it's being plotted (for aesthetic reasons)
                [X4q,Y4q] = meshgrid(linspace(lonmin, lonmax, WN*Nq),linspace(latmin, latmax, WN*Nq))
                W4 = zeros((WN*Nq,WN*Nq))
                for i in range(shape(X4q)[0]):
                    for j in range(shape(X4q)[1]):
                        (u,v,w,_,_,_,_,_,xm,ym) = self.eval_field(X4q[i,j], Y4q[i,j], t)
                        W4[i,j] = w
                # mask so the plot looks nice
                W4 = ma.array(W4,mask=np.isnan(W4))
                xpt,ypt = m(X4q,Y4q)
                m.pcolormesh(xpt,ypt,W4,cmap='RdBu_r')
                clim(climits)
                cb = m.colorbar()
                cb.set_label('Vertical wind [m/s]', ha='center')


            # Plot Wind Field
            xpt,ypt = m(Xq,Yq)
            Q1 = m.quiver(xpt,ypt, U,V, angles='uv', scale_units='inches', scale=1, 
                          pivot=pivot, headwidth=headwidth, headlength=headlength,
                          minshaft=minshaft, headaxislength=headaxislength, color=qcolor)
            xq,yq = m(lons.min(),lats.min()-1.5)
            qk1 = plt.quiverkey(Q1, xq, yq, sc*scale, r'$%i\,\frac{m}{s}$'%scale, color='k', coordinates='data')

            #xpt,ypt = m(xm, ym)
            #m.scatter(xpt,ypt, c='w', s= 20, edgecolor='none')

            title('%s UT' % t.strftime('%Y/%m/%d %H:%M:%S'))
            m.drawparallels(np.arange(-45,46,5.),labels=[1,0,0,0],color='black',linewidth=0) # draw parallels
            m.drawmeridians(np.arange(0,360,5.),labels=[0,0,0,1],color='black',linewidth=0) # draw meridians

        #tight_layout()
        savefn = None
        if save:
            savefn = '%s%s_windfield_summary.png' % (self.dir, self.get_stub())
            savefig(savefn)
        return fig, savefn
      
      
      
      
      
      
    def make_quiver_movie(self, timestep=None, framerate=10, show_vert_wind=True):
        '''
        Make a movie showing how the wind field evolves over the night, using quiver for horizontal wind
        and color for the vertical wind.
        OPTIONAL INPUTS:
            timestep: sec. Cadence with which to plot wind field. If None, use time step from the fit. 
                     (default None)
            framerate: Hz. Frame rate for the output movie
            show_vert_wind: If true, use color to display the vertical wind estimate (default True)
        OUTPUTS:
            movie_fn: full path to the created file. Returns None if it failed.
        '''
        if timestep is None:
            timestep = self.timestep_fit
        all_t = self.all_t
        distthresh = self.distthresh
        lons = self.lons
        lats = self.lats
        
        Nq = 10 # approximately how many arrows across and down you want
        WN = 8 # upsample factor for vertical wind color plot
        figsize=(6,6) # 6,6
        sc = figsize[1]/(Nq*200.) # 1/Nq^th of the plot width <==> 250 m/s, for quiver
        pivot = 'middle' # for quiver
        headwidth = 4 # 3
        headlength = 4 # 3
        headaxislength= headlength-1
        minshaft = 2
        qcolor = 'k'

        # Figure out the color limits for vertical wind based on the data
        wmax = 30.
        for wmat in self.windfield['w']:
            mx = nanmax(abs(wmat))
            if mx > wmax:
                wmax = mx
        climits = (-wmax, wmax)
        
        # Make directory to store pngs
        pngdir = '%s%s' % (self.dir,'movie_quiver_pngs') 
        try: # make directory
            os.mkdir(pngdir)
        except Exception as e:
            pass
        try: # delete all files in that directory
            pngs = glob.glob('%s/*' % pngdir)
            for png in pngs:
                os.remove(png)
        except Exception as e:
            pass
            
        # Define output movie fn
        moviefn = '%s%s_windfield.mov' % (self.dir, self.get_stub())

        lonmin = lons.min() - 2.5
        lonmax = lons.max() + 2.5
        latmin = lats.min() - 2.5
        latmax = lats.max() + 2.5
        xs = linspace(lonmin, lonmax, Nq) # TODO: don't rely on lonmin rolling over from before
        ys = linspace(latmin, latmax, Nq)
        [Xq,Yq] = meshgrid(xs,ys)

        t0, t1 = self.get_start_stop_time(timestep)
        timestepvec = range(0,int((t1-t0).total_seconds()+1), timestep)
        start_t = t0
        start_t = start_t.astimezone(pytz.utc).replace(tzinfo=None)
        
        for timeindex in range(len(timestepvec)):

            ####### PLOT #######
            t = start_t + timedelta(seconds = timestepvec[timeindex])

            # Load up the necessary values
            Uest = zeros((Nq,Nq))
            Vest = zeros((Nq,Nq))
            Uerr = zeros((Nq,Nq))
            Verr = zeros((Nq,Nq))
            Werr = zeros((Nq,Nq))
            for i in range(shape(Xq)[0]):
                for j in range(shape(Xq)[1]):
                    (u,v,w,ue,ve,we,lam0,lam1,xm,ym) = self.eval_field(Xq[i,j], Yq[i,j], t)
                    Uest[i,j] = u
                    Vest[i,j] = v
                    Uerr[i,j] = ue
                    Verr[i,j] = ve
                    Werr[i,j] = we
            # Determine rms uncertainty
            uerms = sqrt(nanmean(Uerr**2))
            verms = sqrt(nanmean(Verr**2))
            werms = sqrt(nanmean(Werr**2))
            # Determine smoothing factor
            smooth_factor = log10(lam0/(1+lam1))+7
            if lam0 == 1e10:
                smooth_factor = inf

            fig = figure(num=None, figsize=figsize, dpi=300, facecolor='w', edgecolor='k')

            # Transform U,V to inches
            Uest = sc*Uest
            Vest = sc*Vest

            # setup Lambert Conformal basemap.
            m = Basemap(llcrnrlon=lons.min()-3,llcrnrlat=lats.min()-3,urcrnrlon=lons.max()+3,urcrnrlat=lats.max()+3,
                    projection='merc', area_thresh=1000,
                    resolution='i')
            m.drawcoastlines()
            #m.fillcontinents(alpha=.5)
            #m.drawstates()
            
            if show_vert_wind:
                # Sample the vertical wind field higher, if it's being plotted (for aesthetic reasons)
                [X4q,Y4q] = meshgrid(linspace(lonmin, lonmax, WN*Nq),linspace(latmin, latmax, WN*Nq))
                W4 = zeros((WN*Nq,WN*Nq))
                for i in range(shape(X4q)[0]):
                    for j in range(shape(X4q)[1]):
                        (u,v,w,_,_,_,_,_,xm,ym) = self.eval_field(X4q[i,j], Y4q[i,j], t)
                        W4[i,j] = w
                # mask so the plot looks nice
                W4 = ma.array(W4,mask=np.isnan(W4))
                xpt,ypt = m(X4q,Y4q)
                m.pcolormesh(xpt,ypt,W4,cmap='RdBu_r')
                clim(climits)
                cb = m.colorbar()
                cb.set_label('Vertical wind [m/s]', ha='center')

            # Plot Wind Field
            xpt,ypt = m(Xq,Yq)
            Q1 = m.quiver(xpt,ypt, Uest,Vest, angles='uv', scale_units='inches', scale=1, 
                          pivot=pivot, headwidth=headwidth, headlength=headlength,
                          minshaft=minshaft, headaxislength=headaxislength, color=qcolor)
            xq,yq = m(lons.min(),lats.min()-1.5)
            qk1 = plt.quiverkey(Q1, xq, yq, sc*200, r'$200\,\frac{m}{s}$', color='k', coordinates='data')

            # Plot Measurement locations
            xpt,ypt = m(xm, ym)
            m.scatter(xpt,ypt, c='k', s= 20, edgecolor='none')

            # Print text for uncertainties and bias (Turned off for now)
            #s = 'Uncertainty: (%2.0f, %2.0f, %2.0f) m/s     Smoothing: %4.1f' % (uerms, verms, werms, smooth_factor)
            #rectprops = dict(boxstyle='round', facecolor='w', linewidth=0.5)
            #fontprops = matplotlib.font_manager.FontProperties(family='monospace', size=6)
            #annotate(s, xy=(0.5,0.03), xycoords='axes fraction', fontproperties=fontprops, ha='center' ,bbox=rectprops)

            title('%s UT' % t.strftime('%Y/%m/%d %H:%M:%S'))
            m.drawparallels(np.arange(-45,46,5.),labels=[1,0,0,0],color='black',linewidth=0) # draw parallels
            m.drawmeridians(np.arange(0,360,5.),labels=[0,0,0,1],color='black',linewidth=0) # draw meridians

            savefig('%s/%05i.png' % (pngdir, timeindex), bbox_inches='tight') 
            # To save memory, close the figure.
            close(fig)
            
        # Make movie
        cmd = 'avconv -framerate %i -i "%s/%%05d.png" -y -vcodec qtrle %s' % (framerate,pngdir,moviefn)
        return_code = subprocess.call(cmd, shell=True)
        if return_code == 0: 
            return moviefn
        else: 
            return None
        
        
        
 
 
        
    def make_quiver_gif(self, timestep=1800, framerate=4, show_vert_wind=False):
        '''
        Make a gif showing how the wind field evolves over the night, using quiver for horizontal wind
        and color for the vertical wind, if desired.
        OPTIONAL INPUTS:
            timestep: sec. Cadence with which to plot wind field. If None, use time step from the fit. 
                     (default 1800)
            framerate: Hz. Frame rate for the output gif
            show_vert_wind: If true, use color to display the vertical wind estimate (default False)
        OUTPUTS:
            gif_fn: full path to the created file. Returns None if it failed.
        '''
        if timestep is None:
            timestep = self.timestep_fit
        all_t = self.all_t
        distthresh = self.distthresh
        lons = self.lons
        lats = self.lats
        
        Nq = 10 # approximately how many arrows across and down you want
        WN = 8 # upsample factor for vertical wind color plot
        figsize=(5,5) # 6,6
        sc = figsize[1]/(Nq*200.) # 1/Nq^th of the plot width <==> 250 m/s, for quiver
        pivot = 'middle' # for quiver
        headwidth = 4 # 3
        headlength = 4 # 3
        headaxislength= headlength-1
        minshaft = 2
        qcolor = 'k'

        # Figure out the color limits for vertical wind based on the data
        wmax = 30.
        for wmat in self.windfield['w']:
            mx = nanmax(abs(wmat))
            if mx > wmax:
                wmax = mx
        climits = (-wmax, wmax)
        
        # Make directory to store pngs
        pngdir = '%s%s' % (self.dir,'gif_quiver_pngs') 
        try: # make directory
            os.mkdir(pngdir)
        except Exception as e:
            pass
        try: # delete all files in that directory
            pngs = glob.glob('%s/*' % pngdir)
            for png in pngs:
                os.remove(png)
        except Exception as e:
            pass
            
        # Define output movie fn
        giffn = '%s%s_windfield.gif' % (self.dir, self.get_stub())

        lonmin = lons.min() - 2.5
        lonmax = lons.max() + 2.5
        latmin = lats.min() - 2.5
        latmax = lats.max() + 2.5
        xs = linspace(lonmin, lonmax, Nq) # TODO: don't rely on lonmin rolling over from before
        ys = linspace(latmin, latmax, Nq)
        [Xq,Yq] = meshgrid(xs,ys)

        t0, t1 = self.get_start_stop_time(timestep)
        timestepvec = range(0,int((t1-t0).total_seconds()+1), timestep)
        start_t = t0
        start_t = start_t.astimezone(pytz.utc).replace(tzinfo=None)
        
        for timeindex in range(len(timestepvec)):

            ####### PLOT #######
            t = start_t + timedelta(seconds = timestepvec[timeindex])

            # Load up the necessary values
            Uest = zeros((Nq,Nq))
            Vest = zeros((Nq,Nq))
            Uerr = zeros((Nq,Nq))
            Verr = zeros((Nq,Nq))
            Werr = zeros((Nq,Nq))
            for i in range(shape(Xq)[0]):
                for j in range(shape(Xq)[1]):
                    (u,v,w,ue,ve,we,lam0,lam1,xm,ym) = self.eval_field(Xq[i,j], Yq[i,j], t)
                    Uest[i,j] = u
                    Vest[i,j] = v
                    Uerr[i,j] = ue
                    Verr[i,j] = ve
                    Werr[i,j] = we
            # Determine rms uncertainty
            uerms = sqrt(nanmean(Uerr**2))
            verms = sqrt(nanmean(Verr**2))
            werms = sqrt(nanmean(Werr**2))
            # Determine smoothing factor
            smooth_factor = log10(lam0/(1+lam1))+7
            if lam0 == 1e10:
                smooth_factor = inf

            fig = figure(num=None, figsize=figsize, dpi=300, facecolor='w', edgecolor='k')

            # Transform U,V to inches
            Uest = sc*Uest
            Vest = sc*Vest

            # setup Lambert Conformal basemap.
            m = Basemap(llcrnrlon=lons.min()-3,llcrnrlat=lats.min()-3,urcrnrlon=lons.max()+3,urcrnrlat=lats.max()+3,
                    projection='merc', area_thresh=1000,
                    resolution='i')
            m.drawcoastlines()
            #m.fillcontinents(alpha=.5)
            #m.drawstates()
            
            if show_vert_wind:
                # Sample the vertical wind field higher, if it's being plotted (for aesthetic reasons)
                [X4q,Y4q] = meshgrid(linspace(lonmin, lonmax, WN*Nq),linspace(latmin, latmax, WN*Nq))
                W4 = zeros((WN*Nq,WN*Nq))
                for i in range(shape(X4q)[0]):
                    for j in range(shape(X4q)[1]):
                        (u,v,w,_,_,_,_,_,xm,ym) = self.eval_field(X4q[i,j], Y4q[i,j], t)
                        W4[i,j] = w
                # mask so the plot looks nice
                W4 = ma.array(W4,mask=np.isnan(W4))
                xpt,ypt = m(X4q,Y4q)
                m.pcolormesh(xpt,ypt,W4,cmap='RdBu_r')
                clim(climits)
                cb = m.colorbar()
                cb.set_label('Vertical wind [m/s]', ha='center')

            # Plot Wind Field
            xpt,ypt = m(Xq,Yq)
            Q1 = m.quiver(xpt,ypt, Uest,Vest, angles='uv', scale_units='inches', scale=1, 
                          pivot=pivot, headwidth=headwidth, headlength=headlength,
                          minshaft=minshaft, headaxislength=headaxislength, color=qcolor)
            xq,yq = m(lons.min()-0.5,lats.min()-1.9)
            qk1 = plt.quiverkey(Q1, xq, yq, sc*200, r'$200\,\frac{m}{s}$', color='k', coordinates='data')

            # Plot Measurement locations
            xpt,ypt = m(xm, ym)
            m.scatter(xpt,ypt, c='k', s= 20, edgecolor='none')

            # Print text for uncertainties and bias
            #s = 'Uncertainty: (%2.0f, %2.0f, %2.0f) m/s     Smoothing: %4.1f' % (uerms, verms, werms, smooth_factor)
            #rectprops = dict(boxstyle='round', facecolor='w', linewidth=0.5)
            #fontprops = matplotlib.font_manager.FontProperties(family='monospace', size=6)
            #annotate(s, xy=(0.5,0.03), xycoords='axes fraction', fontproperties=fontprops, ha='center' ,bbox=rectprops)

            title('%s UT' % t.strftime('%Y/%m/%d %H:%M:%S'))
            m.drawparallels(np.arange(-45,46,5.),labels=[1,0,0,0],color='black',linewidth=0) # draw parallels
            m.drawmeridians(np.arange(0,360,5.),labels=[0,0,0,1],color='black',linewidth=0) # draw meridians

            savefig('%s/%05i.png' % (pngdir, timeindex), bbox_inches='tight') 
            # To save memory, close the figure.
            close(fig)
            
        # Make gif
        delay = int(100./framerate)
        cmd = 'convert -delay %i %s/*.png %s' % (delay, pngdir, giffn)
        #print cmd
        return_code = subprocess.call(cmd, shell=True)
        if return_code == 0:
            return giffn
        else:
            return None
 
 
 
 
 
 
 
        
    def make_tracer_movie(self, timestep=90, framerate=10, show_vert_wind=False):
        '''
        Make a movie showing how the wind field evolves over the night, using tracers to indicate
        how air parcels are advected.
        OPTIONAL INPUTS:
            timestep: sec. Cadence with which to propagate tracers. (default 90, recommended)
            framerate: Hz. Frame rate for the output movie
            show_vert_wind: If true, use color to display the vertical wind estimate (default False)
        TODO:
            Automatically scale "res" and "npts" when a different timestep is used.
        '''

        # Define 3 helper functions:
        def trace_flow(lon,lat,scale=10.,dn=datetime(2011,3,23,9,30)):
            '''
            Trace a particle from a specified point through one step in a wind field.
            INPUTS:
                lon - longitude of point (deg)
                lat - latitude of point (deg)
            OPTIONAL INPUTS:
                scale - a scaling factor for following the point through the wind field.
                dn - datetime of interest.  Defaults to March 23, 2011 at 9:30 UT for no particular reason.
            OUTPUTS:
                ([lon, lat, speed, lon_next, lat_next) - the current location (lon, lat), the speed of the wind at the 
                                          requested location and time (speed), and the next location of
                                          the particle given the wind field and requested scaling factor
                                          (lon_next, lat_next).
            HISTORY:
                Written by Jonathan J. Makela (jmakela@illinois.edu) on 17 March 2014.
            '''
            
            # We will save the information in a deque structure
            a = deque()
            
            # Evaluate the wind field at the requested location and time
            # TODO: should probably pass in windfield instead of making it a global variable. 
            (U,V,W,_,_,_,_,_,_,_) = self.eval_field(lon,lat,dn)
            
            # Generate the output
            a.append([lon,lat,np.sqrt(U*U+V*V),lon+(U/111.0e3/np.cos(lat*np.pi/180.))/scale, lat+(V/111.0e3)/scale])
            return a

        def plot_flow(A,m):
            '''
            Plots a line varying the color based on a third parameter.  Code modified from
            http://matplotlib.org/examples/pylab_examples/multicolored_line.html
            INPUTS:
                A - Matrix containing (x,y) pairs of points in the line in (A[0,:], A[1,:] corresponding to lon,lat)
                    as well as a value to scale the intensity by (A[2,:], in our case the speed at the x,y location).
                m - the Basemap that we are using to plot data on
            OUTPUTS:
                Plots on the currently-defined plot axis
            HISTORY:
                Written by Jonathan J. Makela (jmakela@illinois.edu) on 17 March 2014.
            '''
            
            # Transform from longitude/latitude to (x,y) given the defined Basemap
            xp,yp = m(A[0,:],A[1,:])
            
            # Create line segments
            p = np.array([xp,yp]).T.reshape(-1,1,2)
            s = np.concatenate([p[:-1],p[1:]],axis=1)
            
            # Set the colormap to gray_r so that small values are white and large values are black.
            # Normalize to speeds between 0 and 100 m/s
            lws = 0.4*ones(shape(s)[0])
            l = LineCollection(s,cmap=plt.get_cmap('gray_r'),norm=plt.Normalize(0,100),linewidths=lws)
            #l.set_array(A[2,:])
            l.set_array(100*ones_like(A[2,:]))
            plt.gca().add_collection(l)

        def seed_flow(scale, dn, lonmin, lonmax, latmin, latmax, step):
            '''
            Seeds a wind field by dropping particles randomly within defined cells on the map
            OPTION INPUTS:
                scale - a scaling factor for following the point through the wind field.
                dn - datetime of interest. (naive, presumably UT)
                lonmin - left edge of map (degrees longitude)
                lonmax - right edge of map (degrees longitude)
                latmin - bottom edge of map (degeres latitude)
                latmax - top edge of map (degrees latitude)
                res    - grid size in degrees
            OUTPUTS:
                out - the result of trace_flow at the radomly chosen point
            HISTORY:
                Written by Jonathan J. Makela (jmakela@illinois.edu) on 17 March 2014.
            '''
            
            # The output queue
            out = deque()
            
            # Step through the grid (TODO: ranges should be input parameters)
            for la in arange(latmin,latmax,res):
                for lo in arange(lonmin,lonmax,res):
            
                    # Chose the lat/lon from a uniform distribution in the current cell
                    lat = np.random.uniform(low=la-res/2., high=la+res/2.)
                    lon = np.random.uniform(low=lo-res/2., high=lo+res/2.)

                    # Trace this particle one step in the wind field
                    out.append(trace_flow(lon,lat,scale=scale,dn=dn))
            
            return out
            
                       
        ## "scale" scales the spatial step sizes for a particle's motion
        # "npts" defines the number of time steps saved for each particle.  
        exagg = 3. # how much to exaggerate the transport of air parcels
        scale = 1./(timestep*exagg)
        npts = 5 # this will need to be adjusted with the time step (5 for timestep=90 seconds seems good)
        res  = 1.3 # deg, size of grid for flow tracers, scaled for size of FoV (prev 1.0)
            
        all_t = self.all_t
        distthresh = self.distthresh
        lons = self.lons
        lats = self.lats 
        Nx = self.Nx
        Ny = self.Ny    
        
        lonmin = lons.min() - 2.5
        lonmax = lons.max() + 2.5
        latmin = lats.min() - 2.5
        latmax = lats.max() + 2.5
                
        # Scale res to account for size of field of view, and for wind speed
        res = res * sqrt((latmax-latmin)*(lonmax-lonmin)/232.0)
        allwind = abs(array([self.windfield[d][i][Nx/2,Ny/2] for d in ['u','v'] for i in range(len(self.windfield[d]))]))
        maxwind = prctile(allwind[~isnan(allwind)],95)
        res = res * sqrt(maxwind/200.0)
        
        # Figure out the color limits for vertical wind based on the data
        wmax = 30
        for wmat in self.windfield['w']:
            mx = nanmax(abs(wmat))
            if mx > wmax:
                wmax = mx
        climits = (-wmax, wmax)


        # Make directory to store pngs
        pngdir = '%s%s' % (self.dir,'movie_tracer_pngs') 
        try: # make directory
            os.mkdir(pngdir)
        except Exception as e:
            pass
        try: # delete all files in that directory
            pngs = glob.glob('%s/*' % pngdir)
            for png in pngs:
                os.remove(png)
        except Exception as e:
            pass
            
        # Define output movie fn
        moviefn = '%s%s_windtracer.mov' % (self.dir, self.get_stub())
        
        # Define start/stop time
        t0, t1 = self.get_start_stop_time(timestep)
        timestepvec = range(0,int((t1-t0).total_seconds()+1), timestep)
        start_t = t0
        start_t = start_t.astimezone(pytz.utc).replace(tzinfo=None)

        # Generate the grow and shrink lists
        grow = seed_flow(scale, start_t, lonmin, lonmax, latmin, latmax, res)
        shrink = deque()
        count = 0

        for timeindex in range(len(timestepvec)):
            
            t = start_t + timedelta(seconds = timestepvec[timeindex])
            #print t

            # Create the plot
            fig = plt.figure(figsize=(6,6))
            fig.hold(True)
             
            # Define the Basemap
            m = Basemap(llcrnrlon=lons.min()-3,llcrnrlat=lats.min()-3,urcrnrlon=lons.max()+3,urcrnrlat=lats.max()+3,
                    projection='merc', area_thresh=1000,
                    resolution='i')#,lat_1=45.,lat_2=55,lat_0=40,lon_0=-85.)
            m.drawcoastlines(linewidth=0.5)
            
            # Find out how many elements are in the shrink list
            ns = shape(shrink)[0]
            
            # Go through each element in the shrink list
            for i in range(ns):
                sl = shrink.popleft()
                
                # This is where we would plot a shrinking list
                if shape(sl)[0] > 1:
                    plot_flow(array(sl).T,m)
                
                # Remove the first value put into the list
                sl.popleft()
                if np.shape(sl)[0] > 0:
                    shrink.append(sl)
            
            # Find out how many elements are in the grow list
            ng = shape(grow)[0]
            
            # Go through each element in the grow list
            for i in range(ng):
                gl = grow.popleft()
                
                # This is where we would plot a growing list
                if shape(gl)[0] > 1:
                    plot_flow(array(gl).T,m)
                
                # See how long this list is
                if np.shape(gl)[0] < npts:
                    # The list needs to grow
                    gl.append(trace_flow(gl[-1][3],gl[-1][4],scale=scale,dn=t)[0])
            
                    # Put this back in the grow list
                    grow.append(gl)
                else:
                    # The list is too long, it needs to go into shrink
                    shrink.append(gl)
                    
            # Add a new initiation to the grow list
            g2 = seed_flow(scale, t, lonmin, lonmax, latmin, latmax, res)
            for i in range(shape(g2)[0]):
                grow.append(g2.popleft())
                
            # Plot vertical wind
            if show_vert_wind:
                # Sample the vertical wind field higher, if it's being plotted (for aesthetic reasons)
                [X4q,Y4q] = meshgrid(linspace(lonmin, lonmax, Nx),linspace(latmin, latmax, Ny))
                W4 = zeros(shape(X4q))
                for i in range(shape(X4q)[0]):
                    for j in range(shape(X4q)[1]):
                        (u,v,w,_,_,_,_,_,xm,ym) = self.eval_field(X4q[i,j], Y4q[i,j], t)
                        W4[i,j] = w
                # mask so the plot looks nice
                W4 = ma.array(W4,mask=np.isnan(W4))
                xpt,ypt = m(X4q,Y4q)
                m.pcolormesh(xpt,ypt,W4,cmap='RdBu_r')
                clim(climits)
                cb = m.colorbar()
                cb.set_label('Vertical wind [m/s]', ha='center')
            
            plt.title('%s UT' % t)
            m.drawparallels(np.arange(-45,46,5.),labels=[1,0,0,0],color='black',linewidth=0) # draw parallels
            m.drawmeridians(np.arange(0,360,5.),labels=[0,0,0,1],color='black',linewidth=0) # draw meridians

            #(_,_,_,_,_,_,xm,ym) = eval_field(windfield,lonmin,latmin,t)
            #xpt,ypt = m(xm, ym)
            #m.scatter(xpt,ypt, c='k', s= 20, edgecolor='none')
            #m.plot(xpt,ypt,'kx')
            

            fname = '%05d.png' % count
            count += 1
            
            savefig('%s/%05i.png' % (pngdir, timeindex), bbox_inches='tight', dpi=300) 
            # To save memory, close the figure.
            close(fig)  
            
        # Make movie
        cmd = 'avconv -framerate %i -i "%s/%%05d.png" -y -vcodec qtrle %s' % (framerate,pngdir,moviefn)
        return_code = subprocess.call(cmd, shell=True)
        
        return moviefn
    
    
    
    
    
    def plot_div_vort(self, save=True, timestep=None, delta_deg=1.0, eval_loc=[None, None] ):
        '''
        Make two figures, one each for plots of the divergence and vorticity over the night.
        These parameters are calculated by taking finite differences from a point
        centered on the middle of the field of view of the network.
        OPTIONAL INPUTS:
            save:      bool. If True, save the resulting figures to self.save_dir
            timestep:  sec. Cadence with which to plot wind field. If None, use time step from the fit. 
                      (default None)
            delta_deg: deg. Size of step to calculate finite differences. Make sure that this
                       is larger than the grid size (default 1.0)
            eval_loc: list of length 2: [eval_lat, eval_lon]. Where the divergence and 
                      vorticity will be evaluated lat. If [None, None], then the center
                      of the field of view will be used. (default [None, None])
        RETURNS:
            (figs, fns)
                figs: list containing both figures
                fns:  list containing both saved filenames. [None, None] if save=False.
        '''
        # A note on error/uncertainty analysis: It would be nice if we could put an 
        # error bar on the divergences and vorticities that we estimate with this method.
        # The calculation of uncertainties below assumes that statistical errors 
        # in adjacent pixels are uncorrelated. This really couldn't be farther from the truth. 
        # Since adjacent pixels are determined from pretty much the same measurement(s), the 
        # errors are almost fully correlated. This means that the true statistical error bar on
        # gradient terms is virtually zero. It is a highly precise but potentially
        # inaccurate measurement. Because this is so confusing, we just won't show
        # error bars on this plot at all.
        if timestep is None:
            timestep = self.timestep_fit
        all_t = self.all_t
        distthresh = self.distthresh
        lons = self.lons
        lats = self.lats
        dlat = delta_deg # step for numerical approximation to gradient
        dlon = delta_deg

        #latseval, lonseval = meshgrid(linspace(mean(lats)-boxsize/2, mean(lats)+boxsize/2, Nbox),\
        #                              linspace(mean(lons)-boxsize/2, mean(lons)+boxsize/2, Nbox))
        #latseval = reshape(latseval, size(latseval))
        #lonseval = reshape(lonseval, size(lonseval))
        
        # These are lists because previously I had calculated the
        # average divergence/vorticity over some region. If you
        # want to calculate it at a single point, then just
        # make a list with one element.
        latseval = [mean(lats)]
        lonseval = [mean(lons)]
        # If the user specified a location, use that instead:
        if eval_loc[0] is not None:
            latseval = [eval_loc[0]] 
        if eval_loc[1] is not None:
            lonseval = [eval_loc[1]]
        

        # Define start/stop time
        t0, t1 = self.get_start_stop_time(timestep)
        timestepvec = range(0,int((t1-t0).total_seconds()+1), timestep)
        start_t = t0
        start_t = start_t.astimezone(pytz.utc).replace(tzinfo=None)
        divvec = nan*zeros(len(timestepvec))
        vortvec = nan*zeros(len(timestepvec))
        divvec_err = nan*zeros(len(timestepvec))
        vortvec_err = nan*zeros(len(timestepvec))
        tvec = []
                    
        for timeindex in range(len(timestepvec)):
            
            t = start_t + timedelta(seconds = timestepvec[timeindex])
            divs = zeros(len(latseval))
            div_errs = zeros(len(latseval))
            vorts = zeros(len(latseval))
            vort_errs = zeros(len(latseval))
            for i in range(len(latseval)):
                latx = latseval[i]
                lonx = lonseval[i]

                (u0,v0,w0,u0e,v0e,w0e,_,_,xm,ym) = self.eval_field(lonx, latx, t)
                (ux,vx,wx,uxe,vxe,wxe,_,_,xm,ym) = self.eval_field(lonx+dlon, latx, t)
                (uy,vy,wy,uye,vye,wye,_,_,xm,ym) = self.eval_field(lonx, latx+dlat, t)

                # Calculate gradients in units of s^-1
                dudx = (ux-u0)/(dlon*111e3*cos(latx*pi/180.))
                dudx_err = sqrt(uxe**2 + u0e**2)/(dlon*111e3*cos(latx*pi/180.))
                dudy = (uy-u0)/(dlat*111e3)
                dudy_err = sqrt(uye**2 + u0e**2)/(dlat*111e3)
                dvdx = (vx-v0)/(dlon*111e3*cos(latx*pi/180.))
                dvdx_err = sqrt(vxe**2 + v0e**2)/(dlon*111e3*cos(latx*pi/180.))
                dvdy = (vy-v0)/(dlat*111e3)
                dvdy_err = sqrt(vye**2 + v0e**2)/(dlat*111e3)
                div = dudx + dvdy
                div_err = sqrt(dudx_err**2 + dvdy_err**2)
                vort = dvdx - dudy
                vort_err = sqrt(dvdx_err**2 + dudy_err**2)
                
                divs[i] = div
                div_errs[i] = div_err
                vorts[i] = vort
                vort_errs[i] = vort_err
            
            divvec[timeindex] = nanmean(div)
            vortvec[timeindex] = nanmean(vort)
            # When calculating errors, don't divide by sqrt(N), because the error is probably correlated.
            # The proper way to do this is to propagate the entire covariance matrix through
            # all of these calculations. However, we don't calculate the entire covariance matrix, only the
            # diagonal. This discussion applies to the above uncertainty calculations as well.
            divvec_err[timeindex] = nanmean(div_errs) 
            vortvec_err[timeindex] = nanmean(vort_errs)            
            tvec.append(t)

        tvec = array(tvec)
        
        # Make plots
        yvals = [divvec, vortvec]
        yerrs = [divvec_err, vortvec_err]
        fnstubs = ['_divergence.png', '_vorticity.png']
        ylabels = ['Divergence [1/s]', 'Vorticity [1/s]']
        figs = []
        savefns = []
        for i in range(len(yvals)):
            yval = array(yvals[i])
            yerr = array(yerrs[i])
            fig = figure(figsize=(8,3), dpi=150)
            #fill_between(tvec, yval-yerr, yval+yerr, linewidth=0, edgecolor='k',\
            #             facecolor=[0.5, 0.5, 0.5], alpha = 1.0)
            plot(tvec,yval,'k')
            plot([tvec[0]-timedelta(hours=0.1), tvec[-1]+timedelta(hours=0.1)],[0,0],'k--')
            gca().xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H'))
            xlabel('UT [hours]')
            ylabel(ylabels[i])
            title('%s -- %s UT' % (tvec[0], tvec[-1]))
            #mx = 1.2*(max(abs(yval)+yerr))
            mx = 1.2*max(abs(yval))
            ylim([min(-0.0005,-mx), max(0.0005,mx)])
            tight_layout()
            savefn = None
            if save:
                savefn = '%s%s%s' % (self.dir, self.get_stub(), fnstubs[i])
                savefig(savefn)
            figs.append(fig)
            savefns.append(savefn)
    
        self.t_divvort = tvec
        self.div = divvec
        self.vort = vortvec
        return figs, savefns
          
# End class definition
