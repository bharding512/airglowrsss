#import pdb
from datetime import datetime,timedelta
import numpy as np
import fpiinfo
import pytz
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as md
from collections import defaultdict as _defaultdict
mpl.rcParams['savefig.bbox']='tight'

def cosd(arg):
    return np.cos(arg*np.pi/180.)
def sind(arg):
    return np.sin(arg*np.pi/180.)

def cloudthreshold():
    return 2.0

def moonuplimit():
    return 0.75

def errorbarlimit():
    return 200.

def synctiming():
    return 16*60.

def dn2utc(dn):
    import pytz
    return dn.astimezone(pytz.timezone('UTC'))
            
def dn2lt(dn):
    import pytz
    return dn.astimezone(pytz.timezone('UTC')) + dn.utcoffset()

def azel2lla(el,az,alt,lla0,horizon=10.):
    '''
    Summary:
        Script to calcuate assumed lla point of airglow origin.
    
    Inputs:
        el - elevation angle in degrees
        az - azimuth angle in degrees
        lla0 - latitude longitude altitude of origin [*,*,m]
        alt - assumed emission peak altitude in km
        horizion - cut off elevation angle

    Outputs:
        lat - latitude
        lon - longitude
        
    History:
        3/30/14 - Written by Daniel J. Fisher (dfisher2@illionis.edu)
    '''
    # Constants
    Ea = 6378137.     # semi-major axis of the earth [m]
    Eb = 6356752.3145    # semi-minor axis of the earth [m]  
      
    # Set variables
    zm = alt*1E3   # transform into m
    Re = Ea**2/np.sqrt(Ea**2*cosd(lla0[0])**2+Eb**2*sind(lla0[0])**2)
    
    # Calculate differencial angles
    B = np.arcsin((Re+lla0[2])*sind(el+90.)/(Re+zm))*180./np.pi
    A = 180.-(90+el+B)
    
    # Transfrom to lat,lon
    lat = np.arcsin(sind(lla0[0])*cosd(A) + cosd(lla0[0])*sind(A)*cosd(az))*180./np.pi
    lon = lla0[1] + np.arctan2(sind(az)*sind(A)*cosd(lla0[0]),cosd(A)-sind(lla0[0])*sind(lat))*180./np.pi

    # Nullify if below horizon
    bad = np.where(el < horizon)
    try:
        lat[bad[0],bad[1]] = np.nan
        lon[bad[0],bad[1]] = np.nan
    except:
        try:
            lat[bad[0]] = np.nan
            lon[bad[0]] = np.nan
        except:
            if el < horizon:
                lat = np.nan
                lon = np.nan

    return(lat,lon,alt)

def GetLocation(SITE,DIRECTION,ALT=250.):
    '''
    Summary:
        Script to calcuate assumed lla point of airglow origin (250km).
    
    Inputs:
        SITE - site of origin e.g.='uao'
        DIRECTION - direction keyword e.g.='CV_EKU_UAO_1'

    Outputs:
        latlonalt - (latitude, longitude, altitude)
        
    History:
        5/20/14 - Written by Daniel J. Fisher (dfisher2@illionis.edu)
    '''
    try:
        az = fpiinfo.get_site_info(SITE)['Directions'][DIRECTION]['az']
        el = 90-fpiinfo.get_site_info(SITE)['Directions'][DIRECTION]['ze']
    except:
        az = fpiinfo.get_site_info(SITE)['Directions'][DIRECTION[4:]]['az']
        el = 90-fpiinfo.get_site_info(SITE)['Directions'][DIRECTION[4:]]['ze']
    lla = fpiinfo.get_site_info(SITE)['Location']
    latlonalt = azel2lla(el,az,ALT,lla)
    
    return latlonalt

def GetLevel1(instr_name, dn):
    stub = '/rdata/airglow/fpi/results'
    site_name = fpiinfo.get_site_of(instr_name, dn)
    f = '%s/%s_%s_%4i%02d%02d.npz' \
            % (stub, instr_name, site_name.lower(), dn.year, dn.month, dn.day)
    return Level1(f, site_name, dn, instr_name)

class Level1:
    def __init__(self, f, site, dn, instr_name):
        import FPI
        import numpy as np
        import ephem
        
        self.f = f
        self.site = site
        self.dn = dn
        self.lla = fpiinfo.get_site_info(site)['Location']
        self.instr = instr_name

        self.error = False
        self.log = ""
        self.log += "%-24s" % "[Level 1, init]" \
                + "instance created on %s.\n" % str(datetime.now().strftime('%m/%d/%Y %H:%M:%S %p'))
        self.log += "%-24s" % "[Level 1, init]" \
                + "trying to open %s.\n" % self.f.split('/')[-1]

        # Check if moon is near full
        moon = ephem.Moon(dn)
        self.moonup = moon.moon_phase > moonuplimit()

        # make sure we have data to work with:
        try:
            npzfile = np.load(self.f,allow_pickle=True)
            self.r = npzfile['FPI_Results'].ravel()[0]
            del npzfile.f
            npzfile.close() # Brian's fix.
            # make sure we have zenith data to work with (djf 2/25/14):
            if 'Zenith' not in self.r['direction']:
                self.error = True
                self.log += "%-24s" % "[Level 1, init]" \
                        + "  => no zenith %s\n" %f.split('/')[-1]
                return
        except:
            self.error = True
            self.log += "%-24s" % "[Level 1, init]" \
                    + "  => can't load %s\n" %f.split('/')[-1]
            return

        
        # initialize variables:
        self.los_wind = {}
        self.los_sigma = {}
        self.t = {}
        self.iw = {}
        self.iwe = {}
        self.directions = list(set(self.r['direction']))
        self.allt = self.r['sky_times']
        self.reference = self.r['reference']
        self.T = {}
        self.Te = {}
        self.zenith = {}
        self.azimuth = {}
        self.alliw = np.array([])
        self.alliwe = np.array([])
        self.w = np.array([])
        self.we = np.array([])
        self.i = {}
        self.ie = {}
        self.b = {}
        self.be = {}
        self.flag_wind = {}
        self.flag_T = {}

        try:
            # Check for latest npz format
            self.allc = self.r['wind_quality_flag']
            self.ind = {}
        except:
            self.error = True
            self.log += "%-24s" % "[Level 1, init]" \
                    + "  => file outdated %s\n" %f.split('/')[-1]
            return

        # initialize dictionaries:
        for direction in self.directions:
            self.los_wind [direction] = []
            self.los_sigma[direction] = []
            self.t        [direction] = []
            self.iw       [direction] = []
            self.iwe      [direction] = []
            self.ind      [direction] = []
            self.T        [direction] = []
            self.Te       [direction] = []
            self.zenith   [direction] = []
            self.azimuth  [direction] = []
            self.i        [direction] = []
            self.ie       [direction] = []
            self.b        [direction] = []
            self.be       [direction] = []
            self.flag_wind[direction] = []
            self.flag_T   [direction] = []

        # Verify new file format
        try:
            test = self.r['sigma_fit_LOSwind']
        except:
            self.error = True
            self.log += "%-24s" % "[Level 1, init]" \
                    + "  => old file format; can't load %s\n" %f.split('/')[-1]
            return

        # fill in dictionaries with directions as keys:
        for (kk, (t, direction, los, sigma, T, Te, zenith, azimuth, i, ie, it, b, be, fw, fT)) \
                in enumerate(zip( \
                self.r['sky_times'], self.r['direction'], \
                self.r['LOSwind'], self.r['sigma_LOSwind'], \
                self.r['T'], self.r['sigma_T'], \
                self.r['ze'], self.r['az'], \
                self.r['skyI'], self.r['sigma_skyI'], self.r['sky_intT'],\
                self.r['ccdB'], self.r['sigma_ccdB'], \
                self.r['wind_quality_flag'], self.r['temp_quality_flag'] \
                )):
            self.los_wind[direction].append( los )
            self.los_sigma[direction].append( sigma )
            self.t[direction].append( t )
            self.ind[direction].append(kk)
            self.T[direction].append( T )
            self.Te[direction].append( Te )
            self.zenith[direction].append( zenith )
            self.azimuth[direction].append( azimuth )
            self.i[direction].append( i )
            self.ie[direction].append( ie )
            self.b[direction].append( b )
            self.be[direction].append( be )
            self.flag_wind[direction].append( fw )
            self.flag_T[direction].append( fT )

        # finally, make into arrays
        for direction in self.directions:
            self.los_wind [direction] = np.array(self.los_wind [direction])
            self.los_sigma[direction] = np.array(self.los_sigma[direction])
            self.t        [direction] = np.array(self.t        [direction])
            self.T        [direction] = np.array(self.T        [direction])
            self.Te       [direction] = np.array(self.Te       [direction])
            self.zenith   [direction] = np.array(self.zenith   [direction])
            self.azimuth  [direction] = np.array(self.azimuth  [direction])
            self.i        [direction] = np.array(self.i        [direction])
            self.ie       [direction] = np.array(self.ie       [direction])
            self.b        [direction] = np.array(self.b        [direction])
            self.be       [direction] = np.array(self.be       [direction])
            self.flag_wind[direction] = np.array(self.flag_wind[direction])
            self.flag_T   [direction] = np.array(self.flag_T   [direction])
            

        # call some functions:
        self._interpolate_w()
        self._cleanup()
        
        # clean L0 to fix memory issue? 
        #del self.r # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
        #self.r.close()

        self.log += "%-24s" % "" \
                + "  => after cleaning, has %03d data\n" % len(self.allt)

            
    def _cleanup(self):
        # remove bad indices corresponding to sigma > 50 or los > 1000.
        for direction in self.directions:
            bad_ind = np.where(\
                    self.los_sigma[direction] > errorbarlimit           
                    #(self.los_fit[direction] > errorbarlimit) \
                    #| (self.los_cal[direction] > errorbarlimit) \
                    #| (np.abs(self.los_wind[direction]) > 1e3) \
                    )  

            # log bad indices:
            if len(bad_ind[0]) > 0:
                self.log += "%-24s" % "[Level 1, cleanup]" + "deleted %03d 'bad' indices (los_fit > 50.) for %s\n" \
                % (len(bad_ind[0]), direction) 

            self.los_sigma[direction] = np.delete( self.los_sigma[direction], bad_ind)
            self.los_wind [direction] = np.delete( self.los_wind [direction], bad_ind)
            self.t        [direction] = np.delete( self.t        [direction], bad_ind)
            self.ind      [direction] = np.delete( self.ind      [direction], bad_ind)
            self.iw       [direction] = np.delete( self.iw       [direction], bad_ind)
            self.iwe      [direction] = np.delete( self.iwe      [direction], bad_ind)
            #self.iwef     [direction] = np.delete( self.iwef     [direction], bad_ind)
            #self.iwec     [direction] = np.delete( self.iwec     [direction], bad_ind)
            self.zenith   [direction] = np.delete( self.zenith   [direction], bad_ind)
            self.azimuth  [direction] = np.delete( self.azimuth  [direction], bad_ind)
            self.i        [direction] = np.delete( self.i        [direction], bad_ind)
            self.ie       [direction] = np.delete( self.ie       [direction], bad_ind)
            self.b        [direction] = np.delete( self.b        [direction], bad_ind)
            self.be       [direction] = np.delete( self.be       [direction], bad_ind)
            self.flag_wind[direction] = np.delete( self.flag_wind[direction], bad_ind)
            self.flag_T   [direction] = np.delete( self.flag_T   [direction], bad_ind)
            # temps?
        return

    def _interpolate_w(self):
        import time
        from scipy import interpolate
        from numpy import array as arr

        dointerp = True
        all_times = [time.mktime(dn.timetuple()) for dn in  self.allt] 
        
        # good indices for interpolating are not cloudy and have small los_fita
        # and reasonable los velocity
        if self.los_wind.has_key('Zenith'): 
            good_ind = np.where(\
                    (arr(self.flag_wind['Zenith']) <= cloudthreshold) \
                    * (self.los_sigma['Zenith'] <= errorbarlimit)\
                    * (np.abs(self.los_wind['Zenith']) < 1e3)\
                    )
        else:
            # if we don't have zenith measurements we can't interpolate vertical winds:
            # TODO set to zero instead for safety??? FLAG
            self.log+= "No Zenith for Interpolation \n"
            self.alliw  = np.ones( len(all_times)) * float('nan')
            self.alliwe = np.ones( len(all_times)) * float('nan')
            dointerp = False

        if dointerp:
            self.w  = arr(self.los_wind['Zenith'])
            self.we = arr(self.los_sigma['Zenith'])

            times =  arr([time.mktime(dn.timetuple()) for dn in  self.t['Zenith']])[good_ind]
            times = times[times.argsort()]

            # if we only have less than 4 data points, we can't do a cubic spline,
            # so just return with nothing
            if len(arr(times)) < 4:
                self.log+=  "%-24s" % "[Level 1, interpolate]" + "less than 4 data points, can't do cubic spline \n"
                iw  = np.ones( len(all_times)) * float('nan')
                iwe = np.ones( len(all_times)) * float('nan')
            else:
                # interpolate vertical winds for all times:
                # remember that for interpolation, we need to key in "good indices",
                # which are not cloudy days
                # Try spline, then linear, then zeroth-order interpolation
                #print self.w, self.w[good_ind]
                try: # spline interpolation
                    sfit = interpolate.UnivariateSpline(times, self.w[good_ind], w=1/self.we[good_ind]) #, s=s
                    iw = sfit(all_times)
                except: # Try linear.
                    print 'BAD NEWS BEARS'
                    tck = interpolate.splrep(arr(times),self.w[good_ind])
                    iw = interpolate.splev(all_times, tck)

                # interpolate vertical errors for all times:
                f = interpolate.interp1d( times, self.we[good_ind], \
                    bounds_error=False, fill_value=0.0)
                iwe = f(all_times)
                    
                # Fix all Cloudy times with 0m/s vertical wind since we cannot trust cloudy velocities
                # Find bad times (cloudy points only)
                cld_ind = np.where((np.array(self.allc) >= cloudthreshold) | (np.abs(iwe) > errorbarlimit))
                    # * (self.los_fit['Zenith'] < 50.)* (np.abs(self.los_wind['Zenith']) < 1e3))
                iw[cld_ind] = 0.0
                #iwe[cld_ind] = 0.0
                
            # apply to object:
            self.alliw  = iw
            self.alliwe = iwe

        # fill in dictionaries directions as keys:
        for (t, direction, iw, iwe) \
                in zip( \
                self.r['sky_times'], self.r['direction'], \
                self.alliw, self.alliwe):
            self.iw[direction] .append( iw )
            self.iwe[direction].append( iwe )

        # finally, make into numpy arrays...
        self.w  = np.array( self.w )
        self.we = np.array( self.we )
        for direction in self.directions:
            self.iw[direction] = np.array( self.iw[direction] )
            self.iwe[direction] = np.array( self.iwe[direction] )

    def __add__(self, other):
        '''
        Add two Level 1 objects together, commonly used because
        we can't rely on npz data containing prescribed date information,
        which is especially important for getting two common volume pairs
        together.

        '''
        import copy
        import sys
        
        # make sure we're adding Level 1 objects that have the same site
        # doesn't make sense to do this otherwise:
        if self.site != other.site:
            print "error: cannot combine %s with %s" % (out.site, other.site)
            sys.exit(-1)

        # if one or the other object has an error,
        # just return the 'good' one
        if self.error:
            other.log += self.log
            return other
        if other.error:
            self.log += other.log
            return self

        out  = copy.deepcopy(self)
        out.log += other.log

        '''
        Here are the items we need to sew together:
         1. 
         2. ind
         3. iw
         4. iwe
         5. los_wind
         6. los_sigma
         9. t
        10. T
        11. Te
        12. zenith
        13. azimuth
        '''

        # TODO: is there a cleaner/better way to do this?
        for direction in out.directions:
            if other.ind.has_key(direction):
                out.ind      [direction] = np.append( out.ind[direction], other.ind[direction] )
            if other.iw.has_key(direction):
                out.iw       [direction] = np.append( out.iw[direction], other.iw[direction] )
            if other.iwe.has_key(direction):
                out.iwe      [direction] = np.append( out.iwe[direction], other.iwe[direction] )
            if other.los_wind.has_key(direction):
                out.los_wind [direction] = np.append( out.los_wind[direction], other.los_wind[direction] )
            if other.los_sigma.has_key(direction):
                out.los_sigma[direction] = np.append( out.los_sigma[direction], other.los_sigma[direction] )
            if other.t.has_key(direction):
                out.t        [direction] = np.append( out.t[direction], other.t[direction] )
            if other.T.has_key(direction):
                out.T        [direction] = np.append( out.T[direction], other.T[direction] )
            if other.Te.has_key(direction):
                out.Te       [direction] = np.append( out.Te[direction], other.Te[direction] )
            if other.i.has_key(direction):
                out.i        [direction] = np.append( out.i[direction], other.i[direction] )
            if other.ie.has_key(direction):
                out.ie       [direction] = np.append( out.ie[direction], other.ie[direction] )
            if other.b.has_key(direction):
                out.b        [direction] = np.append( out.b[direction], other.b[direction] )
            if other.Te.has_key(direction):
                out.be       [direction] = np.append( out.be[direction], other.be[direction] )
            if other.zenith.has_key(direction):
                out.zenith   [direction] = np.append( out.zenith[direction], other.zenith[direction] )
            if other.azimuth.has_key(direction):
                out.azimuth  [direction] = np.append( out.azimuth[direction], other.azimuth[direction] )
            if other.flag_wind.has_key(direction):
                out.flag_wind[direction] = np.append( out.flag_wind[direction], other.flag_wind[direction] )
            if other.flag_T.has_key(direction):
                out.flag_T   [direction] = np.append( out.flag_T[direction], other.flag_T[direction] )

        # if by chance there is a direction in 'other' but not 'self':
        for direction in other.directions:
            if direction not in self.directions:
                out.ind       [direction] = other.ind      [direction]
                out.iw        [direction] = other.iw       [direction]
                out.iwe       [direction] = other.iwe      [direction]
                out.los_wind  [direction] = other.los_wind [direction]
                out.los_sigma [direction] = other.los_sigma[direction]
                out.t         [direction] = other.t        [direction]
                out.T         [direction] = other.T        [direction]
                out.Te        [direction] = other.Te       [direction]
                out.i         [direction] = other.i        [direction]
                out.ie        [direction] = other.ie       [direction]
                out.b         [direction] = other.b        [direction]
                out.be        [direction] = other.be       [direction]
                out.zenith    [direction] = other.zenith   [direction]
                out.azimuth   [direction] = other.azimuth  [direction]
                out.flag_wind [direction] = other.flag_wind[direction]
                out.flag_T    [direction] = other.flag_T   [direction]

                out.directions.append(direction)

        # these variables don't have a dictionary of directions, 
        # so it's simple to sew together:
        out.alliw   = np.append( out.alliw,  other.alliw  )
        out.alliwe  = np.append( out.alliwe, other.alliwe )
        out.allt    = np.append( out.allt,   other.allt   )
        out.w       = np.append( out.w,      other.w      )
        out.we      = np.append( out.we,     other.we     )

        out.f += "\n" + other.f
        out.error = out.error and other.error

        # for safety, 
        # it doesn't make much sense to add the original
        # data product
        out.r = None

        return out


    def plot_diagnostics(self, ax=None):
        '''
        TODO: documentation
        '''
        # Temporary note 22 Jul 2013:
        # Brian removed wind plot, in light of FPIDisplay.PlotDay().
        # What is the ultimate plan for plotting single-station data?
        

        if ax is None:
            fig = figure(1); clf()
            ax = fig.add_subplot(111)
        
        fontP = mpl.font_manager.FontProperties() 
        fontP.set_size(6)

        w = -self.los_wind['Zenith'] # los is towards observer
        werr = self.los_sigma['Zenith']

        ax.errorbar(self.allt, self.alliw, yerr=self.alliwe,\
                color='b',marker='.', label='interpolated',\
                )
        ax.errorbar(self.t['Zenith'], w,\
                yerr=werr,\
                marker='.', color='red', label='zenith'\
                )
        # TODO: Why does alliw come up nan sometimes? (e.g., (uao, 2013, 191))
        ax.set_ylim(  [np.nanmin([-70.,1.1*self.alliw.min(), 1.1*w.min()]),
                np.nanmax([70., 1.1*self.alliw.max(), 1.1*w.max()])  ])
        #ax.set_ylim([-70., 70.])
        ax.set_ylabel("Interpolated vertical wind")
        ax.get_xaxis().set_major_formatter(md.DateFormatter('%H'))
        ax.set_xlabel('Universal Time')
        td = timedelta(hours=0.5) # offset for plotting
        t0 = self.r['sky_times'][0]
        t1 = self.r['sky_times'][-1]
        ax.set_xlim([t0 - td, t1 + td])
        ax.legend(loc='best', prop=fontP, numpoints=1)
        ax.grid(True)

        #print self.log
        if ax is None: draw(); show()

        

# Note iw,iwe,iwef,iwec are never used in L2       
class Level2:
    def __init__(self, dn):
        import numpy as np
        import ephem
        
        self.dn = dn
        self.key = ""
        self.lla = np.array([])
        self.f = ""
        self.u  = np.array([]) 
        self.ue = np.array([])
        self.v  = np.array([])
        self.ve = np.array([])
        self.w  = np.array([])
        self.we = np.array([])
        self.wi = np.array([])
        self.wie= np.array([])
        
        self.i  = np.array([])
        self.ie = np.array([])
        self.b  = np.array([])
        self.be = np.array([])

        '''
        self.los_sigma1 = np.array([])
        self.los_sigma2 = np.array([])
        '''

        # I only use these 2 for the CardFinder,
        # specifically for Zenith... (?)
        self.allt = None

        self.flag_wind= np.array([])
        self.flag_T   = np.array([])

        '''
        # interpolated stuff:
        self.it  = np.array([])
        self.iw  = np.array([])
        self.iwe = np.array([])
        '''

        self.T  = np.array([])
        self.Te = np.array([])
        self.t1 = np.array([])
        self.t2 = np.array([])
        self.log = ""
        self.notes = ""
        self.length = 0
        self.error = False
        self.errorT = False
        self.parent = None

        # Check if moon is near full
        moon = ephem.Moon(dn)
        self.moonup = moon.moon_phase > moonuplimit()

    def __str__(self):
        # a print statement that looks nice:
        string = ""
        string += "%11s" % "dn = " + self.dn.strftime("%Y-%m-%d") + "\n"
        string += "%11s" % "f = " + self.f + "\n"
        string += "%11s" % "key = " + self.key + "\n"
        string += "%11s" % "log = " + self.log + "\n"
        string += "%11s" % "notes = " + self.notes + "\n"
        string += "%11s" % "length = " + "%3i" % self.length + "\n"

        return string

    def cut(self, dn1, dn2, inds=None):
        '''
        Cuts the Level 2 object between dn1 and dn2,
        or just selects the inds (if passed in).
        Dan uses the inds in his routines
        '''

        # dn1 and dn2 are LT
        t1 = np.array([dn.replace(tzinfo=None) for dn in self.t1])

        # if inds isn't passed in, then cut between dn1 and dn2:
        if inds is None:
            inds = np.where( (t1 > dn1) * (t1 < dn2) )

        self.length = len(inds)
        
        # cut them:
        '''
        if len(self.it) > 0:
            self.it = self.it[inds]

        if len(self.iw) > 0:
            self.iw = self.iw[inds]

        if len(self.iwe) > 0:
            self.iwe = self.iwe[inds]

        '''
        if len(self.T) > 0:
            self.T = self.T[inds]

        if len(self.t1) > 0:
            self.t1 = self.t1[inds]

        if len(self.t2) > 0:
            self.t2 = self.t2[inds]

        if len(self.Te) > 0:
            self.Te = self.Te[inds]

        if len(self.u) > 0:
            self.u = self.u[inds]

        if len(self.ue) > 0:
            self.ue = self.ue[inds]

        if len(self.v) > 0:
            self.v = self.v[inds]

        if len(self.ve) > 0:
            self.ve = self.ve[inds]

        if len(self.w) > 0:
            self.w = self.w[inds]

        if len(self.we) > 0:
            self.we = self.we[inds]

        if len(self.wi) > 0:
            self.wi = self.wi[inds]
            
        if len(self.wie) > 0:
            self.wie = self.wie[inds]

        if len(self.i) > 0:
            self.i = self.i[inds]

        if len(self.ie) > 0:
            self.ie = self.ie[inds]
            
        if len(self.b) > 0:
            self.b = self.b[inds]

        if len(self.be) > 0:
            self.be = self.be[inds]
        '''
        if len(self.los_sigma1) > 0:
            self.los_sigma1 = self.los_sigma1[inds]

        if len(self.los_sigma2) > 0:
            self.los_sigma2 = self.los_sigma2[inds]
        '''
        if len(self.flag_wind) > 0:
            self.flag_wind = self.flag_wind[inds]

        if len(self.flag_T) > 0:
            self.flag_T = self.flag_T[inds]

        self.log += "%-24s" % "[Level 2, cut]" \
                + "%03d data between %s and %s.\n" % (len(self.t1),dn1, dn2)
        return

    def plot(self, switch_onefig=False):
        mpl.rcParams['font.family'] = 'monospace'

        switch_plot_u = True
        switch_plot_v = True
        switch_plot_w = True
        switch_plot_iw = True

        if len(self.u) <3:
            switch_plot_u = False
        if len(self.v) <3:
            switch_plot_v = False
        if len(self.w) <3:
            switch_plot_w = False
        if len(self.wi) <3:
            switch_plot_iw = False

        if (not (switch_plot_u or switch_plot_v or switch_plot_w)):
            return None
        if self.error:
            return None

        if switch_onefig:
            fig = plt.figure(1); 
        else:
            fig = plt.figure();
        plt.clf()

        ax = fig.add_axes((.1,.2,.8,.7)) # left, bottom, width, height
        
        if switch_plot_u:
            plt.errorbar(self.t1, self.u, yerr=self.ue, \
                    color='b', marker='o', label='u')
        if switch_plot_v:
            plt.errorbar(self.t1, self.v, yerr=self.ve, \
                    color='g', marker='+', label='v')
        if switch_plot_w:
            plt.errorbar(self.t1, self.w, yerr=self.we, \
                    color='r', marker='*', label='w')

        if switch_plot_wi:
            plt.errorbar(self.it, self.wi, yerr=self.wie, \
                color='k', label='iw')

        dnp1 = self.dn + timedelta(days=1)
        plt.ylim([-200.,200.]) 
        ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        plt.legend()
        plt.grid()
        fig.text(.1,.05,self.notes)
        datestr = self.dn.strftime("%Y-%m-%d")
        fig.text(.1,.92,"%12s, %10s" % (self.key, datestr))
        plt.xlim( [datetime(self.dn.year, self.dn.month, self.dn.day, 20), datetime(dnp1.year, dnp1.month, dnp1.day, 12)] )

        plt.plot([datetime(self.dn.year, self.dn.month, self.dn.day, 20), datetime(dnp1.year, dnp1.month, dnp1.day, 12)],[0,0],'k--') 


        fig.text(.7,.030, self.log)
        plt.draw();
        plt.show()
        return 0


def convert2dict(obj):
    '''
    Summary:
        data = convert2dict(obj)
        Converts Level Object to Dictionary

    Inputs:
        obj = Level1 object / list of Level2 objects

    Outputs:
        data = Level1/2 dictonary (for easier useage)

    History:
        6/9/16 -- Written by DJF (dfisher2@illionis.edu)
    '''        
    if not obj:
        return([])
    # Setup dictionary
    d = _defaultdict(dict)
    # Convert L1 to Dictionary
    if isinstance(obj,Level1):
        for l1key in self.t.keys():
            d[l1key]['t']           = obj.t        [l1key]
            d[l1key]['los_wind']    = obj.los_wind [l1key]
            d[l1key]['los_sigma']   = obj.los_sigma[l1key]
            d[l1key]['T']           = obj.T        [l1key]
            d[l1key]['Te']          = obj.Te       [l1key]
            d[l1key]['i']           = obj.i        [l1key]
            d[l1key]['ie']          = obj.ie       [l1key]
            d[l1key]['flag_T']      = obj.flag_T   [l1key]
            d[l1key]['flag_wind']   = obj.flag_wind[l1key]
            d[l1key]['instr']       = obj.instr
            d[l1key]['site']        = obj.site
            d[l1key]['dn']          = obj.dn
            d[l1key]['type']        = 'Level1Dictionary'

    # Convert L2 to Dictionary
    elif type(obj) is list and isinstance(obj[0],Level2):
        for l2 in obj:
            if l2.error and l2.errorT:
                continue # Break if no data
            d[l2.key]['t']          = l2.t1
            d[l2.key]['u']          = l2.u
            d[l2.key]['ue']         = l2.ue
            d[l2.key]['v']          = l2.v
            d[l2.key]['ve']         = l2.ve
            d[l2.key]['w']          = l2.w
            d[l2.key]['we']         = l2.we
            d[l2.key]['T']          = l2.T
            d[l2.key]['Te']         = l2.Te
            d[l2.key]['i']          = l2.i
            d[l2.key]['ie']         = l2.ie
            d[l2.key]['flag_T']     = l2.flag_T
            d[l2.key]['flag_wind']  = l2.flag_wind
            #d[l2.key]['instr']      = l2.instr
            d[l2.key]['dn']         = l2.dn
            d[l2.key]['type']       = 'Level2Dictionary'

    return(d)


def CardFinder(dn, instr1, w_is_0=False):
    '''
    Summary:
        data = CardFinder(dn, instr1, w_is_0)
        Returns data for cardinal mode points,

    Inputs:
        dn = datetime day
        instr = instrument name
        w_is_0 = Flag to set w to 0, default=False

    Outputs:
        data = Level2 object

    History:
        3/26/13 -- Written by DJF (dfisher2@illionis.edu),
                        & TMD (duly2@illinois.edu)
    '''
    import os, sys
    import copy
    from datetime import timedelta

    #print "CardFinder: dn=",dn,"instr1=",instr1

    # Ouput variable is an instances of class Data()
    d = Level2(dn)
    site1 = fpiinfo.get_site_of(instr1, dn)
    d.key = site1.upper()
    d.instr = instr1
    d.log += "%-24s" % "[CardFinder]" \
            + "created on %s.\n" % str(datetime.now().strftime('%m/%d/%Y %H:%M:%S %p'))
    d.log += "%-24s" % "" \
            + "input: CardFinder(%s, '%s')\n" % (dn.strftime('datetime(%Y,%m,%d)'),instr1)

    l1 = GetLevel1(instr1, dn)
    d.moonup = l1.moonup
    d.log += l1.log
    if l1.error:
        # sometimes the file may not have anything in it:
        d.log += l1.log
        d.error = True
        d.errorT = True
        return [d]

    # log the parent Level 1 object:
    d.parent = [l1]
    
    # check to see if we have some vertical wind measurements
    # as we need them to back out u and v for Card measurements
    if len(l1.w)==0:
        d.log += "no vertical wind measurements found\n"
        d.error = True
        d.errorT = True
        return [d]

    d.allt = l1.allt

    # keep cardinal looks only and make it unique:
    dirs = ['Zenith', 'North', 'East', 'West', 'South']
    looks = list(set([l for x in dirs for l in l1.directions if x in l])) 

    # ------------------------------------------------
    # loop thru for different cardinal directions
    # ------------------------------------------------
    ds = []
    for look in looks:
        
        # copy the data instance with 
        # information we have so far:
        d_loop = copy.deepcopy(d)

        # reset output:
        u = np.array([]); ue = np.array([])
        v = np.array([]); ve = np.array([])
        w = np.array([]); we = np.array([])
        #wi= np.array([]); wie= np.array([])

        ind1 = l1.ind[look]
        
        # Record look times
        t1 = l1.t[look]

        # get temperatures
        d_loop.T  = l1.T [look]
        d_loop.Te = l1.Te[look]
        
        # get intensity and background
        d_loop.i  = l1.i [look]
        d_loop.ie = l1.ie[look]
        d_loop.b  = l1.b [look]
        d_loop.be = l1.be[look]

        # Save out interpolated wind
        wi  = l1.iw [look]
        wie = l1.iwe[look]

        # record parent Level 1 LOS errors:
        d_loop.los_sigma1 = l1.los_sigma[look]
        #d_loop.los_fit1 = l1.los_fit[look]
        #d_loop.los_cal1 = l1.los_cal[look]

        # fill in cloud information
        d_loop.flag_wind = l1.flag_wind[look]
        d_loop.flag_T    = l1.flag_T   [look]
    
        # Calculations using w
        if not(w_is_0):
            if 'Zenith' in look:
                # ------------------
                # Zenith measurement
                # ------------------
                w = l1.w
                we = l1.we
                d_loop.notes += 'Vertical wind is measurement\n'           

            elif 'East' in look:
                # ------------------
                # Eastern Zonal measurement
                # ------------------
                u = (l1.los_wind[look]-l1.iw[look]*cosd(l1.zenith[look]))/ \
                        sind(l1.zenith[look])
                ue = np.sqrt( l1.los_sigma[look]**2+l1.iwe[look]**2*cosd(l1.zenith[look])**2)/ \
                        sind(l1.zenith[look])
                d_loop.notes += 'Vertical wind is interpolated\n'

            elif 'West' in look:
                # ------------------
                # Western Zonal measurement
                # ------------------
                # Calculated Horizontal Winds
                u = (l1.los_wind[look]-l1.iw[look]*cosd(l1.zenith[look]))/ \
                        -sind(l1.zenith[look])
                ue = np.sqrt( l1.los_sigma[look]**2+l1.iwe[look]**2*cosd(l1.zenith[look])**2)/ \
                        sind(l1.zenith[look])
                d_loop.notes += 'Vertical wind is interpolated\n'
                
            elif 'North' in look:
                # ----------------------
                # North Meridional measurement
                # ----------------------
                v = (l1.los_wind[look]-l1.iw[look]*cosd(l1.zenith[look])) / \
                        sind(l1.zenith[look])
                ve = np.sqrt( l1.los_sigma[look]**2+l1.iwe[look]**2*cosd(l1.zenith[look])**2)/ \
                        sind(l1.zenith[look])
                d_loop.notes += 'Vertical wind is interpolated\n'

            elif 'South' in look:
                # ----------------------
                # South Meridional measurement
                # ----------------------
                v = (l1.los_wind[look]-l1.iw[look]*cosd(l1.zenith[look])) / \
                        -sind(l1.zenith[look])
                ve = np.sqrt( l1.los_sigma[look]**2+l1.iwe[look]**2*cosd(l1.zenith[look])**2)/ \
                        sind(l1.zenith[look])
                d_loop.notes += 'Vertical wind is interpolated\n'

            else:
                print "very bad ERROR"
                sys.exit(0)

        # Calculations for  w=0 case:
        else:
            if 'Zenith' in look:
                # ------------------
                # Zenith measurement
                # ------------------
                w = l1.w
                we = l1.we
                d_loop.notes += 'w is measurement\n'           
                wi = wi*0
                
            elif 'East' in look:
                # ------------------
                # Eastern Zonal measurement
                # ------------------
                u = l1.los_wind[look] / sind(l1.zenith[look])
                ue = l1.los_sigma[look] / sind(l1.zenith[look])
                d_loop.notes += 'Assumed w = 0\n'

            elif 'West' in look:
                # ------------------
                # Western Zonal measurement
                # ------------------
                u = l1.los_wind[look] / -sind(l1.zenith[look])
                ue = l1.los_sigma[look] / sind(l1.zenith[look])
                d_loop.notes += 'Assumed w = 0\n'
                
            elif 'North' in look:
                # ----------------------
                # North Meridional measurement
                # ----------------------
                v = l1.los_wind[look] / sind(l1.zenith[look])
                ve = l1.los_sigma[look] / sind(l1.zenith[look])
                d_loop.notes += 'Assumed w = 0\n'

            elif 'South' in look:
                # ----------------------
                # South Meridional measurement
                # ----------------------
                v = l1.los_wind[look] / -sind(l1.zenith[look])
                ve = l1.los_sigma[look] / sind(l1.zenith[look])
                d_loop.notes += 'Assumed w = 0\n'

            else:
                print "very bad ERROR"
                sys.exit(0)


        # ------------
        #  end look loop
        # ------------

        # Save information
        d_loop.key = "%s_%s" % (site1.upper(), look)
        d_loop.instr = instr1
        d_loop.lla = GetLocation(d_loop.parent[0].site,d_loop.key)

        d_loop.u = u
        d_loop.ue = ue
        d_loop.v = v
        d_loop.ve = ve
        d_loop.w = w
        d_loop.wi = wi 
        d_loop.wie = wie
        d_loop.we = we
        d_loop.t1 = t1
        d_loop.length = len(t1)

        ds.append(d_loop)
        d_loop.log += "%-24s" % "[CardFinder]" \
                + "%03d data for %s.\n" % (len(t1), d_loop.key)

    return ds

def CVFinder(dn, instr1, instr2, w_is_0=False):
    '''
    Summary:
        data = CVFinder(dn,instr1,instr2,w_is_0)
        Returns data for common value points,

    Inputs:
        dn = datetime day
        instr1 = instrument 1 name
        instr2 = instrument 2 name
        w_is_0 = Flag to set w to 0, default=False

    Outputs:
        data = Level2 object

    History:
    3/11/13 -- Written by DJF (dfisher2@illionis.edu),
                        & TMD (duly2@illinois.edu)
    '''
    #print "dn=",dn,'site1=',site1,'site2=',site2
    import os, sys 
    import copy
    import FPI
    from datetime import timedelta

    stub = '/mnt/FPIData/Results/'

    site1 = fpiinfo.get_site_of(instr1, dn).upper()
    site2 = fpiinfo.get_site_of(instr2, dn).upper()

    # Ouput variable is an instances of class Data()
    d = Level2(dn)

    d.log += "%-24s" % "[CVFinder]" \
            + "created on %s.\n" % str(datetime.now().strftime('%m/%d/%Y %H:%M:%S %p'))
    d.log += "%-24s" % "" \
            + "input: CVFinder(%s, '%s','%s')\n" % \
            (dn.strftime('datetime(%Y,%m,%d)'),instr1,instr2)
    d.log += "%-24s" % "" + "site1 = "+site1+", site2 = "+site2 + "\n"

    d.key = "%s_%s" % (site1, site2)

    # make sure we have different sites:
    if site1==site2:
        d.log += "%-24s" % "[CVFinder]" + "site1 and site2 are the same, let's not do CV" + "\n"
        d.error = True
        d.errorT = True
        return [d]

    # load in data:
    l1_1 = GetLevel1(instr1, dn)
    l1_2 = GetLevel1(instr2, dn)

    d.log += l1_1.log
    d.log += l1_2.log

    # if we have nothing, then just return:
    if l1_1.error and l1_2.error:
        d.error = True
        d.errorT = True
        return [d]


    # this block ensures that if we have 1 good and 1 bad site,
    # we'll use the temperature data from only the good one:
    # -------------------------------------------------------
    if l1_1.error:
        d.error = True
        #d.log += l1_1.log
        # just use l1_2 temps then:
        dirs = ['CV_', 'IN_']
        common_pair = list(set([l for x in dirs for l in l1_2.directions if x in l and site1 in l])) 
        ds = []
        for cv in common_pair:
            d_loop = copy.deepcopy(d)
            d_loop.t1 = l1_2.t[cv]
            d_loop.T  = l1_2.T[cv]
            d_loop.Te = l1_2.Te[cv]
            d_loop.i  = l1_2.i[cv]
            d_loop.ie = l1_2.ie[cv]
            d_loop.b  = l1_2.b[cv]
            d_loop.be = l1_2.be[cv]
            d_loop.u  = -999.*np.ones(len(l1_2.t[cv]))
            d_loop.ue = 999.*np.ones(len(l1_2.t[cv]))
            d_loop.v  = -999.*np.ones(len(l1_2.t[cv]))
            d_loop.ve = 999.*np.ones(len(l1_2.t[cv]))
            d_loop.w  = -999.*np.ones(len(l1_2.t[cv]))
            d_loop.we = 999.*np.ones(len(l1_2.t[cv]))
            d_loop.flag_wind = 3*np.ones(len(l1_2.t[cv]))
            d_loop.flag_T = l1_2.flag_T[cv]
            d_loop.key = cv
            ds.append(d_loop)
        return ds

    if l1_2.error:
        d.error = True
        #d.log += l1_2.log
        # just use l1_1 temps then:# just use l1_2 temps then:
        dirs = ['CV_', 'IN_']
        common_pair = list(set([l for x in dirs for l in l1_1.directions if x in l and site2 in l])) 
        ds = []
        for cv in common_pair:
	    d_loop = copy.deepcopy(d)
            d_loop.t1 = l1_1.t[cv]
            d_loop.T  = l1_1.T[cv]
            d_loop.Te = l1_1.Te[cv]
            d_loop.i  = l1_1.i[cv]
            d_loop.ie = l1_1.ie[cv]
            d_loop.b  = l1_1.b[cv]
            d_loop.be = l1_1.be[cv]
            d_loop.u  = -999.*np.ones(len(l1_1.t[cv]))
            d_loop.ue = 999.*np.ones(len(l1_1.t[cv]))
            d_loop.v  = -999.*np.ones(len(l1_1.t[cv]))
            d_loop.ve = 999.*np.ones(len(l1_1.t[cv]))
            d_loop.w  = -999.*np.ones(len(l1_1.t[cv]))
            d_loop.we = 999.*np.ones(len(l1_1.t[cv]))
            d_loop.flag_wind = 3*np.ones(len(l1_1.t[cv]))
            d_loop.flag_T = l1_1.flag_T[cv]
            d_loop.key = cv
            ds.append(d_loop)
        return ds
    # -------------------------------------------------------


    # log the parents as level 1 data for site 1 and 2:
    d.parent = [l1_1, l1_2]

    # check to see if we have some vertical wind measurements
    # our CV winds depend on these
    if ((len(l1_1.w)==0) or (len(l1_2.w)==0)):
        d.log += "%-24s" % "[CVFinder]" + "no vertical wind measurements found\n"
        d.error = True
        d.errorT = True
        return [d]
    
    # get common value locations between the 2 sites:
    common_pair = [val for val in l1_1.directions if val in l1_2.directions]
    # Keep CV and make unique pairs
    dirs = ['CV_', 'IN_']
    common_pair = list(set([l for x in dirs for l in common_pair if x in l]))
    #print "CV =", common_pair

    # check to make sure we have some common pairs. 
    # if not, exit
    if len(common_pair) == 0:
        d.log += "%-24s" % "[CVFinder]" + "no common pairs found, exiting with error flag \n"
        d.error = True
        d.errorT = True
        return [d]

    # ------------------------------------------------
    # loop thru for different common volume directions
    # ------------------------------------------------
    ds = []
    for cv in common_pair:
        #print cv

        # copy the data instance with 
        # information we have so far:
        d_loop = copy.deepcopy(d)
        
        # reset output:
        u = np.array([]); ue = np.array([]); uef = np.array([]); uec = np.array([])
        v = np.array([]); ve = np.array([]); vef = np.array([]); vec = np.array([])
        w = np.array([]); we = np.array([]); wef = np.array([]); wec = np.array([])
        iw= np.array([]); iwe= np.array([])

        # -----------------------------------------------------
        # this section of code gets rid of times that do 
        # not match their partner CV times
        t1_list = l1_1.t[cv]
        t2_list = l1_2.t[cv]

        if len(t1_list)==0 or len(t2_list)==0: continue
        # (tmd) I used to have this as 'break',
        #       but 'continue' skips to the next iterator

        good_ind1 = []
        good_ind2 = []

        for kk, t1 in enumerate(t1_list):
            min_list = []
            for t2 in t2_list:
                min_list.append( abs(t1-t2).seconds)
            jj = np.argmin(min_list)
            if abs(t1-t2_list[jj]).seconds < synctiming:
                good_ind1.append(kk)
                good_ind2.append(jj)

        t1         = t1_list           [good_ind1]
        los_sigma1 = l1_1.los_sigma[cv][good_ind1]
        iw1        = l1_1.iw[cv]       [good_ind1]
        iwe1       = l1_1.iwe[cv]      [good_ind1]
        T1         = l1_1.T[cv]        [good_ind1]
        Te1        = l1_1.Te[cv]       [good_ind1]
        ind1       = l1_1.ind[cv]      [good_ind1]
        ze1        = l1_1.zenith[cv]   [good_ind1]
        az1        = l1_1.azimuth[cv]  [good_ind1]
        los_wind1  = l1_1.los_wind[cv] [good_ind1]
        i1         = l1_1.i[cv]        [good_ind1]
        ie1        = l1_1.ie[cv]       [good_ind1]
        b1         = l1_1.b[cv]        [good_ind1]
        be1        = l1_1.be[cv]       [good_ind1]
        flag_wind1 = l1_1.flag_wind[cv][good_ind1]
        flag_T1    = l1_1.flag_T[cv]   [good_ind1]

        t2         = t2_list           [good_ind2]
        los_sigma2 = l1_2.los_sigma[cv][good_ind2]
        iw2        = l1_2.iw[cv]       [good_ind2]
        iwe2       = l1_2.iwe[cv]      [good_ind2]
        T2         = l1_2.T[cv]        [good_ind2]
        Te2        = l1_2.Te[cv]       [good_ind2]
        ind2       = l1_2.ind[cv]      [good_ind2]
        ze2        = l1_2.zenith[cv]   [good_ind2]
        az2        = l1_2.azimuth[cv]  [good_ind2]
        los_wind2  = l1_2.los_wind[cv] [good_ind2]
        i2         = l1_2.i[cv]        [good_ind2]
        ie2        = l1_2.ie[cv]       [good_ind2]
        b2         = l1_2.b[cv]        [good_ind2]
        be2        = l1_2.be[cv]       [good_ind2]
        flag_wind2 = l1_2.flag_wind[cv][good_ind2]
        flag_T2    = l1_2.flag_T[cv]   [good_ind2]
        # -----------------------------------------------------

        d_loop.log += "%-24s" % "[CVFinder]" +\
                "'%s' and '%s' have %03d indices that time-match\n\t\t\t(from a possible (%03d,%03d), respectively)\n" % (instr1,instr2, len(good_ind1), len(t1_list), len(t2_list))

        # Get Average Temperature/Error at pair pt
        T = np.average(np.vstack((T1, T2)), \
                axis=0, \
                weights=1./np.vstack((Te1,Te2)))
        Te = np.average(np.vstack((Te1, Te2)), \
                axis=0, \
                weights=1./np.vstack((Te1,Te2)))

        # Get Average Instesity/Background & Error at pair pt
        si = np.average(np.vstack((i1, i2)), \
                axis=0, \
                weights=1./np.vstack((ie1,ie2)))
        sie = np.average(np.vstack((ie1, ie2)), \
                axis=0, \
                weights=1./np.vstack((ie1,ie2)))
        sb = np.average(np.vstack((b1, b2)), \
                axis=0, \
                weights=1./np.vstack((be1,be2)))
        sbe = np.average(np.vstack((be1, be2)), \
                axis=0, \
                weights=1./np.vstack((be1,be2)))
                
        d_loop.flag_wind = np.maximum(flag_wind1,flag_wind2)
        d_loop.flag_T    = np.maximum(flag_T1,flag_T2)

        if 'IN' in cv[:2]:
            # ------------------
            # INline measurement
            # ------------------
    
            # FIX cos[ze] is from look direction, not at the point in the sky where they are equal...
            w = (los_wind1 + los_wind2) / \
                    (cosd(ze1) + cosd(ze2))
            we = np.sqrt( los_sigma1**2 + los_sigma2**2) / \
                    (cosd(ze1) + cosd(ze2))
            den=np.sqrt( (cosd(ze1) + cosd(ze2))**2)
            d_loop.notes += 'Vertical wind is inline measurement\n'
            
        elif 'CV' in cv[:2]:
            # ------------------
            # CV measurement
            # ------------------
            
            # Calculated Average vertical wind at cv points
            iw = np.average(np.vstack((iw1,iw2)), \
                    axis=0, \
                    weights=(1./np.vstack((iwe1,iwe2))**2),\
                    )
            iwe = np.sqrt((iwe1**-2 + iwe2**-2)/(iwe1**-4+ iwe2**-4))


            if not(w_is_0):
                # interpolated winds used.
                d_loop.notes += 'Vertical wind is interpolated\n'
                
                # FIX cos[ze] is from look direction, not at the point in the sky where they are equal...
                # vh = velocity horizontal
                vh1 = (los_wind1+iw*cosd(ze1)) / \
                        sind(ze1)
                vh2 = (los_wind2+iw*cosd(ze2)) / \
                        sind(ze2)
                
                vh1e = np.sqrt( los_sigma1**2+(iwe*cosd(ze1))**2) / \
                        sind(ze1)
                vh2e = np.sqrt( los_sigma2**2+(iwe*cosd(ze2))**2) / \
                        sind(ze2)
            
            # Case with w=0
            else:
                # Don't use iw
                iw = iw*0
                d_loop.notes += 'Assumed w = 0\n'

                # FIX cos[ze] is from look direction, not at the point in the sky where they are equal...
                # vh = velocity horizontal
                vh1 = los_wind1 / sind(ze1)
                vh2 = los_wind2 / sind(ze2)
                
                vh1e = los_sigma1 / sind(ze1)
                vh2e = los_sigma2 / sind(ze2)


            # Calculate winds
            #M = np.array([[sind(az1),cosd(az1)],[sind(az2),cosd(az2)]])
            u = [] ; ue = [] 
            v = [] ; ve = []
            for kk, (myt, a) in enumerate(zip(t1,vh1)):
                M = np.array([[sind(az1[kk]),cosd(az1[kk])],[sind(az2[kk]),cosd(az2[kk])]])
                M = np.matrix(M)# important to make into matrix

                # value calculation:
                i = min(range(len(t2)), key=lambda i: abs(t2[i]-myt))
                b = vh2[i]

                temp = np.linalg.solve(M,np.array([[a],[b]]))
                u.append(temp[0][0])
                v.append(temp[1][0])

                # error calculation:
                S = np.zeros((2,2))
                S[0,0] = vh1e[kk]**2
                S[1,1] = vh2e[i ]**2
                S = np.matrix(S) # important to make into matrix
                temp = np.array(M*S*M.T) # be sure that M and S are matricies, 
                                         # or else it won't multiply properly!
                ue.append(np.sqrt(temp[0][0]))
                ve.append(np.sqrt(temp[1][1]))

        # Save information
        d_loop.key = cv
        d_loop.lla = GetLocation(d_loop.parent[0].site,d_loop.key)

        d_loop.u = np.array(u)
        d_loop.ue = np.array(ue)
        d_loop.v = np.array(v)
        d_loop.ve = np.array(ve)
        d_loop.w = w
        d_loop.wi = iw
        d_loop.we = we
        d_loop.wec = wec
        d_loop.T = T
        d_loop.Te = Te
        d_loop.i = si
        d_loop.ie = sie
        d_loop.b = sb
        d_loop.be = sbe
        d_loop.t1 = t1
        d_loop.t2 = t2
        d_loop.length = len(t1)

        d_loop.log += "%-24s" % "[CVFinder]" \
                + "%03d data for %s.\n" % (len(t1), d_loop.key)

        ds.append(d_loop)

    return ds


def TempFinder(dn, instr1):
    '''
    Summary
    -------

    data = TempFinder(dn, instr1)
    
    Returns temp data for non-card/CV mode points,

    Inputs
    ------
        dn = datetime day

    Outputs
    -------


    History
    -------
    10/19/15 -- Written by DJF (dfisher2@illionis.edu)

    '''
    import os, sys
    import copy
    from datetime import timedelta

    # Ouput variable is an instances of class Data()
    d = Level2(dn)
    site1 = fpiinfo.get_site_of(instr1, dn)
    d.key = site1.upper()
    d.instr = instr1
    d.log += "%-24s" % "[TempFinder]" \
            + "created on %s.\n" % str(datetime.now().strftime('%m/%d/%Y %H:%M:%S %p'))
    d.log += "%-24s" % "" \
            + "input: TempFinder(%s, '%s')\n" % (dn.strftime('datetime(%Y,%m,%d)'),instr1)

    l1 = GetLevel1(instr1, dn)

    d.moonup = l1.moonup
    d.log += l1.log
    if l1.error:
        # sometimes the file may not have anything in it:
        d.log += l1.log
        d.error = True
        d.errorT = True
        return([])

    # log the parent Level 1 object:
    d.parent = [l1]
    d.allt = l1.allt

    # For Temp Only files (ADD MORE and make unique pairs
    dirs = ['MTM_Search', 'Windfield', 'Along_B', 'Aux']
    looks = list(set([l for x in dirs for l in l1.directions if x in l]))

    # ------------------------------------------------
    # loop thru for different temp directions
    # ------------------------------------------------
    ds = []
    for look in looks:
        
        # copy the data instance with 
        # information we have so far:
        d_loop = copy.deepcopy(d)

        
        # Record look times
        t1 = l1.t[look]

        # get temperatures
        d_loop.T  = l1.T [look]
        d_loop.Te = l1.Te[look]
        
        # get intensity and background
        d_loop.i  = l1.i [look]
        d_loop.ie = l1.ie[look]
        d_loop.b  = l1.b [look]
        d_loop.be = l1.be[look]

        # Save out interpolated wind
        wi  = l1.iw [look]
        wie = l1.iwe[look]

        # record parent Level 1 LOS wind info:
        d_loop.los_wind1 = l1.los_wind[look]
        d_loop.los_sigma1 = l1.los_sigma[look]

        # fill in cloud information
        d_loop.flag_wind = l1.flag_wind[look]
        d_loop.flag_T    = l1.flag_T   [look]
    
        d_loop.notes += 'TEMPS ONLY'

        # Save information
        d_loop.key = "%s_%s" % (site1.upper(), look)
        d_loop.instr = instr1
        d_loop.lla = GetLocation(d_loop.parent[0].site,d_loop.key)

        d_loop.wi = wi 
        d_loop.wie = wie
        d_loop.t1 = t1
        d_loop.length = len(t1)

        ds.append(d_loop)
        d_loop.log += "%-24s" % "[TempFinder]" \
                + "%03d data for %s.\n" % (len(t1), d_loop.key)

    return ds


'''
def PlotLatSLT(cvs, dn1, dn2, switch_interpolate_T=False ):
   
    #input: a list of level 2 objects
    
    from scipy.interpolate import interp2d
    from scipy.interpolate import LinearNDInterpolator
    figure(1, figsize=(18,10.5)); clf()
    #figure(1, figsize=(18/1.5,10.5/1.5)); clf() # laptop testing
    matplotlib.rcParams.update({'font.size': 24})
    from datetime import timedelta
    import pytz
    from FPIResults import BinMonthlyData
    #Sites = MasterDictionary()

    # OK FOR NOW!
    project = 'nation'
    

    ALPHA_VALUE = 0.3

    # quiver plot
    # -----------
    qscale = 2.2
    qunits = 'dots'
    qwidth = 3.
    dns = []
    plot_text = {}
    for mm, cv in enumerate(cvs):
        # CV measurement:
        if 'CV' == cv.key[:2] and cv.error is False:
            (lat,lon,alt) = GetLocation(cv.parent[0].site,cv.key)
            slt = [dn2utc(a)+timedelta(hours=lon/15.) for a in cv.t1]
            times = matplotlib.dates.date2num(slt) 
            for t, u, v, cloudy in zip(times, cv.u, cv.v, cv.cloudy):
                Q = quiver(t, lat, u, v, \
                        scale=qscale, \
                        units=qunits, \
                        width=qwidth, \
                        zorder=100, \
                        #alpha = ALPHA_VALUE if cloudy else None, \
                        visible = False if cloudy else True,\
                        )
            if len(slt) > 0:
                if plot_text.has_key(cv.key) is False:
                    plot_text[cv.key] = (dn1 + timedelta(hours=.1), lat, cv.key)
            dns += list(cv.t1)
            
        # Card measurment, use North/East and South/West combos
        elif cv.key[:3] in [a.upper() for a in fpiinfo.get_network_info(project).keys()] and cv.error is False:
            # card points
            # use north/east and south/west pairs to plot wind
            lat = SiteLocations[cv.key]['lat']
            lon = SiteLocations[cv.key]['lon']
            slt = [dn2utc(a)+timedelta(hours=lon/15.) for a in cv.t1]
            times = matplotlib.dates.date2num(slt) 
            site = cv.key.split("_")[0]
            first_loc = cv.key.split("_")[1]
            if first_loc in ['East','West',]:
                Matches = { 'East': 'North', \
                            'North': 'East', \
                            'South': 'West', \
                            'West': 'South', \
                            }
                matching_key = site + "_" + Matches[first_loc]
                print matching_key
                lat = SiteLocations[matching_key]['lat'] 
                lon = SiteLocations[matching_key]['lon'] 
                cv_match = None
                for a_cv in cvs:
                    if matching_key == a_cv.key:
                        cv_match = a_cv
                if cv_match is not None:
                    for u, t1, cloudy in zip(cv.u, cv.t1, cv.cloudy):
                        times = []
                        for t1match in cv_match.t1:
                            times.append(abs(t1match-t1).seconds)
                        if len(times) > 0:
                            kk = np.argmin( np.array(times))
                            t1_slt = dn2utc(t1)+timedelta(hours=lon/15.)
    
                            diff = abs(cv_match.t1[kk]-t1)
                            if diff.seconds < 30*60. and diff.days==0 \
                                    and abs(u) < 500. and abs(cv_match.v[kk]) < 500.:

                                print "matching times:",\
                                        t1.strftime("%H:%M"), cv_match.t1[kk].strftime("%H:%M")
                                print u, cv_match.v[kk]
                    
                                t1_match = cv_match.t1[kk]
                                time = matplotlib.dates.date2num(t1_slt)
                                Q = quiver(time, lat, u, cv_match.v[kk], \
                                        scale=qscale, \
                                        units=qunits, \
                                        width=qwidth, \
                                        color='r', \
                                        zorder=100, \
                                        #alpha = ALPHA_VALUE if cloudy else None, \
                                        visible = False if cloudy else True,\
                                        )
                            if plot_text.has_key(cv.key) is False:
                                plot_text[cv.key] = (dn1 + timedelta(hours=.1), lat, cv.key+"/"+Matches[first_loc])
            dns += list(cv.t1)


        
    # put monthly averages down, using Dan's routine:
    # -----------------------------------------------
    for site in ['PAR', 'ANN', 'UAO', 'EKU']:
        print "-----------"
        print "inputs to BinMonthlyData:"
        print "site, year, month = ", site, dn1.year, dn1.month
        o = BinMonthlyData(site, dn1.year, dn1.month)
        t1 = [ datetime(dn1.year, dn1.month, dn1.day, dn.hour, dn.minute, dn.second, tzinfo= pytz.utc) for dn in o.t1]
        t2 = [ datetime(dn2.year, dn2.month, dn2.day, dn.hour, dn.minute, dn.second, tzinfo= pytz.utc) for dn in o.t1]
        lat = SiteLocations[site+"_Zenith"]['lat']
        lon = SiteLocations[site+"_Zenith"]['lon']
        # copy and paste data for both dn1 and dn2, ensuring that when 
        # plotted, the information is visible:
        t1_slt = [dn2utc(a) + timedelta(hours=lon/15.) for a in t1] 
        t2_slt = [dn2utc(a) + timedelta(hours=lon/15.) for a in t2] 
        for ttt1, ttt2, u, v in zip(t1_slt, t2_slt, o.u, o.v):
            nan = float('nan')
            if (not np.isnan(u)) and (not np.isnan(v))\
                    and np.sqrt(u**2+v**2) < 400.:
                Q = quiver(matplotlib.dates.date2num(ttt1), lat, u, v,\
                        scale=qscale,\
                        units=qunits,\
                        width=qwidth,\
                        edgecolor='k',\
                        facecolor='w',\
                        linewidth=0.5,
                        #linestyle='dashed',\
                        )
                Q = quiver(matplotlib.dates.date2num(ttt2), lat, u, v,\
                        scale=qscale,\
                        units=qunits,\
                        width=qwidth,\
                        edgecolor='k',\
                        facecolor='w',\
                        linewidth=0.5,
                        #linestyle='dashed',\
                        )
    # -------- thanks Dan -------------------------------

    Qbogus = quiver(0,0,0,0,\
            scale=qscale,\
            units=qunits,\
            width=qwidth,\
            color='k',\
            zorder=100,\
            )
    qk = quiverkey(Qbogus, 0.92, 0.92, 100.,\
            r'100 $\mathrm{m/s}$',\
            fontproperties={'size':20},\
            )
    for key in plot_text.keys():
        text( plot_text[key][0], plot_text[key][1], plot_text[key][2],\
                color='r' if "/" in plot_text[key][2] else 'k',\
                fontsize=12,\
                )
    xlabel('Solar Local Time')
    ylabel('Latitude')
    title('%s Thermospheric Winds\n' % (project.upper(),) + \
            "Night of UT " +
            dn1.strftime("%Y-%m-%d") \
              )
    ax = gca()
    # taken from 
    #  http://stackoverflow.com/questions/13950053/quiver-or-barb-with-a-date-axis
    ax.get_xaxis().set_major_formatter(matplotlib.dates.DateFormatter('%H'))
    xlim([dn1, dn2])

    my_yticks = {\
            'nation': range(30,50,2),\
            'RENOIR': range(-40,0,1),\
            }
    my_ylims = {\
            'nation': [31., 45.],\
            'RENOIR': [-12., -2.],\
            }
    yticks(my_yticks[project])
    ylims = my_ylims[project]
    ylim(ylims)
    grid()
    # (end) quiver plot
    # -----------------
    

    ## T pcolor plot
    ## -----------------
    #time_grid = np.linspace( date2num(dn1), date2num(dn2), 101)
    #lat_grid  = np.linspace(ylims[0], ylims[1], 100)
    #x=[];y=[];z=[];
    #
    #flag_colorbar = False
    #for cv in cvs:
    #    if not cv.errorT and cv.key != 'IN_PAR_UAO_EKU':
    #        lat = SiteLocations[cv.key]['lat']
    #        lon = SiteLocations[cv.key]['lon']
    #        slts = [dn2utc(a)+timedelta(hours=lon/15.) for a in cv.t1]
    #        lat_index = np.argmin( np.abs( lat_grid - lat))
    #        for kk, (slt, cloudy) in enumerate(zip(slts, cv.cloudy)):
    #            grid = np.zeros( (len(lat_grid), len(time_grid)) ) * float('nan')
    #            time_index = np.argmin( np.abs( time_grid - date2num(slt) ) )
    #            grid[lat_index, time_index] = cv.T[kk]

    #            TIME_GRID, LAT_GRID = np.meshgrid(time_grid, lat_grid)
    #            pcolor(TIME_GRID, LAT_GRID, np.ma.array(grid, mask=np.isnan(grid)), \
    #                    vmin=600.,\
    #                    vmax=1200.,\
    #                    zorder=10, \
    #                    edgecolors='k',\
    #                    alpha = ALPHA_VALUE if cloudy else None, \
    #                    )
    #            flag_colorbar = True
    #            if cv.T[kk] < 5000.: # gets rid of really high T values (messes up interpolation)
    #                x.append( time_grid[time_index] )
    #                y.append( lat_grid[lat_index] )
    #                z.append( cv.T[kk] )
    #if switch_interpolate_T:
    #    #f = interp2d(x,y,z, bounds_error=False)
    #    #f = np.vectorize(f)
    #    f = LinearNDInterpolator((x,y),z) # this one works better
    #    zi = f(TIME_GRID, LAT_GRID)

    #    pcolormesh(TIME_GRID, LAT_GRID, np.ma.array(zi, mask=np.isnan(zi)), \
    #            vmin=600.,\
    #            vmax=1200.,\
    #            zorder=5, \
    #            )

    #plt.grid();
    #if flag_colorbar:
    #    c = colorbar(alpha=None);
    #    c.set_label(r'Temperature [$^\circ K$]')
    # (end) T pcolor plot
    # -----------------

    draw(); show()
    #show()
'''

def GetLevel2(project,dn,dnstart='noon',dnend='noon',w_is_0=False):
    '''
    Summary:
        Returns all Level 2 data for an entire project for a single night.

    Inputs:
        project = name of poject to collect data from, project = 'NATION'
        dn = datetime of data desired, dn = datetime(2013,2,23)
        dnstart = starting point of time range desired, defaults to noon of dn, dnend(2013,2,23,14,0,0)
        dnend = end point of time range desired, defaults to noon of dn +1 day, dnend(2013,2,24, 6,0,0)
        w_is_0 = Flag to set w to 0, default=False

    Outputs:
        cvs -- a list of Level 2 instances

    History:
        05/23/13 -- Written by DJF (dfisher2@illinois.edu)
        08/07/13 -- Updated to be instrument-based (Timothy Duly, duly2@illinois.edu)
        10/21/14 -- Updated to allow sites & projects (DJF dfisher2@illinois.edu)
    '''
    from datetime import timedelta
    from datetime import datetime

    # look up the sites for project:
    try:
        sites = fpiinfo.get_network_info(project).keys()
        psites = sites
        editlist = False
    except:
        sites = [project]
        psites = fpiinfo.get_network_info(fpiinfo.get_site_info(project)['Network']).keys()
        editlist = True

    #print editlist

    # with these sites, look up the corresponding sites.
    # these are the instruments used for a given day & project.
    instrs = []
    combos = []
    for site in sites:
        instrs += fpiinfo.get_instr_at(site,dn)
    for site in psites:
        combos += fpiinfo.get_instr_at(site,dn)

    # this loops thru instruments for CVFinder,
    # withouth repeating. (Thanks Dan)
    cvs = []
    pairs = ['Unknown']
    for instr in instrs:
        cvs += CardFinder(dn, instr, w_is_0)
        pairs.append(instr)
        for combo in instrs:
            if combo not in pairs:
                cvs += CVFinder(dn, instr, combo, w_is_0)
        cvs += TempFinder(dn,instr)
        
    # remove error parts in cvs:
    #cvs = []
    if editlist:
        cvs = [x for x in cvs if project in x.key.lower()]
    
    if dnstart == 'noon':
        dnstart = dn.replace(hour=12)
    if dnend == 'noon':
        dnend = (dn + timedelta(days=1)).replace(hour=12)
        
    # cut so that we only have data between dnstart and dnend
    for kk, cv in enumerate(cvs):
        cv.cut(dnstart, dnend)
    
    return cvs


def PlotLevel2(site,dn,ut=True,w_is_0=False):
    '''
    Summary
    -------
    Plots u,v,w,T from Level 2 data for a site for a single night.

    Inputs
    ------
    site = name o fsite to collect data from, e.g.: 'car'
    dn = datetime of data desired, e.g.: datetime(2013,2,23)
    UTC = plot time in UTC? [default = True]
    w_is_0 = Flag to set w to 0, default=False

    History
    -------
    11/02/15 -- Written by DJF (dfisher2@illionis.edu)
    '''

    project = fpiinfo.get_site_info(site,dn)['Network']
    slt_shift = 0.
    utc = pytz.UTC
    slt_shift = 24 - fpiinfo.get_site_info(site,dn)['Location'][1]%360/360.*24

    if ut:
        tunit = 'UTC'
        tlim = [dn+timedelta(hours=18+slt_shift),dn+timedelta(hours=30+slt_shift)]
        slt_shift = 0
    else:
        tunit = 'SLT'
        tlim = [dn+timedelta(hours=18),dn+timedelta(hours=30)]
    D = GetLevel2(project,dn,w_is_0) 
    
    # Set colors for different flags
    lalpha = 0.1
    ascale = 2.3
    ax={}
    f,(ax[0],ax[1],ax[2],ax[3]) = plt.subplots(4,sharex=True,sharey=False,figsize=(14,10))
    for d in D:
        for flag in [2,1,0]:
            indw = [x for x in range(len(d.t1)) if d.flag_wind[x] <= flag]
            indt = [x for x in range(len(d.t1)) if d.flag_T[x] <= flag]
            if flag == 2:
                line = '-'
                alpha = 0.2
            elif flag == 1:
                alpha = 0.45
                line = ''
            else:
                alpha = 1
                line = ''
                
            # Time Travel Math
            timew = [x.astimezone(utc).replace(tzinfo=None) for x in d.t1[indw]-timedelta(hours=slt_shift)]
            timet = [x.astimezone(utc).replace(tzinfo=None) for x in d.t1[indt]-timedelta(hours=slt_shift)]

            if not(timew) or not(timet):
                continue

            # Zonal
            if ('East' in d.key and site in d.key.lower()):
                ax[0].errorbar(timew,d.u[indw],yerr=d.ue[indw],fmt='g.'+line,alpha=alpha)
                if flag==0:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='g.'+line,alpha=alpha,label=d.key)
                else:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='g.'+line,alpha=alpha)
            # Zonal
            if ('West' in d.key and site in d.key.lower()):
                ax[0].errorbar(timew,d.u[indw],yerr=d.ue[indw],fmt='y.'+line,alpha=alpha)
                if flag==0:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='y.'+line,alpha=alpha,label=d.key)
                else:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='y.'+line,alpha=alpha)
            # Meridional
            if ('North' in d.key and site in d.key.lower()):
                ax[1].errorbar(timew,d.v[indw],yerr=d.ve[indw],fmt='b.'+line,alpha=alpha)
                if flag==0:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='b.'+line,alpha=alpha,label=d.key)
                else:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='b.'+line,alpha=alpha)
            # Meridional
            if ('South' in d.key and site in d.key.lower()):
                ax[1].errorbar(timew,d.v[indw],yerr=d.ve[indw],fmt='r.'+line,alpha=alpha)
                if flag==0:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='r.'+line,alpha=alpha,label=d.key)
                else:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='r.'+line,alpha=alpha)
            # CV
            if ('CV1' in d.key and site in d.key.lower()):
                ax[0].errorbar(timew,d.u[indw],yerr=d.ue[indw],fmt='m.'+line,alpha=alpha)
                ax[1].errorbar(timew,d.v[indw],yerr=d.ve[indw],fmt='m.'+line,alpha=alpha)
                if flag==0:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='m.'+line,alpha=alpha,label='CV1')
                else:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='m.'+line,alpha=alpha)
            # CV
            if ('CV2' in d.key and site in d.key.lower()):
                ax[0].errorbar(timew,d.u[indw],yerr=d.ue[indw],fmt='c.'+line,alpha=alpha)
                ax[1].errorbar(timew,d.v[indw],yerr=d.ve[indw],fmt='c.'+line,alpha=alpha)
                if flag==0:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='c.'+line,alpha=alpha,label='CV2')
                else:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='c.'+line,alpha=alpha)
            # Vertical
            if ('Zenith' in d.key and site in d.key.lower()):
                ax[2].errorbar(timew,d.w[indw],yerr=d.we[indw],fmt='ko'+line,alpha=alpha)
                if flag==0:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='ko'+line,alpha=alpha,label=d.key)
                else:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='ko'+line,alpha=alpha)
            # Inline
            if ('IN' in d.key and site in d.key.lower()):
                ax[2].errorbar(timew,d.w[indw],yerr=d.we[indw],fmt='k.'+line,alpha=alpha)
                if flag==0:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='k.'+line,alpha=alpha,label=d.key)
                else:
                    ax[3].errorbar(timet,d.T[indt],yerr=d.Te[indt],fmt='k.'+line,alpha=alpha)
    
    for k in range(3):
        ax[k].plot([dn,dn+timedelta(hours=48)],[0,0],'k--')
        ax[k].grid()
    ax[3].grid()
    ax[3].set_xlim(tlim)
    ax[3].xaxis.set_major_formatter(md.DateFormatter('%H'))
    ax[0].set_ylim(-175,175)
    ax[1].set_ylim(-150,150)
    ax[2].set_ylim(-75,75)
    ax[3].set_ylim(600,1200)
    ax[0].set_ylabel('u [m/s]')
    ax[1].set_ylabel('v [m/s]')
    ax[2].set_ylabel('w [m/s]')
    ax[3].set_ylabel('T [K]')
    ax[0].set_title('%s Level2 Data for %02i-%02i-%4i in %s'%(site.upper(),dn.month,dn.day,dn.year,tunit))
    ax[3].legend(bbox_to_anchor=(0,-.45,1,.3),ncol=5,loc=8,mode="expand",borderaxespad=0.)
    f.subplots_adjust(hspace=0.001)
    f.subplots_adjust(wspace=0.003)
    f.show()


def FindKey(cvs, key):
    out = []
    for cv in cvs:
        if cv.key == key:
            out.append(cv)
    return out


def PrintLevel2Summary(cvs, print_log=False):
    print "----------------------------------------------------------------------"
    print "index |","%-15s|" % ("key"), " error |", "len(t1)|"
    print "----------------------------------------------------------------------"
    for kk, cv in enumerate(cvs):
        print "%-6i|" % (kk), "%-15s| " % cv.key, "%5s | " % cv.error, "%-5i | " % len(cv.t1)

    if print_log:
        # prints the log for the indidual Level 2 instances.
        print "\n\n\n----------------"
        print "  log files  "
        print "----------------\n\n\n"
        for kk, cv in enumerate(cvs):
            print kk, cv.key
            print "===================================================================="
            print cv.log
            print "\n\n"

def FindClosestTime(t, t_list):
    '''
    this simple function finds the index and absolute minutes 
    difference of searching (the time-zone aware) t through 
    the list t_list

    ind, minutes_different = FindClosestTime(t, t_list)

    '''

    # some pytz examples:
    # pytz.utc
    # pytz.timezone('US/Eastern')

    # convert to utc:
    t_list = [tt.astimezone(pytz.utc) for tt in t_list]
    t = t.astimezone(pytz.utc)

    diff_list = [(tt-t) for tt in t_list]
    minutes_difference = \
            np.abs([tt.total_seconds()/60. for tt in diff_list])
    ind     = np.argmin(minutes_difference)
    minutes = np.abs(np.min(minutes_difference))
    return ind, minutes

if __name__=="__main__":
    #import matplotlib.pyplot as plt
    #from matplotlib.pyplot import *
    import time
    #import matplotlib
    #from matplotlib.dates import date2num
    #from datetime import timedelta
    import pytz
    import os


    project = 'nation'
    dn = datetime(2013,3,16)
    cvs = GetLevel2(project, dn, w_is_0=False)

    PrintLevel2Summary(cvs, print_log=False)


    #par_east  = cvs[1].t1
    #par_north = cvs[2].t1

    #ind, minutes = FindClosestTime(par_east[10], par_north)


    #dn1 = dn.replace(hour=14)
    #dn2 = (dn + timedelta(days=1)).replace(hour=8)



    #PlotLatSLT(cvs, dn1, dn2)


    cvs = GetLevel2('renoir', datetime(2012,11,2), w_is_0=False)




#    # PlotLatSLT testing
#    # ----------------------
#    savepath = "/home/duly/plt_nation/"
#    close('all')
#    print "removing pngs..."; os.system("rm %s/*.png" % savepath)
#    dn_start = datetime(2013,1,1)
#    dn_end   = datetime(2013,6,1)
#    dns = [dn_start + timedelta(days=kk) for kk in range((dn_end-dn_start).days)]
#    #dns = [datetime(2013,4,25)]
#    #dns = [datetime(2013,3,2), datetime(2013,3,13), datetime(2013,3,16), datetime(2013,4,19)]
#    #dns = [datetime(2013,3,2)]
#    dns = [datetime(2013,3,13)]
#
#    los_e = { 'EKU_PAR': [], 'UAO_EKU': [] }
#    cv_e  = { 'EKU_PAR': [], 'UAO_EKU': [] }
#
#    for dn in dns:
#        print dn
#        #project = "RENOIR"
#        project = "NATION"
#        reference = 'Laser'
#    
#        Sites = MasterDictionary()
#    
#        # Get Data for a Day:
#        # get data for every site location:
#        # -----------------------------------
#        cvs = []
#        pairs = []
#        for site1 in Sites[project]['all']:
#            cvs += CardFinder(dn, site1, project, reference)
#            for combo in Sites[project][site1]['Combos']:
#                if combo not in pairs:
#                    cvs += CVFinder(dn, site1, combo, project, reference)
#    
#        # remove error parts in cvs:
#        cvs = [cv for cv in cvs if cv.errorT is False]
#    
#        #dn1 = datetime(2013,4,25,14)
#        dn1 = dn.replace(hour=14)
#        dn2 = (dn + timedelta(days=1)).replace(hour=8)
#        #datetime(2013,4,26,8)
#        
#        for kk, cv in enumerate(cvs):
#            cv.cut(dn1, dn2)
#            if cv.error is False and len(cv.t1)==0: cv.error=True
#    
#        # plot lat vs slt:
#        # -------------------------------------
#        PlotLatSLT(cvs, dn1, dn2, switch_interpolate_T=False)
#        savefig(savepath + "%s_lat_slt_" % (project,) + dn.strftime("%Y-%m-%d") + ".png")


    #    # testing:
    #    # --------
    #    cvs = [cv for cv in cvs if cv.error is False]

    #    #cvs = [cv for cv in cvs if cv.key=="CV_EKU_UAO_1"]

    #    #for kk, cv in enumerate(cvs):
    #    #    if cv.flag_cloud:
    #    #        print kk, "%3i of %3i are cloudy" % (sum(cv.cloudy), len(cv.cloudy)), cv.key

    #    print " ====================================="
    #    print "                lose            ue/ve"
    #    for cv in cvs:
    #        if 'EKU' in cv.key and 'PAR' in cv.key \
    #                or 'UAO' in cv.key and 'EKU' in cv.key:
    #                    if 'EKU' in cv.key and 'PAR' in cv.key:
    #                        key = 'EKU_PAR'
    #                    elif 'UAO' in cv.key and 'EKU' in cv.key:
    #                        key = 'UAO_EKU'
    #                    else:
    #                        print "ERROR"
    #                        sys.exit(-999)
    #                    if 'IN' == cv.key[:2]:
    #                        a=2
    #                        if len(cv.we) != len(cv.los_sigma1):
    #                           print "ERROR: we, los_sigma1 do not line up"
    #                           sys.exit(-1)
    #                        for we, lose1 in zip(cv.we, cv.los_sigma1):
    #                            cv_e[key].append( we )
    #                            los_e[key].append( lose1 )
    #                    elif 'CV' == cv.key[:2]:
    #                        a = 2
    #                        #if len(cv.ve) != len(cv.ue) != len(cv.los_sigma1) != len(cv.los_sigma2):
    #                        #    print "ERROR: [u,v]e, los_sigma[1,2] do not line up"
    #                        #    sys.exit(-1)
    #                        #for ve, ue, lose1, lose2 in zip(cv.ve, cv.ue, cv.los_sigma1, cv.los_sigma2):
    #                        #    #print cv.key, "%10.2f, %10.2f" %(lose, ue)
    #                        #    cv_e[key].append( ue )
    #                        #    los_e[key].append( lose1 )

    #                        #    cv_e[key].append( ue )
    #                        #    los_e[key].append( lose2 )

    #                        #    cv_e[key].append( ve )
    #                        #    los_e[key].append( lose1 )

    #                        #    cv_e[key].append( ve )
    #                        #    los_e[key].append( lose2 )
    #                    else:
    #                        print "ERROR2"
    #                        sys.exit(-999)
    #                    #print key, cv.key
    #                    a = cv
    #figure(1); clf()
    #plot( los_e['EKU_PAR'], cv_e['EKU_PAR'], 'rx')
    #plot( los_e['UAO_EKU'], cv_e['UAO_EKU'], 'bo')
    #plot( [0., 30.,], [0., 30.,])
    #grid();
    #axis('equal')
    #draw(); show()



'''
        # map plots
        # -----------------
        verbose = False
        close('all')
        os.system('rm -rf ./tim_test/*.png')
        from mpl_toolkits.basemap import Basemap
        import matplotlib.cm as cm
        ALPHA_VALUE = 0.3
        vmin, vmax = 600, 1200.

        utc = pytz.utc

        from ParseSuperDarn import ParseSuperDarn
        print "parsing superdarn..."
        superdarn_fhe = ParseSuperDarn( '20130214.fhe.f2t' )
        superdarn_bks = ParseSuperDarn( '20130214.bks.f2t' )
        print "done."


        bunch_of_times = []
        for cv in cvs:
            for dn in cv.t1:
                bunch_of_times.append( dn.astimezone(utc) )
        dn_min, dn_max = min(bunch_of_times), max(bunch_of_times)
        dn_min = dn_min - timedelta(minutes=10)
        dn_max = dn_max + timedelta(minutes=10)

        #dn_min = datetime(2013,4,26,2,0,tzinfo=utc )
        #dn_max = datetime(2013,4,26,3,10,tzinfo=utc )

        minutes = (dn_max-dn_min).seconds/60.
        plot_times = [dn_min + timedelta(minutes=kk) for kk in range(minutes)][::5]
        #a = []
        #for cv in cvs:
        #    if cv.key=="PAR_West": a.append(cv)
        #    if cv.key=="PAR_South": a.append(cv)
        #cvs = a


        tec_bin = []
        mlims = {\
                'NATION': {'llcrnrlon': -100,\
                           'llcrnrlat': 30,\
                           'urcrnrlon': -75,\
                           'urcrnrlat': 55,\
                           },\
                'RENOIR': {'llcrnrlon': -48,\
                           'llcrnrlat': -12,\
                           'urcrnrlon': -32,\
                           'urcrnrlat': 2,\
                           }\
                }
                    
        for kk, t in enumerate(plot_times):
            print "%04d of %04d" % (kk, len(plot_times))
            fig = figure(1); clf()
            m = Basemap(\
                    llcrnrlon= mlims[project]['llcrnrlon'],\
                    llcrnrlat= mlims[project]['llcrnrlat'],\
                    urcrnrlon= mlims[project]['urcrnrlon'],\
                    urcrnrlat= mlims[project]['urcrnrlat'],\
                    projection='merc', area_thresh=1000,\
                    resolution='c',\
                    #lat_1=45,lat_2=55,lat_0=40,lon_0=-85,\
                    )
            m.drawcoastlines()
            m.drawstates()


            # plot superdarn:
            # --------------------
            def cosd(x): return np.cos(x*180./np.pi)
            def sind(x): return np.sin(x*180./np.pi)
            for superdarn in [superdarn_fhe, superdarn_bks]:
                min_list = [abs(dn-t.astimezone(utc)).seconds for dn in np.sort(superdarn.keys())]
                jj = np.argmin(min_list)
                diff = abs( np.sort(superdarn.keys())[jj] - t)
                if diff.seconds < 20*60. and diff.days==0:
                    dn = np.sort(superdarn.keys())[jj]
                    print "--------------------------------------------"
                    print "for t = ",t
                    print "we have selected"
                    print "dn = ",dn
                    print "a difference of ", diff.seconds/60., "  mins"
                    for LOS, az, sza, gsf, glat, glon in zip(\
                            superdarn[dn]['velocity'],\
                            superdarn[dn]['gazm'],\
                            superdarn[dn]['sza'],\
                            superdarn[dn]['GroundScatterFlag'],\
                            superdarn[dn]['glat'],\
                            superdarn[dn]['glon'],\
                            ):
                        if gsf==0:
                            el = 90.-sza
                            u = LOS*cosd(el)*cosd(az)
                            v = LOS*cosd(el)*sind(az)
                            mlon, mlat = m(glon, glat)
                            Q = m.quiver(mlon,mlat,\
                                    0.,v,\
                                    scale = 1000.,\
                                    units = 'width',\
                                    color = 'g',\
                                    angles = 'uv',\
                                    zorder = 50,\
                                    )

            # plot GPS
            # --------------------
            Nlon = 25
            Nlat = 20
            lons = np.linspace(m.lonmin-1., m.lonmax+1, Nlon)
            lats = np.linspace(m.latmin-1, m.latmax+1, Nlat)
            LON, LAT = np.meshgrid( lons, lats )

            tec_grid = np.zeros((Nlat, Nlon, 50,)) * float('nan')
            counter = np.zeros((Nlat,Nlon))
            
            from ParseGPS import gps_data, gps_dns
            min_list = [abs(dn - t.astimezone(utc)).seconds for dn in gps_dns]
            jj = np.argmin(min_list)
            diff = abs(gps_dns[jj] - t)
            if diff.seconds < 20*60. and diff.days==0:
                print "----------------------"
                print "for t=",t
                print "select gps_dns[jj]=",gps_dns[jj]
                print "a difference of"
                print abs(gps_dns[jj]-t.astimezone(utc)).seconds / 60. , " mins"
                for lat, lon, tec in zip(\
                        gps_data[gps_dns[jj]]['lat'],\
                        gps_data[gps_dns[jj]]['lon'],\
                        gps_data[gps_dns[jj]]['tec'],\
                        ):

                    lon_index = np.argmin(np.abs(lon-lons))
                    lat_index = np.argmin(np.abs(lat-lats))

                    # this is to prevent far away tec values from entering in
                    # on the edges:
                    if 0 < lon_index < Nlon-1 and 0 < lat_index < Nlat-1:
                        tec_grid[ lat_index, lon_index,\
                                counter[lat_index,lon_index]\
                                ] = tec
                        counter[lat_index, lon_index] += 1

                    tec_bin.append(tec)
            MLON, MLAT = m(LON, LAT)
            # collapse down tec_grid:
            tec_grid_collapsed = np.mean( np.ma.masked_array(tec_grid, np.isnan(tec_grid)),\
                    axis=2)
            m.pcolor(MLON,MLAT,\
                    np.log10(tec_grid_collapsed),\
                    vmin=0.0, vmax=1.8,\
                    #cmap = cm.Greys_r,\
                    )
            # --------------------

            # for plotting cloud info @ site location:
            for site_zenith in ['UAO_Zenith', 'EKU_Zenith', 'ANN_Zenith', 'PAR_Zenith']:
                match_cv = None
                for cv in cvs:
                    if cv.key == site_zenith:
                        match_cv = cv
                        break
                if match_cv is not None:
                    (lat,lon) = GetLocation(SITE,match_cv.key)
                    mlon, mlat = m(lon, lat)
                    c = 'bo'
                    ts = [ttt.astimezone(utc) for ttt in match_cv.allt]
                    min_list = [ abs(t-tss).seconds for tss in ts]
                    jj = np.argmin(min_list)
                    diff = abs(t-ts[jj])
                    if diff.seconds < 20.*60 and diff.days==0:
                        cloud = match_cv.allcloud[jj]
                        if cloud < cloudthreshold:
                            c = 'bo'
                        else:
                            c = 'ro'
                        if cloud < -998: c = 'yo'
                    else:
                        c = 'ko'
                    
                    m.plot(mlon,mlat,c)
                    lon = SiteLocations[match_cv.key]['lon'] +0.2
                    lat = SiteLocations[match_cv.key]['lat'] +0.2
                    mlon, mlat = m(lon, lat)
                    text( mlon, mlat,match_cv.key[:3])

            for cv in cvs:
                if cv.error is False and len(cv.t1) > 0:

                    t1 = np.array([dn.astimezone(utc) for dn in cv.t1])

                    min_list = [ (t - t11).seconds for t11 in t1]
                    jj = np.argmin(min_list)
                    #print "l.t. 30 minutes; ", (t-t1[jj]).seconds/60., jj
                    diff = t-t1[jj]
                    if 0. < diff.seconds < 30*60. and diff.days==0:
                        mins_away = abs(t-t1[jj]).seconds / 60.
                        #print "mins_away=",mins_away

                        # average last 20 minutes:
                        ind = np.where(np.array([b.seconds/60. for b in t1[jj]-t1]) < 20.)

                        if verbose:
                            print "------------------------"
                            print "for t=",t
                            print cv.key
                            print "from the times of "
                            for kkk, ttt in enumerate(t1):
                                print "   ",kkk, ttt
                            print "we have selected jj=",jj, t1[jj]
                            print "and within 20 minutes: "
                            print ind
                            print t1[ind]

                        flag_cv = False
                        flag_quiver = False
                        flag_card = False

                        def FindMatch(cvs, t, key):
                            error = False
                            cv_match = None
                            for cv in cvs:
                                if cv.key == key:
                                    cv_match = cv
                                    break
                            ind = None
                            if cv_match is not None:
                                #print "cv_match.key=", cv_match.key
                                ts = np.array([ttt.astimezone(utc) for ttt in cv_match.t1])
                                min_list = [ (t -tss).seconds for tss in ts]
                                jj = np.argmin(min_list)
                                if 0. < (t-ts[jj]).seconds < 30 * 60.:
                                    ind = np.where(np.array([b.seconds/60. for b in ts[jj]-ts]) < 20.)
                            else:
                                error = True
                            if ind == None: error = True
                            #print "ind=",ind
                            return cv_match, ind, error

                        if 'CV' == cv.key[:2]:
                            #u = np.average( cv.u[ind], weights=1./cv.ue[ind] )
                            v = np.average( cv.v[ind], weights=1./cv.ve[ind] )
                            flag_cv = True
                            flag_quiver = True

                            u = 0.
                            if np.abs(v) > 400: flag_quiver = False
                        #elif 'East' in cv.key:
                            #u = np.average( cv.u[ind], weights=1./cv.ue[ind] )
                            #cv_match, ind_match, error = FindMatch(cvs, t, cv.key[:3] + "_" + 'North')
                            #if error is False:
                            #    v = np.average( cv_match.v[ind_match], weights=1./cv_match.ve[ind_match])
                            #    flag_quiver = True
                            #    flag_card = True
                            #else:
                            #    flag_quiver = False

                        #elif 'West' in cv.key:
                            #u = np.average( cv.u[ind], weights=1./cv.ue[ind] )
                            #cv_match, ind_match, error = FindMatch(cvs, t, cv.key[:3] + "_" + 'South')
                            #if error is False:
                            #    v = np.average( cv_match.v[ind_match], weights=1./cv_match.ve[ind_match])
                            #    flag_quiver = True
                            #    flag_card = True
                            #else:
                            #    flag_quiver = False

                        #elif 'East' in cv.key or 'West' in cv.key:

                        #    u = np.average( cv.u[ind], weights=1./cv.ue[ind] )
                        #    v = 0.
                        #    flag_quiver = True
                        #    if np.abs(u) > 400: flag_quiver = False
                        elif 'North' in cv.key or 'South' in cv.key:
                            u = 0.
                            v = np.average( cv.v[ind], weights=1./cv.ve[ind] )
                            flag_quiver = True
                            if np.abs(v) > 400.: flag_quiver = False


                        # jonathan (originally) had restrictions on plotting T and u,v
                        # if Te and Dopplere were less than 100, 25... should I include?
                        lon = SiteLocations[cv.key]['lon']
                        lat = SiteLocations[cv.key]['lat']
                        mlon, mlat = m(lon, lat)
                        #m.scatter(mlon,mlat,\
                        #        #c=cv.T[jj],\
                        #        c = np.average( cv.T[ind], weights=1./cv.Te[ind] ),\
                        #        vmin=vmin, vmax=vmax,\
                        #        s = 250. + (50.-250.)/30. * mins_away,\
                        #        alpha = ALPHA_VALUE if cv.cloudy[jj] else None,\
                        #        edgecolor='None',\
                        #        )
                        if flag_card:
                            # average lat, lon locations for map:
                            if 'East' in cv.key:
                                match_key = cv.key[:3] + "_" + 'North'
                            elif 'West' in cv.key:
                                match_key = cv.key[:3] + "_" + 'South'

                            lon = np.mean([\
                                    SiteLocations[cv.key]['lon'] ,\
                                    SiteLocations[match_key]['lon']\
                                    ])
                            lat = np.mean([\
                                    SiteLocations[cv.key]['lat'] ,\
                                    SiteLocations[match_key]['lat']\
                                    ])
                            mlon, mlat = m(lon, lat)
                        if flag_quiver:
                            Q = m.quiver(mlon,mlat,\
                                    u, v, \
                                    scale=1000.,\
                                    units='width',\
                                    alpha = 1. + (0-1.)/30. * mins_away,\
                                    color = 'k' if flag_cv else 'r',\
                                    angles = 'uv',\
                                    visible = False if cv.cloudy[jj] else True,\
                                    )
                            Qbogus = m.quiver(m(0,0)[0],m(0,0)[1],\
                                    0., 0.,\
                                    scale=1000.,\
                                    width=0.005,\
                                    angles='uv',\
                                    )
                            qk = quiverkey(Qbogus, 0.9, 0.9, 100.,\
                                    r'$100 \mathrm{m/s}$',\
                                    color='k', \
                                    )
                            Qbogus.set_visible(False) # just to make sure...
                    title("%s Thermospheric Winds and Temperatures\n" % (project.upper(),) +
                            t.strftime('%Y-%m-%d') + 
                            " UTC " + t.strftime('%H:%M'))
            #m.scatter(m(0,0)[0],m(0,0)[0],c=0,s=0,vmin=vmin,vmax=vmax) # only to get a colorbar each and every plot...
            m.scatter(m(0,0)[0],m(0,0)[0],c=0,s=0,vmin=0.0,vmax=1.8) # only to get a colorbar each and every plot...
            colorbar()
            #fig.tight_layout()
            draw(); savefig("./tim_test/" + "NATION_%04d.png" % kk)
        draw(); show()
'''







'''
    # put monthly averages down, using Dan's routine:
    # -----------------------------------------------
    for site in Sites[cvs[0].project]['all']:
        print "-----------"
        print "inputs to BinMonthlyData:"
        print "site, year, month = ", site, dn1.year, dn1.month
        o = BinMonthlyData(site, dn1.year, dn1.month)
        t1 = [ datetime(dn1.year, dn1.month, dn1.day, dn.hour, dn.minute, dn.second, tzinfo= pytz.utc) for dn in o.t1]
        t2 = [ datetime(dn2.year, dn2.month, dn2.day, dn.hour, dn.minute, dn.second, tzinfo= pytz.utc) for dn in o.t1]
        lat = SiteLocations[site+"_Zenith"]['lat']
        lon = SiteLocations[site+"_Zenith"]['lon']
        # copy and paste data for both dn1 and dn2, ensuring that when 
        # plotted, the information is visible:
        t1_slt = [dn2utc(a) + timedelta(hours=lon/15.) for a in t1] 
        t2_slt = [dn2utc(a) + timedelta(hours=lon/15.) for a in t2] 
        for ttt1, ttt2, u, v in zip(t1_slt, t2_slt, o.u, o.v):
            nan = float('nan')
            if (not np.isnan(u)) and (not np.isnan(v))\
                    and np.sqrt(u**2+v**2) < 400.:
                Q = quiver(matplotlib.dates.date2num(ttt1), lat, u, v,\
                        scale=qscale,\
                        units=qunits,\
                        width=qwidth,\
                        edgecolor='k',\
                        facecolor='w',\
                        linewidth=0.5,
                        #linestyle='dashed',\
                        )
                Q = quiver(matplotlib.dates.date2num(ttt2), lat, u, v,\
                        scale=qscale,\
                        units=qunits,\
                        width=qwidth,\
                        edgecolor='k',\
                        facecolor='w',\
                        linewidth=0.5,
                        #linestyle='dashed',\
                        )
'''
