# A module containing functions useful for the GITM/FPI comparison paper
# Brian Harding 2018 Mar 7

import pytz
from pyglow import pyglow
from datetime import datetime, timedelta
import FPIprocess
import fpiinfo
import numpy as np
from numpy import array, zeros, ones, arange, sqrt, mean, median, nan, pi, sum, random
import pandas as pd
import glob


def get_max_kp(t, prevhrs = 24.):
    kp = []
    for h in arange(-prevhrs, 1., 3.):
        ti = t + timedelta(hours=h)
        pt = pyglow.Point(ti, 0, 0, 250)
        kp.append(pt.kp)
    return max(kp)


    
def get_raw_fpi(instr_name, year, doy_start, doy_stop, SIGMA_THRESH = 100. , CLOUD_THRESH = np.inf, DROP_NO_CLOUD = True,
                SKYI_THRESH = None, SKYB_THRESH = np.inf, T_THRESH = 150., verbose = True):
    '''
    Load raw line-of-sight FPI data and perform quality control.
    Return a pandas DataFrame with all trustworthy samples. 
    Note that this function may return some samples before and after
    the specified time window, because the UT/LT crossover can be
    confusing, and because if the bin size is large, the previous and
    next day's data may be important.
    
    INPUTS:
    
    instr_name      - e.g., 'minime05'
    year            - e.g., 2013
    doy_start       - UT day-of-year. This is UT, not the FPI "night of" convention.
    doy_stop        - ". Inclusive.
    SIGMA_THRESH    - [m/s, K] threshold for filtering LOS samples based on their uncertainty
    CLOUD_THRESH    - [C] threshold for filtering LOS samples based on cloud indicator. Can be
                      filtered in this step or during the cardinalization step.
    DROP_NO_CLOUD   - If False, keep data with no cloud sensor reading. If True, drop it.
    SKYI_THRESH     - [cnt/s] threshold for filtering LOS samples based on brightness (if None, 
                      use default value from fpiinfo [the lower of the two])
    SKYB_THRESH     - [cnt/s] threshold for filtering LOS samples based on spectral background
    T_THRESH        - [K] temperatures smaller than this are indicative of a bad fit
    verbose         - whether to print quality control information
    
    OUTPUTS:
    
    df   - pandas DataFrame
    '''
    
    t0 = datetime(year, 1, 1) + timedelta(days=doy_start-1)

    if SKYI_THRESH is None:
        instr = fpiinfo.get_instr_info(instr_name, t0)
        SKYI_THRESH = instr['skyI_quality_thresh'][1] # Use the liberal threshold, and assume zenith reference will help

    ############################### Load data ################################

    # Figure out which FPI data to load in order to fill in the right UT bins
    fpi_years = []
    fpi_doys = []
    

    dfs = []
    for doyoff in range(-2, doy_stop-doy_start+2): # That should be wide enough to handle any time zone
        t = pd.to_datetime(t0 + timedelta(days=doyoff))
        try:
            r = FPIprocess.load_level0(instr_name, t.year, t.dayofyear)['FPI_Results']
        except IOError: # No data
            continue

        # Convert night of data to a DataFrame
        d = {} # FPI data to transform to DataFrame
        d['t'] = pd.to_datetime(r['sky_times']).tz_convert(pytz.utc)
        # Variables to pass through directly
        v = ['skyI', 'sigma_skyI', 'LOSwind', 'sigma_LOSwind', 'sky_intT', 'T', 'sigma_T',\
             'sigma_fit_LOSwind', 'sigma_cal_LOSwind','az','ze','skyB','direction','sky_ccd_temperature']
        for key in v:
            d[key] = r[key]
        # Note using max instead of mean. If it's cloudy at all during the exposure, it's bad.
        d['Clouds'] = r['Clouds']['max'] 

        df = pd.DataFrame(d)
        dfs.append(df)

    # if no data, return None
    if len(dfs)==0:
        return None
        
    # Make one big DataFrame for all the FPI data
    dfraw = pd.concat(dfs)
    dfraw.index = arange(len(dfraw)) # make sure indices are unique
    

    #################### Quality control on samples #######################
    bmiss  = dfraw['Clouds'].isnull()
    dfraw.loc[bmiss,'Clouds'] = -999. # Treat as clear
    bcloud = dfraw['Clouds'] > CLOUD_THRESH
    bsigma = (dfraw['sigma_LOSwind'] > SIGMA_THRESH) | (dfraw['sigma_T'] > SIGMA_THRESH)
    bskyI  = dfraw['skyI'] < SKYI_THRESH
    bskyB  = dfraw['skyB'] > SKYB_THRESH
    bT     = dfraw['T'] < T_THRESH
    bunk   = dfraw['direction'] == 'Unknown'

    dfraw[ bcloud | bsigma | bskyI | bskyB | bT | bunk ] = np.nan
    if DROP_NO_CLOUD:
        dfraw[bmiss] = np.nan

    if verbose:
        print ''
        if DROP_NO_CLOUD:
            print('%.2f%% samples dropped for missing cloud sensor' % (100.*sum(bmiss)/len(bmiss)))
        print('%.2f%% samples cloudy' % (100.*sum(bcloud)/len(bcloud)))
        print('%.2f%% samples with large uncertainty' % (100.*sum(bsigma)/len(bsigma)))
        print('%.2f%% samples with low brightness' % (100.*sum(bskyI)/len(bskyI)))
        print('%.2f%% samples with large background' % (100.*sum(bskyB)/len(bskyB)))
        print('%.2f%% samples with low temperature' % (100.*sum(bT)/len(bT)))
        print('%.2f%% samples with unknown direction' % (100.*sum(bunk)/len(bunk)))
        bad = dfraw.isnull().any(axis=1)
        print('  --> %.2f%% samples deemed bad' % (100.*sum(bad)/len(bad)))

    dfraw = dfraw.dropna()
    
    ################### Append metadata #######################
    #dfraw.instr_name = instr_name
    
    return dfraw
    

    


def cardinal_fpi(dfraw, year, doy_start, doy_stop, dtb=1.0, KP_THRESH = 3., DISAGREE_THRESH = 25., 
                 NBIN_THRESH = 3, COND_THRESH = 1e2, SIGMA_W_THRESH=100., CLOUD_THRESH = np.inf, verbose = True):
    '''
    Take line-of-sight winds, and transform them to cardinal winds on 
    a constant time cadence. Time binning is performed simultaneously
    with cardinalization. This works for a single site. Assume that
    vertical wind is zero and use zenith reference as a function
    of time.
    
    Return a pandas DataFrame at the specified time cadence.
    
    INPUTS:
    
    dfraw           - pandas.DataFrame with line-of-sight samples (see get_raw_fpi)
    year            - e.g., 2013
    doy_start       - UT day-of-year, 0 UT, start of time bins. This is UT not the FPI "night of" convention.
    doy_stop        - ", end of time bins. Inclusive.
    dtb             - [hours] time resolution for binning
    KP_THRESH       - threshold for filtering LOS samples based on Kp
    DISAGREE_THRESH - [m/s] threshold for filtering daily bins based on how much the values in 
                      the bin disagree when converting LOS to cardinal winds (RMS of
                      least squares fit)
    NBIN_THRESH     - minimum number of samples per bin (this input needs to be at least 3)
    COND_THRESH     - threshold for filtering daily bins based on the condition number of the 
                      matrix used to convert LOS to cardinal winds. This ensures that enough 
                      varied lines of sight are included.
    SIGMA_THRESH    - [m/s, K] maximum uncertainty of u, w, and T in each bin, others will be filtered.
    CLOUD_THRESH    - [K] maximum cloud indicator. If any samples in bin are larger, this bin will be filtered.
                      Cloud filtering can be done at the raw data level or at this level. This input can either
                      be a scalar, or a function which takes doy and returns a scalar.
    verbose         - whether to print quality control information
                      
    OUTPUTS:
    
    df   - pandas DataFrame
    
    '''

    ################### Binning and Cardinal Winds #####################
    # Two simultaneous steps: time binning, and cardinalizing (i.e., converting line of sight 
    # winds to u,v,w). The assumption is that the vertical wind is zero, which helps us 
    # better handle Doppler reference errors due to laser drift and OH contamination.
    # 
    # Additional quality control steps are included here: number of samples per bin, whether
    # samples within a bin agree with each other (measured by cardinalization residual), and 
    # Kp during the previous 24 hours

    # Define time bins
    start_date = datetime(year, 1, 1) + timedelta(days=doy_start-1) # 0 UT start of day on doy_start
    stop_date  = datetime(year, 1, 1) + timedelta(days=doy_stop) # 0 UT end of day on doy_stop
    timedelta_bin = timedelta(hours=dtb)
    timebin_edges = pd.date_range(start_date-timedelta_bin/2, stop_date-timedelta_bin/2, freq=timedelta_bin)
    #timebin = timebin_edges[:-1] + timedelta_bin/2

    # Bin the data into the desired time bins
    cuts = pd.cut(dfraw['t'], timebin_edges)
    grouped = dfraw.groupby(cuts)

    # Loop through the bins and do the cardinalization
    dfrows = [] # rows which will be constructed into an output dataframe
    for interval, g in grouped:            
        # Formulate least squares problem to estimate cardinal wind:
        # loswind = A * [u,v,doppref]
        A = zeros((len(g),3))
        A[:,0] = np.sin(g['az']*pi/180)*np.sin(g['ze']*pi/180) # zonal contribution
        A[:,1] = np.cos(g['az']*pi/180)*np.sin(g['ze']*pi/180) # meridional contribution
        A[:,2] = ones(len(g)) # Doppler reference contribution
        
        if len(g) >= NBIN_THRESH and np.linalg.matrix_rank(A) >= 3:
        
            # Convert to weighted least squares and solve
            Ap = np.diag(1/g['sigma_LOSwind']).dot(A)
            yp = g['LOSwind']/g['sigma_LOSwind']
            M = np.linalg.inv(Ap.T.dot(Ap)).dot(Ap.T) # least squares matrix, explicitly
            uvd = M.dot(yp)
            sigma_uvd = sqrt(np.diag(M.dot(M.T)))
            resid = g['LOSwind'] - A.dot(uvd)
            if any(np.isnan(resid)):
                print g['LOSwind']
                print A
                print uvd
                raise Exception()
            rms_resid = sqrt(mean(resid**2))
            
            # Also propagate errors using purely sigma_cal and sigma_fit
            Ap = np.diag(1/g['sigma_fit_LOSwind']).dot(A)
            M = np.linalg.inv(Ap.T.dot(Ap)).dot(Ap.T) # least squares matrix, explicitly
            sigma_fit_uvd = sqrt(np.diag(M.dot(M.T)))  
            
            Ap = np.diag(1/g['sigma_cal_LOSwind']).dot(A)
            M = np.linalg.inv(Ap.T.dot(Ap)).dot(Ap.T) # least squares matrix, explicitly
            sigma_cal_uvd = sqrt(np.diag(M.dot(M.T)))      
            
            # Compute airglow gradient (variability within bin, normalized by mean)
            Ivert = g['skyI']/np.cos(g['ze']*np.pi/180.) # converted from slant to vertical
            Ivariab = Ivert.std()/Ivert.mean()
            
            # Compute cloud threshold for this day
            if callable(CLOUD_THRESH): # It's a function of doy
                doy = interval.left.dayofyear
                cthresh = CLOUD_THRESH(doy)
            else: # It's a scalar
                cthresh = CLOUD_THRESH
            
            

            # Create a row for this interval of time
            dfrows.append([uvd[0], # zonal wind
                           uvd[1], # meridional wind
                           uvd[2], # Doppler reference
                           np.nansum(g['T']/g['sigma_T']**2)/np.nansum(1/g['sigma_T']**2), # Temperature
                           np.nanstd(g['T']), # Variability of temperature in the bin (not the uncertainty of the mean)
                           sigma_uvd[0], # uncertainty in u
                           sigma_uvd[1], # uncertainty in v
                           sigma_uvd[2], # uncertainty in Dopp ref
                           sigma_fit_uvd[0], # fit uncertainty in u
                           sigma_fit_uvd[1], # fit uncertainty in v
                           sigma_fit_uvd[2], # fit uncertainty in Dopp ref
                           sigma_cal_uvd[0], # cal uncertainty in u
                           sigma_cal_uvd[1], # cal uncertainty in v
                           sigma_cal_uvd[2], # cal uncertainty in Dopp ref
                           rms_resid,    # disagreement in each 
                           len(g),       # number of samples
                           np.linalg.cond(A), # condition number of cardinalization matrix
                           g['Clouds'].max(), # max cloud indicator
                           cthresh, # cloud threshold for this day
                           g['sky_ccd_temperature'].mean(), # mean CCD temperature
                           Ivariab, # Fractional variability with each bin (proxy for airglow gradient)
                          ]
                         )
        else:
            dfrows.append([np.nan]*21)

    # Create DataFrame object to store the binned data
    df = pd.DataFrame(data = dfrows, 
                      index  = timebin_edges[:-1] + timedelta_bin/2,
                      columns = ['u','v','doppref','T','variab_T', 'sigma_u','sigma_v','sigma_doppref',\
                                 'sigma_fit_u', 'sigma_fit_v', 'sigma_fit_doppref',
                                 'sigma_cal_u', 'sigma_cal_v', 'sigma_cal_doppref',
                                 'residrms','Nsamp','cond','cloud','cloud_thresh', 'ccdtemp','variab_I'])

    # Fill in "previous 24 hour's max kp"
    kp = []
    for t, row in df.iterrows():
        if not row.isnull().any(): # this is a good row, compute Kp
            kp.append(get_max_kp(t, prevhrs=24.))
        else:
            kp.append(np.nan)
    df['kp24'] = kp

    # Quality control
    bkp    = df['kp24'] > KP_THRESH
    bresid = df['residrms'] > DISAGREE_THRESH
    bcond  = df['cond'] > COND_THRESH
    bsigma = (df['sigma_u'] > SIGMA_W_THRESH) | (df['sigma_v'] > SIGMA_W_THRESH)
    bcloud = df['cloud'] > df['cloud_thresh']
    vars_to_nan = ['u','v','T','sigma_u','sigma_v','variab_T','doppref','sigma_doppref',
                                 'sigma_fit_u', 'sigma_fit_v', 'sigma_fit_doppref',
                                 'sigma_cal_u', 'sigma_cal_v', 'sigma_cal_doppref',]
    df.loc[ bkp | bresid | bcond | bsigma | bcloud , vars_to_nan ] = np.nan

    if verbose:
        # print how many bins are being removed for too few samples
        Ndata = np.sum(grouped.size() > 0) # how many bins started off with some data
        if Ndata==0:
            print 'No data'
        else:
            small = (grouped.size() < NBIN_THRESH) & (grouped.size() > 0)
            print '%.2f%% of bins filtered (<%i samples per bin or not enough variety)' % (100.*sum(small)/Ndata, NBIN_THRESH)
            print '%.2f%% of bins filtered (prev 24 hrs Kp > %.2f)' % (100.*sum(bkp)/Ndata, KP_THRESH)
            print '%.2f%% of wind bins filtered (condition number > %.1f)' % (100.*sum(bcond)/Ndata, COND_THRESH)
            print '%.2f%% of wind bins filtered (cardinal residual > %.2f m/s)' % (100.*sum(bresid)/Ndata, DISAGREE_THRESH)
            print '%.2f%% of wind bins filtered (propagated uncertainty > %.2f m/s)' % (100.*sum(bsigma)/Ndata, SIGMA_W_THRESH)
            print '%.2f%% of wind bins filtered (cloud indicator > threshold)' % (100.*sum(bcloud)/Ndata)

    # Metadata
    #df.instr_name = dfraw.instr_name

    return df
    


def reindex_date_slt(df, lon):
    '''
    Take a dataframe which is (singly-)indexed by UT time, and reindex it using a multiindex of date and SLT.
    '''
    lon = np.mod(lon + 180., 360.) - 180.
    hr = df.index.hour + df.index.minute/60.
    date = df.index.date
    slt_hr = hr.values + 24./360*lon
    for i in range(len(slt_hr)):
        if slt_hr[i] > 12:
            slt_hr[i] -= 24
            date[i] += timedelta(days=1)
    df['slt_hr'] = slt_hr
    df['slt_date'] = date
    return df.set_index(['slt_date','slt_hr'])
    
    
    
    
def subtract_running_mean(p, pe=None,w=30, nthresh=4, extrap=True, verbose=False):
    '''
    Remove centered w-day mean from the dataframe with index SLT and columns dates.
    After this, a correction is applied to make the whole dataframe have zero
    mean. This is usually a small correction (3 m/s)
    
    INPUT:
        p - DataFrame, probably created by something like:
            p = df.pivot_table(index='slt_hr',columns = 'slt_date', values='u', dropna=True)
        pe - uncertainty of each element of p. OPTIONAL: if None, a dummy DataFrame full of
             nans is returned for the propagated uncertainty
        w - full size of running window in days
        nthresh - If running mean has less than nthresh entries, void the result
        extrap - If True, extrapolation the running mean to the first and last half-window
    OUTPUT:
        p - Same as input but with running mean removed
        pe - uncertainty of p if pe is provided, otherwise it will be filled with nans
        m - The means that were subtracted (same format as p)
    '''

    # Remove the running mean from each sample
    dates = []
    hours = []
    dvals = []
    dfrows = []
    meanrows = []
    nrows = []
    sigmarows = [] # uncertainty of variability (dfrows)
    
    # If pe was not specified, replace it with dummy DataFrame
    # so calculations can progress.
    sigmavalid = True
    if pe is None:
        sigmavalid = False
        pe = p.copy()
        
    
    tstart = p.columns[0]
    tstop  = p.columns[-1]
    for date in pd.date_range(tstart, tstop):
        if date in p.columns: # there's data on this date
            p2 = p.loc[:,date]
            pe2 = pe.loc[:,date]
            # Compute mean over time window and remove
            t0 = pd.to_datetime(date - timedelta(days=w/2))
            t1 = pd.to_datetime(date + timedelta(days=w/2))
            # Handle extrapolation
            if t0 < tstart:
                if extrap:
                    t1 = t1 + timedelta(days=(tstart-t0).days)
                    t0 = tstart
                else:
                    continue
            if t1 > tstop:
                if extrap:
                    t0 = t0 - timedelta(days=(t1-tstop).days)
                    t1 = tstop
                else:
                    continue
            ave = p.loc[:,t0:t1].mean(axis=1)
            count = p.loc[:,t0:t1].notnull().sum(axis=1)
            
            # Two ways to calculate uncertainty of mean: propagated uncertainty, and
            # scatter within bin. Use the larger of the two.
            avesigma0 = 1./count * np.sqrt((pe.loc[:,t0:t1]**2).sum(axis=1))
            avesigma1 = 1./count * (p.loc[:,t0:t1].std(axis=1))
            avesigma = np.array([avesigma0,avesigma1]).T.max(axis=1)
            
            # Add to series
            dates.append(pd.to_datetime(date))
            dfrows.append(p2.values - ave.values)
            sigmarows.append(np.sqrt(pe2**2 + avesigma**2))
            nrows.append(count)
            meanrows.append(ave.values)

    dfd = pd.DataFrame(data=dfrows, index=dates, columns=p.index).T
    n   = pd.DataFrame(data=nrows,  index=dates, columns=p.index).T
    dfdsigma = pd.DataFrame(data=sigmarows, index=dates, columns=p.index).T
    dfdmeans = pd.DataFrame(data=meanrows,  index=dates, columns=p.index).T

    # Void entries with few samples
    dfd[n < nthresh] = np.nan
    if verbose:
        idx = (n < nthresh) & (n > 0)
        print '%i/%i entries removed (running mean < %i entries)' % (idx.sum().sum(), (n>0).sum().sum(), nthresh)
        
    # Void sigma where value doesn't exist
    dfdsigma[dfd.isnull()] = np.nan
    
    # Remove mean again for mathematical purposes (e.g., SVD)
    # (And adjust the means calculated above)
    dfdmean = dfd.mean(axis=1)
    for date in dfd.columns:
        dfd.loc[:,date] = dfd.loc[:,date] - dfdmean
        dfdmeans.loc[:,date] += dfdmean
    
    if not sigmavalid:
        dfdsigma.loc[:] = np.nan
    
    return dfd, dfdsigma, dfdmeans    
    
    
    
def remove_outliers(df, stddevs=3, verbose=False):
    '''
    Null all values in the dataframe that differ from the mean by more than 3 stddevs 
    (or some other threshold, as set by stddevs).
    '''
    # Remove outliers
    df = df.copy()
    s = df.stack().std()
    idx = []
    for col in df:
        i = df.loc[:,col].abs() > stddevs*s
        idx.append(i)
        df.loc[i,col] = np.nan
    ii = pd.concat(idx)
    if verbose:
        print '%i/%i outliers removed' % (sum(ii),df.notnull().sum().sum())    
    return df



def remove_partial_days(df, samples=5, verbose=False):
    '''
    Remove dates that only have a few times in them
    '''
    # Remove dates with low counts
    fill = (~df.isnull()).sum()
    idx0 = fill == 0
    idx1 = (~idx0) & (fill < samples)
    idx = idx0 | idx1
    if verbose:
        print '%i/%i days removed (%i are empty and %i have <%i samples)' % (sum(idx), len(idx), sum(idx0), sum(idx1), samples)
    df.loc[:,idx] = np.nan
    df2 = df.dropna(axis=1, how='all')
    return df2

    
def my_fill_na(df, method='0'):
    '''
    Fill null values in the dataframe with something
    method = '0': zeros
    method = 'noise': Gaussian random variable with sample stddev from df
    '''
    s = df.stack().std()
    dffill = df.copy()
    if method=='0':
        dffill.loc[:,:] = 0.0
    elif method=='noise':
        dffill.loc[:,:] = s*random.randn(*np.shape(dffill))
    else:
        raise Exception('Input method="%s" not recognized.' % method)
    df2 = df.combine_first(dffill)
    
    return df2
        
    
    
    
"""
Taylor diagram (Taylor, 2001) implementation.
"""

#__version__ = "Time-stamp: <2017-11-24 18:01:03 ycopin>"
#__author__ = "Yannick Copin <yannick.copin@laposte.net>"

import numpy as NP
import matplotlib.pyplot as PLT


class TaylorDiagram(object):
    """
    Taylor diagram.
    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd, fig=None, rect=111, label='_', srange=(0, 1.5)):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.
        Parameters:
        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        """

        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF

        self.refstd = refstd            # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = NP.concatenate((NP.arange(10)/10., [0.95, 0.99]))
        tlocs = NP.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0]*self.refstd
        self.smax = srange[1]*self.refstd

        ghelper = FA.GridHelperCurveLinear(tr,
                                           extremes=(0, NP.pi/2,  # 1st quadrant
                                                     self.smin, self.smax),
                                           grid_locator1=gl1,
                                           tick_formatter1=tf1)

        if fig is None:
            fig = PLT.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("Standard deviation")

        ax.axis["right"].set_axis_direction("top")   # "Y axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("left")
        #ax.axis["right"].label.set_text("Standard deviation")

        ax.axis["bottom"].set_visible(False)         # Useless

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'k*',
                          ls='', ms=10, label=label)
        t = NP.linspace(0, NP.pi/2)
        r = NP.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]

    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """
        sc = self.ax.scatter(NP.arccos(corrcoef), stddev,
                          *args, **kwargs) # (theta,radius)
        self.samplePoints.append(sc)
        return sc

        #l, = self.ax.plot(NP.arccos(corrcoef), stddev,
        #                  *args, **kwargs)  # (theta,radius)
        #self.samplePoints.append(l)#
        #return l

    def add_grid(self, *args, **kwargs):
        """Add a grid."""

        self.ax.grid(*args, **kwargs)

    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """

        rs, ts = NP.meshgrid(NP.linspace(self.smin, self.smax),
                             NP.linspace(0, NP.pi/2))
        # Compute centered RMS difference
        rms = NP.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*NP.cos(ts))

        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours
    
    
    
class TaylorDiagram2(object):
    """
    Taylor diagram, modified to include 2nd quadrant
    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd, fig=None, rect=111, label='_', srange=(0, 1.5)):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.
        Parameters:
        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        """

        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF

        self.refstd = refstd            # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = NP.concatenate((NP.arange(10)/10., [0.95, 0.99]))
        tlocs = NP.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0]*self.refstd
        self.smax = srange[1]*self.refstd

        ghelper = FA.GridHelperCurveLinear(tr,
                                           extremes=(0, NP.pi,  # 1st quadrant
                                                     self.smin, self.smax),
                                           grid_locator1=gl1,
                                           tick_formatter1=tf1)

        if fig is None:
            fig = PLT.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("Standard deviation")

        ax.axis["right"].set_axis_direction("top")   # "Y axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("left")
        #ax.axis["right"].label.set_text("Standard deviation")

        ax.axis["bottom"].set_visible(False)         # Useless

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'k*',
                          ls='', ms=10, label=label)
        t = NP.linspace(0, NP.pi/2)
        r = NP.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]

    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """
        sc = self.ax.scatter(NP.arccos(corrcoef), stddev,
                          *args, **kwargs) # (theta,radius)
        self.samplePoints.append(sc)
        return sc

        #l, = self.ax.plot(NP.arccos(corrcoef), stddev,
        #                  *args, **kwargs)  # (theta,radius)
        #self.samplePoints.append(l)#
        #return l

    def add_grid(self, *args, **kwargs):
        """Add a grid."""

        self.ax.grid(*args, **kwargs)

    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """

        rs, ts = NP.meshgrid(NP.linspace(self.smin, self.smax),
                             NP.linspace(0, NP.pi/2))
        # Compute centered RMS difference
        rms = NP.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*NP.cos(ts))

        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours
    
    
    
    
    
    
def load_gitm(instr_name, month=None, new=False, hires=False):
    '''
    Load a year (2013) of GITM data from Aaron's files, and return it as a pandas 
    DataFrame.
    
    If month (int) is provided, only load the file from that month.
    
    new: boolean, whether to load new files (07/2018) or old files (2017ish?)
    hires: boolean, only used if new=True. If True, use partial files sent by
        Aaron Aug 20, 2018, which are the same as new=True except with higher
        resolution.
    '''
    
    direc = '/home/bhardin2/GITM_2013/'
    names = ['tdelta','T','v','u','w'] # note swap u,v
    if new:
        direc = '/home/bhardin2/GITM_2013_new/'
        names = ['tdelta','T','u','v','w'] # Aaron seems to have fixed u,v
        if hires:
            direc = '/home/bhardin2/GITM_2013_new2/'

    
    site_gitm = {'minime01':'Cariri',
                 'minime02':'Cajazeiras',
                 'minime03':'Morocco',
                 'minime05':'Urbana',
                 'minime06':'Pari',
                 'minime07':'Kentucky',
                 'minime08':'Ann',
                 }
    year = 2013
    
    if month is None:
        months = arange(1,13)
    else:
        months = [month]
    dfs = []
    for month in months:
        globstr = '%s%s_%i%02i*.*' % (direc, site_gitm[instr_name], year, month)
        fns = glob.glob(globstr)
        if len(fns)==0:
            print('No file found: %s' % (globstr))
            continue
        if len(fns)>1:
            raise Exception('Multiple files found: %s' % (globstr))
        fn = fns[0]
        with open(fn,'r') as f:
            f.readline() #  dummy 
            s = f.readline() # t0
            tv = [int(x) for x in s.split()]
        t0 = datetime(tv[0],tv[1],tv[2],tv[3],tv[4])
        
        df = pd.read_table(fn, skiprows=123, names=names,sep='\s+') # Note swap u,v
        tdelt = np.around(2*df.tdelta.values, 1)/2.0 # round to nearest half-hour (Why is this necessary?)
        df['t'] = t0 + pd.to_timedelta(tdelt, unit='h')
        df.set_index('t', inplace=True)
        df.drop('tdelta',axis=1, inplace=True)
        dfs.append(df)
        
    df = pd.concat(dfs)
    
    if len(months)>1:
        # Fix up indexing issues
        print 'Manually fixing month overlaps...'
        df = df[~df.index.duplicated(keep='first')]
        df = df.reindex(pd.date_range('2013-01-01 00:00','2013-12-31 23:59', freq='0.5H'))
        
    # Metadata
    df.instr_name = instr_name
    df.site_name = site_gitm[instr_name]
    
    return df