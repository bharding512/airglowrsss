'''
Summary
-------
FPIResults contains functions to do both CV and binning/filtering daily/monthly/yearly analysis
This script is meant to contain Level3 routines such as data averaging and plotting.

Included functions are:
    BinDailyData
    BinMonthlyData
    **CreateL2ASCII
    **CreateL2ASCII_Legacy
    CreateMonthlyASCII
    GetModels
    PlotAverages
    PlotClimatology
    PlotClimatologyF107
    PlotLocalTime
    PlotGridMonth
    PlotSpaghetti
    SetBinTime
    SetFilter
        
'''
import matplotlib as _mpl
import matplotlib.dates as _md
import matplotlib.pyplot as _plt
from pyglow import pyglow as _pyglow
import datetime as _dt
import calendar as _cal
import numpy as _np
from scipy import stats as _stats
import pytz as _pytz
import FPIprocessLevel2 as L2
import fpiinfo as _fpiinfo
from multiprocessing import Pool as _pool
from functools import partial as _partial

# I have no idea where these ipython.html.widgets and matplotlib.use are called, 
# but these warings occur in BinMonthlyData... TOBERECTIFIED DJF
import warnings
warnings.filterwarnings("ignore")
# This one removes divided by nan warnings in Filtering OK - DJF
_np.seterr(divide='ignore',invalid='ignore')

def SetBinTime(MIN):
    '''
    Summary:
        Sets Bin time for Binning functions
        
    Input:
        MIN = Amount of time per bin in minutes
        
    History:
        10/6/14 -- Written by DJF (dfisher2@illinois.edu)
    '''
    global b_len
    global btime
    global times
    b_len = 24*60/MIN
    btime = _np.array([(arbdate + _dt.timedelta(minutes=_x)).time() for _x in range(0,b_len*MIN,MIN)])
    times = _np.array([(arbdate + _dt.timedelta(minutes=_x)) for _x in range(0,b_len*MIN,MIN)])

    print 'Bin length is set to %i minutes.'%MIN


def SetFilter(UMIN=-250.,UMAX=250.,VMIN=-250.,VMAX=250.,WMIN=-75.,WMAX=75.,TMIN=600.,TMAX=1400., \
        TERR=50.,WERR=50.,VERBOSE=True):
    '''
    Summary:
        Returns filter limits for averages. This is used because the Level0 automated processing \
            quality flags don't always capture 100%. 

    Inputs:
        UMIN = Zonal wind minimum
        UMAX = Zonal wind maximum
        VMIN = Meridional wind minimum
        VMAX = Meridional wind maximum
        WMIN = Vertical wind minimum
        WMAX = Vertical wind maximum
        TMIN = Temperature minimum
        TMAX = Temperature maximum
        TERR = Temperature uncertainity (error) limit
        WERR = Wind uncertainty (error) limit

    History:
        4/18/16 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # physical limits for data
    global lim
    lim = {'Tmax':TMAX,'Tmin':TMIN,'Te':TERR, \
           'umax':UMAX,'umin':UMIN,'ue':WERR, \
           'vmax':VMAX,'vmin':VMIN,'ve':WERR, \
           'wmax':WMAX,'wmin':WMIN,'we':WERR}
   
    if VERBOSE:
        print 'The following limits are set:'
        print sorted(lim.items())


# Set up default parameters
_mpl.rcParams.update({'font.size': 11})
dirout = '/rdata/airglow/database/L2/plots/'
_n_cores = 16 

# Binning Time
arbdate = _dt.datetime(1970,1,1)
_utc = _pytz.UTC
SetBinTime(30)
SetFilter()
loader = _pyglow.Point(arbdate,0,0,250.) # First run is slow, so get it out of the way now
print 'FPIResults Ready\n'


def FilterData(DATA,QF=1):
    '''
    Summary:
        Returns filtered single input day. This is because automated filtering is not perfect.

    Inputs:
        DATA = The data object
        QF   = Quality Flag max [default = 1] 

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)
        4/18/16 -- Modified by DJF
    '''
    # Set default flag limits
    qualitylimit = QF+1 
    bcnt = 0
    
    # For each direction in DATA
    for r1 in DATA:
        if len(r1.t1)>0:
            #r1 = data[link]
            #ind1 = range(0,len(r1.t1))
            ind1 = range(len(r1.t1))
            alld = len(ind1)

            # Filter using flags 
            if len(r1.flag_T) >= 1:
                ind1 = _np.delete(ind1, _np.where(r1.flag_T[ind1] >= qualitylimit))
                ind1 = _np.delete(ind1, _np.where(r1.flag_wind[ind1] >= qualitylimit))

            # Bad Data Filtering Temps
            if len(r1.T) >= 1:
                ind1 = _np.delete(ind1, _np.where(lim['Tmin'] > r1.T[ind1]))
                ind1 = _np.delete(ind1, _np.where(r1.T[ind1] > lim['Tmax']))
                ind1 = _np.delete(ind1, _np.where(_np.abs(r1.Te[ind1]) > lim['Te']))

            # Bad Data Filtering Winds -both limits, fits, & spikes(cal)
            if len(r1.u) >= 1:
                ind1 = _np.delete(ind1, _np.where(lim['umin'] > r1.u[ind1]))
                ind1 = _np.delete(ind1, _np.where(r1.u[ind1] > lim['umax']))
                ind1 = _np.delete(ind1, _np.where(r1.ue[ind1] > lim['ue']))
            if len(r1.v) >= 1:
                ind1 = _np.delete(ind1, _np.where(lim['vmin'] > r1.v[ind1]))
                ind1 = _np.delete(ind1, _np.where(r1.v[ind1] > lim['vmax']))
                ind1 = _np.delete(ind1, _np.where(r1.ve[ind1] > lim['ve']))
            if len(r1.w) >= 1:
                ind1 = _np.delete(ind1, _np.where(lim['wmin'] > r1.w[ind1]))
                ind1 = _np.delete(ind1, _np.where(r1.w[ind1] > lim['wmax']))
                ind1 = _np.delete(ind1, _np.where(r1.we[ind1] > lim['we']))

            r1.cut(arbdate,arbdate+_dt.timedelta(days=1),ind1)
            bcnt += alld - len(ind1)
    return(bcnt)
    


def WeightedAverage(VAL,STD,CNT=None,AXIS=0,test=False):
    '''
    Summary:
        Returns weighted mean and std of FPI data
        Returns weighted mean/std, and variability/std of FPI data if CNT is given

    Inputs:
        VAL  = data to average
        STD  = standard deviation of data
        CNT  = count of data
        AXIS = axis to average across?

    Outputs:
        WM  = weighted mean
        WE  = weighted std
    Outputs2:
        SS  = sample std (variability)
        SE  = sample std uncertainty
        WSS = weighted sampel std
        AE  = average uncertainty (average error/std)
        P16 = 16th percentile
        P25 = 25th percentile
        P50 = 50th percentile
        P75 = 75th percentile
        P84 = 84th percentile

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)
        9/23/14 -- Redid as separate function w/ corrected errors: DJF
        10/14/15-- Added Percentile to function: DJF
    '''
    # Rotate averaging axis if deisred
    if AXIS==1:
        VAL = VAL.T
        STD = STD.T
        if CNT is not None:
            CNT = CNT.T
    
    # weighted mean and weighted std
    wt = (VAL/VAL)/STD**2
    V1 = _np.nansum(wt,axis=1)
    V2 = _np.nansum(wt**2,axis=1)
    WM = _np.nansum(VAL*wt,axis=1)/V1
    WE = _np.sqrt(_np.nansum(STD**2*wt**2,axis=1)/V1**2)
    if CNT == None:
        # Return wt mean and wt std only
        return(WM,WE)
    elif test:
        # return wt mean and std and monthly variability and percentiles
        SS = _np.sqrt(_np.nansum(_np.subtract(VAL.T,WM)**2,axis=0)/(_np.array(CNT)-1.))
        SE = 2.*WE**4/(_np.array(CNT)-1.)
        WSS= _np.sqrt(_np.nansum(wt.T*_np.subtract(VAL.T,WM)**2,axis=0)/(V1-V2/V1))
        AE = _np.sqrt(_np.nansum(STD**2,axis=1)/CNT) 
        P16 = _np.nanpercentile(VAL,16,axis=1)
        P25 = _np.nanpercentile(VAL,25,axis=1)
        P50 = _np.nanpercentile(VAL,50,axis=1)
        P75 = _np.nanpercentile(VAL,75,axis=1)
        P84 = _np.nanpercentile(VAL,84,axis=1)
        return(WM,WE,SS,SE,WSS,AE,P16,P25,P50,P75,P84)
    else:
        # return wt mean and std and monthly variability/std
        SS = _np.sqrt(_np.nansum(_np.subtract(VAL.T,WM)**2,axis=0)/(_np.array(CNT)-1.))
        SE = 2.*WE**4/(_np.array(CNT)-1.)
        WSS = _np.sqrt(_np.nansum(wt.T*_np.subtract(VAL.T,WM)**2,axis=0)/(V1-V2/V1))
        AE = _np.sqrt(_np.nansum(STD**2,axis=1)/CNT) 
        return(WM,WE,SS,SE)
        #return(WM,WE,SS,SE,WSS,AE)
    
    

def BinDailyData(SITE,YEAR,DOY,SPLIT=False,KP=[0,10],CV=True,QF=1,WIS0=False):
    '''
    Summary:
        Returns filted and binned single data for a single instrument over one night.

    Inputs:
        SITE  = site, e.g. uao
        YEAR  = year, e.g. 2013
        DOY   = doy of year, e.g. 47
        SPLIT = Split look directions in binning [default = False]
        KP    = Filter days by KP [default = [0,10] - all kp]
        CV    = Include CV directions [default = True]
        QF    = Quality Flag Limit Allowed [default = 1]
        WIS0  = Flag to use w=0 in processing [default = False]
        
    Outputs:
        DATA  = Object with winds, Temps, and more

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)

    '''
    
    # Create the YYYYMMDD date format
    dn = _dt.datetime(YEAR,1,1) + _dt.timedelta(days = DOY-1)
    year = dn.strftime('%Y')
    date = dn.strftime('%Y%m%d')
    
    # Load in Day's Data
    if SITE in ['renoir','peru','nation']:
        nets = _fpiinfo.get_network_info(SITE)
        tots = _fpiinfo.get_all_sites_info()
        sites = [x for x in nets.keys() if x in tots]
        lla = nets['mean_location']
        project = SITE.lower()
    else:
        sites = [SITE.lower()]
        siteinfo = _fpiinfo.get_site_info(sites[0])
        project = siteinfo['Network']
        lla = siteinfo['Location']
    data = L2.GetLevel2(project,dn,w_is_0=WIS0)
    bc = FilterData(data,QF)
    
    # Get Empty
    d = _BinnedData(dn,SITE)
    d.key = "Daily"
    d.lla = lla
    d.t = times
    count_len = 500
    uData   = _np.empty((b_len,count_len))*_np.nan
    ueData  = _np.empty((b_len,count_len))*_np.nan
    vData   = _np.empty((b_len,count_len))*_np.nan
    veData  = _np.empty((b_len,count_len))*_np.nan
    u2Data  = _np.empty((b_len,count_len))*_np.nan
    u2eData = _np.empty((b_len,count_len))*_np.nan
    v2Data  = _np.empty((b_len,count_len))*_np.nan
    v2eData = _np.empty((b_len,count_len))*_np.nan
    wData   = _np.empty((b_len,count_len))*_np.nan
    weData  = _np.empty((b_len,count_len))*_np.nan
    TData   = _np.empty((b_len,count_len))*_np.nan
    TeData  = _np.empty((b_len,count_len))*_np.nan
    iData   = _np.empty((b_len,count_len))*_np.nan
    ieData  = _np.empty((b_len,count_len))*_np.nan
    # For harmonic testing
    uefData = _np.empty((b_len,count_len))*_np.nan
    vefData = _np.empty((b_len,count_len))*_np.nan
    wefData = _np.empty((b_len,count_len))*_np.nan
    count   = b_len*[0]
    cvc = 0
    cc  = 0

    # Bin all data  TODO: More Pythonic
    for r1 in data:
        d.moonup = r1.moonup
        # for all sites in list
        for s in sites:
            # TODO: Verify site is good 
            try:
                instr = _fpiinfo.get_instr_at(s,dn)[0]
            except:
                break
            # for each location with that site
            if s in r1.key.lower():
                # Do cv or card counts
                if 'cv_' in r1.key.lower() or 'in_' in r1.key.lower():
                    if CV:
                        cvc += len(r1.t1)
                    else:
                        continue
                else:
                    cc += len(r1.t1)
                # For the total number of exposures (time)
                for zelda in range(len(r1.t1)):
                    # If the kp is in the specified range
                    pt = _pyglow.Point(r1.t1[zelda],0,0,0)
                    kpi = pt.kp
                    if KP[0]<= kpi <=KP[1] or _np.isnan(kpi):
                        bin = int(_np.floor(_np.mod((r1.t1[zelda].astimezone(_utc).replace(tzinfo=None)-arbdate).total_seconds(),60*60*24)/(60*24*60/b_len)))
                        # If the data exists
                        if len(r1.u) > 0:
                            if SPLIT and ('west' in r1.key.lower() or '_2' in r1.key.lower()):
                            # Split E&W
                                u2Data[bin,count[bin]] = r1.u[zelda]
                                u2eData[bin,count[bin]] = r1.ue[zelda]
                            else:
                                uData[bin,count[bin]] = r1.u[zelda]
                                ueData[bin,count[bin]] = r1.ue[zelda]
                        if len(r1.v) > 0:
                            if SPLIT and ('south' in r1.key.lower() or '_2' in r1.key.lower()):
                            # Split N&S
                                v2Data[bin,count[bin]] = r1.v[zelda]
                                v2eData[bin,count[bin]] = r1.ve[zelda]
                            else:
                                vData[bin,count[bin]] = r1.v[zelda]
                                veData[bin,count[bin]] = r1.ve[zelda]
                        #if len(r1.w) > 0 and ('zenith' in r1.key.lower() or 'in' in r1.key.lower()) and r1.parent[0].reference == 'laser':
                        if len(r1.w) > 0 and r1.parent[0].reference == 'laser':
                            wData[bin,count[bin]] = r1.w[zelda]
                            weData[bin,count[bin]] = r1.we[zelda]
                        if len(r1.T) > 0:
                            TData[bin,count[bin]] = r1.T[zelda]
                            TeData[bin,count[bin]] = r1.Te[zelda]
                        if len(r1.i) > 0:
                            iData[bin,count[bin]] = r1.i[zelda]
                            ieData[bin,count[bin]] = r1.ie[zelda]
                        count[bin] += 1
            
    
    ## Weighted Mean of Winds
    # Zonal
    uD,ueD = WeightedAverage(uData,ueData)
    # Meridional
    vD,veD = WeightedAverage(vData,veData)
    # Vert            
    wD,weD = WeightedAverage(wData,weData)
    # Temps
    TD,TeD = WeightedAverage(TData,TeData)
    # Intensities
    iD,ieD = WeightedAverage(iData,ieData)

    if SPLIT:
        # Zonal2 - West
        u2D,u2eD = WeightedAverage(u2Data,u2eData)
        # Meridional2 - South
        v2D,v2eD = WeightedAverage(v2Data,v2eData) 

    # Save Averages
    d.u = uD
    d.ue = ueD
    d.uc = count_len-_np.nansum(_np.isnan(uData),1)
    d.v = vD
    d.ve = veD
    d.vc = count_len-_np.nansum(_np.isnan(vData),1)
    d.w = wD
    d.we = weD
    d.wc = count_len-_np.nansum(_np.isnan(wData),1)
    d.T = TD
    d.Te = TeD
    d.Tc = count_len-_np.nansum(_np.isnan(TData),1)
    d.i = iD
    d.ie = ieD
    d.ic = count_len-_np.nansum(_np.isnan(iData),1)
    if SPLIT:
        d.u2 = u2D
        d.u2e = u2eD
        d.u2c = count_len-_np.nansum(_np.isnan(u2Data),1)
        d.v2 = v2D
        d.v2e = v2eD
        d.v2c = count_len-_np.nansum(_np.isnan(v2Data),1)
    d.cards = cc
    d.cvs = cvc
    d.bads = bc
    d.doabarrelroll()
    return d
    
    
def GetModels(SITELLA,YEAR,DOY,WMODEL,TMODEL='msis',ALT=250.,WEIGHTED=False,QUIET=False,MULTICORE=True):
    '''
    Summary:
        Returns HWM for a single instrument over one night.

    Inputs:
        SITELLA   = site latitude, longitude, altitude
        YEAR      = year, e.g. 2013
        DOY       = doy of year, e.g. 47
        WMODEL    = name of wind model, e.g. 'hwm93'
        TMODEL    = name of temp model [default = 'msis']
        ALT       = altitude of desired profile in km [default = 250]
        WEIGHTED  = flag to intensity weight winds [default = True]
        QUIET     = flag to set Kp=0 and Ap=0.0 [default = False]
        MULTICORE = flag to allow multicore processing [default = True]

    Outputs:
        DATA      = Object with winds, Temps, and more

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)
        4/21/16 -- Mods by DJF (QUIET & REDLINE)
    '''

    # Create the YYYYMMDD date format
    dn = _dt.datetime(YEAR,1,1) + _dt.timedelta(days = DOY-1)
    year = dn.strftime('%Y')
    date = dn.strftime('%Y%m%d')
    
    # allocate Arrays
    uData  = _np.zeros((b_len,1))
    ueData = _np.zeros((b_len,1))
    vData  = _np.zeros((b_len,1))
    veData = _np.zeros((b_len,1))
    wData  = _np.zeros((b_len,1))
    weData = _np.zeros((b_len,1))
    TData  = _np.zeros((b_len,1))
    TeData = _np.zeros((b_len,1))
    iData  = _np.zeros((b_len,1))
    ieData = _np.zeros((b_len,1))
    aData  = _np.zeros((b_len,1))
    
    # Get Empty
    d = _BinnedData(dn,WMODEL)
    d.key = 'DailyModel'
    d.t = times
    d.lla = SITELLA
    if _np.isnan(SITELLA).any():
        print 'Bad LLA'
        return

    # Altitudes
    if WEIGHTED:
        alts = range(100,450,25) # This has been tested and shows it is accurate as 1 km spacing
    else:
        alts = [ALT]
        d.notes+= 'AG peak is set to '+ str(ALT) +'km'

    # Lets see if we can use multicores:
    if WEIGHTED and MULTICORE:
        # Loop through all times
        t_list = [t.replace(year=dn.year,month=dn.month,day=dn.day) for t in times]

        # Prep multicores
        singlemodel = _partial(_MPsinglemodel,ALTS=alts,SITELLA=SITELLA, \
                WMODEL=WMODEL,TMODEL=TMODEL,QUIET=QUIET)
        pool = _pool(processes=_n_cores)
        results = pool.map(singlemodel,t_list)

        # Unwrap results
        for tind in range(b_len):
            uData[tind]  = results[tind][0]
            ueData[tind] = results[tind][1]
            vData[tind]  = results[tind][2]
            veData[tind] = results[tind][3]
            TData[tind]  = results[tind][4]
            TeData[tind] = results[tind][5]
            iData[tind]  = results[tind][6]
            ieData[tind] = results[tind][7]
            aData[tind]  = results[tind][8]

        pool.close()
        pool.join()

    else:
        # Faux model
        #ag6300 = Fmodel('high') #high 
        
        # Fill Data
        for tind,t in enumerate(times):
            for aind,alt in enumerate(alts):
                pt = _pyglow.Point(t.replace(year=dn.year,month=dn.month,day=dn.day),SITELLA[0],SITELLA[1],alt)
                if QUIET:
                    pt.kp = 0.0
                    pt.ap = 0.0
        
                # Intensity
                pt.run_airglow() #FIX THIS AFTER ANALYSIS
                #pt.ag6300 = ag6300[aind]
                iData[tind] += pt.ag6300
                ieData[tind] += 1.

                # Wind
                if WMODEL.lower() == 'hwm93':
                    pt.run_hwm(version=1993)
                elif WMODEL.lower() == 'hwm07':
                    pt.run_hwm(version=2007)
                elif WMODEL.lower() == 'hwm14':
                    pt.run_hwm(version=2014)
                else:
                    print 'Bad Wind Model'
                uData[tind] += pt.u*pt.ag6300
                ueData[tind] += 1.
                vData[tind] += pt.v*pt.ag6300
                veData[tind] += 1.
                #wData[tind] += pt.w*pt.ag6300
                #weData[tind] += 1.
                
                # Temp
                if TMODEL.lower() == 'msis':
                    pt.run_msis()
                    TData[tind] += pt.Tn_msis*pt.ag6300
                elif WMODEL.lower() == 'iri':
                    pt.run_iri()
                    TData[tind] += pt.Tn_iri*pt.ag6300
                else:
                    print 'Bad Temp Model'
                TeData[tind] += 1.
                
                # Altitude
                aData[tind] += alt*pt.ag6300

            # Get profile-weighted average
            uData[tind] = uData[tind]/iData[tind]
            vData[tind] = vData[tind]/iData[tind]
            #wData[tind] = wData[tind]/iData[tind]
            TData[tind] = TData[tind]/iData[tind]
            aData[tind] = aData[tind]/iData[tind]

    # Save Averages
    d.u = uData[:,0]
    d.ue = ueData[:,0]
    d.uc = _np.ones([len(times)])
    d.v = vData[:,0]
    d.ve = veData[:,0]
    d.vc = _np.ones([len(times)])
    d.w = wData[:,0]
    d.we = weData[:,0]
    d.wc = _np.zeros([len(times)])
    d.T = TData[:,0]
    d.Te = TeData[:,0]
    d.Tc = _np.ones([len(times)])
    d.i = iData[:,0]
    d.ie = ieData[:,0]
    d.ic = _np.ones([len(times)])
    d.alts = aData[:,0]
    d.doabarrelroll()
    
    return d
        
        
def _MPsinglemodel(T,ALTS,SITELLA,WMODEL,TMODEL,QUIET):
    '''
    Summary:
        Returns Model results for a single time. (Multicore Code)

    Inputs:
        T = datetime to use
        ALTS = list of altitudes to use
        SITELLA = site latitude, longitude, altitude
        WMODEL = name of wind model, e.g. 'hwm93'
        TMODEL = name of temp model [default = 'msis']
        QUIET = flag to set Kp=0 and Ap=Ap_daily [default = False]

    Outputs:
        uData,ueData,vData,veData,TData,TeData,iData,ieData,aData
            = u,v,T,i,and alt values and uncertainty

    History:
        5/03/16 -- Created by DJF (dfisher2@illinois.edu)
    '''
    # allocate incrementers
    uData  = 0
    ueData = 0
    vData  = 0
    veData = 0
    wData  = 0
    weData = 0
    TData  = 0
    TeData = 0
    iData  = 0
    ieData = 0
    aData  = 0

    for aind,alt in enumerate(ALTS):
        pt = _pyglow.Point(T,SITELLA[0],SITELLA[1],alt)
        if QUIET:
            pt.kp = 0.0
            pt.ap = 0.0
    
        # Intensity
        pt.run_airglow() #FIX THIS AFTER ANALYSIS
        #pt.ag6300 = ag6300[aind]
        iData += pt.ag6300
        ieData += 1.

        # Wind
        if WMODEL.lower() == 'hwm93':
            pt.run_hwm(version=1993)
        elif WMODEL.lower() == 'hwm07':
            pt.run_hwm(version=2007)
        elif WMODEL.lower() == 'hwm14':
            pt.run_hwm(version=2014)
        else:
            print 'Bad Wind Model'
        uData  += pt.u*pt.ag6300
        ueData += 1.
        vData  += pt.v*pt.ag6300
        veData += 1.
        #wData  += pt.w*pt.ag6300
        #weData += 1.
        
        # Temp
        if TMODEL.lower() == 'msis':
            pt.run_msis()
            TData += pt.Tn_msis*pt.ag6300
        elif WMODEL.lower() == 'iri':
            pt.run_iri()
            TData += pt.Tn_iri*pt.ag6300
        else:
            print 'Bad Temp Model'
        TeData += 1.
        
        # alt
        aData += alt*pt.ag6300

    # Get profile-weighted average
    uData = uData/iData
    vData = vData/iData
    #wData = wData/iData
    TData = TData/iData
    aData = aData/iData

    return(uData,ueData,vData,veData,TData,TeData,iData,ieData,aData)
    
    

def BinMonthlyData(SITE,YEAR,MONTH,DLIST=[],YLIST=[],SPLIT=False,KP=[0,10],CV=True,QF=1, \
        SITELLA=[],TMODEL='msis',ALT=250.,WEIGHTED=False,QUIET=False,VERBOSE=True):
    '''
    Summary:
        Returns filted and binned data/models over one month.

    Inputs:
        SITE     = site of interest, e.g. 'UAO' or model, e.g. 'hwm14'
        YEAR     = year, e.g. 2013
        MONTH    = month of year, e.g. 2 (February)
        DLIST    = list of doys in year  [default = [] - all doys in MONTH,YEAR used]
                   or [doy,spread] for storms
        YLIST    = list of years in year [default = [] - only YEAR used]
        SPLIT    = [data param] split look directions in binning [default = False]
        KP       = [data param] limits of kp for filtering days [default = [0,10] - all kp used]
        CV       = [data param] use CV modes [default = True]
        QF       = [data param] Quality Flag Limit Allowed [default = 1]
        SITELLA  = [model param] site location array for model run [lat,lon,alt_km]
        TMODEL   = [model param] name of temp model [default = 'msis']
        ALT      = [model param] altitude of desired profile in km [default = 250]
        WEIGHTED = [model param] flag to intensity weight winds [default = False]
        QUIET    = [model param] flag to set Kp=0 and Ap=0.0 [default = False]
        VERBOSE  = Print information to stdout [default = True]

    Outputs:
        DATA = dictionary of data whose keys are Zonal, Meridional, or Temp 
               Temp contains Temp and Temp_Error (averaged from all directions)
               Zonal/Meridional contains Wind and Wind_Error

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)
        10/14/15-- Added 2 unit DLIST input, (doy,spread) to allow doy+/-spread averages.
        12/02/15-- Added VERBOSE optional input (BJH)

    '''
    # Define Outputs...
    dn = _dt.datetime(YEAR,MONTH,1)
    mon = _cal.month_name[MONTH]

    # Set Output Variable
    dimset = 3  # Require x days in each month for an variability set
    d = _BinnedData(dn,SITE)
    d.t = times
    d.key = '{0:%B} {0:%Y}'.format(dn)
    
    # Get Empty 
    dim = 31*5*10
    uData   = _np.empty((b_len,dim))*_np.nan
    ueData  = _np.empty((b_len,dim))*_np.nan
    uCount  = _np.zeros((b_len))
    vData   = _np.empty((b_len,dim))*_np.nan
    veData  = _np.empty((b_len,dim))*_np.nan
    vCount  = _np.zeros((b_len))
    u2Data  = _np.empty((b_len,dim))*_np.nan
    u2eData = _np.empty((b_len,dim))*_np.nan
    u2Count = _np.zeros((b_len))
    v2Data  = _np.empty((b_len,dim))*_np.nan
    v2eData = _np.empty((b_len,dim))*_np.nan
    v2Count = _np.zeros((b_len))
    wData   = _np.empty((b_len,dim))*_np.nan
    weData  = _np.empty((b_len,dim))*_np.nan
    wCount  = _np.zeros((b_len))
    TData   = _np.empty((b_len,dim))*_np.nan
    TeData  = _np.empty((b_len,dim))*_np.nan
    TCount  = _np.zeros((b_len))
    iData   = _np.empty((b_len,dim))*_np.nan
    ieData  = _np.empty((b_len,dim))*_np.nan
    iCount  = _np.zeros((b_len))
    F107    = _np.empty((dim))*_np.nan
    mflag = False
    oscar = []
    lat = []
    lon = []
    alt = []
    cards = 0
    cvs = 0
    bads = 0
    
    # Collect month's Data:
    if DLIST == []:
        doystart = (dn - _dt.datetime(YEAR,1,1)).days+1
        doyend = doystart + _cal.monthrange(YEAR,MONTH)[1]
        dl = range(doystart,doyend)
        yl = list(_np.array(dl)/_np.array(dl)*YEAR)
    elif len(DLIST) == 2 and len(YLIST) == 0:
        # Use listed DOY +/- SPREAD disregarding month.
        dl = []; yl = []
        d.key = str(DLIST[0])+'-'+str(YEAR)+' +/-'+str(DLIST[1])
        dn0 = _dt.datetime(YEAR,1,1)+_dt.timedelta(days=DLIST[0]-1)
        for x in range(-DLIST[1],DLIST[1]+1):
            dn = dn0 + _dt.timedelta(days=x)
            dl.append((dn-_dt.datetime(dn.year,1,1)).days+1)
            yl.append(dn.year)
    else:
        # Use listed DOYS and make sure doys actually exist in the month
        d.key = "Listed days in {0:%B}".format(dn)
        dl = []; yl = []
        for doy,yr in zip(DLIST,YLIST):
            dl.append(doy)
            yl.append(yr)
        if not(dl):
            print 'Doy List contains no days in desired month.'

    # Make a site list
    if SITE in ['renoir','peru','nation']:
        nets = _fpiinfo.get_network_info(SITE)
        tots = _fpiinfo.get_all_sites_info()
        sites = [x for x in nets.keys() if x in tots]
    else:
        sites = [SITE]
    if 'hwm' in SITE:
        mflag = True   

    # Get inputs for multicores
    doy_arg = [a for s in sites for a in dl]
    yr_arg = [y for s in sites for y in yl]
    s_arg = [s for s in sites for a in dl]

    # Process all days using multicores
    singleday = _partial(_MPsingleday,SPLIT=SPLIT,KP=KP,CV=CV,QF=QF, \
            SITELLA=SITELLA,TMODEL=TMODEL,ALT=ALT,WEIGHTED=WEIGHTED,QUIET=QUIET)
    pool = _pool(processes=_n_cores)
    results = pool.map(singleday,zip(s_arg,yr_arg,doy_arg))

    # Unwrap results
    for ind,doy in enumerate(doy_arg):
        DD = results[ind][0]
        F107[ind] = results[ind][1]

        # for printing stats
        if not(mflag):
            cards += DD.cards
            cvs += DD.cvs
            bads += DD.bads

        # Undo shift for easy averaging
        DD.doabarrelroll()
        
        # Count days total used
        if sum(_np.isfinite(DD.T)):
            oscar.append(doy)
        # Average Location 
        if len(DD.lla) == 3:
            lat.append(DD.lla[0])
            lon.append(DD.lla[1])
            alt.append(DD.lla[2])
        # Add data
        if len(DD.u) > 0:
            uData[:,ind] = DD.u
            ueData[:,ind] = DD.ue
            uCount += DD.uc
        if len(DD.v) > 0:
            vData[:,ind] = DD.v
            veData[:,ind] = DD.ve
            vCount += DD.vc
        if len(DD.w) > 0:
            wData[:,ind] = DD.w
            weData[:,ind] = DD.we
            wCount += DD.wc
        if len(DD.T) > 0:
            TData[:,ind] = DD.T
            TeData[:,ind] = DD.Te
            TCount += DD.Tc
        if len(DD.i) > 0:
            iData[:,ind] = DD.i
            ieData[:,ind] = DD.ie
            iCount += DD.ic
        if SPLIT and not(mflag) and len(DD.u2) > 0:
            u2Data[:,ind] = DD.u2
            u2eData[:,ind] = DD.u2e
            u2Count += DD.u2c
        if SPLIT and not(mflag) and len(DD.v2) > 0:
            v2Data[:,ind] = DD.v2
            v2eData[:,ind] = DD.v2e
            v2Count += DD.v2c
    pool.close()
    pool.join()
        
    # Get count of days used in each bin
    uDays = dim - sum(_np.isnan(uData.T))
    uDays = [_np.nan if x<dimset else x for x in uDays]
    vDays = dim - sum(_np.isnan(vData.T))
    vDays = [_np.nan if x<dimset else x for x in vDays]
    wDays = dim - sum(_np.isnan(wData.T))
    wDays = [_np.nan if x<dimset else x for x in wDays]
    TDays = dim - sum(_np.isnan(TData.T))
    TDays = [_np.nan if x<dimset else x for x in TDays]
    iDays = dim - sum(_np.isnan(iData.T))
    iDays = [_np.nan if x<dimset else x for x in iDays]
    if SPLIT:
        u2Days = dim - sum(_np.isnan(u2Data.T))
        u2Days = [_np.nan if x<dimset else x for x in u2Days]
        v2Days = dim - sum(_np.isnan(v2Data.T))
        v2Days = [_np.nan if x<dimset else x for x in v2Days]
    
    # Zonal
    uD,uDe,uV,uVe,uV2,uU,u16,u25,u50,u75,u84 = WeightedAverage(uData,ueData,uDays,test=True)
    # Meridional
    vD,vDe,vV,vVe,vV2,vU,v16,v25,v50,v75,v84 = WeightedAverage(vData,veData,vDays,test=True)
    # Vert            
    wD,wDe,wV,wVe,wV2,wU,w16,w25,w50,w75,w84 = WeightedAverage(wData,weData,wDays,test=True)
    #  Temps
    TD,TDe,TV,TVe,TV2,TU,T16,T25,T50,T75,T84 = WeightedAverage(TData,TeData,TDays,test=True)
    # Intensity
    iD,iDe,iV,iVe,iV2,iU,i16,i25,i50,i75,i84 = WeightedAverage(iData,ieData,iDays,test=True)
    if SPLIT and not(mflag):
        # Zonal2 - West
        u2D,u2De,u2V,u2Ve,u2V2,u2U,u216,u225,u250,u275,u284 = WeightedAverage(u2Data,u2eData,u2Days,test=True)
        # Meridional2 - South
        v2D,v2De,v2V,v2Ve,v2V2,v2U,v216,v225,v250,v275,v284 = WeightedAverage(v2Data,v2eData,v2Days,test=True)
    # F107
    if mflag:
        d.f107 = _np.nanmean(F107)
    else:
        d.f107 = _np.nansum(F107)/(cards+cvs)
    
    # Save Averages
    d.lla = _np.array([_np.nanmean(lat),_np.nanmean(lon),_np.nanmean(alt)])
    d.u  = uD
    d.ue = uDe
    d.uv = uV
    d.uve= uVe
    d.uv2= uV2
    d.uu = uU
    d.uc = uCount
    d.ud = uDays
    d.u16= u16
    d.u25= u25
    d.u50= u50
    d.u75= u75
    d.u84= u84
    d.v  = vD
    d.ve = vDe
    d.vv = vV
    d.vve= vVe
    d.vv2= vV2
    d.vu = vU
    d.vc = vCount
    d.vd = vDays
    d.v16= v16
    d.v25= v25
    d.v50= v50
    d.v75= v75
    d.v84= v84
    d.w  = wD
    d.we = wDe
    d.wv = wV
    d.wve= wVe
    d.wv2= wV2
    d.wu = wU
    d.wc = wCount
    d.wd = wDays
    d.w16= w16
    d.w25= w25
    d.w50= w50
    d.w75= w75
    d.w84= w84
    d.T  = TD
    d.Te = TDe
    d.Tv = TV
    d.Tve= TVe
    d.Tu = TU
    d.Tc = TCount
    d.Td = TDays
    d.Tv2= TV2
    d.T16= T16
    d.T25= T25
    d.T50= T50
    d.T75= T75
    d.T84= T84
    d.i  = iD
    d.ie = iDe
    d.iv = iV
    d.ive= iVe
    d.iv2= iV2
    d.iu = iU
    d.ic = iCount
    d.id = iDays
    d.i16= i16
    d.i25= i25
    d.i50= i50
    d.i75= i75
    d.i84= i84
    
    if SPLIT and not(mflag):
        d.u2  = u2D
        d.u2e = u2De
        d.u2v = u2V
        d.u2ve= u2Ve
        d.u2v2= u2V2
        d.u2u = u2U
        d.u2c = u2Count
        d.u2d = u2Days
        d.u216= u216
        d.u225= u225
        d.u250= u250
        d.u275= u275
        d.u284= u284
        d.v2  = v2D
        d.v2e = v2De
        d.v2v = v2V
        d.v2ve= v2Ve
        d.v2v2= v2V2
        d.v2u = v2U
        d.v2c = v2Count
        d.v2d = v2Days
        d.v216= v216
        d.v225= v225
        d.v250= v250
        d.v275= v275
        d.v284= v284
    d.cards = cards
    d.cvs = cvs
    d.daysused = _np.unique(oscar)
    d.doabarrelroll()

    #print '%s %s %04d'%(SITE, mon,YEAR),'Days used:',len(d.daysused),' Ave pts/time: %02.2f'%(_np.nanmin(_np.array([_np.ma.masked_array(ucount,_np.isnan(ucount)).mean(), _np.ma.masked_array(vcount,_np.isnan(vcount)).mean(), _np.ma.masked_array(Tcount,_np.isnan(Tcount)).mean()])))
    if VERBOSE:
        print '%s %s %04d'%(SITE, mon,YEAR),' Days used:',len(d.daysused)
        if SPLIT and not(mflag):
            print 'Ave pts/bin: u-%02.2f  u2-%02.2f  v-%02.2f  v2-%02.2f  w-%02.2f  T-%02.2f'%(_np.nanmean(uCount),_np.nanmean(u2Count),_np.nanmean(vCount),_np.nanmean(v2Count),_np.nanmean(wCount),_np.nanmean(TCount))
            print 'Ave day/bin: u-%02.2f  u2-%02.2f  v-%02.2f  v2-%02.2f  w-%02.2f  T-%02.2f'%(_np.nanmean(uDays),_np.nanmean(u2Days),_np.nanmean(vDays),_np.nanmean(v2Days),_np.nanmean(wDays),_np.nanmean(TDays))

        else:
            print 'Ave pts/bin: u-%02.2f  v-%02.2f  w-%02.2f  T-%02.2f'%(_np.nanmean(uCount),_np.nanmean(vCount),_np.nanmean(wCount),_np.nanmean(TCount))
            print 'Ave day/bin: u-%02.2f  v-%02.2f  w-%02.2f  T-%02.2f'%(_np.nanmean(uDays),_np.nanmean(vDays),_np.nanmean(wDays),_np.nanmean(TDays))
        print 'Cards:',cards,' CVs:',cvs,' %%CV: %02.2f'%(100.*cvs/(cvs+cards+.000001))
        print 'Bads:',bads, ' %%Good: %02.2f'%(100.*(cards+cvs)/(cards+cvs+bads+.000001))
        print 'Days:',d.daysused,'\n'

    return d


def _MPsingleday(SITE_YEAR_DOY,SPLIT,KP,CV,QF,SITELLA,TMODEL,ALT,WEIGHTED,QUIET):
    '''
    Summary:
        Returns Daily results for a single time. (Multicore Code)

    Inputs:
        SITE_YEAR_DOY = tuple of the following variable inputs:
            SITE      = sites of interest, e.g. 'UAO'
            YEAR      = year, e.g. 2013
            DOY       = day of year, e.g. 47 
        SPLIT         = [data param] Split look directions in binning [default = False]
        KP            = [data param] Filter days by KP [default = [0,10] - all kp]
        CV            = [data param] use CV modes [default = True]
        QF            = [data param] Quality Flag Limit Allowed [default = 1]
        SITELLA       = [model param] site location array for model run [lat,lon,alt_km]
        TMODEL        = [model param] name of temp model [default = 'msis']
        ALT           = [model param] altitude of desired profile in km [default = 250]
        WEIGHTED      = [model param] flag to intensity weight winds [default = True]
        QUIET         = [model param] flag to set Kp=0 and Ap=0.0 [default = False]

    Outputs:
        DD = Daily Data binned object
        F107 = Weighted average F107 value

    History:
        5/03/16 -- Written by DJF (dfisher2@illinois.edu)
    '''
    
    site = SITE_YEAR_DOY[0]
    year = SITE_YEAR_DOY[1]
    doy  = SITE_YEAR_DOY[2]

    if 'hwm' in site:
        if len(SITELLA) == 0:
            raise ValueError('Need location for model')
        DD = GetModels(SITELLA,year,doy,site,TMODEL=TMODEL,ALT=ALT, \
                WEIGHTED=WEIGHTED,QUIET=QUIET,MULTICORE=False)

        # get F107 weighted at midnight of data (assume constant for night)
        point = _pyglow.Point(DD.dn,0,0,250)
        F107  = (point.f107 + point.f107a)/2.

    else:
        DD = BinDailyData(site,year,doy,SPLIT=SPLIT,KP=KP,CV=CV,QF=QF)

        # get F107 weighted at midnight of data (assume constant for night)
        point = _pyglow.Point(DD.dn,0,0,250)
        F107 = (point.f107 + point.f107a)/2.*(DD.cards+DD.cvs)

    return(DD,F107)


def PlotClimatology(SITE,YEAR,MONTHSTART=1,NMONTHS=12,SPLIT=False,KP=[0,10],UT=True,QF=1,VERBOSE=True):
    '''
    Summary:
        Plots monthly averages in a 2x6 month plot, binning by year
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTHSTART = month to start yearly climatology [default = 1 - Jan]
        NMONTHS = number of months desired from MONTHSTART [default = 12]
        SPLIT = Split look directions in binning [default = False]
        KP = Filter days by KP [default = [0,10] - all kp]
        UT = Plot in UT or SLT [default = True]
        QF = Quality Flags [default = 1]
        VERBOSE = Print monthly average info [default = True]

    Outputs:

    History:
        6/13/13 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    
    # Set up Figures
    axlim = 0
    _mpl.rcParams.update({'font.size':8})
    _mpl.rcParams['savefig.dpi'] = 300
    _mpl.rcParams['figure.figsize'] = (6,4)
    _plt.figure(0)
    axz={};fz,((axz[0],axz[1]),(axz[2],axz[3]),(axz[4],axz[5]), \
               (axz[6],axz[7]),(axz[8],axz[9]),(axz[10],axz[11]))  \
               = _plt.subplots(6,2, sharex=True, sharey=False)
    _plt.figure(1)
    axm={};fm,((axm[0],axm[1]),(axm[2],axm[3]),(axm[4],axm[5]), \
               (axm[6],axm[7]),(axm[8],axm[9]),(axm[10],axm[11])) \
               = _plt.subplots(6,2, sharex=True, sharey=False)
    _plt.figure(2)
    axv={};fv,((axv[0],axv[1]),(axv[2],axv[3]),(axv[4],axv[5]), \
               (axv[6],axv[7]),(axv[8],axv[9]),(axv[10],axv[11])) \
               = _plt.subplots(6,2, sharex=True, sharey=False)
    _plt.figure(3)
    axt={};ft,((axt[0],axt[1]),(axt[2],axt[3]),(axt[4],axt[5]), \
               (axt[6],axt[7]),(axt[8],axt[9]),(axt[10],axt[11])) \
               = _plt.subplots(6,2, sharex=True, sharey=False)
    csize = 0.0001
    lwide = 2
    msize = 3
    
    # Get UT offset if SLT needed
    ut_offset = 0
    if not(UT):
        if SITE == 'renoir':
            ut_offset = - _dt.timedelta(hours=2,minutes=12).total_seconds()
        else:
            try:
                ut_offset = _fpiinfo.get_site_info(SITE)['Location'][1]/360.*(24*60*60)
            except:
                print 'No UT/LT conversion... Plots lie! Actually in UT'

    # Get Monthly Average for basis
    dn = _dt.datetime(YEAR,MONTHSTART,1)
    if SPLIT and NMONTHS > 12:
        print "Only 1 year allowed in SPLIT Mode"
        NMONTHS = 12
    for nsp,mon in enumerate(range(0,NMONTHS)):
        nsp = ((nsp*2)%24 -11*((nsp*2)%24>11))%12  # fill left column then right.
        month = MONTHSTART + mon
        year = YEAR
        if month > 12:
            year = YEAR + (month-1)/12
            month = month - 12*((month-1)/12)
        #print nsp,year,month

        MD = BinMonthlyData(SITE,year,month,SPLIT=SPLIT,KP=KP,QF=QF,VERBOSE=VERBOSE)
        tlim = [MD.t[len(MD.t)/5],MD.t[-len(MD.t)/5]]
        MD.t  = MD.t + _dt.timedelta(seconds=ut_offset)
#	tlim = [MD.t[len(MD.t)/5],MD.t[-len(MD.t)/5]]
        MD.t2 = MD.t + _dt.timedelta(minutes=3)
        gflag = 'on'
 
        try:
	    tstart = _np.where(_np.isnan(MD.T[:-1])*_np.isfinite(MD.T[1:]))[0][0]-1
    	    tend   = _np.where(_np.isfinite(MD.T[:-1])*_np.isnan(MD.T[1:]))[0][0]+1
        except:
            continue
       
        ## Subplot data
        tits = ['Zonal Wind','Meridional Wind','Veritcal Wind','Temperature']
        xlabs = ['Wind Speed [m/s]','Wind Speed [m/s]','Wind Speed [m/s]','Temperature [K]']
        title = ['Clima-Z','Clima-M','Clima-V','Clima-T']
        #years  = [2009, 2010, 2011, 2012, 2013, 2014, 2015]
        colors = ['r.-','b.-','g.-','m.-','c.-','y.-','k.-','rs-','bs-','gs-','ms-','cs-','ys-','ks-']
        colors2= ['ro-','bo-','go-','mo-']
        markerwide = 0
        
        # Zonal
        # ylabel("%s %4.0f" % (_cal.month_abbr[month], year))
        axz[nsp].set_ylabel("%s" % (_cal.month_abbr[month]),labelpad=2)
        axz[nsp].set_ylim(-50, 150)
	if UT:
	    axz[nsp].set_xlim(tlim)
	else:
#        print MD.t[tstart], MD.t[tend]
#        axz[nsp].set_xlim(MD.t[tstart],MD.t[tend])
            axz[nsp].set_xlim([_dt.datetime(1969,12,31,17,00,00),_dt.datetime(1970,1,1,7,0,0)])
        axz[nsp].set_yticks([0, 50, 100])  #[-50, 0, 50, 100]
        if SPLIT:
            axz[nsp].errorbar(MD.t,MD.u,yerr=MD.uv,fmt=colors2[2],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize,label='East') 
            axz[nsp].errorbar(MD.t2,MD.u2,yerr=MD.u2v,fmt=colors2[3],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize,label='West')
        else:
            axz[nsp].errorbar(MD.t+_dt.timedelta(minutes=3*(year-YEAR)),MD.u,yerr=MD.uv,fmt=colors[year-YEAR],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize,label='Data')
        axz[nsp].grid(gflag)
        axz[nsp].plot([MD.t[0],MD.t[-1]],[0,0],'k--',label=None)

        # Meridional
        axm[nsp].set_ylabel("%s" % (_cal.month_abbr[month]),labelpad=2) #"%s %4.0f" % (_cal.month_abbr[month], year)
        axm[nsp].set_ylim(-100, 100)
#        axm[nsp].set_xlim(tlim)
        if UT:
            axm[nsp].set_xlim(tlim)
        else:
#        print MD.t[tstart], MD.t[tend]
#        axz[nsp].set_xlim(MD.t[tstart],MD.t[tend])
            axm[nsp].set_xlim([_dt.datetime(1969,12,31,17,00,00),_dt.datetime(1970,1,1,7,0,0)])
        axm[nsp].set_yticks([-50, 0, 50])
        if SPLIT:
            axm[nsp].errorbar(MD.t,MD.v,yerr=MD.vv,fmt=colors2[0],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize,label='North') 
            axm[nsp].errorbar(MD.t2,MD.v2,yerr=MD.v2v,fmt=colors2[1],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize,label='South')
        else:
            axm[nsp].errorbar(MD.t+_dt.timedelta(minutes=3*(year-YEAR)),MD.v,yerr=MD.vv,fmt=colors[year-YEAR],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize,label='Data')
        axm[nsp].grid(gflag)
        axm[nsp].plot([MD.t[0],MD.t[-1]],[0,0],'k--',label=None)
        
        # Vertical
        axv[nsp].set_ylabel("%s" % (_cal.month_abbr[month]),labelpad=2)
        axv[nsp].set_ylim(-75, 75)
#        axv[nsp].set_xlim(tlim)
        if UT:
            axv[nsp].set_xlim(tlim)
        else:
#        print MD.t[tstart], MD.t[tend]
#        axz[nsp].set_xlim(MD.t[tstart],MD.t[tend])
            axv[nsp].set_xlim([_dt.datetime(1969,12,31,17,00,00),_dt.datetime(1970,1,1,7,0,0)])
        axv[nsp].set_yticks([-25, 0, 25])
        axv[nsp].errorbar(MD.t+_dt.timedelta(minutes=3*(year-YEAR)),MD.w,yerr=MD.wv,fmt=colors[year-YEAR],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize,label='Data')
        axv[nsp].grid(gflag)
        axv[nsp].plot([MD.t[0],MD.t[-1]],[0,0],'k--',label=None)
        
        # Temps
        axt[nsp].set_ylabel("%s" % (_cal.month_abbr[month]),labelpad=2)
        axt[nsp].set_ylim(600, 1200)
#        axt[nsp].set_xlim(tlim)
        if UT:
            axt[nsp].set_xlim(tlim)
        else:
#        print MD.t[tstart], MD.t[tend]
#        axz[nsp].set_xlim(MD.t[tstart],MD.t[tend])
            axt[nsp].set_xlim([_dt.datetime(1969,12,31,17,00,00),_dt.datetime(1970,1,1,7,0,0)])
        axt[nsp].set_yticks([700, 800, 900, 1000, 1100])
        axt[nsp].grid(gflag)
        axt[nsp].errorbar(MD.t+_dt.timedelta(minutes=3*(year-YEAR)),MD.T,yerr=MD.Tv,fmt=colors[year-YEAR],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize)
    
    # Finalize Plots
        axz[0].xaxis.set_major_formatter(_md.DateFormatter('%H')) # :%M
        axm[0].xaxis.set_major_formatter(_md.DateFormatter('%H')) 
        axv[0].xaxis.set_major_formatter(_md.DateFormatter('%H')) 
        axt[0].xaxis.set_major_formatter(_md.DateFormatter('%H')) 
    for i,f in enumerate([fz,fm,fv,ft]):
        f.suptitle("Average %s %s" % (SITE.upper(), tits[i]), fontsize=16, fontweight='bold')
        f.subplots_adjust(hspace = 0.001)
        f.subplots_adjust(wspace = 0.320)
        f.subplots_adjust(left   = 0.140)
        if UT:
            f.text(0.5,0.05,'Hour [UTC]',ha='center',va='center', fontsize=13)
        else:
            f.text(0.5,0.05,'Hour [SLT]',ha='center',va='center', fontsize=13)
        f.text(0.05,0.5,xlabs[i],ha='center',va='center',rotation='vertical', fontsize=14)
        #fig.savefig("%s%s-%4.0f-%s.eps" % (dirout,title[i],YEAR,SITE))
        

def PlotClimatologyF107(SITE,DNSTART,DNEND,SPLIT=False,KP=[0,10],UT=True,QF=1,F_VAL=[50,100,150,200,250],VERBOSE=True):
    '''
    Summary:
        Plots monthly averages in a 2x6 month plot binning by average F10.7
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        DNSTART = Datetime to start
        DNEND = Datetime to end
        SPLIT = Split look directions in binning [default = False]
        KP = Filter days by KP [default = [0,10] - all kp]
        UT = Plot in UT [default = True]
        QF = Quality Flags [default = 1]
        F_VALS = List of F10.7 Cutoff values
        VERBOSE = Print monthly average information [default = True]

    Outputs:

    History:
        6/13/13 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    _mpl.rcParams.update({'font.size':8})
    _mpl.rcParams['savefig.dpi'] = 200
    _mpl.rcParams['figure.figsize'] = (6,4)   
    
    # Set up Figures
    if len(F_VAL) == 2:
        colors= ['#0000FF','#FF0000']
    elif len(F_VAL) == 3:
        colors= ['#0000FF','#008000','#FF0000']
    else:
        colors= ['#000000','#AA0000','#FF0000','#FFAA00','#FFFF00','#AAFF00','#00FF00','#00FFAA','#00FFFF','#00AAFF','#0000FF','#0000AA','#000000','#AA0000','#FF0000','#FFAA00','#FFFF00','#AAFF00','#00FF00','#00FFAA','#00FFFF','#00AAFF','#0000FF','#0000AA'][::-1]
    markerwide = 0.0001
    linewide = 1
    ms = 3
    
    _plt.figure(0)
    axz={};fz,((axz[0],axz[1]),(axz[2],axz[3]),(axz[4],axz[5]), \
               (axz[6],axz[7]),(axz[8],axz[9]),(axz[10],axz[11]))  = _plt.subplots(6,2, sharex=True, sharey=False)
    _plt.figure(1)
    axm={};fm,((axm[0],axm[1]),(axm[2],axm[3]),(axm[4],axm[5]), \
               (axm[6],axm[7]),(axm[8],axm[9]),(axm[10],axm[11]))  = _plt.subplots(6,2, sharex=True, sharey=False)
    _plt.figure(2)
    axv={};fv,((axv[0],axv[1]),(axv[2],axv[3]),(axv[4],axv[5]), \
               (axv[6],axv[7]),(axv[8],axv[9]),(axv[10],axv[11]))  = _plt.subplots(6,2, sharex=True, sharey=False)
    _plt.figure(3)
    axt={};ft,((axt[0],axt[1]),(axt[2],axt[3]),(axt[4],axt[5]), \
               (axt[6],axt[7]),(axt[8],axt[9]),(axt[10],axt[11]))  = _plt.subplots(6,2, sharex=True, sharey=False)
    
    # First Find F10.7 Range Bins
    f_doy_index = []
    f_yr_index = []
    for a in F_VAL:
        f_doy_index.append([])
        f_yr_index.append([])

    for step in range((DNEND-DNSTART).days):
        dn = DNSTART+_dt.timedelta(days=step)
        
        pt = _pyglow.Point(dn,0,0,250)
        try:
            Fb = (pt.f107a+pt.f107)/2.
            fbin = [k for k,i in enumerate(F_VAL) if Fb>=i and Fb-(F_VAL[1]-F_VAL[0])<i][0]
        except:
            fbin = [k for k,i in enumerate(F_VAL) if pt.f107a>=i and pt.f107a-(F_VAL[1]-F_VAL[0])<i][0]
        f_doy_index[fbin].append(int((dn-_dt.datetime(dn.year,1,1)).total_seconds()/(60*60*24)+1))
        f_yr_index[fbin].append(dn.year)
    
    # Get Monthly Average for basis
    for nsp,mon in enumerate(range(1,13)):
        nsp = ((nsp*2)%24 - 11*((nsp*2)%24>11))%12
        for k,(f_doy,f_yr) in enumerate(zip(f_doy_index,f_yr_index)):
            print mon,'-',F_VAL[k]
            MD = BinMonthlyData(SITE,arbdate.year,mon,DLIST=f_doy,YLIST=f_yr,SPLIT=SPLIT,KP=KP,QF=QF,VERBOSE=VERBOSE)
            #print '|_ F107b:',MD.f107,'\n\n'
            MD.t2 = MD.t + _dt.timedelta(minutes=3)
            '''
            # Remove certain garbage data points:
            if SITE == 'renoir' and mon == 8:
                MD.w = MD.w*_np.nan
            if SITE == 'renoir' and mon == 12*5+7:
                MD.T[b_len/2:b_len/2+10] = _np.ones(10)*_np.nan
                MD.v[b_len/2:b_len/2+10] = _np.ones(10)*_np.nan
                MD.u[b_len/2:b_len/2+10] = _np.ones(10)*_np.nan
                MD.w[b_len/2:b_len/2+10] = _np.ones(10)*_np.nan
            if SITE == 'par' and (mon == 12*2+5 or mon == 12*2+7 or mon == 12*5+7):
                MD.T = MD.T*_np.nan
                MD.v = MD.v*_np.nan
                MD.u = MD.u*_np.nan
                MD.w = MD.w*_np.nan
            if SITE == 'par' and mon == 12*3+5:
                MD.v = MD.v*_np.nan
            '''
            for l,w in enumerate(MD.w):
                if abs(w) > 40:
                    MD.w[l] = _np.nan
            gflag = 'on'
            if not(UT):
                # FIX THIS TO BE UNIVERSAL OR DELETE LATER...
                if SITE == 'mor':
                    sltoffset = - _dt.timedelta(minutes=31)
                elif SITE == 'par':
                    sltoffset = - _dt.timedelta(hours=5,minutes=31)
                elif SITE == 'renoir':
                    sltoffset = - _dt.timedelta(hours=2,minutes=12)
                else:
                    sltoffset = _dt.timedelta(minutes=0)
                    print 'No UT/LT conversion... in UT'
            else:
                sltoffset = _dt.timedelta(minutes=0)
            MD.t = MD.t + sltoffset
            #tlim = [MD.t[len(MD.t)/5] + sltoffset,MD.t[-len(MD.t)/5] + sltoffset]
            tlim = [arbdate - _dt.timedelta(hours=7),arbdate + _dt.timedelta(hours=7)]
            tshift = 10

            # Zonal
            axz[nsp].set_ylabel("%s" % (_cal.month_abbr[mon]),labelpad=1)
            axz[nsp].set_ylim(-100, 200)
            axz[nsp].set_xlim(tlim)
            axz[nsp].set_yticks([0, 100])
            if SPLIT:
                axz[nsp].errorbar(MD.t,MD.u,yerr=MD.uv,fmt='.-',color=colors[2],linewidth=linewide,elinewidth=linewide,capsize=markerwide,label='East',markersize=ms) 
                axz[nsp].errorbar(MD.t2,MD.u2,yerr=MD.u2v,fmt='.-',color=colors[3],linewidth=linewide,elinewidth=linewide,capsize=markerwide,label='West',markersize=ms)
            else:
                axz[nsp].errorbar(MD.t+_dt.timedelta(minutes=tshift*k),MD.u,yerr=MD.uv,fmt='.-',color=colors[k],linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms,label='Data')
            axz[nsp].grid(gflag)
            axz[nsp].plot([MD.t[0],MD.t[-1]],[0,0],'k--')

            # Meridional
            axm[nsp].set_ylabel("%s" % (_cal.month_abbr[mon]),labelpad=1)
            axm[nsp].set_ylim(-150, 150)
            axm[nsp].set_xlim(tlim)
            axm[nsp].set_yticks([ -75, 0, 75])
            if SPLIT:
                axm[nsp].errorbar(MD.t,MD.v,yerr=MD.vv,fmt='.-',color=colors[0],linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms,label='East') 
                axm[nsp].errorbar(MD.t2,MD.v2,yerr=MD.v2v,fmt='.-',color=colors[1],linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms,label='West')
            else:
                axm[nsp].errorbar(MD.t+_dt.timedelta(minutes=tshift*k),MD.v,yerr=MD.vv,fmt='.-',color=colors[k],linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms,label='Data')
            axm[nsp].grid(gflag)
            axm[nsp].plot([MD.t[0],MD.t[-1]],[0,0],'k--',label=None)

            # Vertical
            axv[nsp].set_ylabel("%s" % (_cal.month_abbr[mon]),labelpad=1)
            axv[nsp].set_ylim(-60, 60)
            axv[nsp].set_xlim(tlim)
            axv[nsp].set_yticks([ -25, 0, 25])
            axv[nsp].grid(gflag)
            axv[nsp].errorbar(MD.t+_dt.timedelta(minutes=tshift*k),MD.w,yerr=MD.wv,fmt='.-',color=colors[k],linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms)
            axv[nsp].plot([MD.t[0],MD.t[-1]],[0,0],'k--',label=None)

            # Temps
            axt[nsp].set_ylabel("%s" % (_cal.month_abbr[mon]),labelpad=1)
            axt[nsp].set_ylim(600, 1200)
            axt[nsp].set_xlim(tlim)
            axt[nsp].set_yticks([700, 900, 1100])
            #axt[nsp].set_yticks([800, 1000])
            axt[nsp].grid(gflag)
            axt[nsp].errorbar(MD.t+_dt.timedelta(minutes=tshift*k),MD.T,yerr=MD.Tv,fmt='.-',color=colors[k],linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms)

    # Finalize legend
    s = _plt.errorbar([-2,-1],[0,0],yerr=[1,1],fmt='.-',color='b',linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms)
    m = _plt.errorbar([-2,-1],[0,0],yerr=[1,1],fmt='.-',color='g',linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms)
    l = _plt.errorbar([-2,-1],[0,0],yerr=[1,1],fmt='.-',color='r',linewidth=linewide,elinewidth=linewide,capsize=markerwide,markersize=ms)
    
    if len(F_VAL) == 2:
        titles = [r'$\overline{F_{10.7}} < 125$',r'$\overline{F_{10.7}} \geq 125$']
        fz.legend([s,l],titles,bbox_to_anchor=(.15, .9, .76, .102), loc=3,ncol=2, mode="expand", borderaxespad=0.,frameon=False)
        fm.legend([s,l],titles,bbox_to_anchor=(.15, .9, .76, .102), loc=3,ncol=2, mode="expand", borderaxespad=0.,frameon=False)
        fv.legend([s,l],titles,bbox_to_anchor=(.15, .9, .76, .102), loc=3,ncol=2, mode="expand", borderaxespad=0.,frameon=False)
        ft.legend([s,l],titles,bbox_to_anchor=(.15, .9, .76, .102), loc=3,ncol=2, mode="expand", borderaxespad=0.,frameon=False)
    elif len(F_VAL) == 3:
        titles = [r'$\overline{F_{10.7}} < 100$',r'$100 \leq \overline{F_{10.7}} < 200$',r'$\overline{F_{10.7}} \geq 200$']
        fz.legend([s,m,l],titles,bbox_to_anchor=(.15, .9, .76, .102), loc=3,ncol=2, mode="expand", borderaxespad=0.,frameon=False)
        fm.legend([s,m,l],titles,bbox_to_anchor=(.15, .9, .76, .102), loc=3,ncol=2, mode="expand", borderaxespad=0.,frameon=False)
        fv.legend([s,m,l],titles,bbox_to_anchor=(.15, .9, .76, .102), loc=3,ncol=2, mode="expand", borderaxespad=0.,frameon=False)
        ft.legend([s,m,l],titles,bbox_to_anchor=(.15, .9, .76, .102), loc=3,ncol=2, mode="expand", borderaxespad=0.,frameon=False)
       
    # Finalize Zonal
    #fz.suptitle("Monthly Averaged Zonal Wind from %s" % name[SITE], fontsize=12, fontweight='bold')
    fz.subplots_adjust(hspace = 0.001)
    fz.subplots_adjust(wspace = .32)
    fz.subplots_adjust(left=.14)
    if UT:
        fz.text(0.5,0.05,'Hour [UTC]',ha='center',va='center', fontsize=11)
    else:
        fz.text(0.5,0.05,'Hour [SLT]',ha='center',va='center', fontsize=11)
    fz.text(0.02,0.5,'Wind Speed [$m/s$]',ha='center',va='center',rotation='vertical', fontsize=11)
    axz[0].xaxis.set_major_formatter(_md.DateFormatter('%H')) # :%M
    fz.savefig("%s%s-F107b-%s.eps" % (dirout,'Clima-Z',SITE),format='eps')
    #fz.savefig("%s%s-F107b-%s.pdf" % (dirout,'Clima-Z',SITE),format='pdf')
    
    # Finalize Meridional
    #fm.suptitle("Monthly Averaged Meridional Wind from %s" % name[SITE], fontsize=12, fontweight='bold')
    fm.subplots_adjust(hspace = 0.001)
    fm.subplots_adjust(wspace = .32)
    fm.subplots_adjust(left=.14)
    if UT:
        fm.text(0.5,0.05,'Hour [UTC]',ha='center',va='center', fontsize=11)
    else:
        fm.text(0.5,0.05,'Hour [SLT]',ha='center',va='center', fontsize=11)
    fm.text(0.02,0.5,'Wind Speed [$m/s$]',ha='center',va='center',rotation='vertical', fontsize=11)
    axm[0].xaxis.set_major_formatter(_md.DateFormatter('%H')) # :%M
    fm.savefig("%s%s-F107b-%s.eps" % (dirout,'Clima-M',SITE),format='eps')
    #fm.savefig("%s%s-F107b-%s.pdf" % (dirout,'Clima-M',SITE),format='pdf')
    
    # Finalize Vert
    #fv.suptitle("Monthly Averaged Vertical Wind from %s" % name[SITE], fontsize=12, fontweight='bold')
    fv.subplots_adjust(hspace = 0.001)
    fv.subplots_adjust(wspace = .32)
    fv.subplots_adjust(left=.14)
    if UT:
        fv.text(0.5,0.05,'Hour [UTC]',ha='center',va='center', fontsize=11)
    else:
        fv.text(0.5,0.05,'Hour [SLT]',ha='center',va='center', fontsize=11)
    fv.text(0.02,0.5,'Wind Speed [$m/s$]',ha='center',va='center',rotation='vertical', fontsize=11)
    axv[0].xaxis.set_major_formatter(_md.DateFormatter('%H')) # :%M

    fv.savefig("%s%s-F107b-%s.eps" % (dirout,'Clima-V',SITE),format='eps')
    #fv.savefig("%s%s-F107b-%s.pdf" % (dirout,'Clima-V',SITE),format='pdf')
    
    # Finalize Temp
    #ft.suptitle("Monthly Averaged Temperature from %s" % name[SITE], fontsize=11, fontweight='bold')
    ft.subplots_adjust(hspace = 0.001)
    ft.subplots_adjust(wspace = .32)
    ft.subplots_adjust(left=.14)
    if UT:
        ft.text(0.5,0.05,'Hour [UTC]',ha='center',va='center', fontsize=11)
    else:
        ft.text(0.5,0.05,'Hour [SLT]',ha='center',va='center', fontsize=11)
    ft.text(0.02,0.5,'Temperature [$K$]',ha='center',va='center',rotation='vertical', fontsize=11)
    axt[0].xaxis.set_major_formatter(_md.DateFormatter('%H')) # :%M
    ft.savefig("%s%s-F107b-%s.eps" % (dirout,'Clima-T',SITE),format='eps')
    #ft.savefig("%s%s-F107b-%s.pdf" % (dirout,'Clima-T',SITE),format='pdf')
    
    
def PlotLocalTime(SITE,DNSTART,DNEND,HOUR,WMODEL=[],TMODEL='msis',KP=[0,10],QF=1,QUIET=False):
    '''
    Summary:
        Plots single site at a local time w/ models
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        DNSTART = datetime of start date
        DNEND = datetime of end date
        HOUR = List of Local Time Hour to use, e.g. 4
        WMODEL = list of wind models to plot with data [default = []]
        TMODEL = temp model to plot with data [default = 'msis']
        KP = Filter days by KP [default = [0,10] - all kp]
        QF = Quality Flag Limit Allowed [default = 1]
        QUIET  = [model param] flag to set Kp=0 and Ap=Ap_daily [default = False]

    History:
        5/5/16 -- Written by DJF (dfisher2@illinois.edu)
    '''
    _mpl.rcParams.update({'font.size':8})
    _mpl.rcParams['savefig.dpi'] = 200
    _mpl.rcParams['figure.figsize'] = (6,4)   
    
    # Set up Figures
    markerwide = 0
    linewide = 1
    linewide2= 3
    ms  = 2
    ms2 = 2
    color = {'hwm93':'y.','hwm07':'r.','hwm14':'g.'}
    gflag = True
    if isinstance(HOUR, int):
        HOUR = [HOUR]
    h_len = len(HOUR)

    _plt.figure(0)
    fz,((axz))  = _plt.subplots(h_len,1, sharex=True, sharey=False)
    _plt.figure(1)
    fm,((axm))  = _plt.subplots(h_len,1, sharex=True, sharey=False)
    _plt.figure(2)
    fv,((axv))  = _plt.subplots(h_len,1, sharex=True, sharey=False)
    _plt.figure(3)
    ft,((axt))  = _plt.subplots(h_len,1, sharex=True, sharey=False)
    
    # Get LT bin
    lt = []
    dh = [t.hour+t.minute/60. for t in times]
    for hr in HOUR:
        lt.append(int(_np.argmin(abs(_np.roll(dh,b_len/2)-hr))))
    yl = []; dl = [];
    t_lim = [_dt.datetime(1969,12,31),_dt.datetime(1971,1,1)]
    
    # allocate arrays for data and models
    t_len = (DNEND-DNSTART).days
    m_len = 12.  # for each month
    t   = list(_np.empty((t_len))*_np.nan)
    Td  = _np.empty((t_len,h_len))*_np.nan
    Ted = _np.empty((t_len,h_len))*_np.nan
    Ud  = _np.empty((t_len,h_len))*_np.nan
    Ued = _np.empty((t_len,h_len))*_np.nan
    Vd  = _np.empty((t_len,h_len))*_np.nan
    Ved = _np.empty((t_len,h_len))*_np.nan
    Wd  = _np.empty((t_len,h_len))*_np.nan
    Wed = _np.empty((t_len,h_len))*_np.nan
    mt  = list(_np.empty((m_len))*_np.nan)
    Ta  = _np.empty((m_len,h_len))*_np.nan
    Tea = _np.empty((m_len,h_len))*_np.nan
    Ua  = _np.empty((m_len,h_len))*_np.nan
    Uea = _np.empty((m_len,h_len))*_np.nan
    Va  = _np.empty((m_len,h_len))*_np.nan
    Vea = _np.empty((m_len,h_len))*_np.nan
    Wa  = _np.empty((m_len,h_len))*_np.nan
    Wea = _np.empty((m_len,h_len))*_np.nan
    Tm = {}; Um = {}; Vm = {}; Wm = {};
    if isinstance(WMODEL, str):
        WMODEL = [WMODEL]        
    for m in WMODEL:
        Tm[m]  = _np.empty((t_len,h_len))*_np.nan
        Um[m]  = _np.empty((t_len,h_len))*_np.nan
        Vm[m]  = _np.empty((t_len,h_len))*_np.nan
        Wm[m]  = _np.empty((t_len,h_len))*_np.nan

    # Loop though days
    for step in range(t_len):
        dn = DNSTART+_dt.timedelta(days=step)
        doy = (dn-_dt.datetime(dn.year,1,1)).total_seconds()/(60.*60*24)+1
        yl.append(dn.year)   # for monthly binning
        dl.append(doy)   # for monthly binning
        t[step] = dn.replace(year=1970)

        # Grab data for day TODO Weekly!
        DD = BinDailyData(SITE,dn.year,doy,KP=KP,QF=QF)
        for h_ind,hr in enumerate(lt):
            Td[step,h_ind]  = DD.T [hr]
            Ted[step,h_ind] = DD.Te[hr]
            Ud[step,h_ind]  = DD.u [hr]
            Ued[step,h_ind] = DD.ue[hr]
            Vd[step,h_ind]  = DD.v [hr]
            Ved[step,h_ind] = DD.ve[hr]
            Wd[step,h_ind]  = DD.w [hr]
            Wed[step,h_ind] = DD.we[hr]

        # Grab model for day TODO Weekly
        for m in WMODEL:
            MM = GetModels(DD.lla,dn.year,doy,m,TMODEL=TMODEL,ALT=250.,WEIGHTED=True,QUIET=QUIET)
            for h_ind,hr in enumerate(lt):
                Tm[m] [step,h_ind] = MM.T [hr]
                Um[m] [step,h_ind] = MM.u [hr]
                Vm[m] [step,h_ind] = MM.v [hr]
                Wm[m] [step,h_ind] = MM.w [hr]

    for m_ind,mon in enumerate(range(1,13)):
        AD = BinMonthlyData(SITE,1970,mon,DLIST=dl,YLIST=yl,KP=KP,QF=QF,VERBOSE=False)
        mt[m_ind] = AD.dn
        for h_ind,hr in enumerate(lt):
            Ta[m_ind,h_ind]  = AD.T [hr]
            Tea[m_ind,h_ind] = AD.Tv[hr]
            Ua[m_ind,h_ind]  = AD.u [hr]
            Uea[m_ind,h_ind] = AD.uv[hr]
            Va[m_ind,h_ind]  = AD.v [hr]
            Vea[m_ind,h_ind] = AD.vv[hr]
            Wa[m_ind,h_ind]  = AD.w [hr]
            Wea[m_ind,h_ind] = AD.wv[hr]
            

    # Plot this stuff
    for h_ind,hr in enumerate(HOUR):
        # Zonal
        axz[h_ind].set_ylabel("%02i LT" % (hr))
        axz[h_ind].set_xlim(t_lim)
        axz[h_ind].set_ylim(-100, 200)
        axz[h_ind].set_yticks([0, 100])
        axz[h_ind].plot(t,Ud[:,h_ind],'k.',linewidth=linewide,markersize=ms,label='Data')
        axz[h_ind].errorbar(mt,Ua[:,h_ind],yerr=Uea[:,h_ind],fmt='.',color='b', \
                linewidth=linewide,elinewidth=linewide2,capsize=markerwide,markersize=ms2,label='Mon_Ave')
        for m in WMODEL:
            axz[h_ind].plot(t,Um[m][:,h_ind],color[m],alpha=0.33,markersize=ms2)
        axz[h_ind].grid(gflag)
        axz[h_ind].plot(t_lim,[0,0],'k--')

        # Meridional 
        axm[h_ind].set_ylabel("%02i LT" % (hr))
        axm[h_ind].set_xlim(t_lim)
        axm[h_ind].set_ylim(-150, 150)
        axm[h_ind].set_yticks([-75,0,75])
        axm[h_ind].plot(t,Vd[:,h_ind],'k.',linewidth=linewide,markersize=ms,label='Data')
        axm[h_ind].errorbar(mt,Va[:,h_ind],yerr=Vea[:,h_ind],fmt='.',color='b', \
                linewidth=linewide,elinewidth=linewide2,capsize=markerwide,markersize=ms2,label='Mon_Ave')
        for m in WMODEL:
            axm[h_ind].plot(t,Vm[m][:,h_ind],color[m],alpha=0.33,markersize=ms2)
        axm[h_ind].grid(gflag)
        axm[h_ind].plot(t_lim,[0,0],'k--')

        # Vertical
        axv[h_ind].set_ylabel("%02i LT" % (hr))
        axv[h_ind].set_xlim(t_lim)
        axv[h_ind].set_ylim(-75,75)
        axv[h_ind].set_yticks([-25,0,25])
        axv[h_ind].plot(t,Wd[:,h_ind],'k.',linewidth=linewide,markersize=ms,label='Data')
        axv[h_ind].errorbar(mt,Wa[:,h_ind],yerr=Wea[:,h_ind],fmt='.',color='b', \
                linewidth=linewide,elinewidth=linewide2,capsize=markerwide,markersize=ms2,label='Mon_Ave')
        for m in WMODEL:
            axv[h_ind].plot(t,Wm[m][:,h_ind],color[m],alpha=0.33,markersize=ms2)
        axv[h_ind].grid(gflag)
        axv[h_ind].plot(t_lim,[0,0],'k--')

        # Temps 
        axt[h_ind].set_ylabel("%02i LT" % (hr))
        axt[h_ind].set_xlim(t_lim)
        axt[h_ind].set_ylim(600, 1150)
        axt[h_ind].set_yticks([700,800,900,1000,1100])
        axt[h_ind].plot(t,Td[:,h_ind],'k.',linewidth=linewide,markersize=ms,label='Data')
        axt[h_ind].errorbar(mt,Ta[:,h_ind],yerr=Tea[:,h_ind],fmt='.',color='b', \
                linewidth=linewide,elinewidth=linewide2,capsize=markerwide,markersize=ms2,label='Mon_Ave')
        for m in WMODEL:
            axt[h_ind].plot(t,Tm[m][:,h_ind],color[m],alpha=0.33,markersize=ms2)
        axt[h_ind].grid(gflag)

    # Finalize Zonal
    fz.suptitle('Zonal winds')
    fz.subplots_adjust(hspace = 0.001)
    fz.subplots_adjust(wspace = .32)
    fz.subplots_adjust(left=.14)
    fz.text(0.5,0.05,'Date',ha='center',va='center', fontsize=11)
    fz.text(0.02,0.5,'Wind Speed [$m/s$]',ha='center',va='center',rotation='vertical', fontsize=11)
    axz[0].xaxis.set_major_formatter(_md.DateFormatter('%b')) # :%M
    fz.savefig("%s%s-LT-%s.eps" % (dirout,'Clima-Z',SITE),format='eps')
    #fz.savefig("%s%s-F107b-%s.pdf" % (dirout,'Clima-Z',SITE),format='pdf')

    # Finalize Meridional
    fm.suptitle('Meridional winds')
    fm.subplots_adjust(hspace = 0.001)
    fm.subplots_adjust(wspace = .32)
    fm.subplots_adjust(left=.14)
    fm.text(0.5,0.05,'Date',ha='center',va='center', fontsize=11)
    fm.text(0.02,0.5,'Wind Speed [$m/s$]',ha='center',va='center',rotation='vertical', fontsize=11)
    axm[0].xaxis.set_major_formatter(_md.DateFormatter('%b')) # :%M
    fm.savefig("%s%s-LT-%s.eps" % (dirout,'Clima-M',SITE),format='eps')
    #fm.savefig("%s%s-F107b-%s.pdf" % (dirout,'Clima-M',SITE),format='pdf')

    # Finalize Vertical
    fv.suptitle('Vertical winds')
    fv.subplots_adjust(hspace = 0.001)
    fv.subplots_adjust(wspace = .32)
    fv.subplots_adjust(left=.14)
    fv.text(0.5,0.05,'Date',ha='center',va='center', fontsize=11)
    fv.text(0.02,0.5,'Wind Speed [$m/s$]',ha='center',va='center',rotation='vertical', fontsize=11)
    axv[0].xaxis.set_major_formatter(_md.DateFormatter('%b')) # :%M
    fv.savefig("%s%s-LT-%s.eps" % (dirout,'Clima-V',SITE),format='eps')
    #fv.savefig("%s%s-F107b-%s.pdf" % (dirout,'Clima-V',SITE),format='pdf')

    # Finalize Temps
    ft.suptitle('Neutral Temperatures')
    ft.subplots_adjust(hspace = 0.001)
    ft.subplots_adjust(wspace = .32)
    ft.subplots_adjust(left=.14)
    ft.text(0.5,0.05,'Date',ha='center',va='center', fontsize=11)
    ft.text(0.02,0.5,'Temperature [$K$]',ha='center',va='center',rotation='vertical', fontsize=11)
    axt[0].xaxis.set_major_formatter(_md.DateFormatter('%b')) # :%M
    ft.savefig("%s%s-LT-%s.eps" % (dirout,'Clima-T',SITE),format='eps')
    #ft.savefig("%s%s-F107b-%s.pdf" % (dirout,'Clima-T',SITE),format='pdf')



def PlotAverages(SITE,YEAR,MONTH,WMODEL=[],TMODEL='msis',SPLIT=False,KP=[0,10],UT=True,QF=1,HIST=False):
    '''
    Summary:
        Plots single site monthly averages w/ models
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTH = month to plot, e.g. 4
        WMODEL = list of wind models to plot with data [default = []]
        TMODEL = temp model to plot with data [default = 'msis']
        SPLIT = Split look directions in binning [default = False]
        KP = Filter days by KP [default = [0,10] - all kp]
        UT = Plot in UT or SLT [default = True]
        QF = Quality Flag Limit Allowed [default = 1]
 	HIST = Plot histogram of number of data in each bin [default = False]

    History:
        5/8/14 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # Get UT offset if SLT needed
    ut_offset = 0
    if not(UT):
        ut_offset = _fpiinfo.get_site_info(SITE)['Location'][1]/360.*(24*60*60)
        
    # Get Monthly Average for basis
    MD = BinMonthlyData(SITE,YEAR,MONTH,SPLIT=SPLIT,KP=KP,QF=QF)
    mon = _cal.month_name[MONTH]
    MD.t  = MD.t + _dt.timedelta(seconds=ut_offset)
    MD.t2 = MD.t + _dt.timedelta(minutes=3) 

    # Get models
    tm = {}; Tm = {}; Tem = {}; Um = {}; Uem = {}; Vm = {}; Vem = {}; Wm = {}; Wem = {};
    if isinstance(WMODEL, str):
        WMODEL = [WMODEL]
    for nm,m in enumerate(WMODEL):
        if _np.nan in MD.lla:
            break
        MM = BinMonthlyData(m,YEAR,MONTH,SPLIT=SPLIT,SITELLA=MD.lla,TMODEL=TMODEL,WEIGHTED=True,QUIET=False)
        tm[m] = MM.t + _dt.timedelta(minutes=3*(nm+2)) + _dt.timedelta(seconds=ut_offset)
#        if not(UT):
#            tm[m] = tm[m] - _dt.timedelta(seconds=ut_offset) # JJM
        Tm[m] = MM.T
        Tem[m] = MM.Tv
        Um[m] = MM.u
        Uem[m] = MM.uv
        Vm[m] = MM.v
        Vem[m] = MM.vv
        Wm[m] = MM.w
        Wem[m] = MM.wv

    # Get start stop times
    tstart = _np.where(_np.isnan(MD.T[:-1])*_np.isfinite(MD.T[1:]))[0][0]-1
    tend   = _np.where(_np.isfinite(MD.T[:-1])*_np.isnan(MD.T[1:]))[0][0]+1

    color = ['r--','g--','m--','y--']
    facecolor = {'MSIS': 'r','HWM93': 'g', 'HWM07': 'y', 'HWM14': 'm'}
    markerwide = 0
    try:
        name = _fpiinfo.get_site_info(SITE)['Name']
    except:
        name = SITE.title()
    # Temp Figure
    fig = _plt.figure(0,figsize=(10,6)); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    (_, caps1, _) = _plt.errorbar(MD.t,MD.T,yerr=MD.Tv,label='Data',color='k',linewidth=2)
    # Hardcoded Tmodel -- Fix
    for nm,m in enumerate(WMODEL):
        if nm == 0:
	    _plt.fill_between(tm[m],Tm[m]-Tem[m],Tm[m]+Tem[m],alpha=0.5,linewidth=0,facecolor=facecolor['MSIS'])
	    _plt.plot(tm[m],Tm[m],'r',label='MSIS')
#           _plt.errorbar(tm[m],Tm[m],yerr=Tem[m],fmt=color[nm],capsize=markerwide,label='MSIS')
    _plt.plot([MD.t[0],MD.t[-1]],[0,0],'k--')
    _plt.ylim([600,1200])
    _plt.xlim([MD.t[tstart],MD.t[tend]])
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.ylabel('Temperature [K]')
    if UT:
        _plt.xlabel('Time [UT]')
    else:
        _plt.xlabel('Time [SLT]')
    _plt.title('%04i %s Average Temperatures at %s' % (YEAR,mon,name))
    _plt.legend()
    _plt.grid()

    if(HIST):
        bin_width = (MD.t[1]-MD.t[0])/4
        ax2 = ax.twinx()
        ax2.bar(MD.t-bin_width,MD.Tc,.015,alpha=0.5,linewidth=0)
        ax2.set_ylabel('N/bin')
        ax2.set_ylim([0,100])
    _plt.xlim([MD.t[tstart],MD.t[tend]])

    _plt.draw();# _plt.show()
    #_plt.savefig('%s%s_%04d-%s_temps.png' % (dirout,SITE,YEAR,mon))
    
    # Winds Figure Zonal
    fig = _plt.figure(1,figsize=(10,6)); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    for nm,m in enumerate(WMODEL):
        _plt.fill_between(tm[m],Um[m]-Uem[m],Um[m]+Uem[m],alpha=0.5,linewidth=0,facecolor=facecolor[m.upper()])
        _plt.plot(tm[m],Um[m],color=facecolor[m.upper()],capsize=markerwide,label=m.upper())
#        _plt.errorbar(tm[m],Um[m],yerr=Uem[m],fmt=color[nm],capsize=markerwide,label=m.upper())
    if SPLIT:
        _plt.errorbar(MD.t,MD.u,yerr=MD.uv,fmt='b-',capsize=markerwide,label='East')
        _plt.errorbar(MD.t2,MD.u2,yerr=MD.u2v,fmt='c-',capsize=markerwide,label='West')
    else:
        _plt.errorbar(MD.t,MD.u,yerr=MD.uv,label='Data',color='k',capsize=markerwide,linewidth=2)

    _plt.plot([MD.t[0],MD.t[-1]],[0,0],'k--')
    _plt.ylim([-200,200])
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.ylabel('Wind Speed [m/s]')
    if UT:
        _plt.xlabel('Time [UT]')
    else:
        _plt.xlabel('Time [SLT]')
    _plt.title('%04i %s Average Zonal Winds at %s' % (YEAR,mon,name))
    _plt.legend()
    _plt.grid()

    if(HIST):
	bin_width = (MD.t[1]-MD.t[0])/4
        ax2 = ax.twinx()
        ax2.bar(MD.t-bin_width,MD.uc,.015,alpha=0.5,linewidth=0)
        ax2.set_ylabel('N/bin',position=(1,.1875))
	ticklabelpad = _mpl.rcParams['xtick.major.pad']
        ax2.set_ylim([0,80])
	ax2.set_yticks([0,10,20,30])

    _plt.xlim([MD.t[tstart],MD.t[tend]])

    _plt.draw(); #_plt.show()
    _plt.savefig('%s%s_%04d-%s_zonal_winds.png' % (dirout,SITE,YEAR,mon))
    
    # Winds Figure Meridional
    fig = _plt.figure(2,figsize=(10,6)); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    for nm,m in enumerate(WMODEL):
        _plt.fill_between(tm[m],Vm[m]-Vem[m],Vm[m]+Vem[m],alpha=0.5,linewidth=0,facecolor=facecolor[m.upper()])
        _plt.plot(tm[m],Vm[m],color=facecolor[m.upper()],capsize=markerwide,label=m.upper())
#        _plt.errorbar(tm[m],Um[m],yerr=Uem[m],fmt=color[nm],capsize=markerwide,label=m.upper())
    if SPLIT:
        _plt.errorbar(MD.t,MD.v,yerr=MD.vv,fmt='b-',capsize=markerwide,label='North')
        _plt.errorbar(MD.t2,MD.v2,yerr=MD.v2v,fmt='c-',capsize=markerwide,label='South')
    else:
        _plt.errorbar(MD.t,MD.v,yerr=MD.vv,label='Data',color='k',capsize=markerwide,linewidth=2)
    _plt.plot([MD.t[0],MD.t[-1]],[0,0],'k--')
    _plt.ylim([-200,200])
    _plt.xlim([MD.t[tstart],MD.t[tend]])
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.ylabel('Wind Speed [m/s]')
    if UT:
        _plt.xlabel('Time [UT]')
    else:
        _plt.xlabel('Time [SLT]')
    _plt.title('%04i %s Average Meridional Winds at %s' % (YEAR,mon,name))
    _plt.legend()
    _plt.grid()

    if(HIST):
	bin_width = (MD.t[1]-MD.t[0])/4
        ax2 = ax.twinx()
        ax2.bar(MD.t-bin_width,MD.vc,.015,alpha=0.5,linewidth=0)
        ax2.set_ylabel('N/bin',position=(1,.1875))
        ax2.set_ylim([0,80])
        ax2.set_yticks([0,10,20,30])

    _plt.xlim([MD.t[tstart],MD.t[tend]])

    _plt.draw(); #_plt.show()
    #_plt.savefig('%s%s_%04d-%s_meridional_winds.png' % (dirout,SITE,YEAR,mon))
    
    # Winds Figure Vertical
    fig = _plt.figure(3,figsize=(10,6)); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    _plt.errorbar(MD.t,MD.w,yerr=MD.wv,label='Data',color='k',capsize=markerwide,linewidth=2)
#    for nm,m in enumerate(WMODEL):
#        _plt.errorbar(tm[m],Wm[m],yerr=Wem[m],fmt=color[nm],capsize=markerwide,label=m.upper())
    _plt.plot([MD.t[0],MD.t[-1]],[0,0],'k--')
    _plt.ylim([-75,75])
    _plt.xlim([MD.t[tstart],MD.t[tend]])
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.ylabel('Wind Speed [m/s]')
    if UT:
        _plt.xlabel('Time [UT]')
    else:
        _plt.xlabel('Time [SLT]')
    _plt.title('%04i %s Average Vertical Winds at %s' % (YEAR,mon,name))
    _plt.legend()
    _plt.grid()

    if(HIST):
	bin_width = (MD.t[1]-MD.t[0])/4
        ax2 = ax.twinx()
        ax2.bar(MD.t-bin_width,MD.wc,.015,alpha=0.5,linewidth=0)
        ax2.set_ylabel('N/bin',position=(1,.1875))
        ax2.set_ylim([0,80])
        ax2.set_yticks([0,10,20,30])

    _plt.xlim([MD.t[tstart],MD.t[tend]])

    _plt.draw(); #_plt.show()
    #_plt.savefig('%s%s_%04d-%s_vertical_winds.png' % (dirout,SITE,YEAR,mon))
    
    return(MD)


def PlotSpaghetti(SITE,YEAR,MONTH,SPLIT=False,LIST=[],CV=False,KP=[0,10],QF=1,WIS0=False):
    '''
    Summary:
        Plots all raw data for one month in spaghetti plot!
        TODO: multicolor lines w/ legend or transparent blue lines...
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTH = month to plot, e.g. 4
        SPLIT = Split look directions in binning [default = False]
        LIST = List of doys to use in averaging
        CV = Include CV directions [default = False]
        KP = Filter days by KP [default = [0,10] - all kp]
        QF = Quality Flag Allowed [default = 1]
        WIS0  = Flag to use w=0 in processing [default = False]

    History:
        10/17/14 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # Get Monthly Average for basis
    D = {}
    dn = _dt.datetime(YEAR,MONTH,1)
    if not LIST:
        for k,d in enumerate(range(_cal.monthrange(YEAR,MONTH)[1])):
            L = L2.GetLevel2(SITE,dn+_dt.timedelta(days=d),w_is_0=WIS0)
            FilterData(L,QF); D[k] = L
    else:
        for k,d in enumerate(LIST):
            L = L2.GetLevel2(SITE,_dt.datetime(YEAR,1,1)+_dt.timedelta(days=d-1),w_is_0=WIS0)
            FilterData(L,QF); D[k] = L
    M = BinMonthlyData(SITE,YEAR,MONTH,SPLIT=SPLIT,CV=CV,DLIST=LIST,KP=KP,QF=QF)
    M.t = M.t + _dt.timedelta(days=(_dt.datetime(YEAR-1,1,2+_cal.isleap(YEAR-1))-_dt.datetime(M.t[0].year,1,1)).days)

    markerwide = 0
    lalpha = 0.1
    lw=3
    
    # Get xlimits w/ data (IMPROVE)
    #tlim = [(D[0][0].allt[0]-_dt.timedelta(days=(dn-_dt.datetime(YEAR,1,1)).days+1,minutes=80)).astimezone(_utc).replace(tzinfo=None), (D[0][0].allt[-1]-_dt.timedelta(days=(dn-_dt.datetime(YEAR,1,1)).days+1,minutes=50)+_dt.timedelta(hours=1)).astimezone(_utc).replace(tzinfo=None)]
    tlim = [M.t[len(M.t)/5],M.t[-len(M.t)/5]]
 
    _plt.close('all');ax={}
    # Figure
    f,(ax[0],ax[1],ax[3])  = _plt.subplots(3, sharex=True, sharey=False)
    # Plot Data
    for d in D:
        try:
            doy = (D[d][0].dn-_dt.datetime(YEAR,1,1)).days+1
        except:
            continue
        # Temps
        for l in [y for y in D[d]]:
            ax[3].fill_between(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy-1)]),l.T-l.Te,l.T+l.Te,alpha=lalpha,linewidth=0,facecolor='k')
            #print [x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy-1)]
            '''
            ax[3].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,linewidth=lw,alpha=lalpha,label='Doy %03i'%doy)
            '''
        # Zonal
        for l in [y for y in D[d] if ('East' in y.key or 'West' in y.key)]:
        #for l in [y for y in D[d] if ('East' in y.key)]:
            ax[0].fill_between(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy-1)]),l.u-l.ue,l.u+l.ue,alpha=lalpha,linewidth=0,facecolor='k')
            '''
            ax[0].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,linewidth=lw,alpha=lalpha,label='Doy %03i'%doy)
            '''
        # Meridional
        for l in [y for y in D[d] if ('South' in y.key or 'North' in y.key)]:
        #for l in [y for y in D[d] if ('South' in y.key)]:
            ax[1].fill_between(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy-1)]),l.v-l.ve,l.v+l.ve,alpha=lalpha,linewidth=0,facecolor='k')
            '''
            ax[1].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,linewidth=lw,alpha=lalpha,label='Doy %03i'%doy)
            '''
    # Overlay Monthly Average
    ax[3].errorbar(M.t,M.T,yerr=M.Tv,fmt='b.-')
    ax[0].errorbar(M.t,M.u,yerr=M.uv,fmt='b.-')
    ax[1].errorbar(M.t,M.v,yerr=M.vv,fmt='b.-')

    # Overlay 16th and 84th percentiles
    '''
    ax[0].plot(M.t,M.u16,'r.')
    ax[0].plot(M.t,M.u84,'r.')
    ax[1].plot(M.t,M.v16,'r.')
    ax[1].plot(M.t,M.v84,'r.')
    ax[3].plot(M.t,M.T16,'r.')
    ax[3].plot(M.t,M.T84,'r.')
    '''
    
    # finalize Temps
    f.subplots_adjust(hspace=0)
    ax[3].set_ylim([600,1200])
    _plt.xlim(tlim)
    ax[3].set_ylabel('Temperature [K]')
    yticks = ax[3].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    # finalize Zonal
    ax[0].plot(tlim,[0,0],'k--')
    f.subplots_adjust(hspace=0)
    ax[0].set_ylim([-100,200])
    _plt.xlim(tlim)
    ax[0].set_ylabel('Zonal Wind [m/s]')
    yticks = ax[0].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    # finalize Meridional
    ax[1].plot(tlim,[0,0],'k--')
    f.subplots_adjust(hspace=0)
    ax[1].set_ylim([-100,100])
    _plt.xlim(tlim)
    ax[1].set_ylabel('Meridional Wind [m/s]')
    yticks = ax[1].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    '''
    ax[2].plot(MR.t,MR.v)
    # finalize Vertical
    f.subplots_adjust(hspace=0)
    ax[3].ylim([-75,75])
    _plt.xlim(tlim)
    ax[3].set_ylabel('Vertical Wind [m/s]')
    yticks = ax[3].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    '''
    #ax[0].annotate('%04d %s'%(YEAR,_cal.month_name[MONTH]), xy=(0.53,0.92), xycoords='axes fraction', fontsize=12)
    ax[3].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[0].set_title('Spaghetti %s %04d'%(_cal.month_name[MONTH],YEAR))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    #ax[0].legend(loc=2, prop={'size':12}, bbox_to_anchor=(1.1, 0.5), fancybox=True)
    _plt.draw();
    _plt.savefig('%s%s_Spaghetti_%04d-%02d.png'%(dirout,SITE,YEAR,MONTH))



def PlotGridMonth(SITE,YEAR,MONTH,SPLIT=True,WIS0=False):
    '''
    Summary:
        Plots all raw cardinal data for one month in 5x7 plot! 
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTH = month to plot, e.g. 4
        SPLIT = Split look directions in binning [default = True]
        WIS0  = Flag to use w=0 in processing [default = False]

    History:
        10/17/14 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # Get Monthly Average for basis
    D = {}
    dn = _dt.datetime(YEAR,MONTH,1)
    for k,d in enumerate(range(_cal.monthrange(YEAR,MONTH)[1])):
        D[k] = L2.GetLevel2(SITE,dn+_dt.timedelta(days=d),w_is_0=WIS0)
    #M = BinMonthlyData(SITE,YEAR,MONTH,SPLIT=SPLIT)
    #M.t = M.t + _dt.timedelta(seconds=(_dt.datetime(YEAR-1,1,1)-_dt.datetime(M.t[0].year,1,1)).total_seconds())

    cs = 0
    lalpha = .3
    lw=2
    
    # Get xlimits w/ data (IMPROVE)
    #tlim = [(D[0][0].allt[0]-_dt.timedelta(days=(dn-_dt.datetime(YEAR,1,1)).days+1,minutes=80)).astimezone(_utc).replace(tzinfo=None), (D[0][0].allt[-1]-_dt.timedelta(days=(dn-_dt.datetime(YEAR,1,1)).days+1,minutes=50)+_dt.timedelta(hours=1)).astimezone(_utc).replace(tzinfo=None)]
    tlim = [_dt.datetime(YEAR,1,1)-_dt.timedelta(hours=4), _dt.datetime(YEAR,1,1)+_dt.timedelta(hours=12)]
    
    _plt.close('all');ax={}
    # Figure 1
    f,((ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]),(ax[7],ax[8],ax[9],ax[10],ax[11],ax[12],ax[13]),(ax[14],ax[15],ax[16],ax[17],ax[18],ax[19],ax[20]),(ax[21],ax[22],ax[23],ax[24],ax[25],ax[26],ax[27]),(ax[28],ax[29],ax[30],ax[31],ax[32],ax[33],ax[34]))  = _plt.subplots(5,7, sharex=True, sharey=True)
    # Plot Data
    for d in D:
        try:
            doy = (D[d][0].dn-_dt.datetime(YEAR,1,1)).days+1
        except:
            continue
        # Temps
        for l in [y for y in D[d] if ('North' in y.key)]:
            #print (l.t1[0]-_dt.timedelta(days=doy)).astimezone(_utc).replace(tzinfo=None)
            l1 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,capsize=cs,fmt='r.')
            m1 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'r-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('South' in y.key)]:
            l2 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,capsize=cs,fmt='g.')
            m2 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'g-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('East' in y.key)]:
            l3 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,capsize=cs,fmt='y.')
            m3 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'y-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('West' in y.key)]:
            l4 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,capsize=cs,fmt='m.')
            m4 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'m-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('Zenith' in y.key)]:
            l5 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,capsize=cs,fmt='b.')
            m5 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'b-',linewidth=lw,alpha=lalpha)
        ax[d].annotate('Doy %03d'%doy, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    ax[0].set_ylim([500,1000])
    _plt.xlim(tlim)
    for d in [0,7,14,21,28]:
        yticks = ax[d].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
    for d in range(21,35):
        xticks = ax[d].xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)
        xticks[2].label1.set_visible(False)
        xticks[4].label1.set_visible(False)
        xticks[-3].label1.set_visible(False)
        xticks[-1].label1.set_visible(False)
    try:
        ax[34].legend([(l1,m1),(l2,m2),(l3,m3),(l4,m4),(l5,m5)],['North','South','East','West','Z'],ncol=3)
    except:
        print 'not enough data for legend...'
    ax[0].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[3].set_title('%s Temperature Grid %s %04d'%(SITE,_cal.month_name[MONTH],YEAR))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    _plt.draw();
    _plt.savefig('%s%s_GridT_%04d-%02d.png'%(dirout,SITE,YEAR,MONTH))
            
    # Figure 2        
    ax={}; f,((ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]),(ax[7],ax[8],ax[9],ax[10],ax[11],ax[12],ax[13]),(ax[14],ax[15],ax[16],ax[17],ax[18],ax[19],ax[20]),(ax[21],ax[22],ax[23],ax[24],ax[25],ax[26],ax[27]),(ax[28],ax[29],ax[30],ax[31],ax[32],ax[33],ax[34]))  = _plt.subplots(5,7, sharex=True, sharey=True)
    # Plot Data
    for d in D:
        try:
            doy = (D[d][0].dn-_dt.datetime(YEAR,1,1)).days+1
        except:
            continue
        # Zonal
        for l in [y for y in D[d] if ('East' in y.key)]:
            l1 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,yerr=l.ue,capsize=cs,fmt='yo')
            m1 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,'y-',linewidth=lw,alpha=lalpha,label='Doy %03i'%doy)
        for l in [y for y in D[d] if ('West' in y.key)]:
            l2 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,yerr=l.ue,capsize=cs,fmt='mo')
            m2 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,'m-',linewidth=lw,alpha=lalpha,label='Doy %03i'%doy)
        ax[d].plot(tlim,[0,0],'k--')
        ax[d].annotate('Doy %03d'%doy, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    ax[0].set_ylim([-100,200])
    _plt.xlim(tlim)
    for d in [0,7,14,21,28]:
        yticks = ax[d].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
    for d in range(21,35):
        xticks = ax[d].xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)
        xticks[2].label1.set_visible(False)
        xticks[4].label1.set_visible(False)
        xticks[-3].label1.set_visible(False)
        xticks[-1].label1.set_visible(False)
    try:
        ax[34].legend([(l1,m1),(l2,m2)],['East','West'])
    except:
        print 'not enough data for legend...'
    ax[0].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[3].set_title('%s Zonal Wind Grid %s %04d'%(SITE,_cal.month_name[MONTH],YEAR))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    _plt.draw();
    _plt.savefig('%s%s_GridU_%04d-%02d.png'%(dirout,SITE,YEAR,MONTH))
            
    # Figure 3
    ax={};f,((ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]),(ax[7],ax[8],ax[9],ax[10],ax[11],ax[12],ax[13]),(ax[14],ax[15],ax[16],ax[17],ax[18],ax[19],ax[20]),(ax[21],ax[22],ax[23],ax[24],ax[25],ax[26],ax[27]),(ax[28],ax[29],ax[30],ax[31],ax[32],ax[33],ax[34]))  = _plt.subplots(5,7, sharex=True, sharey=True)
    # Plot Data
    for d in D:
        try:
            doy = (D[d][0].dn-_dt.datetime(YEAR,1,1)).days+1
        except:
            continue
        # Meridional
        for l in [y for y in D[d] if ('North' in y.key)]:
            l1 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,yerr=l.ve,capsize=cs,fmt='ro')
            m1 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'r-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('South' in y.key)]:
            l2 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,yerr=l.ve,capsize=cs,fmt='go')
            m2 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'g-',linewidth=lw,alpha=lalpha)
        ax[d].plot(tlim,[0,0],'k--')
        ax[d].annotate('Doy %03d'%doy, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    ax[0].set_ylim([-150,150])
    _plt.xlim(tlim)
    for d in [0,7,14,21,28]:
        yticks = ax[d].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
    for d in range(21,35):
        xticks = ax[d].xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)
        xticks[2].label1.set_visible(False)
        xticks[4].label1.set_visible(False)
        xticks[-3].label1.set_visible(False)
        xticks[-1].label1.set_visible(False)
    try:
        ax[34].legend([(l1,m1),(l2,m2)],['North','South'])
    except:
        print 'not enough data for legend...'
    ax[0].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[3].set_title('%s Meridional Wind Grid %s %04d'%(SITE,_cal.month_name[MONTH],YEAR))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    _plt.draw();
    _plt.savefig('%s%s_GridV_%04d-%02d.png'%(dirout,SITE,YEAR,MONTH))
    
    '''
    # Figure 4 
    ax={};f,((ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]),(ax[7],ax[8],ax[9],ax[10],ax[11],ax[12],ax[13]),(ax[14],ax[15],ax[16],ax[17],ax[18],ax[19],ax[20]),(ax[21],ax[22],ax[23],ax[24],ax[25],ax[26],ax[27]),(ax[28],ax[29],ax[30],ax[31],ax[32],ax[33],ax[34]))  = _plt.subplots(5,7, sharex=True, sharey=True)
    # Plot Data
    for d in D:
        try:
            doy = (D[d][0].dn-_dt.datetime(YEAR,1,1)).days+1
        except:
            continue
        # Vertical
        for l in [y for y in D[d] if ('Zenith' in y.key)]:
            l1 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.we,capsize=cs,fmt='bo')
            m1 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,'b-',linewidth=lw,alpha=lalpha)
        ax[d].plot(tlim,[0,0],'k--')
        ax[d].annotate('Doy %03d'%doy, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    ax[0].set_ylim([-75,75])
    _plt.xlim(tlim)
    for d in [0,7,14,21,28]:
        yticks = ax[d].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
    for d in range(21,35):
        xticks = ax[d].xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)
        xticks[2].label1.set_visible(False)
        xticks[4].label1.set_visible(False)
        xticks[-3].label1.set_visible(False)
        xticks[-1].label1.set_visible(False)
    ax[34].legend([(l1,m1)],['Zenith'])
    ax[0].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[3].set_title('%s Vertical Wind Grid %s %04d'%(SITE,_cal.month_name[MONTH],YEAR))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    _plt.draw();
    _plt.savefig('%s%s_GridW_%04d-%02d.png'%(dirout,SITE,YEAR,MONTH))
    ''' 
    
    # Figure 5
    ax={};f,((ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]),(ax[7],ax[8],ax[9],ax[10],ax[11],ax[12],ax[13]),(ax[14],ax[15],ax[16],ax[17],ax[18],ax[19],ax[20]),(ax[21],ax[22],ax[23],ax[24],ax[25],ax[26],ax[27]),(ax[28],ax[29],ax[30],ax[31],ax[32],ax[33],ax[34]))  = _plt.subplots(5,7, sharex=True, sharey=True)
    # Plot Data
    for d in D:
        try:
            doy = (D[d][0].dn-_dt.datetime(YEAR,1,1)).days+1
        except:
            continue
        # Meridional
        for l in [y for y in D[d] if ('CV' in y.key and '_1' in y.key)]:
            l1 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.ve,capsize=cs,fmt='ro')
            m1 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'r-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('CV' in y.key and '_2' in y.key)]:
            l2 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.ve,capsize=cs,fmt='go')
            m2 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'g-',linewidth=lw,alpha=lalpha)
        '''
        for l in [y for y in D[d] if ('IN' in y.key)]:
            l3 = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,yerr=l.ve,capsize=cs,fmt='bo')
            m3 = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'b-',linewidth=lw,alpha=lalpha)
        '''
        ax[d].plot(tlim,[0,0],'k--')
        ax[d].annotate('Doy %03d'%doy, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    ax[0].set_ylim([-150,150])
    _plt.xlim(tlim)
    for d in [0,7,14,21,28]:
        yticks = ax[d].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
    for d in range(21,35):
        xticks = ax[d].xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)
        xticks[2].label1.set_visible(False)
        xticks[4].label1.set_visible(False)
        xticks[-3].label1.set_visible(False)
        xticks[-1].label1.set_visible(False)
    try:
        ax[34].legend([(l1,m1),(l2,m2)],['CV1','CV2'])
    except:
        print 'not enough data for legend...'
    ax[0].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[3].set_title('%s CV Wind Grid %s %04d'%(SITE,_cal.month_name[MONTH],YEAR))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    _plt.draw();
    _plt.savefig('%s%s_GridX_%04d-%02d.png'%(dirout,SITE,YEAR,MONTH))   


       


def CreateL2ASCII(PROJECT,YEAR,DOY,WIS0=False):
    '''
    Summary:
        Script to save out ASCII of all L2 winds for a Project.  Filtered, w=0 assumed.
    
    Inputs:
        PROJECT - Name of project for sites, ex: 'NATION'
        YEAR = year, e.g. 2013
        DOY = day of year, e.g. 64
        WIS0  = Flag to use w=0 in processing [default = False]

    Outputs:
        ASCII text file of data in /mnt/FPIData/Results/ASCII

    History:
        5/17/13 -- Written by Daniel J. Fisher (dfisher2@illinois.edu)
    '''
    # Debug statment... if needed
    #print "Year-",YEAR,"Doy-",DOY
    
    # Create the YYYYMMDD date format
    process_dn = _dt.datetime(YEAR,1,1) + _dt.timedelta(days = DOY-1)
    year = process_dn.strftime('%Y')
    date = process_dn.strftime('%Y%m%d')

    # Create Folder/file labels
    results_stub = '/rdata/airglow/database/L2/'
    notename = results_stub + PROJECT + '_' + date + '.txt'
    D = L2.GetLevel2(PROJECT,process_dn,w_is_0=WIS0)

    # Write out ASCII
    note = open(notename,'w')
    
    # Write out level 2 Data
    note.write('\nLEVEL 2 DATA PRODUCT:\n---------------------\n')

    note.write('Direction  Time[UTC]  Lat  Lon  Alt[km]  u[m/s]  uError[m/s]  v[m/s]  vError[m/s]  CloudDT[C]\n')
    
    for x in D:
	    for i in range(len(x.t1)):
	        if not(x.error): # Remove bad look combo
	            dn = x.t1[i].astimezone(_utc)
	            utctime = dn.strftime("%Y-%m-%d %H:%M:%S")
	            if('North' in x.key or 'South' in x.key):
	                line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  250.0  -999.00  -999.00  {:7.2f}  {:7.2f}  {:1.0f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.v[i], x.ve[i], x.flag_wind[i])
	            elif('East' in x.key or 'West' in x.key):
	                line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  250.0  {:7.2f}  {:7.2f}  -999.00  -999.00  {:1.0f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.u[i], x.ue[i], x.flag_wind[i])
	            elif('CV_VTI_EKU_PAR' in x.key):
	                line = ""
	            elif('CV' in x.key):
	                line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  250.0  {:7.2f}  {:7.2f}  {:7.2f}  {:7.2f}  {:1.0f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.u[i], x.ue[i], x.v[i], x.ve[i], x.flag_wind[i])
	            # For debugg
	            elif('Zenith' in x.key or 'IN' in x.key):
	                line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  250.0  VERTICAL  {:7.2f}  {:7.2f}  {:1.0f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.w[i], x.we[i], x.flag_wind[i])
	            else:
	                line = ""
	            note.write(line)
		    
    note.close()
    
    print 'Results saved to: ' + notename
    print 'This is a test case, not for official L2 use'
    
    

def CreateL2ASCII_Legacy(PROJECT,YEAR,DOY,QF=1,WIS0=False):
    '''
    Summary:
        Script to save out Legacy ASCII of all L2 date for a Project. [for Meriwether]
    
    Inputs:
        PROJECT - Name of project for sites, ex: 'NATION'
        YEAR = year, e.g. 2013
        DOY = day of year, e.g. 64
        QF = Quality Flag Limit Allowed [default = 1]
        WIS0  = Flag to use w=0 in processing [default = False]

    Outputs:
        ASCII text file of data in /mnt/FPIData/Results/ASCII

    History:
        5/17/13 -- Written by Daniel J. Fisher (dfisher2@illinois.edu)
    '''
    # Debug statment... if needed
    #print "Year-",YEAR,"Doy-",DOY
    
    # Create the YYYYMMDD date format
    process_dn = _dt.datetime(YEAR,1,1) + _dt.timedelta(days = DOY-1)
    year = process_dn.strftime('%Y')
    date = process_dn.strftime('%Y%m%d')

    # Create Folder/file labels
    results_stub = '/rdata/airglow/database/L2/'
    notename = results_stub + PROJECT + '_' + date + 'L.txt'
    D = L2.GetLevel2(PROJECT,process_dn,w_is_0=WIS0)
    FilterData(D,QF)

    # Write out ASCII
    note = open(notename,'w')
    
    # Write out level 2 Data
    note.write('\nLEVEL 2 DATA PRODUCT:\n---------------------\n')

    note.write('Direction  Time[UTC]  Lat  Lon  T[K]  TError[K]  u[m/s]  uError[m/s]  v[m/s]  vError[m/s]  w[m/s]  wError[m/s]  Intensity  Intensity_Error  Background  Background_Error  Note\n')
    
    for x in D:
        for i in range(len(x.t1)):
            if not(x.error): # Remove bad look combo
		        dn = x.t1[i].astimezone(_utc)
		        utctime = dn.strftime("%Y-%m-%d %H:%M:%S")
		        if('Zenith' in x.key or 'IN' in x.key):
		            line = "{:16s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  -------  ------  -------  ------  {:7.2f}  {:6.2f}  {:6.4f}  {:4.4f}  {:6.1f}  {:4.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.w[i], x.we[i], x.i[i], x.ie[i], x.b[i], x.be[i])
		        elif('North' in x.key or 'South' in x.key):
		            line = "{:16s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  -------  ------  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:6.4f}  {:4.4f}  {:6.1f}  {:4.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.v[i], x.ve[i], x.wi[i], x.wie[i], x.i[i], x.ie[i], x.b[i], x.be[i])
		        elif('East' in x.key or 'West' in x.key):
		            line = "{:16s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  -------  ------  {:7.2f}  {:6.2f}  {:6.4f}  {:4.4f}  {:6.1f}  {:4.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.u[i], x.ue[i], x.wi[i], x.wie[i], x.i[i], x.ie[i], x.b[i], x.be[i])
		        elif('CV_VTI_EKU_PAR' in x.key):
		            line = ""
		        elif('Unknown' in x.key):
		            line = ""
                        # Code to add MTM Temps (REMOVE SOON)
		        elif('MTM_Search' in x.key):
		            line = "{:16s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  -------  ------  -------  ------  -------  ------  {:6.4f}  {:4.4f}  {:6.1f}  {:4.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.i[i], x.ie[i], x.b[i], x.be[i])
                        elif('CV' in x.key):
		            line = "{:16s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:6.4f}  {:4.4f}  {:6.1f}  {:4.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.u[i], x.ue[i], x.v[i], x.ve[i], x.wi[i], x.wie[i], x.i[i], x.ie[i], x.b[i], x.be[i])
		        #line = "%14s  %19s  %3.1f  %3.1f  %4.2f  %2.2f  %3.2f  %2.2f  %3.2f  %2.2f  %3.2f  %2.2f  %1.3f  %1.3f  %5s  %30s" % (x.key, utctime, lat, lon, x.T[i], x.Te[i], x.u[i], x.ue[i], x.v[i], x.ve[i], x.w[i], x.we[i], x.I[i], x.Ie[i], x.cloudy[i], x.notes)
		        note.write(line)
		    
    note.close()
    

def CreateMonthlyASCII(PROJECT,YEAR,MONTH,QF=1):
    '''
    Summary:
        Script to output ASCII Monthly averaged results for a project.
    
    Inputs:
        PROJECT - Name of project for sites, ex: 'NATION'
        YEAR = year, e.g. 2013
        MONTH = month, e.g. 6
        QF = Quality Flag Limit Allowed [default = 1]

    Outputs:
        ASCII text file of data in /mnt/FPIData/Results/ASCII

    History:
        5/17/13 -- Written by Daniel J. Fisher (dfisher2@illinois.edu)
    '''
    # Debug if needed
    #print "Year-",YEAR,"Month-",MONTH
    
    # Create Folder/file labels
    results_stub = '/rdata/airglow/database/L2/'
    S = _fpiinfo.get_network_info(PROJECT)
    SITES = S.keys()
    
    notename = results_stub + PROJECT + '_' + str(YEAR) + 'M' + str(MONTH).zfill(2) + '.txt'

    # Write out ASCII
    note = open(notename,'w')
    
    # Write out level 3 Data
    note.write('\nLEVEL 3 MONTHLY PRODUCT:\n---------------------\n')
    note.write('Hour[UTC]  Lat  Lon  Temp[K]  Temp_Error[K]  Zonal_Wind[m/s]  Zonal_Wind_Error[m/s] Meridional_Wind[m/s]  Meridional_Wind_Error[m/s] Vertical_Wind[m/s]  Vertical_Wind_Error[m/s] \n')
        
    for site in SITES:
        MD = BinMonthlyData(site,YEAR,MONTH,QF=QF,VERBOSE=False)
        
        for ind,dn in enumerate(MD.t):
            #dn = MD.t[ind] #.astimezone(_utc)
            utctime = dn.strftime("%H:%M")
            line = "%5s  %3.1f  %3.1f  %4.2f  %2.2f  %3.2f  %2.2f  %3.2f  %2.2f  %3.2f  %2.2f\n" % (utctime, S[site]['Location'][0], S[site]['Location'][1], MD.T[ind], MD.Tv[ind], MD.u[ind], MD.uv[ind], MD.v[ind], MD.vv[ind], MD.w[ind], MD.wv[ind])
            note.write(line)

    note.close
    
    

class _BinnedData:
    '''
    This class holds all variables containing binned date.
    Note: BinnedData has T,u,v,w (u2,v2) data
        Each also has c-count of data points per time, e-weighted std, v-monthly variablility
    '''
    def __init__(self, dn, site):

        self.dn = dn
        self.site = site
        self.moonup = False
        self.lla = _np.array([])
        self.key = ""
        self.f = ""
        self.u  = _np.array([]) 
        self.ue = _np.array([])
        self.v  = _np.array([])
        self.ve = _np.array([])
        self.w  = _np.array([])
        self.we = _np.array([])
        self.cards = _np.array([])
        self.cvs = _np.array([])
        self.alts = _np.array([])
        self.barrelroll = False
        try:
            self.project = _fpiinfo.get_site_info(site)['Network']
        except:
            self.project = site

        # interpolated stuff:
        self.it  = _np.array([])
        self.iw  = _np.array([])
        self.iwe = _np.array([])

        self.T  = _np.array([])
        self.Te = _np.array([])
        self.t = _np.array([])
        self.ut_time = _np.array([])
        self.log = ""
        self.notes = ""
        self.rev = "??"
        self.error = False
        self.errorT = False
        
    def __str__(self):
        string = ""
        string += "%11s" % "dn = "     + self.dn.strftime("%Y-%m-%d") + "\n"
        string += "%11s" % "project = " + self.project + "\n"
        string += "%11s" % "f = " + self.f + "\n"
        string += "%11s" % "key = " + self.key + "\n"
        string += "%11s" % "log = " + self.log + "\n"
        string += "%11s" % "notes = " + self.notes + "\n"
        return string
        
    def cut(self, dn1, dn2, inds=None):
        # dn1 and dn2 are LT

        if inds is None:
            t1 = _np.array([dn.replace(tzinfo=None) for dn in self.t])
            inds = _np.where( (t1 > dn1) * (t1 < dn2) )

        self.length = len(inds)
        
        if len(self.flag_wind) > 0:
            self.flag_wind = self.flag_wind[inds]
            
        if len(self.flag_T) > 0:
            self.flag_T = self.flag_T[inds]
            
        if len(self.it) > 0:
            self.it = self.it[inds]

        if len(self.iw) > 0:
            self.iw = self.iw[inds]

        if len(self.iwe) > 0:
            self.iwe = self.iwe[inds]

        if len(self.T) > 0:
            self.T = self.T[inds]

        if len(self.t) > 0:
            self.t = self.t[inds]

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

        if len(self.alts) > 0:
            self.alts = self.alts[inds]
        return

    def plot(self, ):
        
        _mpl.rcParams['font.family'] = 'monospace'
        xlim = [self.t[len(times)/5], self.t[-len(times)/5]]

        if self.error:
            return None

        # Plot Winds
        fig = _plt.figure();
        _plt.clf()

        ax = fig.add_axes((.1,.2,.8,.7)) # left, bottom, width, height
        
        if self.key is 'Daily':
            _plt.errorbar(self.t, self.u, yerr=self.ue, \
                    color='b', marker='o', label='u')
            _plt.errorbar(self.t, self.v, yerr=self.ve, \
                    color='g', marker='.', label='v')
            _plt.errorbar(self.t, self.w, yerr=self.we, \
                    color='r', marker='*', label='w')
        else:
            _plt.errorbar(self.t, self.u, yerr=self.uv, \
                    color='b', marker='o', label='u')
            _plt.errorbar(self.t, self.v, yerr=self.vv, \
                    color='g', marker='.', label='v')
            _plt.errorbar(self.t, self.w, yerr=self.wv, \
                    color='r', marker='*', label='w')
        _plt.plot(xlim,[0,0],'k--')

        _plt.xlim( xlim )
        _plt.ylim([-200.,200.]) 
        ax.xaxis.set_major_formatter(_md.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        _plt.legend()
        _plt.grid()
        fig.text(.1,.05,self.notes)
        datestr = self.dn.strftime("%Y-%m-%d")
        fig.text(.1,.92,"%10s, %12s, %10s" % (self.site, self.key, datestr))
        
        fig.text(.7,.030, self.log)
        _plt.draw();
        _plt.show()

        if self.errorT:
            return None

        # Plot Temps
        fig = _plt.figure();
        _plt.clf()

        ax = fig.add_axes((.1,.2,.8,.7)) # left, bottom, width, height
        
        if 'Daily' in self.key:
            _plt.errorbar(self.t, self.T, yerr=self.Te, \
                    color='r', marker='.', label='T')
        else:
            _plt.errorbar(self.t, self.T, yerr=self.Tv, \
                    color='r', marker='.', label='T')

        _plt.xlim( [self.t[len(times)/5], self.t[-len(times)/5]] )
        _plt.ylim([500.,1200.]) 
        ax.xaxis.set_major_formatter(_md.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        _plt.legend()
        _plt.grid()
        fig.text(.1,.05,self.notes)
        datestr = self.dn.strftime("%Y-%m-%d")
        fig.text(.1,.92,"%10s, %12s, %10s" % (self.site, self.key, datestr))
        
        fig.text(.7,.030, self.log)
        _plt.draw();
        _plt.show()
        return 0
        
    def doabarrelroll(self):
        '''
        Summary:
            Since we are binning in UT, some sites start before or after 00:00. This shifts the data to correct for this and make plots look nicer...

        History:
            9/20/14 -- Written by DJF (dfisher2@illinois.edu)

        '''
        # Get lenth of 24hr span
        fox = len(self.t)
        # Find longitude as fraction of 24hour span
        peppy = self.lla[1]%360/360.
        if _np.isnan(peppy):
            try:
                peppy = _fpiinfo.get_site_info(self.site)['Location'][1]%360/360.
            except:
                if self.site is 'renoir':
                    peppy = -37.5%360/360.
                elif self.site is 'peru':
                    peppy = -76%360/360.
                else:
                    peppy = 0
                    print 'No Location Specified, assume UTC'

        # To force data in middle of 00:00-23:59, need +/-12 UTC offset == 0.5 loc
        roll = int(round((0.5+peppy)%1*fox))
        # Check to see if already done & revert
        if self.barrelroll:
            for ti in range(roll):
                self.t[ti] = self.t[ti]+_dt.timedelta(days=1)
            roll = -1*roll
        # Get variables in _BinnedData
        ship = dir(self)
        
        # Do a barrel roll!
        self.t  = _np.roll(self.t,roll)
        for ti in range(roll):
            self.t[ti] = self.t[ti]-_dt.timedelta(days=1)
        self.ut_time = self.dn + (self.t-arbdate)
        self.T   = _np.roll(self.T ,roll)
        self.Te  = _np.roll(self.Te,roll)
        self.u   = _np.roll(self.u ,roll)
        self.ue  = _np.roll(self.ue,roll)
        self.v   = _np.roll(self.v ,roll)
        self.ve  = _np.roll(self.ve,roll)
        self.w   = _np.roll(self.w ,roll)
        self.we  = _np.roll(self.we,roll)
        self.i   = _np.roll(self.i ,roll)
        self.ie  = _np.roll(self.ie,roll)
        if 'uc' in ship:
            self.Tc   = _np.roll(self.Tc ,roll)
            self.uc   = _np.roll(self.uc ,roll)
            self.vc   = _np.roll(self.vc ,roll)
            self.wc   = _np.roll(self.wc ,roll)
            self.ic   = _np.roll(self.ic ,roll)
        if 'ud' in ship:
            self.Td   = _np.roll(self.Td ,roll)
            self.ud   = _np.roll(self.ud ,roll)
            self.vd   = _np.roll(self.vd ,roll)
            self.wd   = _np.roll(self.wd ,roll)
            self.id   = _np.roll(self.id ,roll)
        if 'uv' in ship:
            self.Tv   = _np.roll(self.Tv ,roll)
            self.uv   = _np.roll(self.uv ,roll)
            self.vv   = _np.roll(self.vv ,roll)
            self.wv   = _np.roll(self.wv ,roll)
            self.iv   = _np.roll(self.iv ,roll)
            self.Tve  = _np.roll(self.Tve,roll)
            self.uve  = _np.roll(self.uve,roll)
            self.vve  = _np.roll(self.vve,roll)
            self.wve  = _np.roll(self.wve,roll)
            self.ive  = _np.roll(self.ive,roll)
        if 'u2' in ship:
            self.u2   = _np.roll(self.u2 ,roll)
            self.u2e  = _np.roll(self.u2e,roll)
            self.v2   = _np.roll(self.v2 ,roll)
            self.v2e  = _np.roll(self.v2e,roll)
        if 'u2c' in ship:
            self.u2c  = _np.roll(self.u2c ,roll)
            self.v2c  = _np.roll(self.v2c ,roll)
        if 'u2d' in ship:
            self.u2d  = _np.roll(self.u2d ,roll)
            self.v2d  = _np.roll(self.v2d ,roll)
        if 'u2v' in ship:
            self.u2v  = _np.roll(self.u2v ,roll)
            self.v2v  = _np.roll(self.v2v ,roll)
            self.u2ve = _np.roll(self.u2ve,roll)
            self.v2ve = _np.roll(self.v2ve,roll)
        if 'uu' in ship:
            self.uu  = _np.roll(self.uu ,roll)
            self.Tu  = _np.roll(self.Tu ,roll)
            self.vu  = _np.roll(self.vu ,roll)
            self.wu  = _np.roll(self.wu ,roll)
            self.iu  = _np.roll(self.iu ,roll)
            self.uv2 = _np.roll(self.uv2,roll)
            self.Tv2 = _np.roll(self.Tv2,roll)
            self.vv2 = _np.roll(self.vv2,roll)
            self.wv2 = _np.roll(self.wv2,roll)
            self.iv2 = _np.roll(self.iv2,roll)
        if 'u2u' in ship:
            self.u2u = _np.roll(self.u2u,roll)
            self.v2u = _np.roll(self.v2u,roll)
            self.u2v2= _np.roll(self.u2v2,roll)
            self.v2v2= _np.roll(self.v2v2,roll)
        if 'u50' in ship:
            self.u16 = _np.roll(self.u16,roll)
            self.v16 = _np.roll(self.v16,roll)
            self.w16 = _np.roll(self.w16,roll)
            self.T16 = _np.roll(self.T16,roll)
            self.i16 = _np.roll(self.i16,roll)
            self.u25 = _np.roll(self.u25,roll)
            self.v25 = _np.roll(self.v25,roll)
            self.w25 = _np.roll(self.w25,roll)
            self.T25 = _np.roll(self.T25,roll)
            self.i25 = _np.roll(self.i25,roll)
            self.u50 = _np.roll(self.u50,roll)
            self.v50 = _np.roll(self.v50,roll)
            self.w50 = _np.roll(self.w50,roll)
            self.T50 = _np.roll(self.T50,roll)
            self.i50 = _np.roll(self.i50,roll)
            self.u75 = _np.roll(self.u75,roll)
            self.v75 = _np.roll(self.v75,roll)
            self.w75 = _np.roll(self.w75,roll)
            self.T75 = _np.roll(self.T75,roll)
            self.i75 = _np.roll(self.i75,roll)
            self.u84 = _np.roll(self.u84,roll)
            self.v84 = _np.roll(self.v84,roll)
            self.w84 = _np.roll(self.w84,roll)
            self.T84 = _np.roll(self.T84,roll)
            self.i84 = _np.roll(self.i84,roll)
        if 'u250' in ship:
            self.u216 = _np.roll(self.u216,roll)
            self.u225 = _np.roll(self.u225,roll)
            self.u250 = _np.roll(self.u250,roll)
            self.u275 = _np.roll(self.u275,roll)
            self.u284 = _np.roll(self.u284,roll)
            self.v216 = _np.roll(self.v216,roll)
            self.v225 = _np.roll(self.v225,roll)
            self.v250 = _np.roll(self.v250,roll)
            self.v275 = _np.roll(self.v275,roll)
            self.v284 = _np.roll(self.v284,roll)
        if 'alts' in ship:
            self.alts = _np.roll(self.alts,roll)

        self.barrelroll = not(self.barrelroll)


def Fmodel(alt):
    if alt == 'low':
        I = _np.array([  5.06418555e-03,   3.79121460e-01,   1.21130929e+01,
         4.33707174e+01,   6.76137501e+01,   6.46582112e+01,
         3.97066263e+01,   1.66496603e+01,   5.56633141e+00,
         1.65258386e+00,   4.65893631e-01,   1.29035828e-01,
         3.56951975e-02,   9.92605253e-03])
    
    else:
        I = _np.array([  1.44081159e-05,   9.37163934e-05,   3.63936919e-04,
         1.46931147e-02,   1.24256301e-01,   7.29099119e-01,
         2.92154760e+00,   5.00024489e+00,   4.64812608e+00,
         2.74066378e+00,   1.17773166e+00,   4.26615701e-01,
         1.43556192e-01,   4.63650114e-02])
        
    return(I)


if __name__=="__main__":

    print "Insert Coin"
