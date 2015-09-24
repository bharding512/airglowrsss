'''
Summary
-------
FPIResults contains functions to do both CV and binning/filtering daily/monthly/yearly analysis

Included functions are:
    BinDailyData
    BinMonthlyData
    CreateL1ASCII
    CreateL2ASCII
    CreateL2ASCII_Legacy
    CreateMonthlyASCII
    FilterData
    GetModels
    PlotAverages
    PlotClimatology
    PlotGridMonth
    PlotSpaghetti
    SetBinTime
    WeightedAverage
    laser_is_drifting
    
History
-------
3/20/13 -- Written by DJF (dfisher2@illinois.edu)
'''

import matplotlib as _mpl
import matplotlib.dates as _md
import matplotlib.pyplot as _plt
#matplotlib.use('AGG')
from pyglow import pyglow as _pyglow
import datetime as _dt
import calendar as _cal
import numpy as _np
from scipy import stats as _stats
import pytz as _pytz
import FPIprocessLevel2_Legacy as L2
import FPIprocessLevel2_Legacy as _L2L
import fpiinfo as _fpiinfo

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


# Set up default parameters
_mpl.rcParams.update({'font.size': 11})
dirout = '/rdata/airglow/database/L2/plots/'

# Binning Time
arbdate = _dt.datetime(1970,1,1)
_utc = _pytz.UTC
SetBinTime(30)



def FilterData(DATA):
    '''
    Summary:
        Returns filters single input day.

    Inputs:
        DATA = The data object

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)

    '''
    #winderrorlimit = 50.
    #temperrorlimit = 100.
    calerrorlimit = 300.
    cloudlimit = 2
    winderrorlimit = 25.
    temperrorlimit = 50.

    lim = {'Tmax':1400.,'Tmin':600.,'Te':temperrorlimit, \
           'umax':250.,'umin':-250.,'uef':winderrorlimit, 'uec':calerrorlimit,\
           'vmax':250.,'vmin':-250.,'vef':winderrorlimit, 'vec':calerrorlimit,\
           'wmax':75., 'wmin':-75., 'wef':winderrorlimit, 'wec':calerrorlimit}
    '''
    lim = {'Tmax':1250.,'Tmin':600.,'Te':temperrorlimit, \
           'umax':200.,'umin':-100.,'uef':winderrorlimit,'uec':calerrorlimit, \
           'vmax':150.,'vmin':-150.,'vef':winderrorlimit,'vec':calerrorlimit, \
           'wmax':75., 'wmin':-75., 'wef':winderrorlimit,'wec':calerrorlimit}
    '''
    bcnt = 0
    
    # For each direction in DATA
    for r1 in DATA:
        if len(r1.t1)>0:
            #r1 = data[link]
            #ind1 = range(0,len(r1.t1))
            ind1 = range(len(r1.t1))
            alld = len(ind1)

            # Filter using flags TODO: Fix with BJH
            if len(r1.flag_T) >= 1:
                ind1 = _np.delete(ind1, _np.where(r1.flag_T[ind1] >= cloudlimit))
                ind1 = _np.delete(ind1, _np.where(r1.flag_wind[ind1] >= cloudlimit))

            # Bad Data Filtering Temps
            if len(r1.T) >= 1:
                ind1 = _np.delete(ind1, _np.where(lim['Tmin'] > r1.T[ind1]))
                ind1 = _np.delete(ind1, _np.where(r1.T[ind1] > lim['Tmax']))
                ind1 = _np.delete(ind1, _np.where(_np.abs(r1.Te[ind1]) > lim['Te']))

            # Bad Data Filtering Winds -both limits, fits, & spikes(cal)
            if len(r1.u) >= 1:
                ind1 = _np.delete(ind1, _np.where(lim['umin'] > r1.u[ind1]))
                ind1 = _np.delete(ind1, _np.where(r1.u[ind1] > lim['umax']))
                ind1 = _np.delete(ind1, _np.where(r1.uef[ind1] > lim['uef']))
                ind1 = _np.delete(ind1, _np.where(r1.uec[ind1] > lim['uec']))
            if len(r1.v) >= 1:
                ind1 = _np.delete(ind1, _np.where(lim['vmin'] > r1.v[ind1]))
                ind1 = _np.delete(ind1, _np.where(r1.v[ind1] > lim['vmax']))
                ind1 = _np.delete(ind1, _np.where(r1.vef[ind1] > lim['vef']))
                ind1 = _np.delete(ind1, _np.where(r1.vec[ind1] > lim['vec']))
            if len(r1.w) >= 1:
                ind1 = _np.delete(ind1, _np.where(lim['wmin'] > r1.w[ind1]))
                ind1 = _np.delete(ind1, _np.where(r1.w[ind1] > lim['wmax']))
                ind1 = _np.delete(ind1, _np.where(r1.wef[ind1] > lim['wef']))
                ind1 = _np.delete(ind1, _np.where(r1.wec[ind1] > lim['wec']))

            r1.cut(arbdate,arbdate+_dt.timedelta(days=1),ind1)
            bcnt += alld - len(ind1)
    return(bcnt)
    


def WeightedAverage(VAL,STD,CNT=None,AXIS=0,test=False):
    '''
    Summary:
        Returns weighted mean and std of FPI data
        Returns weighted mean/std, and variability/std of FPI data if CNT is given

    Inputs:
        VAL = data to average
        STD = standard deviation of data
        CNT = count of data
        AXIS = axis to average across?

    Outputs:
        WM = weighted mean
        WE = weighted std
        SV = sample variance (variability)
        SE = std of sample variance (std of variability)
        AE = average uncertainty (average error/std)

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)
        9/23/14 -- Redid as separate function w/ corrected errors: DJF
    '''
    # Rotate averaging axis if deisred
    if AXIS==1:
        VAL = VAL.T
        STD = STD.T
        if CNT is not None:
            CNT = CNT.T
    
    # weighted mean and weighted std
    #mV = _np.ma.masked_array(VAL,_np.isnan(VAL))
    #mE = _np.ma.masked_array(STD,_np.isnan(STD))
    wt = (VAL/VAL)/STD**2
    V1 = _np.nansum(wt,axis=1)
    V2 = _np.nansum(wt**2,axis=1)
    WM = _np.nansum(VAL*wt,axis=1)/V1
    WE = _np.sqrt(_np.nansum(STD**2*wt**2,axis=1)/V1**2)
    if CNT == None:
        # Return wt mean and wt std only
        return(WM,WE)
    elif test:
        # return wt mean and std and monthly variability/std
        SV = _np.sqrt(_np.nansum(_np.subtract(VAL.T,WM)**2,axis=0)/(_np.array(CNT)-1.))
        
        WSV = _np.sqrt(_np.nansum(wt.T*_np.subtract(VAL.T,WM)**2,axis=0)/(V1-V2/V1))
        SE = 2.*WE**4/(_np.array(CNT)-1.)
        AE = _np.sqrt(_np.nansum(STD**2,axis=1)/CNT) 
        return(WM,WE,SV,SE,WSV,AE)
    else:
        # return wt mean and std and monthly variability/std
        SV = _np.sqrt(_np.nansum(_np.subtract(VAL.T,WM)**2,axis=0)/(_np.array(CNT)-1.))
        
        WSV = _np.sqrt(_np.nansum(wt.T*_np.subtract(VAL.T,WM)**2,axis=0)/(V1-V2/V1))
        SE = 2.*WE**4/(_np.array(CNT)-1.)
        return(WM,WE,SV,SE)
    
    

def BinDailyData(SITE,YEAR,DOY,SPLIT=False,KP=[0,10],CV=True):
    '''
    Summary:
        Returns filted and binned single data for a single instrument over one night.

    Inputs:
        SITE = site, e.g. uao
        YEAR = year, e.g. 2013
        DOY = doy of year, e.g. 47
        SPLIT = Split look directions in binning [default = False]
        KP = Filter days by KP [default = [0,10] - all kp]
        CV = Include CV directions [default = True]

    Outputs:
        DATA = Object with winds, Temps, and more

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)

    '''
    
    # Create the YYYYMMDD date format
    dn = _dt.datetime(YEAR,1,1) + _dt.timedelta(days = DOY-1)
    year = dn.strftime('%Y')
    date = dn.strftime('%Y%m%d')
    
    # Load in Day's Data
    if SITE in ['renoir','peru','nation']:
        nets = _fpiinfo.get_network_info(SITE).keys()
        tots = _fpiinfo.get_all_sites_info()
        sites = [x for x in nets if x in tots]
        project = SITE.lower()
    else:
        sites = [SITE.lower()]
        project = _fpiinfo.get_site_info(sites[0])['Network']
    data = L2.GetLevel2(project,dn)
    bc = FilterData(data)
    
    # Get Empty
    d = _BinnedData(dn,SITE)
    d.key = "Daily"
    d.t = times
    uData = _np.empty((b_len,500))*_np.nan
    ueData = _np.empty((b_len,500))*_np.nan
    vData = _np.empty((b_len,500))*_np.nan
    veData = _np.empty((b_len,500))*_np.nan
    u2Data = _np.empty((b_len,500))*_np.nan
    u2eData = _np.empty((b_len,500))*_np.nan
    v2Data = _np.empty((b_len,500))*_np.nan
    v2eData = _np.empty((b_len,500))*_np.nan
    wData = _np.empty((b_len,500))*_np.nan
    weData = _np.empty((b_len,500))*_np.nan
    TData = _np.empty((b_len,500))*_np.nan
    TeData = _np.empty((b_len,500))*_np.nan
    iData = _np.empty((b_len,500))*_np.nan
    ieData = _np.empty((b_len,500))*_np.nan
    count = _np.zeros((b_len))
    lat = []
    lon = []
    alt = []
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
                # Accumulate AIPP locations
                if len(r1.lla) == 3:
                    lat.append(r1.lla[0])
                    lon.append(r1.lla[1])
                    alt.append(r1.lla[2])
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
                    try:
                        pt = _pyglow.Point(r1.t1[zelda],0,0,0)
                        kpi = pt.kp
                    except:
                        kpi = 999.9 
                    if KP[0]<= kpi <=KP[1]:
                        bin = _np.floor(_np.mod((r1.t1[zelda].astimezone(_utc).replace(tzinfo=None)-arbdate).total_seconds(),60*60*24)/(60*24*60/b_len))
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
            
    # Mean LLA
    d.lla = _np.array([_np.nanmean(lat),_np.nanmean(lon),_np.nanmean(alt)])
    
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
    d.v = vD
    d.ve = veD
    d.w = wD
    d.we = weD
    d.T = TD
    d.Te = TeD
    d.i = iD
    d.ie = ieD
    if SPLIT:
        d.u2 = u2D
        d.u2e = u2eD
        d.v2 = v2D
        d.v2e = v2eD
    d.cards = cc
    d.cvs = cvc
    d.bads = bc
    d.doabarrelroll()
    
    return d
    
    
    
def GetModels(SITELLA,YEAR,DOY,WMODEL,TMODEL='msis'):
    '''
    Summary:
        Returns HWMfor a single instrument over one night.

    Inputs:
        SITELLA = site latitude, longitude, altitude
        YEAR = year, e.g. 2013
        DOY = doy of year, e.g. 47
        WMODEL = name of wind model, e.g. hwm93
        TMODEL = name of temp model, e.g. msis

    Outputs:
        DATA = Object with winds, Temps, and more

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)
    '''
    
    # Create the YYYYMMDD date format
    dn = _dt.datetime(YEAR,1,1) + _dt.timedelta(days = DOY-1)
    year = dn.strftime('%Y')
    date = dn.strftime('%Y%m%d')
    
    uData = _np.empty((b_len,1))*_np.nan
    ueData = _np.empty((b_len,1))*_np.nan
    vData = _np.empty((b_len,1))*_np.nan
    veData = _np.empty((b_len,1))*_np.nan
    wData = _np.empty((b_len,1))*_np.nan
    weData = _np.empty((b_len,1))*_np.nan
    TData = _np.empty((b_len,1))*_np.nan
    TeData = _np.empty((b_len,1))*_np.nan
    iData = _np.empty((b_len,1))*_np.nan
    ieData = _np.empty((b_len,1))*_np.nan
    
    # Get Empty
    d = _BinnedData(dn,WMODEL)
    d.key = 'DailyModel'
    d.t = times
    d.lla = SITELLA
    if _np.nan in SITELLA:
        print 'Bad LLA'
        return
    
    # Fill Data
    for tind,t in enumerate(times):
        pt = _pyglow.Point(t.replace(year=dn.year,month=dn.month,day=dn.day),SITELLA[0],SITELLA[1],SITELLA[2])
        # Wind
        if WMODEL.lower() == 'hwm93':
            pt.run_hwm93()
        elif WMODEL.lower() == 'hwm07':
            pt.run_hwm07()
        elif WMODEL.lower() == 'hwm14':
            pt.run_hwm14()
        else:
            print 'Bad Wind Model'
        uData[tind] = pt.u
        ueData[tind] = 1.
        vData[tind] = pt.v
        veData[tind] = 1.
        #wData[tind] = pt.w
        #weData[tind] = 1.
        
        # Temp
        if TMODEL.lower() == 'msis':
            pt.run_msis()
            TData[tind] = pt.Tn_msis
        elif WMODEL.lower() == 'iri':
            pt.run_iri()
            TData[tind] = pt.Tn_iri
        else:
            print 'Bad Temp Model'
        TeData[tind] = 1.
    
        # Intensity
        pt.run_airglow()
        iData[tind] = pt.ag6300
        ieData[tind] = 1.

    # Save Averages
    d.u = uData[:,0]
    d.ue = ueData[:,0]
    d.v = vData[:,0]
    d.ve = veData[:,0]
    d.w = wData[:,0]
    d.we = weData[:,0]
    d.T = TData[:,0]
    d.Te = TeData[:,0]
    d.i = iData[:,0]
    d.ie = ieData[:,0]
    d.doabarrelroll()
    
    return d
        
        

def BinMonthlyData(SITE,YEAR,MONTH,SPLIT=False,SITELLA=[_np.nan,_np.nan,_np.nan],DLIST=[],YLIST=[],KP=[0,10],CV=True):
    '''
    Summary:
        Returns filted and binned data over one month.

    Inputs:
        SITE = site of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTH = month of year, e.g. 2 (February)
        SPLIT = Split look directions in binning [default = False]
        SITELLA = Lat, Lon, Alt for HWM [default = nan's]
        DLIST = list of doys in year  [default = [] - all doys in MONTH,YEAR used]
        YLIST = list of years in year [default = [] - only YEAR used]
        KP = limits of kp for filtering days [default = [0,10] - all kp used]
        CV = use CV modes [default = True]

    Outputs:
        DATA = dictionary of data whose keys are Zonal, Meridional, or Temp 
               Temp contains Temp and Temp_Error (averaged from all directions)
               Zonal/Meridional contains Wind and Wind_Error

    History:
        3/20/13 -- Written by DJF (dfisher2@illinois.edu)

    '''
    # Define Outputs...
    dn = _dt.datetime(YEAR,MONTH,1)
    mon = _cal.month_name[MONTH]

    # Set Output Variable
    dimset = 3  # Require x days in each month for an variability set
    d = _BinnedData(dn,SITE)
    d.t = times
    d.key = '{0:%B}'.format(dn, "month")
    
    # Get Empty 
    dim = 31*5*10
    uData = _np.empty((b_len,dim))*_np.nan
    ueData = _np.empty((b_len,dim))*_np.nan
    vData = _np.empty((b_len,dim))*_np.nan
    veData = _np.empty((b_len,dim))*_np.nan
    u2Data = _np.empty((b_len,dim))*_np.nan
    u2eData = _np.empty((b_len,dim))*_np.nan
    v2Data = _np.empty((b_len,dim))*_np.nan
    v2eData = _np.empty((b_len,dim))*_np.nan
    wData = _np.empty((b_len,dim))*_np.nan
    weData = _np.empty((b_len,dim))*_np.nan
    TData = _np.empty((b_len,dim))*_np.nan
    TeData = _np.empty((b_len,dim))*_np.nan
    iData = _np.empty((b_len,dim))*_np.nan
    ieData = _np.empty((b_len,dim))*_np.nan
    F107 = _np.empty((dim))*_np.nan
    count = 0
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
    else:
        # make sure doys actually exist in the month
        dl = []; yl = []
        for doy,yr in zip(DLIST,YLIST):
            doystart = (_dt.datetime(yr,MONTH,1) - _dt.datetime(yr,1,1)).days+1
            doyend = doystart + _cal.monthrange(yr,MONTH)[1]
            if doystart <= doy and doy <=doyend:
                dl.append(doy)
                yl.append(yr)
        if not(dl):
            print 'Doy List contains no days in desired month.'

    for doy,yr in zip(dl,yl):
        #print doy
        
        if SITE in ['renoir','peru','nation']:
            nets = _fpiinfo.get_network_info(SITE).keys()
            tots = _fpiinfo.get_all_sites_info()
            sites = [x for x in nets if x in tots]
        else:
            sites = [SITE]

        for s in sites:
            if 'hwm' in s:
                DD = GetModels(SITELLA,YEAR,doy,s)
                mflag = True

                # get F107 weighted at midnight of data (assume constant for night)
                point = _pyglow.Point(DD.dn,0,0,250)
                F107[count] = (point.f107 + point.f107a)/2.

            else:
                DD = BinDailyData(s,yr,doy,SPLIT,KP,CV)
                mflag = False
                cards += DD.cards
                cvs += DD.cvs
                bads += DD.bads

                # get F107 weighted at midnight of data (assume constant for night)
                point = _pyglow.Point(DD.dn,0,0,250)
                F107[count] = (point.f107 + point.f107a)/2.*(DD.cards+DD.cvs)

            # Undo shift for easy averaging
            DD.doabarrelroll()
            # Debug, Overplot all Horizontal winds.
            #_plt.plot(DD.t,DD.u)
            
            # Count days total used
            if sum(_np.isfinite(DD.T)):
                oscar.append(doy)
            
            if len(DD.lla) == 3:
                lat.append(DD.lla[0])
                lon.append(DD.lla[1])
                alt.append(DD.lla[2])
            # Add data
            if len(DD.u) > 0:
                uData[:,count] = DD.u
                ueData[:,count] = DD.ue
            if len(DD.v) > 0:
                vData[:,count] = DD.v
                veData[:,count] = DD.ve
            if len(DD.w) > 0:
                wData[:,count] = DD.w
                weData[:,count] = DD.we
            if len(DD.T) > 0:
                TData[:,count] = DD.T
                TeData[:,count] = DD.Te
            if len(DD.i) > 0:
                iData[:,count] = DD.i
                ieData[:,count] = DD.ie
            if SPLIT and not(mflag) and len(DD.u2) > 0:
                u2Data[:,count] = DD.u2
                u2eData[:,count] = DD.u2e
            if SPLIT and not(mflag) and len(DD.v2) > 0:
                v2Data[:,count] = DD.v2
                v2eData[:,count] = DD.v2e

            count += 1
    
    # Get count of days used in each bin
    ucount = dim - sum(_np.isnan(uData.T))
    ucount = [_np.nan if x<dimset else x for x in ucount]
    vcount = dim - sum(_np.isnan(vData.T))
    vcount = [_np.nan if x<dimset else x for x in vcount]
    wcount = dim - sum(_np.isnan(wData.T))
    wcount = [_np.nan if x<dimset else x for x in wcount]
    Tcount = dim - sum(_np.isnan(TData.T))
    Tcount = [_np.nan if x<dimset else x for x in Tcount]
    icount = dim - sum(_np.isnan(iData.T))
    icount = [_np.nan if x<dimset else x for x in icount]
    u2count = dim - sum(_np.isnan(u2Data.T))
    u2count = [_np.nan if x<dimset else x for x in u2count]
    v2count = dim - sum(_np.isnan(v2Data.T))
    v2count = [_np.nan if x<dimset else x for x in v2count]
    '''
    ## Weighted Mean & Statistical Variance of Winds
    # Zonal
    uD,uDe,uV,uVe = WeightedAverage(uData,ueData,ucount)
    # Meridional
    vD,vDe,vV,vVe = WeightedAverage(vData,veData,vcount)
    # Vert            
    wD,wDe,wV,wVe = WeightedAverage(wData,weData,wcount)
    #  Temps
    TD,TDe,TV,TVe = WeightedAverage(TData,TeData,Tcount)
    if SPLIT and not(mflag):
        # Zonal2 - West
        u2D,u2De,u2V,u2Ve = WeightedAverage(u2Data,u2eData,u2count)
        # Meridional2 - South
        v2D,v2De,v2V,v2Ve = WeightedAverage(v2Data,v2eData,v2count)
        
    # Save Averages
    d.lla = _np.array([_np.nanmean(lat),_np.nanmean(lon),_np.nanmean(alt)])
    d.u  = uD
    d.ue = uDe
    d.uv = uV
    d.uve= uVe
    d.uc = ucount
    d.v  = vD
    d.ve = vDe
    d.vv = vV
    d.vve= vVe
    d.vc = vcount
    d.w  = wD
    d.we = wDe
    d.wv = wV
    d.wve= wVe
    d.wc = wcount
    d.T  = TDdoabarr
    d.Te = TDe
    d.Tv = TV
    d.Tve= TVe
    d.Tc = Tcount
    if SPLIT and not(mflag):
        d.u2  = u2D
        d.u2e = u2De
        d.u2v = u2V
        d.u2ve= u2Ve
        d.u2c = u2count
        d.v2  = v2D
        d.v2e = v2De
        d.v2v = v2V
        d.v2ve= v2Ve
        d.v2c = v2count
    d.cards = cards
    d.cvs = cvs
    d.daysused = _np.unique(oscar)
    '''
    
    # Zonal
    uD,uDe,uV,uVe,uV2,uU = WeightedAverage(uData,ueData,ucount,test=True)
    # Meridional
    vD,vDe,vV,vVe,vV2,vU = WeightedAverage(vData,veData,vcount,test=True)
    # Vert            
    wD,wDe,wV,wVe,wV2,wU = WeightedAverage(wData,weData,wcount,test=True)
    #  Temps
    TD,TDe,TV,TVe,TV2,TU = WeightedAverage(TData,TeData,Tcount,test=True)
    # Intensity
    iD,iDe,iV,iVe,iV2,iU = WeightedAverage(iData,ieData,icount,test=True)
    if SPLIT and not(mflag):
        # Zonal2 - West
        u2D,u2De,u2V,u2Ve,u2V2,u2U = WeightedAverage(u2Data,u2eData,u2count,test=True)
        # Meridional2 - South
        v2D,v2De,v2V,v2Ve,v2V2,v2U = WeightedAverage(v2Data,v2eData,v2count,test=True)
    # F107
    if 'hwm' in SITE:
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
    d.uc = ucount
    d.v  = vD
    d.ve = vDe
    d.vv = vV
    d.vve= vVe
    d.vv2= vV2
    d.vu = vU
    d.vc = vcount
    d.w  = wD
    d.we = wDe
    d.wv = wV
    d.wve= wVe
    d.wv2= wV2
    d.wu = wU
    d.wc = wcount
    d.T  = TD
    d.Te = TDe
    d.Tv = TV
    d.Tve= TVe
    d.Tu = TU
    d.Tc = Tcount
    d.Tv2= TV2
    d.i  = iD
    d.ie = iDe
    d.iv = iV
    d.ive= iVe
    d.iv2= iV2
    d.iu = iU
    d.ic = icount
    
    if SPLIT and not(mflag):
        d.u2  = u2D
        d.u2e = u2De
        d.u2v = u2V
        d.u2ve= u2Ve
        d.u2v2= u2V2
        d.u2u = u2U
        d.u2c = u2count
        d.v2  = v2D
        d.v2e = v2De
        d.v2v = v2V
        d.v2ve= v2Ve
        d.v2v2= v2V2
        d.v2u = v2U
        d.v2c = v2count
    d.cards = cards
    d.cvs = cvs
    d.daysused = _np.unique(oscar)
    d.doabarrelroll()

    #print '%s %s %04d'%(SITE, mon,YEAR),'Days used:',len(d.daysused),' Ave pts/time: %02.2f'%(_np.nanmin(_np.array([_np.ma.masked_array(ucount,_np.isnan(ucount)).mean(), _np.ma.masked_array(vcount,_np.isnan(vcount)).mean(), _np.ma.masked_array(Tcount,_np.isnan(Tcount)).mean()])))
    print '%s %s %04d'%(SITE, mon,YEAR),' Days used:',len(d.daysused)
    if SPLIT and not(mflag):
        print 'Ave pts/time: u-%02.2f  u2-%02.2f  v-%02.2f  v2-%02.2f  w-%02.2f  T-%02.2f'%(_np.nanmean(ucount),_np.nanmean(u2count),_np.nanmean(vcount),_np.nanmean(v2count),_np.nanmean(wcount),_np.nanmean(Tcount))

    else:
        print 'Ave pts/time: u-%02.2f  v-%02.2f  w-%02.2f  T-%02.2f'%(_np.nanmean(ucount),_np.nanmean(vcount),_np.nanmean(wcount),_np.nanmean(Tcount))
    print 'Cards:',cards,' CVs:',cvs,' %%CV: %02.2f'%(100.*cvs/(cvs+cards+.000001))
    print 'Bads:',bads, ' %%Good: %02.2f'%(100.*(cards+cvs)/(cards+cvs+bads+.000001))
    print 'Days:',d.daysused,'\n'

    return d



def PlotClimatology(SITE,YEAR,MONTHSTART=1,NMONTHS=12,SPLIT=False,KP=[0,10],UT=True):
    '''
    Summary:
        Plots monthly averages in a 2x6 month plot
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTHSTART = month to start yearly climatology [default = 1 - Jan]
        NMONTHS = number of months desired from MONTHSTART [default = 12]
        SPLIT = Split look directions in binning [default = False]
        KP = Filter days by KP [default = [0,10] - all kp]
        UT = Plot in UT or SLT [default = True]

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
        print nsp,year,month
            
        MD = BinMonthlyData(SITE,year,month,SPLIT,KP)
        tlim = [MD.t[len(MD.t)/5],MD.t[-len(MD.t)/5]]
        MD.t  = MD.t + _dt.timedelta(seconds=ut_offset)
        MD.t2 = MD.t + _dt.timedelta(minutes=3)
        gflag = 'on'
        
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
        axz[nsp].set_xlim(tlim)
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
        axm[nsp].set_xlim(tlim)
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
        axv[nsp].set_xlim(tlim)
        axv[nsp].set_yticks([-25, 0, 25])
        axv[nsp].errorbar(MD.t+_dt.timedelta(minutes=3*(year-YEAR)),MD.w,yerr=MD.wv,fmt=colors[year-YEAR],linewidth=lwide,elinewidth=lwide,capsize=csize,markersize=msize,label='Data')
        axv[nsp].grid(gflag)
        axv[nsp].plot([MD.t[0],MD.t[-1]],[0,0],'k--',label=None)
        
        # Temps
        axt[nsp].set_ylabel("%s" % (_cal.month_abbr[month]),labelpad=2)
        axt[nsp].set_ylim(600, 1200)
        axt[nsp].set_xlim(tlim)
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
        

    
def PlotAverages(SITE,YEAR,MONTH,MODEL=[],SPLIT=False,KP=[0,10],UT=True):
    '''
    Summary:
        Plots single site monthly averages w/ models
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTH = month to plot, e.g. 4
        MODEL = list of model keys to plot with data [default is empty]
        SPLIT = Split look directions in binning [default = False]
        KP = Filter days by KP [default = [0,10] - all kp]
        UT = Plot in UT or SLT [default = True]

    History:
        5/8/14 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # Get UT offset if SLT needed
    ut_offset = 0
    if not(UT):
        ut_offset = _fpiinfo.get_site_info(SITE)['Location'][1]/360.*(24*60*60)
        
    # Get Monthly Average for basis
    MD = BinMonthlyData(SITE,YEAR,MONTH,SPLIT,KP)
    mon = _cal.month_name[MONTH]
    MD.t  = MD.t + _dt.timedelta(seconds=ut_offset)
    MD.t2 = MD.t + _dt.timedelta(minutes=3) 

    # Get models
    tm = {}; Tm = {}; Tem = {}; Um = {}; Uem = {}; Vm = {}; Vem = {}; Wm = {}; Wem = {};
    if isinstance(MODEL, str):
        MODEL = [MODEL]
    for nm,m in enumerate(MODEL):
        if _np.nan in MD.lla:
            break
        MM = BinMonthlyData(m,YEAR,MONTH,SPLIT,SITELLA=MD.lla)
        tm[m] = MM.t + _dt.timedelta(minutes=3*(nm+2)) + _dt.timedelta(seconds=ut_offset)
        if not(UT):
            tm[m] = tm[m] + _dt.timedelta
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
    markerwide = 0
    try:
        name = _fpiinfo.get_site_info(SITE)['Name']
    except:
        name = SITE.title()
    # Temp Figure
    fig = _plt.figure(0,figsize=(10,6)); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    (_, caps1, _) = _plt.errorbar(MD.t,MD.T,yerr=MD.Tv,label='Data')
    # Hardcoded Tmodel -- Fix
    for nm,m in enumerate(MODEL):
        if nm == 0:
            (_, caps2, _) = _plt.errorbar(tm[m],Tm[m],yerr=Tem[m],fmt=color[nm],label='MSIS')
            for cap in caps2:
                cap.set_markeredgewidth(markerwide)
    for cap in caps1:
        cap.set_markeredgewidth(markerwide)
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
    _plt.draw();# _plt.show()
    #_plt.savefig('%s%s_%04d-%s_temps.png' % (dirout,SITE,YEAR,mon))
    
    # Winds Figure Zonal
    fig = _plt.figure(1,figsize=(10,6)); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    if SPLIT:
        (_, caps1, _) = _plt.errorbar(MD.t,MD.u,yerr=MD.uv,fmt='b-',label='East')
        (_, caps3, _) = _plt.errorbar(MD.t2,MD.u2,yerr=MD.u2v,fmt='c-',label='West')
        for cap in caps3:
            cap.set_markeredgewidth(markerwide)
    else:
        (_, caps1, _) = _plt.errorbar(MD.t,MD.u,yerr=MD.uv,label='Data')
    for nm,m in enumerate(MODEL):
        (_, caps2, _) = _plt.errorbar(tm[m],Um[m],yerr=Uem[m],fmt=color[nm],label=m.upper())
        for cap in caps2:
            cap.set_markeredgewidth(markerwide)
    for cap in caps1:
        cap.set_markeredgewidth(markerwide)
    _plt.plot([MD.t[0],MD.t[-1]],[0,0],'k--')
    _plt.ylim([-200,200])
    _plt.xlim([MD.t[tstart],MD.t[tend]])
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.ylabel('Wind Speed [m/s]')
    if UT:
        _plt.xlabel('Time [UT]')
    else:
        _plt.xlabel('Time [SLT]')
    _plt.title('%04i %s Average Zonal Winds at %s' % (YEAR,mon,name))
    _plt.legend()
    _plt.grid()
    _plt.draw(); #_plt.show()
    _plt.savefig('%s%s_%04d-%s_zonal_winds.png' % (dirout,SITE,YEAR,mon))
    
    # Winds Figure Meridional
    fig = _plt.figure(2,figsize=(10,6)); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    if SPLIT:
        (_, caps1, _) = _plt.errorbar(MD.t,MD.v,yerr=MD.vv,fmt='b-',label='North')
        (_, caps3, _) = _plt.errorbar(MD.t2,MD.v2,yerr=MD.v2v,fmt='c-',label='South')
        for cap in caps3:
            cap.set_markeredgewidth(markerwide)
    else:
        (_, caps1, _) = _plt.errorbar(MD.t,MD.v,yerr=MD.vv,label='Data')
    for nm,m in enumerate(MODEL):
        (_, caps2, _) = _plt.errorbar(tm[m],Vm[m],yerr=Vem[m],fmt=color[nm],label=m.upper())
        for cap in caps2:
            cap.set_markeredgewidth(markerwide)
    for cap in caps1:
        cap.set_markeredgewidth(markerwide)
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
    _plt.draw(); #_plt.show()
    #_plt.savefig('%s%s_%04d-%s_meridional_winds.png' % (dirout,SITE,YEAR,mon))
    
    # Winds Figure Vertical
    fig = _plt.figure(3,figsize=(10,6)); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    (_, caps1, _) = _plt.errorbar(MD.t,MD.w,yerr=MD.wv,label='Data')
    for nm,m in enumerate(MODEL):
        (_, caps2, _) = _plt.errorbar(tm[m],Wm[m],yerr=Wem[m],fmt=color[nm],label=m.upper())
        for cap in caps2:
            cap.set_markeredgewidth(markerwide)
    for cap in caps1:
        cap.set_markeredgewidth(markerwide)
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
    _plt.draw(); #_plt.show()
    #_plt.savefig('%s%s_%04d-%s_vertical_winds.png' % (dirout,SITE,YEAR,mon))
    


def PlotSpaghetti(SITE,YEAR,MONTH,SPLIT=False,LIST=[],CV=False,KP=[0,10]):
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

    History:
        10/17/14 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # Get Monthly Average for basis
    D = {}
    dn = _dt.datetime(YEAR,MONTH,1)
    if not LIST:
        for k,d in enumerate(range(_cal.monthrange(YEAR,MONTH)[1])):
            L = L2.GetLevel2(SITE,dn+_dt.timedelta(days=d))
            FilterData(L); D[k] = L
    else:
        for k,d in enumerate(LIST):
            L = L2.GetLevel2(SITE,_dt.datetime(YEAR,1,1)+_dt.timedelta(days=d-1))
            FilterData(L); D[k] = L
    M = BinMonthlyData(SITE,YEAR,MONTH,SPLIT,CV,LIST,KP)
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
    
    # finalize Temps
    f.subplots_adjust(hspace=0)
    ax[3].set_ylim([500,1000])
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



def PlotGridMonth(SITE,YEAR,MONTH,SPLIT=True):
    '''
    Summary:
        Plots all raw cardinal data for one month in 5x7 plot! 
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTH = month to plot, e.g. 4
        SPLIT = Split look directions in binning [default = True]

    History:
        10/17/14 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # Get Monthly Average for basis
    D = {}
    dn = _dt.datetime(YEAR,MONTH,1)
    for k,d in enumerate(range(_cal.monthrange(YEAR,MONTH)[1])):
        D[k] = L2.GetLevel2(SITE,dn+_dt.timedelta(days=d))
    #M = BinMonthlyData(SITE,YEAR,MONTH,SPLIT)
    #M.t = M.t + _dt.timedelta(seconds=(_dt.datetime(YEAR-1,1,1)-_dt.datetime(M.t[0].year,1,1)).total_seconds())

    markerwide = 0
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
            (l1, caps1, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,fmt='r.')
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
            m1, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'r-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('South' in y.key)]:
            (l2, caps2, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,fmt='g.')
            for cap in caps2:
                cap.set_markeredgewidth(markerwide)
            m2, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'g-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('East' in y.key)]:
            (l3, caps3, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,fmt='y.')
            for cap in caps3:
                cap.set_markeredgewidth(markerwide)
            m3, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'y-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('West' in y.key)]:
            (l4, caps4, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,fmt='m.')
            for cap in caps4:
                cap.set_markeredgewidth(markerwide)
            m4, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'m-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('Zenith' in y.key)]:
            (l5, caps5, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,yerr=l.Te,fmt='b.')
            for cap in caps5:
                cap.set_markeredgewidth(markerwide)
            m5, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.T,'b-',linewidth=lw,alpha=lalpha)
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
            (l1, caps1, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,yerr=l.ue,fmt='yo')
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
            m1, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,'y-',linewidth=lw,alpha=lalpha,label='Doy %03i'%doy)
        for l in [y for y in D[d] if ('West' in y.key)]:
            (l2, caps2, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,yerr=l.ue,fmt='mo')
            for cap in caps2:
                cap.set_markeredgewidth(markerwide)
            m2, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.u,'m-',linewidth=lw,alpha=lalpha,label='Doy %03i'%doy)
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
            (l1, caps1, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,yerr=l.ve,fmt='ro')
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
            m1, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'r-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('South' in y.key)]:
            (l2, caps2, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,yerr=l.ve,fmt='go')
            for cap in caps2:
                cap.set_markeredgewidth(markerwide)
            m2, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'g-',linewidth=lw,alpha=lalpha)
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
            (l1, caps1, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.we,fmt='bo')
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
            m1, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,'b-',linewidth=lw,alpha=lalpha)
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
            (l1, caps1, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.ve,fmt='ro')
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
            m1, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'r-',linewidth=lw,alpha=lalpha)
        for l in [y for y in D[d] if ('CV' in y.key and '_2' in y.key)]:
            (l2, caps2, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.ve,fmt='go')
            for cap in caps2:
                cap.set_markeredgewidth(markerwide)
            m2, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'g-',linewidth=lw,alpha=lalpha)
        '''
        for l in [y for y in D[d] if ('IN' in y.key)]:
            (l3, caps3, _) = ax[d].errorbar(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,yerr=l.ve,fmt='bo')
            for cap in caps3:
                cap.set_markeredgewidth(markerwide)
            m3, = ax[d].plot(_np.array([x.astimezone(_utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.v,'b-',linewidth=lw,alpha=lalpha)
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


       
def CreateL1ASCII(NPZ,OUT):
    '''
    Summary:
        Script to read in processed FPI npz file and output a text file of results.
    
    Inputs:
        NPZ = full file name and path of npz file 
        OUT = full file name and path of ASCII results output

    History:
        2/21/13 -- Written by Daniel J. Fisher (dfisher2@illinois.edu)
        2/16/15 -- Modified for complete npz file (no 
    '''

    # Load in FPI processed Data
    data = _np.load(NPZ)
    FPI_Results = data['FPI_Results'].ravel()[0]
    reference = data['FPI_Results'].ravel()[0]['reference']
    instr = data['instrument'].ravel()[0]['Abbreviation']
    direction = data['FPI_Results'].ravel()[0]['direction']
    temps = data['FPI_Results'].ravel()[0]['T']
    e_temps = data['FPI_Results'].ravel()[0]['sigma_T']
    winds = data['FPI_Results'].ravel()[0]['LOSwind']
    e_winds = data['FPI_Results'].ravel()[0]['sigma_LOSwind']
    e_fit = data['FPI_Results'].ravel()[0]['sigma_fit_LOSwind']
    e_cal = data['FPI_Results'].ravel()[0]['sigma_cal_LOSwind']
    i = data['FPI_Results'].ravel()[0]['skyI']
    e_i = data['FPI_Results'].ravel()[0]['sigma_skyI']
    b = data['FPI_Results'].ravel()[0]['ccdB']
    e_b = data['FPI_Results'].ravel()[0]['sigma_ccdB']
    az = data['FPI_Results'].ravel()[0]['az']
    ze = data['FPI_Results'].ravel()[0]['ze']
    #laser_t = data['FPI_Results'].ravel()[0]['laser_times']
    #laser_v = data['FPI_Results'].ravel()[0]['laser_value']
    #laser_chi = data['FPI_Results'].ravel()[0]['laser_chisqr']
    chisqr = data['FPI_Results'].ravel()[0]['sky_chisqr']
    wavelength = data['FPI_Results'].ravel()[0]['lam0']
    version = data['FPI_Results'].ravel()[0]['SVNRevision']
    inttime = data['FPI_Results'].ravel()[0]['sky_intT']
    timeywimey = data['FPI_Results'].ravel()[0]['sky_times']
    if(data['FPI_Results'].ravel()[0]['Clouds']):
        sky_temp = data['FPI_Results'].ravel()[0]['Clouds']['mean']
    else:
        sky_temp = _np.ones(len(timeywimey))*-999.
    wind_flag = data['FPI_Results'].ravel()[0]['wind_quality_flag']
    temp_flag = data['FPI_Results'].ravel()[0]['temp_quality_flag']

    del data.f # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
    data.close()

    # Write out ASCII
    with open(OUT,'w') as note:

        note.write('LEVEL 1 DATA PRODUCT:\n---------------------\n')
        '''
        note.write('VARIABLES:\n
            UTCTime - Start time of image in UTC
            Az - Azimuth angle in degrees (compass coordinates)
            Ze - Zenith angle in degrees (0 is zenith)
            Temp - Estimated temperature of neutral layer in K. Biases in temperature may exist for different instruments
            Temp_Sig - Estimated uncertainty of temperature estimate in K
            LOS_Wind - Estimated line-of-sight winds of measurement in m/s (+ away from instrument). This is NOT projected
            LOS_Wind_Sig - Estimated uncertainty of wind estimate in m/s
            Fit_Sig - Estimated uncertainty of wind due to LM fit of sky data in m/s
            Cal_Sig - Estimated uncertainty of wind due to laser calibration accuracy in m/s
            I - Estimated airglow intensity in arbitrary units
            I_Sig - Estimated uncertainty of intensity estimate
            Bkgd - Estimated Background intensity of CCD 
            Bkgd_Sig - Estimated uncertainty of background estimate
            Int_Time - Length of exposure for the measurement in s
            Chisqr - Chi-squared of the model fit
            Cld_Ind - Cloudiness indicator: ambient minus sky/cloud temperature in C. If less than -25 skies are assumed clear (cloudy skies are warmer than clear skies).  -999.9 indicates no data available and thus assumes good skies
            T_Flag - Temperature error flag: 2 is bad, 1 known issues/ iffy data, 0 is good
            W_Flag - Wind error flag: 2 is bad, 1 known issues/ iffy data, 0 is good
            Ref - How Doppler reference was calculated: Laser uses laser images to calibrate the doppler offset assuming that on average nighttime veritcal winds are zero, Zenith assumes zenith winds are zero and uses this as a Doppler zero
            Wl - Wavelength of measured emission line in m
            Vers - Current version of python analysis code (to verify up-to-date product)\n')
        note.write('NOTES:\nAssumed emission altitude of 250 km\n')
        
        note.write('\n---------------------\nData:\n')
        '''
        note.write('UTCTime  Az  Ze  Temp  Temp_Sig  LOS_Wind  LOS_Wind_Sig  Fit_Sig  Cal_Sig  I  I_Sig  Bkgd  Bkgd_Sig  Int_Time  Chisqr  Cld_Ind  T_Flag  W_Flag  Ref  Wl  Vers\n')

        for a_tw, a_az, a_ze, a_t, a_e_t, a_w, a_e_w, a_ef, a_ec, a_i, a_e_i, a_b, a_e_b, a_it, a_cs, a_dt,a_tf,a_wf in zip(timeywimey, az, ze, temps, e_temps, winds, e_winds, e_fit, e_cal, i, e_i, b, e_b, inttime, chisqr, sky_temp,temp_flag,wind_flag):
            dn = a_tw.astimezone(_utc)
            utctime = dn.strftime("%Y-%m-%d %H:%M:%S")
            line = "{:19s}  {:6.1f}  {:5.1f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:6.2f}  {:6.2f}  {:6.1f}  {:4.2f}  -999.00  -999  {:3.0f}  {:6.2f}  {:6.1f}  {:1.0f}  {:1.0f}  {:6s}  {:7.1e}  {:5s}\n".format(utctime, a_az, a_ze, a_t, a_e_t, a_w, a_e_w, a_ef, a_ec, a_i, a_e_i, a_it, a_cs, a_dt, a_tf, a_wf, reference, wavelength, version)
            #line = "{:19s}  {:6.1f}  {:5.1f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:6.1f}  {:4.2f}  {:6.1f}  {:4.2f}  {:3.0f}  {:6.2f}  {:6.1f}  {:1d}  {:1d}  {:6s}  {:7.1e}  {:5s}\n".format(utctime, a_az, a_ze, a_t, a_e_t, a_w, aax.xaxis_date()_e_w, a_i, a_e_i, a_b, a_e_b, a_it, a_cs, a_dt, t_flag, w_flag, reference, wavelength, version)
            note.write(line)
	        
    note.closed



def CreateL2ASCII(PROJECT,YEAR,DOY):
    '''
    Summary:
        Script to save out ASCII of all L2 date for a Project.  Filtered, w=0 assumed.
    
    Inputs:
        PROJECT - Name of project for sites, ex: 'NATION'
        YEAR = year, e.g. 2013
        DOY = day of year, e.g. 64

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
    D = L2.GetLevel2(PROJECT,process_dn)
    #FilterData(D)

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
	                line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  250.0  -999.00  -999.00  {:7.2f}  {:7.2f}  {:7.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.v[i], x.ve[i], x.cloudy[i])
	            elif('East' in x.key or 'West' in x.key):
	                line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  250.0  {:7.2f}  {:7.2f}  -999.00  -999.00  {:7.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.u[i], x.ue[i], x.cloudy[i])
	            elif('CV_VTI_EKU_PAR' in x.key):
	                line = ""
	            elif('CV' in x.key):
	                line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  250.0  {:7.2f}  {:7.2f}  {:7.2f}  {:7.2f}  {:7.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.u[i], x.ue[i], x.v[i], x.ve[i], x.cloudy[i])
	            # For debugg
	            elif('Zenith' in x.key or 'IN' in x.key):
	                line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  250.0  VERTICAL  {:7.2f}  {:7.2f}  {:7.2f}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.w[i], x.we[i], x.cloudy[i])
	            else:
	                line = ""
	            note.write(line)
		    
    note.close()
    
    print 'Results saved to: ' + notename
    print 'This is a test case, not for official L2 use'
    
    

def CreateL2ASCII_Legacy(PROJECT,YEAR,DOY):
    '''
    Summary:
        Script to save out Legacy ASCII of all L2 date for a Project.
    
    Inputs:
        PROJECT - Name of project for sites, ex: 'NATION'
        YEAR = year, e.g. 2013
        DOY = day of year, e.g. 64

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
    D = _L2L.GetLevel2(PROJECT,process_dn)
    FilterData(D)

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
		        line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  -------  ------  -------  ------  {:7.2f}  {:6.2f}  {:6.1f}  {:4.2f}  {:6.1f}  {:4.2f}  {:30s}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.w[i], x.we[i], x.i[i], x.ie[i], x.b[i], x.be[i], x.notes)
		    elif('North' in x.key or 'South' in x.key):
		        line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  -------  ------  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:6.1f}  {:4.2f}  {:6.1f}  {:4.2f}  {:30s}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.v[i], x.ve[i], x.w[i], x.we[i], x.i[i], x.ie[i], x.b[i], x.be[i], x.notes)
		    elif('East' in x.key or 'West' in x.key):
		        line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  -------  ------  {:7.2f}  {:6.2f}  {:6.1f}  {:4.2f}  {:6.1f}  {:4.2f}  {:30s}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.u[i], x.ue[i], x.w[i], x.we[i], x.i[i], x.ie[i], x.b[i], x.be[i], x.notes)
		    elif('CV_VTI_EKU_PAR' in x.key):
		        line = ""
		    elif('Unknown' in x.key):
		        line = ""
		    else:
		        line = "{:14s}  {:19s}  {:5.1f}  {:5.1f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:6.1f}  {:4.2f}  {:6.1f}  {:4.2f}  {:30s}\n".format(x.key, utctime, x.lla[0], x.lla[1], x.T[i], x.Te[i], x.u[i], x.ue[i], x.v[i], x.ve[i], x.w[i], x.we[i], x.i[i], x.ie[i], x.b[i], x.be[i], x.notes)
		    #line = "%14s  %19s  %3.1f  %3.1f  %4.2f  %2.2f  %3.2f  %2.2f  %3.2f  %2.2f  %3.2f  %2.2f  %1.3f  %1.3f  %5s  %30s" % (x.key, utctime, lat, lon, x.T[i], x.Te[i], x.u[i], x.ue[i], x.v[i], x.ve[i], x.w[i], x.we[i], x.I[i], x.Ie[i], x.cloudy[i], x.notes)
		    note.write(line)
		    
    note.close()
    

def CreateMonthlyASCII(PROJECT,YEAR,MONTH):
    '''
    Summary:
        Script to output ASCII Monthly averaged results for a project.
    
    Inputs:
        PROJECT - Name of project for sites, ex: 'NATION'
        YEAR = year, e.g. 2013
        MONTH = month, e.g. 6

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
        MD = BinMonthlyData(site,YEAR,MONTH)
        
        for ind,dn in enumerate(MD.t):
            #dn = MD.t[ind] #.astimezone(_utc)
            utctime = dn.strftime("%H:%M")
            line = "%5s  %3.1f  %3.1f  %4.2f  %2.2f  %3.2f  %2.2f  %3.2f  %2.2f  %3.2f  %2.2f\n" % (utctime, S[site]['Location'][0], S[site]['Location'][1], MD.T[ind], MD.Tv[ind], MD.u[ind], MD.uv[ind], MD.v[ind], MD.vv[ind], MD.w[ind], MD.wv[ind])
            note.write(line)

    note.close
    
    

class _BinnedData:
    '''
    Note: BinnedData has T,u,v,w (u2,v2) data
        Each also has c-count of data points per time, e-weighted std, v-monthly variablility
    '''
    def __init__(self, dn, site):
        import datetime as dt
        import ephem

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
        if 'Tv' in ship:
            self.Tv   = _np.roll(self.Tv ,roll)
            self.Tc   = _np.roll(self.Tc ,roll)
            self.Tve  = _np.roll(self.Tve,roll)
            self.uv   = _np.roll(self.uv ,roll)
            self.uve  = _np.roll(self.uve,roll)
            self.uc   = _np.roll(self.uc ,roll)
            self.vv   = _np.roll(self.vv ,roll)
            self.vve  = _np.roll(self.vve,roll)
            self.vc   = _np.roll(self.vc ,roll)
            self.wv   = _np.roll(self.wv ,roll)
            self.wve  = _np.roll(self.wve,roll)
            self.wc   = _np.roll(self.wc ,roll)
            self.iv   = _np.roll(self.iv ,roll)
            self.ive  = _np.roll(self.ive,roll)
            self.ic   = _np.roll(self.ic ,roll)
        if 'u2' in ship:
            self.u2   = _np.roll(self.u2 ,roll)
            self.u2e  = _np.roll(self.u2e,roll)
            self.v2   = _np.roll(self.v2 ,roll)
            self.v2e  = _np.roll(self.v2e,roll)
        if 'u2v' in ship:
            self.u2v  = _np.roll(self.u2v ,roll)
            self.u2ve = _np.roll(self.u2ve,roll)
            self.u2c  = _np.roll(self.u2c ,roll)
            self.v2v  = _np.roll(self.v2v ,roll)
            self.v2ve = _np.roll(self.v2ve,roll)
            self.v2c  = _np.roll(self.v2c ,roll)
        if 'Tu' in ship:
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
        self.barrelroll = not(self.barrelroll)



def laser_is_drifting(instr_name, year, doy):
    '''
    Return True if the laser is suspected to be drifting.
    Return False if zenith reference was used.
    Raise an exception if there is no data.
    
    The laser is "suspected to be drifting" if the following 2 criteria are satisfied:
    1) The vertical wind at the beginning of the night and the end of the night are different by > 30 m/s.
    2) The laser intensity varies by more than 10% from the median across the night.
    '''
    import FPIprocess as _FPIprocess

    try:
        r = _FPIprocess.load_level0(instr_name, year, doy)
    except:
        raise Exception('No data found: No npz file for %s_%i_%i' % (instr_name, year, doy))
    
    fpir  = r['FPI_Results']
    if fpir['reference']=='zenith':
        return False
    
    direction = fpir['direction']
    LOSwind = fpir['LOSwind']
        
    direc = 'Zenith'
    w = _np.array([si for (si,d) in zip(LOSwind, direction) if d == direc])

    lasI = fpir['laser_value']['I']
    lasIe = fpir['laser_stderr']['I']

    # Check if laser varies by more than 10 %
    las_flag = sum(abs(lasI - _np.median(lasI))/_np.median(lasI) > 0.2) > 2
    # Check if vertical wind drifts by more than 30 m/s
    wstart = _np.median(w[:5])
    wend   = _np.median(w[-5:])
    w_flag = abs(wend-wstart) > 30.
    
    return w_flag and las_flag



if __name__=="__main__":

    print "Insert Coin"
