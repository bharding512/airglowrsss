'''
Summary:
    Contains codes to plot various things for papers (standard plots are in FPIResults)
    
Included are:
    
    
History:
    3/20/13 -- Written by DJF (dfisher2@illinois.edu)
    10/17/14 -- Split from FPIResults by DJF
'''

import matplotlib as _mpl
import matplotlib.dates as _md
import matplotlib.pyplot as _plt
#matplotlib.use('AGG')
from pyglow import pyglow as _pyglow
import datetime as _dt
import calendar as _cal
import numpy as _np
import FPIprocessLevel2 as _L2
import FPIprocessLevel2_Legacy as _L2L
import fpiinfo as _fpiinfo
import FPIResults as FPIResults


# Set up default parameters
_mpl.rcParams.update({'font.size': 11}) 
dirout = '/rdata/airglow/database/L2/plots/'
FPIResults.SetBinTime(15)


def PlotStorm(SITE,YEAR,DOYSTART):
    '''
    Summary:
        DATA = PlotStorm(SITE,YEAR,DOYSTART)
    
    Inputs
    ------
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        DOYSTART = day of year for storm day, e.g. 4

    Outputs
    -------
        DATA = dictionary of data whose keys are Zonal, Meridional, or Temp 
               Temp contains Temp and Temp_Error (averaged from all directions)
               Zonal/Meridional contains Wind and Wind_Error

    History
    -------
    3/22/13 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    
    # Get Monthly Average for basis
    dn = _dt.datetime(YEAR,1,1)+_dt.timedelta(days = DOYSTART-1)
    month = dn.month
    SD = []
    MD = FPIResults.BinMonthlyData(SITE,YEAR,month)
    
    # Get Daily Data for comparison
    for doy in range(DOYSTART,DOYSTART+4):
        dd = FPIResults.BinDailyData(SITE,YEAR,doy)
        SD.append(dd)
        ######## WHAT DO I WANT TO DO HERE WITH TIMS DATA STRUCT
    
    # Temp Figure
    fig = _plt.figure(1); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    _plt.errorbar(doy+times,MD.T,yerr=MD.Tv,label='M')
    _plt.errorbar(doy+times,SD.T,yerr=SD.Tv,marker='o',label='D')
    _plt.plot([doy,doy+5*15*b_len],[0,0],'k--')
    _plt.ylim([600,1200])
    _plt.xlim([doy-3.5,doy+1.5])
    _plt.ylabel('Temp [K]')
    _plt.xlabel('Time')
    _plt.legend()
    _plt.draw(); _plt.show()
    
    # Winds Figure Zonal
    fig = _plt.figure(1); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    _plt.errorbar(doy+times,MD.u,yerr=MD.uv,label='M')
    _plt.errorbar(doy+times,SD.u,yerr=SD.uv,marker='o',label='D')
    _plt.plot([doy,doy+5*15*b_len],[0,0],'k--')
    _plt.ylim([-200,200])
    _plt.xlim([doy-3.5,doy+1.5])
    _plt.ylabel('Zonal Wind Speed [m/s]')
    _plt.xlabel('Time')
    _plt.legend()
    _plt.draw(); _plt.show()
    
    # Winds Figure Meridional
    fig = figure(1); clf()
    ax = fig.add_subplot(1,1,1)
    _plt.errorbar(doy+times,MD.v,yerr=MD.vv,label='M')
    _plt.errorbar(doy+times,SD.v,yerr=SD.vv,marker='o',label='D')
    _plt.plot([doy,doy+5*15*b_len],[0,0],'k--')
    _plt.ylim([-200,200])
    _plt.xlim([doy-3.5,doy+1.5])
    _plt.ylabel('Meridional Wind Speed [m/s]')
    _plt.xlabel('Time')
    _plt.legend()
    _plt.draw(); _plt.show()
    #title('Storm from %s beginning on %s' % SITE, dn.strftime('%Y-%m-%d'))
    
    '''
    fig = figure(1); clf()
    ax = fig.add_subplot(1,1,1)
    _plt.errorbar(times,MD['Zonal']['Wind'],yerr=MD['Zonal']['Wind_Error'],marker='o',label='M')
    _plt.errorbar(times,SD['East']['Wind'],yerr=SD['East']['Wind_Error'],marker='o',label='E')
    _plt.errorbar(times,SD['West']['Wind'],yerr=SD['West']['Wind_Error'],marker='o',label='S')
    _plt.plot([0,4*15*b_len],[0,0],'k--')
    _plt.ylim([-200,200])
    _plt.ylabel('Wind Speed [m/s]')
    _plt.xlabel('Time')
    _plt.title('Storm from %s beginning on %s' % SITE, dn.strftime('%Y-%m-%d'))
    #ax.xaxis.seTmajor_formatter(_md.DateFormatter('%H:%M'))
    #fig.autofmt_xdate()
    _plt.legend()
    _plt.draw(); _plt.show()
    '''
    
    
def PlotCompSADay(YEAR,MONTH,DAY):
    '''
    Summary:
        Plot Comparisions of RENOIR vs PERU for a day
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTH = month to plot, e.g. 4

    History:
    5/8/14 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    
    # Get Monthly Average for basis
    DOY = (_dt.datetime(YEAR,MONTH,DAY) - _dt.datetime(YEAR,1,1)).days+1
    D1 = FPIResults.BinDailyData('renoir',YEAR,DOY)
    D2 = FPIResults.BinDailyData('peru',YEAR,DOY)
    
    # Get variables to plot
    markerwide = 0
    D2.t = D2.t + _dt.timedelta(minutes=5)
    
    # Get xlimits w/ data
    tstart = _np.where(_np.isnan(D1.T[:-1])*_np.isfinite(D1.T[1:]))[0][0]-1
    tend   = _np.where(_np.isfinite(D2.T[:-1])*_np.isnan(D2.T[1:]))[0][0]+1
    tlim = [D1.t[tstart],D2.t[tend]]
        
    # Temp Figure
    fig = _plt.figure(0); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    (_, caps1, _) = _plt.errorbar(D1.t,D1.T,yerr=D1.Tv,fmt='-ro',label='RENOIR')
    (_, caps2, _) = _plt.errorbar(D2.t,D2.T,yerr=D2.Tv,fmt='-bs',label='PERU')
    _plt.plot([D1.t[tstart],D2.t[tend]],[0,0],'k--')
    _plt.ylim([500,1000])
    _plt.xlim(tlim)
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for cap in caps1:
        #cap.set_color('red')
        cap.set_markeredgewidth(markerwide)
    for cap in caps2:
        cap.set_markeredgewidth(markerwide)
    _plt.ylabel('Temperature [K]')
    _plt.xlabel('Time [UT]')
    _plt.title('%04d-%02d-%02d Average Temperatures' % (YEAR,MONTH,DAY))
    _plt.legend()
    _plt.draw(); _plt.show()
    _plt.savefig('%sComp_%04d-%02d-%02d_temps.png' % (dirout,YEAR,MONTH,DAY))
    
    # Winds Figure Zonal
    fig = _plt.figure(1); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    (_, caps1, _) = _plt.errorbar(D1.t,D1.u,yerr=D1.uv,fmt='-ro',label='RENOIR')
    (_, caps2, _) = _plt.errorbar(D2.t,D2.u,yerr=D2.uv,fmt='-bs',label='PERU')
    _plt.plot([D1.t[tstart],D2.t[tend]],[0,0],'k--')
    _plt.ylim([-200,200])
    _plt.xlim(tlim)
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for cap in caps1:
        #cap.set_color('red')
        cap.set_markeredgewidth(markerwide)
    for cap in caps2:
        cap.set_markeredgewidth(markerwide)
    _plt.ylabel('Wind Speed [m/s]')
    _plt.xlabel('Time [UT]')
    _plt.title('%04d-%02d-%02d Average Zonal Winds' % (YEAR,MONTH,DAY))
    _plt.legend()
    _plt.draw(); _plt.show()
    _plt.savefig('%sComp_%04d-%02d-%02d_zonal_winds.png' % (dirout,YEAR,MONTH,DAY))
    
    # Winds Figure Meridional
    fig = _plt.figure(2); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    (_, caps1, _) = _plt.errorbar(D1.t,D1.v,yerr=D1.vv,fmt='-ro',label='RENOIR')
    (_, caps2, _) = _plt.errorbar(D2.t,D2.v,yerr=D2.vv,fmt='-bs',label='PERU')
    _plt.plot([D1.t[tstart],D2.t[tend]],[0,0],'k--')
    _plt.ylim([-200,200])
    _plt.xlim(tlim)
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for cap in caps1:
        #cap.set_color('red')
        cap.set_markeredgewidth(markerwide)
    for cap in caps2:
        cap.set_markeredgewidth(markerwide)
    _plt.ylabel('Wind Speed [m/s]')
    _plt.xlabel('Time [UT]')
    _plt.title('%04d-%02d-%02d Average Meridional Winds' % (YEAR,MONTH,DAY))
    _plt.legend()
    _plt.draw(); _plt.show()
    _plt.savefig('%sComp_%04d-%02d-%02d_meridional_winds.png' % (dirout,YEAR,MONTH,DAY))
    
    '''
    # Winds Figure Vertical
    fig = _plt.figure(3); _plt.clf()
    ax = fig.add_subplot(1,1,1)
    (_, caps1, _) = _plt.errorbar(D1.t,D1.w,yerr=D1.wv,fmt='-ro',label='RENOIR')
    (_, caps2, _) = _plt.errorbar(D2.t,D2.w,yerr=D2.wv,fmt='-bs',label='PERU')
    _plt.plot([D1.t[tstart],D2.t[tend]],[0,0],'k--')
    _plt.ylim([-200,200])
    _plt.xlim(tlim)
    ax.xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for cap in caps1:
        #cap.set_color('red')
        cap.set_markeredgewidth(markerwide)
    for cap in caps2:
        cap.set_markeredgewidth(markerwide)
    _plt.ylabel('Wind Speed [m/s]')
    _plt.xlabel('Time [UT]')
    _plt.title('%04d-%02d-%02d Average Vertical Winds' % (YEAR,MONTH,DAY))
    _plt.legend()
    _plt.draw(); _plt.show()
    _plt.savefig('%sComp_%04d-%02d-%02d_vertical_winds.png' % (dirout,YEAR,MONTH,DAY))
    '''


def PlotCompSAStacked(setnum):
    '''
    Summary:
        Plot Comparisons of RENOIR vs PERU for a few hardcoded days using binned averages
    
    History:
        5/30/14 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    SetBinTime(30)
    
    if setnum==1:
        doyset = [246,247,248]
        YEAR = 2010
    elif setnum==2:
        doyset = [251,252,255]
        YEAR = 2010
    elif setnum==3:
        doyset = [199,200,201]
        YEAR = 2011
    else:
        return
    
    # Get Monthly Average for basis
    R={};P={}; ax={}
    Rt={}; RT={}; RTv={}; RU={}; RUv={}; RV={}; RVv={}; RW={}; RWv={}; Pt={}; PT={}; PTv={}; PU={}; PUv={}; PV={}; PVv={}; PW={}; PWv={};
    for k,d in enumerate(doyset):
        DR = FPIResults.BinDailyData('renoir',YEAR,d)
        R[k] = DR
        DP = FPIResults.BinDailyData('peru',YEAR,d)
        DP.t = DP.t + _dt.timedelta(minutes=5)
        P[k] = DP
        MR = FPIResults.GetModels(DR.lla,YEAR,d,WMODEL='hwm14',TMODEL='msis')
        MP = FPIResults.GetModels(DP.lla,YEAR,d,WMODEL='hwm14',TMODEL='msis')
        Rt[k] = MR.t
        RT[k] = MR.T
        RTe[k] = MR.Te
        Pt[k] = MP.t + _dt.timedelta(minutes=5)
        PT[k] = MP.T-75
        PTe[k] = MP.Te

    markerwide = 0
    lalpha = .3

    # Get xlimits w/ data
    tstart = _np.where(_np.isnan(R[0].T[:-1])*_np.isfinite(R[0].T[1:]))[0][0]-1
    tend   = _np.where(_np.isfinite(P[0].T[:-1])*_np.isnan(P[0].T[1:]))[0][0]+1
    tlim = [R[0].t[tstart],P[0].t[tend]]
    
    _plt.close('all')
    # Temp Figure
    f, (ax[0], ax[1], ax[2]) = _plt.subplots(3, sharex=True, sharey=True)
    for k,d in enumerate(doyset):
        #MODEL
        l3, = ax[k].plot(Rt[k],RT[k],'m--',linewidth=2.,label='MSIS-R')
        l4, = ax[k].plot(Pt[k],PT[k],'c--',linewidth=2.,label='MSIS-P')
        #DATA
        (m1, caps1, _) = ax[k].errorbar(R[k].t,R[k].T,yerr=R[k].Te,fmt='ro',label='RENOIR')
        (m2, caps2, _) = ax[k].errorbar(P[k].t,P[k].T,yerr=P[k].Te,fmt='bs',label='PERU')
        l1, = ax[k].plot(R[k].t,R[k].T,'r-',alpha=lalpha)
        l2, = ax[k].plot(P[k].t,P[k].T,'b-',alpha=lalpha)
        ax[k].plot(tlim,[0,0],'k--')
        for cap in caps1:
            cap.set_markeredgewidth(markerwide)
        for cap in caps2:
            cap.set_markeredgewidth(markerwide)
    f.subplots_adjust(hspace=0)
    _plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    _plt.ylim([500,1000])
    for k,d in enumerate(doyset):
        yticks = ax[k].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    ax[0].set_title('Average Temperatures')
    ax[1].set_ylabel('Temperature [K]')
    ax[2].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.xlim(tlim)
    _plt.xlabel('Time [UT]')
    ax[0].legend([(l1,m1),(l2,m2),l3,l4],['RENOIR','PERU','MSIS-RENOIR','MSIS-PERU'],loc='upper right', ncol=4, prop={'size':12}, fancybox=True)
    _plt.draw();
    _plt.savefig('%sComp_set_%03d_temps.png' % (dirout,setnum))
    
    # Winds Figure Zonal
    f, (ax[0], ax[1], ax[2]) = _plt.subplots(3, sharex=True, sharey=True)
    for k,d in enumerate(doyset):
        (m1, caps1, _) = ax[k].errorbar(R[k].t,R[k].u,yerr=R[k].ue,fmt='ro',label='RENOIR')
        (m2, caps2, _) = ax[k].errorbar(P[k].t,P[k].u,yerr=P[k].ue,fmt='bs',label='PERU')
        l1, = ax[k].plot(R[k].t,R[k].u,'r-',alpha=lalpha)
        l2, = ax[k].plot(P[k].t,P[k].u,'b-',alpha=lalpha)
        ax[k].plot(tlim,[0,0],'k--')
        for cap in caps1:
            cap.set_markeredgewidth(markerwide)
        for cap in caps2:
            cap.set_markeredgewidth(markerwide)
    f.subplots_adjust(hspace=0)
    _plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    _plt.ylim([-200,200])
    for k,d in enumerate(doyset):
        yticks = ax[k].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    ax[0].set_title('Average Zonal Winds')
    ax[1].set_ylabel('Wind Speed [m/s]')
    ax[2].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.xlim(tlim)
    _plt.xlabel('Time [UT]')
    ax[0].legend([(m1,l1),(m2,_L2)],['RENOIR','PERU'],loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    _plt.draw();
    _plt.savefig('%sComp_set_%03d_zonal_winds.png' % (dirout,setnum))
    
    # Winds Figure Meridional
    f, (ax[0], ax[1], ax[2]) = _plt.subplots(3, sharex=True, sharey=True)
    for k,d in enumerate(doyset):
        (m1, caps1, _) = ax[k].errorbar(R[k].t,R[k].v,yerr=R[k].ve,fmt='ro',label='RENOIR')
        (m2, caps2, _) = ax[k].errorbar(P[k].t,P[k].v,yerr=P[k].ve,fmt='bs',label='PERU')
        l1, = ax[k].plot(R[k].t,R[k].v,'r-',alpha=lalpha)
        l2, = ax[k].plot(P[k].t,P[k].v,'b-',alpha=lalpha)
        ax[k].plot(tlim,[0,0],'k--')
        for cap in caps1:
            cap.set_markeredgewidth(markerwide)
        for cap in caps2:
            cap.set_markeredgewidth(markerwide)
    f.subplots_adjust(hspace=0)
    _plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    _plt.ylim([-100,100])
    for k,d in enumerate(doyset):
        yticks = ax[k].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    ax[0].set_title('Average Meridional Winds')
    ax[1].set_ylabel('Wind Speed [m/s]')
    ax[2].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.xlim(tlim)
    _plt.xlabel('Time [UT]')
    ax[0].legend([(m1,l1),(m2,_L2)],['RENOIR','PERU'],loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    _plt.draw();
    _plt.savefig('%sComp_set_%03d_meridional_winds.png' % (dirout,setnum))

    '''
    # Winds Figure Vertical
    fig = _plt.figure(3); _plt.clf(); ax = {}
    f, (ax[0], ax[1], ax[2]) = _plt.subplots(3, sharex=True, sharey=True)
    for k,d in enumerate(doyset):
        (m1, caps1, _) = ax[k].errorbar(R[k].t,R[k].w,yerr=R[k].we,fmt='ro',label='RENOIR')
        (m2, caps2, _) = ax[k].errorbar(P[k].t,P[k].w,yerr=P[k].we,fmt='bs',label='PERU')
        l1, = ax[k].plot(R[k].t,R[k].w,'r-',alpha=lalpha)
        l2, = ax[k].plot(P[k].t,P[k].w,'b-',alpha=lalpha)
        ax[k].plot(tlim,[0,0],'k--')
        for cap in caps1:
            cap.set_markeredgewidth(markerwide)
        for cap in caps2:
            cap.set_markeredgewidth(markerwide)
    f.subplots_adjust(hspace=0)
    _plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    _plt.ylim([-100,100])
    for k,d in enumerate(doyset):
        yticks = ax[k].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    ax[0].set_title('Average Vertical Winds')
    ax[1].set_ylabel('Wind Speed [m/s]')
    ax[2].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.xlim(tlim)
    _plt.xlabel('Time [UT]')
    ax[0].legend([(m1,l1),(m2,_L2)],['RENOIR','PERU'],loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    _plt.draw(); _plt.show()
    _plt.savefig('%sComp_set_%04d_vertical_winds.png' % (dirout,setnum))
    '''
    SetBinTime(15)


def PlotCompSAStackedv2(setnum):
    '''
    Summary:
        Plot Comparisons of RENOIR vs PERU for a few hardcoded days using raw east & south w/ MSIS model for temps
    
    History:
        5/30/14 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    if setnum==1:
        doyset = [246,247,248]
        YEAR = 2010
        RToff = 25
    elif setnum==2:
        doyset = [251,252,255]
        YEAR = 2010
        RToff = 25
    elif setnum==3:
        doyset = [199,200,201]
        YEAR = 2011
        RToff = 0
    else:
        return
    
    # Get Monthly Average for basis
    R={};P={}; ax={}
    Rt={}; RT={}; RTv={}; RU={}; RUv={}; RV={}; RVv={}; RW={}; RWv={}; Pt={}; PT={}; PTv={}; PU={}; PUv={}; PV={}; PVv={}; PW={}; PWv={};
    for k,d in enumerate(doyset):
        #DR = _L2.GetLevel('minime01',_dt.datetime(YEAR,1,1)+_dt.timedelta(days=d-1))
        DR = _L2.GetLevel2('renoir',_dt.datetime(YEAR,1,1)+_dt.timedelta(days=d-1))
        FPIResults.FilterData(DR); R[k] = DR
        #DP = _L2.GetLevel1('minime90',_dt.datetime(YEAR,1,1)+_dt.timedelta(days=d-1))
        DP = _L2.GetLevel2('peru',_dt.datetime(YEAR,1,1)+_dt.timedelta(days=d-1))
        FPIResults.FilterData(DP); P[k] = DP
        indr = [x for x,y in enumerate(DR) if 'Zenith' in y.key][0]
        indp = [x for x,y in enumerate(DP) if 'Zenith' in y.key][0]
        MR = FPIResults.GetModels(DR[indr].lla,YEAR,d,WMODEL='hwm14',TMODEL='msis')
        MP = FPIResults.GetModels(DP[indp].lla,YEAR,d,WMODEL='hwm14',TMODEL='msis')
        Rt[k] = MR.t + _dt.timedelta(seconds=(_dt.datetime(YEAR-1,1,1)-_dt.datetime(MP.t[0].year,1,1)).total_seconds())
        RT[k] = MR.T-RToff
        Pt[k] = MP.t + _dt.timedelta(seconds=(_dt.datetime(YEAR-1,1,1)-_dt.datetime(MP.t[0].year,1,1)).total_seconds()) + _dt.timedelta(minutes=5)
        PT[k] = MP.T-75
        RU[k] = MR.u
        RV[k] = MR.v 
        RW[k] = MR.w
        PU[k] = MP.u
        PV[k] = MP.v
        PW[k] = MP.w
        
    Rmidnight = [Rt[0][-1].replace(hour=2,minute=12),Rt[0][-1].replace(hour=2,minute=12)]
    Pmidnight = [Pt[0][-1].replace(hour=4,minute=53),Pt[0][-1].replace(hour=4,minute=53)]
    markerwide = 0
    lalpha = .3

    # Get xlimits w/ data
    tlim = [(R[0][0].allt[0]-_dt.timedelta(days=doyset[0],minutes=80)).astimezone(FPIResults._utc).replace(tzinfo=None), (P[0][0].allt[-1]-_dt.timedelta(days=doyset[0],minutes=50)+_dt.timedelta(hours=1)).astimezone(FPIResults._utc).replace(tzinfo=None)]
    
    _plt.close('all')
    # Temp Figure
    f, ((ax[0],ax[1]), (ax[2],ax[3]), (ax[4],ax[5])) = _plt.subplots(3,2, sharex=True, sharey=True)
    # Plot Data
    for k,d in enumerate(doyset):
        #MODEL
        l3, = ax[k*2].plot(Rt[k],RT[k],'m--',linewidth=2.,label='MSIS-R')
        l4, = ax[k*2+1].plot(Pt[k],PT[k],'c--',linewidth=2.,label='MSIS-P')
        #DATA
        #for l in [y for y in R[k] if ('CAR_West' in y.key) or ('CAR_South' in y.key)]:
        for l in [y for y in R[k] if ('CAR_West' in y.key)]:
            (m1, caps1, _) = ax[k*2].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.T,yerr=l.Te,fmt='ro',label='RENOIR')
            l1, = ax[k*2].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.T,'r-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        for l in [y for y in R[k] if ('CAR_South' in y.key)]:
            (m1, caps1, _) = ax[k*2].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.T,yerr=l.Te,fmt='go',label='RENOIR')
            l1, = ax[k*2].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.T,'g-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        #for l in [y for y in P[k] if ('MRH_East' in y.key) or ('MRH_North' in y.key)]:
        for l in [y for y in P[k] if ('MRH_East' in y.key)]:
            (m2, caps1, _) = ax[k*2+1].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.T,yerr=l.Te,fmt='bs',label='PERU')
            l2, = ax[k*2+1].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.T,'b-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        for l in [y for y in P[k] if ('MRH_North' in y.key)]:
            (m2, caps1, _) = ax[k*2+1].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.T,yerr=l.Te,fmt='ys',label='PERU')
            l2, = ax[k*2+1].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.T,'y-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        ax[k*2].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(Pmidnight,[-2000,2000],'c--')
        ax[k*2].plot(Rmidnight,[-2000,2000],'m--')
    #ax[0].legend([(l1,m1),(l2,m2),l3,l4],['RENOIR','PERU','MSIS-RENOIR','MSIS-PERU'],loc='upper right', ncol=4, prop={'size':12}, fancybox=True)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    #_plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    _plt.ylim([500,1000])
    _plt.xlim(tlim)
    ax[4].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[5].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for k,d in enumerate(doyset):
        yticks = ax[k*2+1].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        yticks = ax[k*2].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k*2+1].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
        ax[k*2].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    ax[1].set_title('Peru')
    ax[0].set_title('Brazil')
    ax[2].set_ylabel('Temperature [K]')
    #_plt.xlabel('Time [UT]')
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    #ax[0].legend(loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    _plt.draw();
    #_plt.savefig('%sComp_set_%03d_temps2.png' % (dirout,setnum))
    _plt.savefig('%sComp_set_%03d_temps2.pdf' % (dirout,setnum))
    
    # Winds Figure Zonal
    f, ((ax[0],ax[1]),(ax[2],ax[3]),(ax[4],ax[5])) = _plt.subplots(3,2, sharex=True, sharey=True)
    # Data
    for k,d in enumerate(doyset):
        for l in [y for y in R[k] if ('CAR_West' in y.key)]:
            #MODEL
            l3, = ax[k*2].plot(Rt[k],RU[k],'m--',linewidth=2.,label='MSIS-R')
            l4, = ax[k*2+1].plot(Pt[k],PU[k],'c--',linewidth=2.,label='MSIS-P')
            (m1, caps1, _) = ax[k*2].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.u,yerr=l.ue,fmt='ro',label='RENOIR')
            l1, = ax[k*2].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.u,'r-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        for l in [y for y in P[k] if ('MRH_East' in y.key)]:
            (m2, caps1, _) = ax[k*2+1].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.u,yerr=l.ue,fmt='bs',label='PERU')
            l2, = ax[k*2+1].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.u,'b-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        ax[k*2].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(Pmidnight,[-2000,2000],'c--')
        ax[k*2].plot(Rmidnight,[-2000,2000],'m--')
    #ax[0].legend([(l1,m1),(l2,m2)],['RENOIR','PERU'],loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    #_plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    _plt.ylim([-200,200])
    _plt.xlim(tlim)
    ax[4].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[5].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for k,d in enumerate(doyset):
        yticks = ax[k*2].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        yticks = ax[k*2+1].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k*2].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
        ax[k*2+1].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    ax[1].set_title('Peru')
    ax[0].set_title('Brazil')
    ax[2].set_ylabel('Zonal Wind Speed [m/s]')
    #_plt.xlabel('Time [UT]')
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    #ax[0].legend(loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    _plt.draw();
    #_plt.savefig('%sComp_set_%03d_zonal_winds2.png' % (dirout,setnum))
    _plt.savefig('%sComp_set_%03d_zonal_winds2.pdf' % (dirout,setnum))
    
    # Winds Figure Meridional
    f, ((ax[0],ax[1]), (ax[2],ax[3]), (ax[4],ax[5])) = _plt.subplots(3,2, sharex=True, sharey=True)
    # Data
    for k,d in enumerate(doyset):
        for l in [y for y in R[k] if ('CAR_South' in y.key)]:
            l3, = ax[k*2].plot(Rt[k],RV[k],'m--',linewidth=2.,label='MSIS-R')
            l4, = ax[k*2+1].plot(Pt[k],PV[k],'c--',linewidth=2.,label='MSIS-P')
            (m1, caps1, _) = ax[k*2].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.v,yerr=l.ve,fmt='go',label='RENOIR')
            l1, = ax[k*2].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.v,'g-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        for l in [y for y in P[k] if ('MRH_South' in y.key)]:
            (m2, caps1, _) = ax[k*2+1].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.v,yerr=l.ve,fmt='ys',label='PERU')
            l2, = ax[k*2+1].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.v,'y-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        ax[k*2].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(Pmidnight,[-2000,2000],'c--')
        ax[k*2].plot(Rmidnight,[-2000,2000],'m--')
    #ax[0].legend([(l1,m1),(l2,m2)],['RENOIR','PERU'],loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    _plt.ylim([-100,100])
    _plt.xlim(tlim)
    ax[4].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[5].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for k,d in enumerate(doyset):
        yticks = ax[k*2].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        yticks = ax[k*2+1].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k*2].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
        ax[k*2+1].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    ax[1].set_title('Peru')
    ax[0].set_title('Brazil')
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    f.text(0.06, 0.5, 'Meridional Wind Speed [m/s]', ha='center', va='center', rotation='vertical')
    #ax[0].legend([(m1,l1),(m2,_L2)],['RENOIR','PERU'],loc='upper left', prop={'size':12}, fancybox=True)[figure()
    _plt.draw();
    #_plt.savefig('%sComp_set_%03d_meridional_winds2.png' % (dirout,setnum))
    _plt.savefig('%sComp_set_%03d_meridional_winds2.pdf' % (dirout,setnum))
    
    '''
    # Winds Figure Vertical
    fig = _plt.figure(3); _plt.clf(); ax = {}
    f, ((ax[0],ax[1]), (ax[2],ax[3]), (ax[4],ax[5])) = _plt.subplots(3,2, sharex=True, sharey=True)
    # Data
    for k,d in enumerate(doyset):
        #MODEL
        l3, = ax[k*2].plot(Rt[k],RW[k],'m--',linewidth=2.,label='MSIS-R')
        l4, = ax[k*2+1].plot(Pt[k],PW[k],'c--',linewidth=2.,label='MSIS-P')
        for l in [y for y in R[k] if ('CAR_Zenith' in y.key)]:
            (m1, caps1, _) = ax[k*2].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.w,yerr=l.we,fmt='ro',label='RENOIR')
            l1, = ax[k*2].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.w,'r-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        for l in [y for y in P[k] if ('MRH_Zenith' in y.key)]:
            (m2, caps1, _) = ax[k*2+1].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.w,yerr=l.we,fmt='bs',label='PERU')
            l2, = ax[k*2+1].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=d)]),l.w,'b-',alpha=lalpha)
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
        ax[k*2+1].plot(tlim,[0,0],'k--')
        ax[k*2].plot(tlim,[0,0],'k--')
    #ax[0].legend([(l1,m1),(l2,m2)],['RENOIR','PERU'],loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    f.subplots_adjust(hspace=0)
    _plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    _plt.ylim([-100,100])
    for k,d in enumerate(doyset):
        yticks = ax[k*2].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        yticks = ax[k*2+1].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k].annotate('Doy %03d'%d, xy=(0.02,0.03), xycoords='axes fraction', fontsize=12)
    ax[0].set_title('PERU')
    ax[1].set_title('RENOIR')
    ax[2].set_ylabel('Vertical Wind Speed [m/s]')
    ax[4].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[5].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    _plt.xlim(tlim)
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    #_plt.xlabel('Time [UT]')
    #ax[0].legend(loc='lower center', ncol=2, prop={'size':12}, fancybox=True)
    _plt.draw(); _plt.show()
    _plt.savefig('%sComp_set_%04d_vertical_winds2.png' % (dirout,setnum))
    '''

 
       
def PlotCompSAStacked2():
    '''
    Summary:
        Plot Comparisons of RENOIR vs PERU for a few hardcoded monthly averages
    
    History:
        9/29/14 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    FPIResults.SetBinTime(20)
    
    monset = _np.array([8,9,10,10])
    yrset  = _np.array([2010,2010,2010,2011])
    RToff = _np.array([50,25,25,50])
    Rl = [[216,229,232,234,235,240,242,243], [245,246,247,248,249,250,251,253,254,255,258,259,261,262,263,264,266,269,270,271,273], [274,275,276,277,279,280,281,282,303,304], [299,276,277,278,279,280,294,298,300]]
    Pl = [[224,225,226,228,230,228,231,232,237,238,239,240,243], [246,247,248,249,251,252,253,254,255,256,257,258,259,260,261,262,263,264,266,267,273,272,271],     [275,276,277,279,280,281,282,283,284,286,287,288,289,290,291,292,293,295,297,298,299,300,302], [274,275,276,277,278,302,303,304,300,274,275,277,278,279,280,282,284,295,300]]
    
    # Get Monthly Average for basis
    R={};P={}
    Rt={}; RT={}; RTv={}; RU={}; RUv={}; RV={}; RVv={}; RW={}; RWv={}; Pt={}; PT={}; PTv={}; PU={}; PUv={}; PV={}; PVv={}; PW={}; PWv={};
    for k,d in enumerate(zip(yrset,monset)):
        #MR = FPIResults.BinMonthlyData('car',d[0],d[1],SPLIT=True,KP=[0,4],LIST=Rl[k])
        MR = FPIResults.BinMonthlyData('car',d[0],d[1],SPLIT=True,KP=[0,4])
        R[k] = MR
        R[k].midnight = [MR.t[-1].replace(hour=2,minute=12),MR.t[-1].replace(hour=2,minute=12)]
        #MP = FPIResults.BinMonthlyData('mrh',d[0],d[1],SPLIT=True,KP=[0,4],LIST=Pl[k],CV=False)
        MP = FPIResults.BinMonthlyData('mrh',d[0],d[1],SPLIT=True,KP=[0,4],CV=False)
        MP.t = MP.t + _dt.timedelta(minutes=5)
        P[k] = MP
        P[k].midnight =[MP.t[-1].replace(hour=4,minute=53),MP.t[-1].replace(hour=4,minute=53)]
        MR = FPIResults.BinMonthlyData('hwm14',d[0],d[1],SITELLA=MR.lla,KP=[0,4])
        MP = FPIResults.BinMonthlyData('hwm14',d[0],d[1],SITELLA=MP.lla,KP=[0,4])
        Rt[k] = MR.t
        RT[k] = MR.T-RToff[k]
        RTv[k] = MR.Tv
        RU[k] = MR.u
        RUv[k] = MR.uv
        RV[k] = MR.v
        RVv[k] = MR.vv 
        RW[k] = MR.w
        RWv[k] = MR.wv
        Pt[k] = MP.t + _dt.timedelta(minutes=3)
        PT[k] = MP.T-75
        PTv[k] = MP.Tv
        PU[k] = MP.u
        PUv[k] = MP.uv
        PV[k] = MP.v
        PVv[k] = MP.vv
        PW[k] = MP.w
        PWv[k] = MP.wv
    
    markerwide = 0
    lalpha = 0.3

    # Get xlimits w/ data (IMPROVE)
    tstart = _np.where(_np.isnan(R[0].T[:-1])*_np.isfinite(R[0].T[1:]))[0][0]-1
    tend   = _np.where(_np.isfinite(P[0].T[:-1])*_np.isnan(P[0].T[1:]))[0][0]+1
    tlim = [R[0].t[tstart],P[0].t[tend]]
    
    _plt.close('all');ax={}
    
    ##### Temp Figure
    f,((ax[1],ax[0]),(ax[3],ax[2]),(ax[5],ax[4]),(ax[7],ax[6]))  = _plt.subplots(4,2, sharex=True, sharey=True)
    for k,d in enumerate(zip(yrset,monset)):
        ax[k*2+1].fill_between(Rt[k],RT[k]-RTv[k],RT[k]+RTv[k],alpha=0.5,linewidth=0,facecolor='m',label='MSIS-R')
        ax[k*2].fill_between(Pt[k],PT[k]-PTv[k],PT[k]+PTv[k],alpha=0.5,linewidth=0,facecolor='c',label='MSIS-P')
        (m1, caps1, _) = ax[k*2+1].errorbar(R[k].t,R[k].T,yerr=R[k].Tv,fmt='ro',label='R')
        l1, = ax[k*2+1].plot(R[k].t,R[k].T,'r-',alpha=lalpha)
        (m2, caps2, _) = ax[k*2].errorbar(P[k].t,P[k].T,yerr=P[k].Tv,fmt='bs',label='P')
        l2, = ax[k*2].plot(P[k].t,P[k].T,'b-',alpha=lalpha)
        ax[k*2].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(tlim,[0,0],'k--')
        ax[k*2].plot(P[k].midnight,[-2000,2000],'c--')
        ax[k*2+1].plot(R[k].midnight,[-2000,2000],'m--')
        for cap in caps1+caps2:
            cap.set_markeredgewidth(markerwide)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    _plt.ylim([500,1000])
    _plt.xlim(tlim)
    ax[6].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[7].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for k,d in enumerate(zip(yrset,monset)):
        xticks = ax[k].xaxis.get_major_ticks()
        #xticks[0].label1.set_visible(False)
        yticks = ax[k*2].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        yticks = ax[k*2+1].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k*2].annotate('%04d %s'%(d[0],_cal.month_name[d[1]]), xy=(0.05,0.05), xycoords='axes fraction', fontsize=12)
        ax[k*2+1].annotate('%04d %s'%(d[0],_cal.month_name[d[1]]), xy=(0.05,0.05), xycoords='axes fraction', fontsize=12)
    ax[0].set_title('Peru')
    ax[1].set_title('Brazil')
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    f.text(0.06, 0.5, 'Temperature [K]', ha='center', va='center', rotation='vertical')
    #ax[0].legend([(m1,l1),(m2,_L2)],['RENOIR','PERU'],loc='upper right',ncol=2, prop={'size':12}, fancybox=True)
    _plt.draw();
    #_plt.savefig('%sCompM_set_temps.png'%(dirout))
    _plt.savefig('%sCompM_set_temps.pdf'%(dirout))
    
    ##### Winds Figure Zonal
    f,((ax[1],ax[0]),(ax[3],ax[2]),(ax[5],ax[4]),(ax[7],ax[6])) = _plt.subplots(4,2, sharex=True, sharey=True)
    for k,d in enumerate(zip(yrset,monset)[0:-1]):
        ax[k*2+1].fill_between(Rt[k],RU[k]-RUv[k],RU[k]+RUv[k],alpha=0.5,linewidth=0,facecolor='m',label='MSIS-R')
        ax[k*2].fill_between(Pt[k],PU[k]-PUv[k],PU[k]+PUv[k],alpha=0.5,linewidth=0,facecolor='c',label='MSIS-P')
        (m1, caps1, _) = ax[k*2+1].errorbar(R[k].t,R[k].u,yerr=R[k].uv,fmt='ro',label='R')
        l1, = ax[k*2+1].plot(R[k].t,R[k].u,'r-',alpha=lalpha)
        (m2, caps2, _) = ax[k*2].errorbar(P[k].t,P[k].u2,yerr=P[k].u2v,fmt='bs',label='P')
        l2, = ax[k*2].plot(P[k].t,P[k].u2,'b-',alpha=lalpha)
        ax[k*2].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(tlim,[0,0],'k--')
        ax[k*2].plot(P[k].midnight,[-2000,2000],'c--')
        ax[k*2+1].plot(R[k].midnight,[-2000,2000],'m--')
        for cap in caps1+caps2:
            cap.set_markeredgewidth(markerwide)
    # Manual fix for weird drive issues with Oct2011
    #fix = FPIResults.BinMonthlyData('peru',2011,10,SPLIT=True)
    #fix.t = fix.t + _dt.timedelta(minutes=5)
    k = 3
    ax[k*2+1].fill_between(Rt[k],RU[k]-RUv[k],RU[k]+RUv[k],alpha=0.5,linewidth=0,facecolor='m',label='MSIS-R')
    ax[k*2].fill_between(Pt[k],PU[k]-PUv[k],PU[k]+PUv[k],alpha=0.5,linewidth=0,facecolor='c',label='MSIS-P')
    (m1, caps1, _) = ax[k*2+1].errorbar(R[k].t,R[k].u,yerr=R[k].uv,fmt='ro',label='RENOIR')
    l1, = ax[k*2+1].plot(R[k].t,R[k].u,'r-',alpha=lalpha)
    (m2, caps2, _) = ax[k*2].errorbar(P[k].t,P[k].u,yerr=P[k].uv,fmt='bs',label='PERU')
    l2, = ax[k*2].plot(P[k].t,P[k].u,'b-',alpha=lalpha)
    ax[k*2].plot(tlim,[0,0],'k--')
    ax[k*2+1].plot(tlim,[0,0],'k--')
    ax[k*2].plot(P[k].midnight,[-2000,2000],'c--')
    ax[k*2+1].plot(R[k].midnight,[-2000,2000],'m--')
    for cap in caps1+caps2:
        cap.set_markeredgewidth(markerwide)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    _plt.ylim([-50,200])
    _plt.xlim(tlim)
    ax[6].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[7].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for k,d in enumerate(zip(yrset,monset)):
        xticks = ax[k].xaxis.get_major_ticks()
        #xticks[0].label1.set_visible(False)
        yticks = ax[k*2].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        yticks = ax[k*2+1].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k*2].annotate('%04d %s'%(d[0],_cal.month_name[d[1]]), xy=(0.05,0.05), xycoords='axes fraction', fontsize=12)
        ax[k*2+1].annotate('%04d %s'%(d[0],_cal.month_name[d[1]]), xy=(0.05,0.05), xycoords='axes fraction', fontsize=12)
    ax[1].set_title('Brazil')
    ax[0].set_title('Peru')
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    f.text(0.06, 0.5, 'Zonal Wind Speed [m/s]', ha='center', va='center', rotation='vertical')
    #ax[0].legend([(m1,l1),(m2,_L2)],['RENOIR','PERU'],loc='upper left', prop={'size':12}, fancybox=True)
    _plt.draw();
    #_plt.savefig('%sCompM_set_zonal_winds.png'%(dirout))
    _plt.savefig('%sCompM_set_zonal_winds.pdf'%(dirout))
    
    ##### Winds Figure Meridional
    f,((ax[1],ax[0]),(ax[3],ax[2]),(ax[5],ax[4]),(ax[7],ax[6])) = _plt.subplots(4,2, sharex=True, sharey=True)
    for k,d in enumerate(zip(yrset,monset)[0:-1]):
        ax[k*2+1].fill_between(Rt[k],RV[k]-RVv[k],RV[k]+RVv[k],alpha=0.5,linewidth=0,facecolor='m',label='MSIS-R')
        ax[k*2].fill_between(Pt[k],PV[k]-PVv[k],PV[k]+PVv[k],alpha=0.5,linewidth=0,facecolor='c',label='MSIS-P')
        (m1, caps1, _) = ax[k*2+1].errorbar(R[k].t,R[k].v2,yerr=R[k].v2v,fmt='ro',label='R')
        l1, = ax[k*2+1].plot(R[k].t,R[k].v2,'r-',alpha=lalpha)
        (m2, caps2, _) = ax[k*2].errorbar(P[k].t,P[k].v,yerr=P[k].vv,fmt='bs',label='P')
        l2, = ax[k*2].plot(P[k].t,P[k].v,'b-',alpha=lalpha)
        ax[k*2].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(tlim,[0,0],'k--')
        ax[k*2].plot(P[k].midnight,[-2000,2000],'c--')
        ax[k*2+1].plot(R[k].midnight,[-2000,2000],'m--')
        for cap in caps1+caps2:
            cap.set_markeredgewidth(markerwide)
    # Manual fix for car - no south data
    k = 3
    ax[k*2+1].fill_between(Rt[k],RV[k]-RVv[k],RV[k]+RVv[k],alpha=0.5,linewidth=0,facecolor='m',label='MSIS-R')
    ax[k*2].fill_between(Pt[k],PV[k]-PVv[k],PV[k]+PVv[k],alpha=0.5,linewidth=0,facecolor='c',label='MSIS-P')
    (m1, caps1, _) = ax[k*2+1].errorbar(R[k].t,R[k].v,yerr=R[k].vv,fmt='ro',label='RENOIR')
    l1, = ax[k*2+1].plot(R[k].t,R[k].v,'r-',alpha=lalpha)
    (m2, caps2, _) = ax[k*2].errorbar(P[k].t,P[k].v,yerr=P[k].vv,fmt='bs',label='P')
    l2, = ax[k*2].plot(P[k].t,P[k].v,'b-',alpha=lalpha)
    ax[k*2].plot(tlim,[0,0],'k--')
    ax[k*2+1].plot(tlim,[0,0],'k--')
    ax[k*2].plot(P[k].midnight,[-2000,2000],'c--')
    ax[k*2+1].plot(R[k].midnight,[-2000,2000],'m--')
    for cap in caps1+caps2:
        cap.set_markeredgewidth(markerwide)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    _plt.ylim([-100,100])
    _plt.xlim(tlim)
    ax[6].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[7].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for k,d in enumerate(zip(yrset,monset)):
        yticks = ax[k*2].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        yticks = ax[k*2+1].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k*2].annotate('%04d %s'%(d[0],_cal.month_name[d[1]]), xy=(0.05,0.05), xycoords='axes fraction', fontsize=12)
        ax[k*2+1].annotate('%04d %s'%(d[0],_cal.month_name[d[1]]), xy=(0.05,0.05), xycoords='axes fraction', fontsize=12)
    ax[1].set_title('Brazil')
    ax[0].set_title('Peru')
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    f.text(0.06, 0.5, 'Meridional Wind Speed [m/s]', ha='center', va='center', rotation='vertical')
    #ax[0].legend([(m1,l1),(m2,_L2)],['RENOIR','PERU'],loc='upper left', prop={'size':12}, fancybox=True)
    _plt.draw();
    #_plt.savefig('%sCompM_set_meridional_winds.png'%(dirout))
    _plt.savefig('%sCompM_set_meridional_winds.pdf'%(dirout))

    '''
    f,((ax[1],ax[0]),(ax[3],ax[2]),(ax[5],ax[4]),(ax[7],ax[6])) = _plt.subplots(4,2, sharex=True, sharey=True)
    for k,d in enumerate(zip(yrset,monset)[0:-1]):
        ax[k*2+1].fill_between(Rt[k],RW[k]-RWv[k],RW[k]+RWv[k],alpha=0.5,linewidth=0,facecolor='m',label='MSIS-R')
        ax[k*2].fill_between(Pt[k],PW[k]-PWv[k],PW[k]+PWv[k],alpha=0.5,linewidth=0,facecolor='c',label='MSIS-P')
        (m1, caps1, _) = ax[k*2+1].errorbar(R[k].t,R[k].w,yerr=R[k].wv,fmt='ro',label='RENOIR')
        l1, = ax[k*2+1].plot(R[k].t,R[k].w,'r-',alpha=lalpha)
        (m2, caps2, _) = ax[k*2].errorbar(P[k].t,P[k].w,yerr=P[k].wv,fmt='bs',label='PERU')
        l2, = ax[k*2].plot(P[k].t,P[k].w,'b-',alpha=lalpha)
        ax[k*2].plot(tlim,[0,0],'k--')
        ax[k*2+1].plot(tlim,[0,0],'k--')
        ax[k*2].plot(P[k].midnight,[-2000,2000],'c--')
        ax[k*2+1].plot(R[k].midnight,[-2000,2000],'m--')
        for cap in caps1+caps2:
            cap.set_markeredgewidth(markerwide)
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    _plt.ylim([-75,75])
    _plt.xlim(tlim)
    ax[6].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[7].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    for k,d in enumerate(zip(yrset,monset)):
        xticks = ax[k].xaxis.get_major_ticks()
        #xticks[0].label1.set_visible(False)
        yticks = ax[k].yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        yticks[-1].label1.set_visible(False)
        ax[k*2].annotate('%04d %s'%(d[0],_cal.month_name[d[1]]), xy=(0.05,0.05), xycoords='axes fraction', fontsize=12)
        ax[k*2+1].annotate('%04d %s'%(d[0],_cal.month_name[d[1]]), xy=(0.05,0.05), xycoords='axes fraction', fontsize=12)
    ax[1].set_title('Brazil')
    ax[0].set_title('Peru')
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    f.text(0.06, 0.5, 'Vertical Wind Speed [m/s]', ha='center', va='center', rotation='vertical')
    #ax[0].legend([(m1,l1),(m2,_L2)],['RENOIR','PERU'],loc='upper left', prop={'size':12}, fancybox=True)
    _plt.draw(); _plt.show()
    _plt.savefig('%sCompM_set_vertical_winds.png'%(dirout))
    '''
    
    
def PlotCompSAStacked3UT(SET):
    '''
    Summary:
        Plot Stacked U,V,T of RENOIR vs PERU for a single month
    
    History:
        10/1/14 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    FPIResults.SetBinTime(20)
    
    monset = _np.array([8,9,10,10])
    yrset  = _np.array([2010,2010,2010,2011])
    RToff = _np.array([50,25,25,50])
    Rl = [[216,229,232,234,235,240,242,243], [245,246,247,248,249,250,251,253,254,255,258,259,261,262,263,264,266,269,270,271,273], [274,275,276,277,279,280,281,282,303,304], [299,276,277,278,279,280,294,298,300]]
    Pl = [[224,225,226,228,230,228,231,232,237,238,239,240,243], [246,247,248,249,251,252,253,254,255,256,257,258,259,260,261,262,263,264,266,267,273,272,271],     [275,276,277,279,280,281,282,283,284,286,287,288,289,290,291,292,293,295,297,298,299,300,302], [274,275,276,277,278,302,303,304,300,274,275,277,278,279,280,282,284,295,300]]
    
    # Get Monthly Average for basis
    MR = FPIResults.BinMonthlyData('car',yrset[SET],monset[SET],LIST=Rl[SET],SPLIT=True,KP=[0,4])
    MP = FPIResults.BinMonthlyData('mrh',yrset[SET],monset[SET],LIST=Rl[SET],SPLIT=True,KP=[0,4])
    MP.t = MP.t + _dt.timedelta(minutes=5)
    MsisR = FPIResults.BinMonthlyData('hwm14',yrset[SET],monset[SET],SITELLA=MR.lla,KP=[0,4])
    MsisP = FPIResults.BinMonthlyData('hwm14',yrset[SET],monset[SET],SITELLA=MP.lla,KP=[0,4])
    Rt = MsisR.t
    RT = MsisR.T-25
    RTv = MsisR.Tv
    Ru = MsisR.u
    Ruv = MsisR.uv
    Rv = MsisR.v
    Rvv = MsisR.vv
    Pt = MsisP.t + _dt.timedelta(minutes=3)
    PT = MsisP.T-75
    PTv = MsisP.Tv
    Pu = MsisP.u
    Puv = MsisP.uv
    Pv = MsisP.v
    Pvv = MsisP.vv

    markerwide = 0

    # Get xlimits w/ data (IMPROVE)
    tstart = _np.where(_np.isnan(MR.T[:-1])*_np.isfinite(MR.T[1:]))[0][0]-1
    tend   = _np.where(_np.isfinite(MP.T[:-1])*_np.isnan(MP.T[1:]))[0][0]+1
    tlim = [MR.t[tstart],MP.t[tend]]
    
    _plt.close('all');ax={}
    # Figure
    f,(ax[0],ax[1],ax[2])  = _plt.subplots(3, sharex=True, sharey=False)
    
    # Temp Figure
    l1, = ax[2].plot(MR.t,MR.T,'-ro',label='Brazil')
    l2, = ax[2].plot(MP.t,MP.T,'-bs',label='Peru')
    ax[2].fill_between(Rt,RT-RTv,RT+RTv,alpha=0.5,linewidth=0,facecolor='m')
    ax[2].fill_between(Pt,PT-PTv,PT+PTv,alpha=0.5,linewidth=0,facecolor='c')
    #l3, = ax[2].plot(Rt,RT,'m--',label='MSIS-RENOIR')
    #l4, = ax[2].plot(Pt,PT,'c--',label='MSIS-PERU')
    f.subplots_adjust(hspace=0)
    ax[2].set_ylim([500,1000])
    _plt.xlim(tlim)
    ax[2].set_ylabel('Temperature [K]')
    yticks = ax[2].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    # Winds Figure Zonal
    ax[0].plot(MR.t,MR.u,'-ro',label='Brazil')
    ax[0].plot(MP.t,MP.u2,'-bs',label='Peru')
    ax[0].fill_between(Rt,Ru-Ruv,Ru+Ruv,alpha=0.5,linewidth=0,facecolor='m')
    ax[0].fill_between(Pt,Pu-Puv,Pu+Puv,alpha=0.5,linewidth=0,facecolor='c')
    ax[0].plot(tlim,[0,0],'k--')
    f.subplots_adjust(hspace=0)
    ax[0].set_ylim([-50,200])
    _plt.xlim(tlim)
    ax[0].set_ylabel('Zonal Wind [m/s]')
    yticks = ax[0].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    # Winds Figure Meridional
    ax[1].plot(MR.t,MR.v2,'-ro',label='Brazil')
    ax[1].plot(MP.t,MP.v,'-bs',label='Peru')
    ax[1].fill_between(Rt,Rv-Rvv,Rv+Rvv,alpha=0.5,linewidth=0,facecolor='m')
    ax[1].fill_between(Pt,Pv-Pvv,Pv+Pvv,alpha=0.5,linewidth=0,facecolor='c')
    ax[1].plot(tlim,[0,0],'k--')
    f.subplots_adjust(hspace=0)
    ax[1].set_ylim([-100,100])
    _plt.xlim(tlim)
    ax[1].set_ylabel('Meridional Wind [m/s]')
    yticks = ax[1].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    '''
    # Winds Figure Vertical
    ax[3].plot(MR.t,MR.w,'-ro',label='Brazil')
    ax[3].plot(MP.t,MP.w,'-bs',label='Peru')
    ax[k].plot(tlim,[0,0],'k--')
    f.subplots_adjust(hspace=0)
    ax[3].ylim([-75,75])
    _plt.xlim(tlim)
    ax[3].set_ylabel('Vertical Wind [m/s]')
    yticks = ax[3].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    '''
    #ax[0].annotate('%04d %s'%(YEAR,_cal.month_name[MONTH]), xy=(0.53,0.92), xycoords='axes fraction', fontsize=12)
    ax[2].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[0].set_title('%s %04d'%(_cal.month_name[monset[SET]],yrset[SET]))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    #f.text(0.06, 0.5, 'Temperature [K]', ha='center', va='center', rotation='vertical')
    ax[0].legend(loc='upper right', prop={'size':12}, fancybox=True)
    _plt.draw();
    #_plt.savefig('%sCompM_%04d-%02d_all.png'%(dirout,yrset[SET],monset[SET]))
    _plt.savefig('%sCompM_%04d-%02d_all_UT.pdf'%(dirout,yrset[SET],monset[SET]))           
         
               
def PlotCompSAStacked3LT(SET):
    '''
    Summary:
        Plot Stacked U,V,T of RENOIR vs PERU for a single month
    
    History:
        10/1/14 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    FPIResults.SetBinTime(20)
    
    dtR = -_dt.timedelta(hours=2,minutes=12)
    dtP = -_dt.timedelta(hours=4,minutes=53)
        
    monset = _np.array([8,9,10,10])
    yrset  = _np.array([2010,2010,2010,2011])
    RToff = _np.array([50,25,25,50])
    Rl = [[216,229,232,234,235,240,242,243], [245,246,247,248,249,250,251,253,254,255,258,259,261,262,263,264,266,269,270,271,273], [274,275,276,277,279,280,281,282,303,304], [299,276,277,278,279,280,294,298,300]]
    Pl = [[224,225,226,228,230,228,231,232,237,238,239,240,243], [246,247,248,249,251,252,253,254,255,256,257,258,259,260,261,262,263,264,266,267,273,272,271],     [275,276,277,279,280,281,282,283,284,286,287,288,289,290,291,292,293,295,297,298,299,300,302], [274,275,276,277,278,302,303,304,300,274,275,277,278,279,280,282,284,295,300]]
    
    # Get Monthly Average for basis
    MR = FPIResults.BinMonthlyData('car',yrset[SET],monset[SET],LIST=Rl[SET],SPLIT=True,KP=[0,4])
    MR.t = MR.t + dtR
    MP = FPIResults.BinMonthlyData('mrh',yrset[SET],monset[SET],LIST=Rl[SET],SPLIT=True,KP=[0,4])
    MP.t = MP.t + dtP
    MsisR = FPIResults.BinMonthlyData('hwm14',yrset[SET],monset[SET],SITELLA=MR.lla,KP=[0,4])
    MsisP = FPIResults.BinMonthlyData('hwm14',yrset[SET],monset[SET],SITELLA=MP.lla,KP=[0,4])
    Rt = MsisR.t + dtR
    RT = MsisR.T-25
    RTv = MsisR.Tv
    Ru = MsisR.u
    Ruv = MsisR.uv
    Rv = MsisR.v
    Rvv = MsisR.vv
    Pt = MsisP.t + dtP # + _dt.timedelta(minutes=3)
    PT = MsisP.T-75
    PTv = MsisP.Tv
    Pu = MsisP.u
    Puv = MsisP.uv
    Pv = MsisP.v
    Pvv = MsisP.vv

    markerwide = 0

    # Get xlimits w/ data (IMPROVE)
    tstart = _np.where(_np.isnan(MR.T[:-1])*_np.isfinite(MR.T[1:]))[0][0]-1
    tend   = _np.where(_np.isfinite(MP.T[:-1])*_np.isnan(MP.T[1:]))[0][0]+1
    tlim = [MR.t[tstart],MP.t[tend]]
    
    _plt.close('all');ax={}
    # Figure
    f,(ax[0],ax[1],ax[2])  = _plt.subplots(3, sharex=True, sharey=False)
    
    # Temp Figure
    l1, = ax[2].plot(MR.t,MR.T,'-ro',label='Brazil')
    l2, = ax[2].plot(MP.t,MP.T,'-bs',label='Peru')
    ax[2].fill_between(Rt,RT-RTv,RT+RTv,alpha=0.5,linewidth=0,facecolor='m')
    ax[2].fill_between(Pt,PT-PTv,PT+PTv,alpha=0.5,linewidth=0,facecolor='c')
    #l3, = ax[2].plot(Rt,RT,'m--',label='MSIS-RENOIR')
    #l4, = ax[2].plot(Pt,PT,'c--',label='MSIS-PERU')
    f.subplots_adjust(hspace=0)
    ax[2].set_ylim([500,1000])
    _plt.xlim(tlim)
    ax[2].set_ylabel('Temperature [K]')
    yticks = ax[2].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    # Winds Figure Zonal
    ax[0].plot(MR.t,MR.u,'-ro',label='Brazil')
    ax[0].plot(MP.t,MP.u2,'-bs',label='Peru')
    ax[0].fill_between(Rt,Ru-Ruv,Ru+Ruv,alpha=0.5,linewidth=0,facecolor='m')
    ax[0].fill_between(Pt,Pu-Puv,Pu+Puv,alpha=0.5,linewidth=0,facecolor='c')
    ax[0].plot(tlim,[0,0],'k--')
    f.subplots_adjust(hspace=0)
    ax[0].set_ylim([-50,200])
    _plt.xlim(tlim)
    ax[0].set_ylabel('Zonal Wind [m/s]')
    yticks = ax[0].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    # Winds Figure Meridional
    ax[1].plot(MR.t,MR.v2,'-ro',label='Brazil')
    ax[1].plot(MP.t,MP.v,'-bs',label='Peru')
    ax[1].fill_between(Rt,Rv-Rvv,Rv+Rvv,alpha=0.5,linewidth=0,facecolor='m')
    ax[1].fill_between(Pt,Pv-Pvv,Pv+Pvv,alpha=0.5,linewidth=0,facecolor='c')
    ax[1].plot(tlim,[0,0],'k--')
    f.subplots_adjust(hspace=0)
    ax[1].set_ylim([-100,100])
    _plt.xlim(tlim)
    ax[1].set_ylabel('Meridional Wind [m/s]')
    yticks = ax[1].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    '''
    # Winds Figure Vertical
    ax[3].plot(MR.t,MR.w,'-ro',label='Brazil')
    ax[3].plot(MP.t,MP.w,'-bs',label='Peru')
    ax[k].plot(tlim,[0,0],'k--')
    f.subplots_adjust(hspace=0)
    ax[3].ylim([-75,75])
    _plt.xlim(tlim)
    ax[3].set_ylabel('Vertical Wind [m/s]')
    yticks = ax[3].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    '''
    #ax[0].annotate('%04d %s'%(YEAR,_cal.month_name[MONTH]), xy=(0.53,0.92), xycoords='axes fraction', fontsize=12)
    ax[2].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[0].set_title('%s %04d'%(_cal.month_name[monset[SET]],yrset[SET]))
    f.text(0.5, 0.04, 'Time [LT]', ha='center', va='center')
    #f.text(0.06, 0.5, 'Temperature [K]', ha='center', va='center', rotation='vertical')
    ax[0].legend(loc='upper right', prop={'size':12}, fancybox=True)
    _plt.draw();
    #_plt.savefig('%sCompM_%04d-%02d_all.png'%(dirout,yrset[SET],monset[SET]))
    _plt.savefig('%sCompM_%04d-%02d_all_LT.pdf'%(dirout,yrset[SET],monset[SET])) 
               
         
def PlotCompSAStacked3UTLT(SET):
    '''
    Summary:
        Plot Stacked U,V,T of RENOIR vs PERU for a single month
    
    History:
        10/1/14 -- Written by DJF (dfisher2@illinois.edu)
    
    '''
    FPIResults.SetBinTime(20)
    
    dtR = -_dt.timedelta(hours=2,minutes=12)
    dtP = -_dt.timedelta(hours=4,minutes=53)
    monset = _np.array([8,9,10,10])
    yrset  = _np.array([2010,2010,2010,2011])
    RToff = _np.array([50,25,25,50])
    Rl = [[216,229,232,234,235,240,242,243], [245,246,247,248,249,250,251,253,254,255,258,259,261,262,263,264,266,269,270,271,273], [274,275,276,277,279,280,281,282,303,304], [299,276,277,278,279,280,294,298,300]]
    Pl = [[224,225,226,228,230,228,231,232,237,238,239,240,243], [246,247,248,249,251,252,253,254,255,256,257,258,259,260,261,262,263,264,266,267,273,272,271],     [275,276,277,279,280,281,282,283,284,286,287,288,289,290,291,292,293,295,297,298,299,300,302], [274,275,276,277,278,302,303,304,300,274,275,277,278,279,280,282,284,295,300]]
    
    # Get Monthly Average for basis
    MR = FPIResults.BinMonthlyData('car',yrset[SET],monset[SET],LIST=Rl[SET],SPLIT=True,KP=[0,4])
    MR2t = MR.t + dtR
    MP = FPIResults.BinMonthlyData('mrh',yrset[SET],monset[SET],LIST=Rl[SET],SPLIT=True,KP=[0,4])
    MP.t = MP.t + _dt.timedelta(minutes=5)
    MP2t = MP.t + dtP
    MsisR = FPIResults.BinMonthlyData('hwm14',yrset[SET],monset[SET],SITELLA=MR.lla,KP=[0,4])
    MsisP = FPIResults.BinMonthlyData('hwm14',yrset[SET],monset[SET],SITELLA=MP.lla,KP=[0,4])
    Rt = MsisR.t
    Rt2 = MsisR.t + dtR
    RT = MsisR.T-25
    RTv = MsisR.Tv
    Ru = MsisR.u
    Ruv = MsisR.uv
    Rv = MsisR.v
    Rvv = MsisR.vv
    Pt = MsisP.t + _dt.timedelta(minutes=3)
    Pt2 = MsisP.t + dtP
    PT = MsisP.T-75
    PTv = MsisP.Tv
    Pu = MsisP.u
    Puv = MsisP.uv
    Pv = MsisP.v
    Pvv = MsisP.vv

    markerwide = 0

    # Get xlimits w/ data (IMPROVE)
    tstart = _np.where(_np.isnan(MR.T[:-1])*_np.isfinite(MR.T[1:]))[0][0]-1
    tend   = _np.where(_np.isfinite(MP.T[:-1])*_np.isnan(MP.T[1:]))[0][0]+1
    tlim = [MR.t[tstart],MP.t[tend]]
    tlim2 = [MR2t[tstart],MP2t[tend]]
    
    _plt.close('all');ax={}
    # Figure
    f,((ax[0],ax[3]),(ax[1],ax[4]),(ax[2],ax[5]))  = _plt.subplots(3,2, sharex='col', sharey='row')
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
        
    # UT  ##################
    # Temp Figure
    l1, = ax[2].plot(MR.t,MR.T,'-ro',label='Brazil')
    l2, = ax[2].plot(MP.t,MP.T,'-bs',label='Peru')
    ax[2].fill_between(Rt,RT-RTv,RT+RTv,alpha=0.5,linewidth=0,facecolor='m')
    ax[2].fill_between(Pt,PT-PTv,PT+PTv,alpha=0.5,linewidth=0,facecolor='c')
    #l3, = ax[2].plot(Rt,RT,'m--',label='MSIS-RENOIR')
    #l4, = ax[2].plot(Pt,PT,'c--',label='MSIS-PERU')
    ax[2].set_ylim([500,1000])
    ax[2].set_ylabel('Temperature [K]')
    yticks = ax[2].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    # Winds Figure Zonal
    ax[0].plot(MR.t,MR.u,'-ro',label='Brazil')
    ax[0].plot(MP.t,MP.u2,'-bs',label='Peru')
    ax[0].fill_between(Rt,Ru-Ruv,Ru+Ruv,alpha=0.5,linewidth=0,facecolor='m')
    ax[0].fill_between(Pt,Pu-Puv,Pu+Puv,alpha=0.5,linewidth=0,facecolor='c')
    ax[0].plot(tlim,[0,0],'k--')
    ax[0].set_ylim([-50,200])
    ax[0].set_ylabel('Zonal Wind [m/s]')
    yticks = ax[0].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    # Winds Figure Meridional
    ax[1].plot(MR.t,MR.v2,'-ro',label='Brazil')
    ax[1].plot(MP.t,MP.v,'-bs',label='Peru')
    ax[1].fill_between(Rt,Rv-Rvv,Rv+Rvv,alpha=0.5,linewidth=0,facecolor='m')
    ax[1].fill_between(Pt,Pv-Pvv,Pv+Pvv,alpha=0.5,linewidth=0,facecolor='c')
    ax[1].plot(tlim,[0,0],'k--')
    ax[1].set_ylim([-100,100])
    ax[1].set_ylabel('Meridional Wind [m/s]')
    ax[1].set_xlim(tlim)
    yticks = ax[1].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)
    
    
    # LT  ##################
    # Temp Figure
    l1, = ax[2+3].plot(MR2t,MR.T,'-ro',label='Brazil')
    l2, = ax[2+3].plot(MP2t,MP.T,'-bs',label='Peru')
    ax[2+3].fill_between(Rt2,RT-RTv,RT+RTv,alpha=0.5,linewidth=0,facecolor='m')
    ax[2+3].fill_between(Pt2,PT-PTv,PT+PTv,alpha=0.5,linewidth=0,facecolor='c')
    
    # Winds Figure Zonal
    ax[0+3].plot(MR2t,MR.u,'-ro',label='Brazil')
    ax[0+3].plot(MP2t,MP.u2,'-bs',label='Peru')
    ax[0+3].fill_between(Rt2,Ru-Ruv,Ru+Ruv,alpha=0.5,linewidth=0,facecolor='m')
    ax[0+3].fill_between(Pt2,Pu-Puv,Pu+Puv,alpha=0.5,linewidth=0,facecolor='c')
    ax[0+3].plot(tlim2,[0,0],'k--')

    
    # Winds Figure Meridional
    ax[1+3].plot(MR2t,MR.v2,'-ro',label='Brazil')
    ax[1+3].plot(MP2t,MP.v,'-bs',label='Peru')
    ax[1+3].fill_between(Rt2,Rv-Rvv,Rv+Rvv,alpha=0.5,linewidth=0,facecolor='m')
    ax[1+3].fill_between(Pt2,Pv-Pvv,Pv+Pvv,alpha=0.5,linewidth=0,facecolor='c')
    ax[1+3].plot(tlim2,[0,0],'k--')
    ax[1+3].set_xlim(tlim2)


    #ax[0].annotate('%04d %s'%(YEAR,_cal.month_name[MONTH]), xy=(0.53,0.92), xycoords='axes fraction', fontsize=12)
    ax[2].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[5].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[0].set_title('%s %04d'%(_cal.month_name[monset[SET]],yrset[SET]))
    f.text(0.3, 0.04, 'Time [UTC]', ha='center', va='center')
    f.text(0.7, 0.04, 'Time [SLT]', ha='center', va='center')
    #f.text(0.06, 0.5, 'Temperature [K]', ha='center', va='center', rotation='vertical')
    ax[3].legend(loc='upper right', prop={'size':12}, fancybox=True)
    _plt.draw();
    #_plt.savefig('%sCompM_%04d-%02d_all.png'%(dirout,yrset[SET],monset[SET]))
    _plt.savefig('%sCompM_%04d-%02d_all.pdf'%(dirout,yrset[SET],monset[SET]))  
         
               
def PlotSpaghettiSA(SITE,SET):
    '''
    Summary:
        Plots all raw data for one month in spaghetti plot!
        TODO: multicolor lines w/ legend or transparent blue lines...
    
    Inputs:
        SITE = sites of interest, e.g. 'UAO'
        YEAR = year, e.g. 2013
        MONTH = month to plot, e.g. 4
        SPLIT = Split look directions in binning [default = False]
        LIST = List of doys to use in averaging.

    History:
        10/17/14 -- Written by DJF (dfisher2@illinois.edu)
    '''
    FPIResults.SetBinTime(30)

    monset = _np.array([8,9,10,10])
    yrset  = _np.array([2010,2010,2010,2011])
    if SITE is 'car':
        Toff = _np.array([50,25,25,50])
        CV = True
        dlist = [[216,229,232,234,235,240,242,243], [245,246,247,248,249,250,251,253,254,255,258,259,261,262,263,264,266,269,270,271,273], [274,275,276,277,279,280,281,282,303,304], [299,276,277,278,279,280,294,298,300]]
    elif SITE is 'mrh':
        Toff = _np.array([75,75,75,75])
        CV = False
        dlist = [[224,225,226,228,230,228,231,232,237,238,239,240,243], [246,247,248,249,251,252,253,254,255,256,257,258,259,260,261,262,263,264,266,267,273,272,271],     [275,276,277,279,280,281,282,283,284,286,287,288,289,290,291,292,293,295,297,298,299,300,302], [274,275,276,277,278,302,303,304,300,274,275,277,278,279,280,282,284,295,300]]
    
    # Get Monthly Average for basis
    D = {}
    dn = _dt.datetime(yrset[SET],monset[SET],1)
    for k,d in enumerate(dlist[SET]):
        L = FPIResults.L2.GetLevel2(SITE,_dt.datetime(yrset[SET],1,1)+_dt.timedelta(days=d-1))
        FPIResults.FilterData(L); D[k] = L
    M = FPIResults.BinMonthlyData(SITE,yrset[SET],monset[SET],SPLIT=True,CV=CV,LIST=dlist[SET],KP=[0,4])
    M.t = M.t + _dt.timedelta(days=(_dt.datetime(yrset[SET]-1,1,2+_cal.isleap(yrset[SET]-1))-_dt.datetime(M.t[0].year,1,1)).days)
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
            doy = (D[d][0].dn-_dt.datetime(yrset[SET],1,1)).days+1
        except:
            continue
        # Temps
        for l in [y for y in D[d]]:
            ax[3].fill_between(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy-1)]),l.T-l.Te,l.T+l.Te,alpha=lalpha,linewidth=0,facecolor='k')
        # Zonal
        for l in [y for y in D[d] if ('East' in y.key and SITE is 'car') or ('West' in y.key and SITE is 'mrh') or ('East' in y.key and SITE is 'mrh' and yrset[SET] == 2011) or ('CV' in y.key and '2' in y.key and SITE is 'car') or ('CV' in y.key and '1' in y.key and SITE is 'car' and yrset[SET] == 2011)]:
            ax[0].fill_between(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy-1)]),l.u-l.ue,l.u+l.ue,alpha=lalpha,linewidth=0,facecolor='k')
        # Meridional
        for l in [y for y in D[d] if ('South' in y.key and SITE is 'car' and yrset[SET] != 2011) or ('North' in y.key and SITE is 'mrh') or ('CV' in y.key and '2' in y.key and SITE is 'car' and yrset[SET] != 2011) or ('CV' in y.key and '1' in y.key and SITE is 'car' and yrset[SET] == 2011) or ('North' in y.key and SITE is 'car' and yrset[SET] == 2011)]:
            ax[1].fill_between(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy-1)]),l.v-l.ve,l.v+l.ve,alpha=lalpha,linewidth=0,facecolor='k')
    # Overlay Monthly Average
    ax[3].errorbar(M.t,M.T,yerr=M.Tv,fmt='b.-')
    if yrset[SET] == 2011:
        ax[0].errorbar(M.t,M.u,yerr=M.uv,fmt='b.-')
        ax[1].errorbar(M.t,M.v,yerr=M.vv,fmt='b.-')
    else:
        if SITE is 'car':
            ax[0].errorbar(M.t,M.u,yerr=M.uv,fmt='b.-')
            ax[1].errorbar(M.t,M.v2,yerr=M.v2v,fmt='b.-')
        elif SITE is 'mrh':
            ax[0].errorbar(M.t,M.u2,yerr=M.u2v,fmt='b.-')
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
    ax[0].set_title('%s Spaghetti %s %04d'%(SITE.upper(),_cal.month_name[monset[SET]],yrset[SET]))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    #ax[0].legend(loc=2, prop={'size':12}, bbox_to_anchor=(1.1, 0.5), fancybox=True)
    _plt.draw();
    #_plt.savefig('%s%s_Spaghetti_%04d-%02d.png'%(dirout,SITE,yrset[SET],monset[SET]))
    _plt.savefig('%s%s_Spaghetti_%04d-%02d.pdf'%(dirout,SITE,yrset[SET],monset[SET]))
    

def PlotGridMonthV(SITE,YEAR,MONTH,SPLIT=True):
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
        D[k] = _L2.GetLevel2(SITE,dn+_dt.timedelta(days=d))
    #M = BinMonthlyData(SITE,YEAR,MONTH,SPLIT)
    #M.t = M.t + _dt.timedelta(seconds=(_dt.datetime(YEAR-1,1,1)-_dt.datetime(M.t[0].year,1,1)).total_seconds())

    markerwide = 0
    lalpha = .3
    lw=2
    
    # Get xlimits w/ data (IMPROVE)
    #tlim = [(D[0][0].allt[0]-_dt.timedelta(days=(dn-_dt.datetime(YEAR,1,1)).days+1,minutes=80)).astimezone(_utc).replace(tzinfo=None), (D[0][0].allt[-1]-_dt.timedelta(days=(dn-_dt.datetime(YEAR,1,1)).days+1,minutes=50)+_dt.timedelta(hours=1)).astimezone(_utc).replace(tzinfo=None)]
    tlim = [_dt.datetime(YEAR,1,1)-_dt.timedelta(hours=4), _dt.datetime(YEAR,1,1)+_dt.timedelta(hours=12)]
    
    _plt.close('all');ax={}
    
    # Figure 4 
    ax={};f,((ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]),(ax[7],ax[8],ax[9],ax[10],ax[11],ax[12],ax[13]),(ax[14],ax[15],ax[16],ax[17],ax[18],ax[19],ax[20]),(ax[21],ax[22],ax[23],ax[24],ax[25],ax[26],ax[27]),(ax[28],ax[29],ax[30],ax[31],ax[32],ax[33],ax[34]))  = _plt.subplots(5,7, sharex=True, sharey=True)
    # Plot Data
    for d in D:
        try:
            doy = (D[d][0].dn-_dt.datetime(YEAR,1,1)).days+1
        except:
            continue
        # Vertical
        for l in [y for y in D[d] if ('CAR_Zenith' in y.key)]:
            (l1, caps1, _) = ax[d].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.we,fmt='go')
            for cap in caps1:
                cap.set_markeredgewidth(markerwide)
            m1, = ax[d].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,'g-',linewidth=lw,alpha=lalpha)
        # Vertical
        for l in [y for y in D[d] if ('CAJ_Zenith' in y.key)]:
            (l2, caps2, _) = ax[d].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.we,fmt='bo')
            for cap in caps2:
                cap.set_markeredgewidth(markerwide)
            m2, = ax[d].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,'b-',linewidth=lw,alpha=lalpha)
        # Inline
        for l in [y for y in D[d] if ('IN' in y.key)]:
            (l3, caps3, _) = ax[d].errorbar(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,yerr=l.we,fmt='ro')
            for cap in caps3:
                cap.set_markeredgewidth(markerwide)
            m3, = ax[d].plot(_np.array([x.astimezone(FPIResults._utc).replace(tzinfo=None) for x in l.t1-_dt.timedelta(days=doy)]),l.w,'r-',linewidth=lw,alpha=lalpha)
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
    try:
        ax[34].legend([(l1,m1)],['Car'])
        ax[34].legend([(l2,m2)],['Caj'])
        ax[34].legend([(l3,m3)],['Inline'])
    except:
        print 'Not enough data for legend'
    ax[0].xaxis.set_major_formatter(_md.DateFormatter('%H'))
    ax[3].set_title('%s Vertical Wind Grid %s %04d'%(SITE,_cal.month_name[MONTH],YEAR))
    f.text(0.5, 0.04, 'Time [UT]', ha='center', va='center')
    _plt.draw();
    _plt.savefig('%s%s_GridW_%04d-%02d.png'%(dirout,SITE,YEAR,MONTH))
    

def PDF2EPS():

    import os
    from glob import glob
    
    print "Transform PDF -> EPS"
    files = glob('/rdata/airglow/database/L2/plots/*.pdf')
    for f in files:
        os.system('pdftops -eps %s %s'%(f,f[:-4]+'.eps'))
           
           
           
if __name__=="__main__":

    print 'convert PDF -> EPS'

    
    
