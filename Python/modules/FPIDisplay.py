import glob
import BoltwoodSensor
import FPI
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.dates as dates
from numpy import ma
from matplotlib.ticker import FuncFormatter
import datetime
from scipy import interpolate
import os as os
import FPIprocessLevel2


from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import pytz
from pytz import timezone
import os.path
from calendar import Calendar, SUNDAY
import ephem
import fpiinfo
import FPI

import sys

def MonthlySummary(site_name, year, month):
#
# Generates PDF pages summarizing a month of FPI data for the requested site.  Pulls Kp values and displays
# data coverage on a plot of the 3-hour Kp.  Displays meridional and zonal winds along with monthly averages.
# Displays vertical winds.  Displays temperatures with montly averages.  Displays whether the moon is up/down.
#
# INPUTS:
#	site_name - site name, matching those defined in fpiiinfo
#	year - year to plot
#	month - month to plot [1-12]
# HISTORY:
#	Written on 29 January 2015 by Jonathan J. Makela (jmakela@illinois.edu)

    from pyglow import pyglow
    
    # Reduce the font size for the x-axis and set the output to be standard paper size
    matplotlib.rc('xtick', labelsize=8) 

    site = fpiinfo.get_site_info(site_name)

    sitetz = timezone(site['Timezone'])

    instrument = fpiinfo.get_instr_at(site_name,datetime.datetime(year,month,1))[0]

    obs = ephem.Observer()
    obs.pressure = 0
    obs.horizon = '0'
    obs.lat, obs.lon = '%.2f' % site['Location'][0], '%.2f' % site['Location'][1]

    # Determines files to be averaged
    f = glob.glob('/rdata/airglow/fpi/results/%s_%s_%04d%02d*.npz' % (instrument,site_name,year,month))
    
    # Calculate sunrise and sunset time
    npzfile = np.load(f[-1],allow_pickle=True)
    d = npzfile['FPI_Results'].ravel()[0]
    del npzfile.f
    npzfile.close()
    psst = ephem.Date(obs.previous_setting(ephem.Sun(), start=d['sky_times'][0].astimezone(pytz.utc))).datetime()
    nsrt = ephem.Date(obs.next_rising(ephem.Sun(), start=d['sky_times'][0].astimezone(pytz.utc))).datetime()
    
#    t0 = psst.hour - 1
#    t1 = nsrt.hour + 1
    t0 = psst.replace(tzinfo=pytz.utc).astimezone(sitetz).hour - 1
    t1= nsrt.replace(tzinfo=pytz.utc).astimezone(sitetz).hour + 1


    if t0 > t1:
        t1 = t1+24
    
    bins = np.arange(t0,t1,0.5)
    if bins[0] < 12:
        bins = bins+24.
    
    center_time, (bin_T, bin_eT), (bin_U, bin_eU), (bin_U2, bin_eU2), (bin_V, bin_eV), (bin_V2, bin_eV2) = FPI.average_data(f, bins=bins)

    rcParams['figure.figsize']= 8.5,11

    # The start date
#    base = sitetz.localize(datetime.datetime(2014,month,1))
    base = pytz.utc.localize(datetime.datetime(year,month,1))
    pdf_pages = PdfPages('/home/jmakela/%s_%s_%04d%02d.pdf' % (instrument,site_name,year,month))

    # Number of dates to consider per page
    span = 7

    # Find the number of days in the requested month
    cal = Calendar(SUNDAY)
    days_in_month = max(cal.itermonthdays(year,month))
    dateList = [base + datetime.timedelta(days = x) for x in range(0,days_in_month,1)]

    grid_size = (8, 4)

    for dn in dateList:
        if mod(dn.day-1,span) == 0:
            # Create a new page
            fig = figure()

            # Files read in
            files = []

            agraph = subplot2grid(grid_size, (0,0), rowspan=1, colspan=4)
            # Plot Kp over the time period
            colors = ['k|-','b|-','r|-','g|-','c|-']
            kpc = ['green','green','green','green','yellow','red','red','red','red']
            dn_end = dn + datetime.timedelta(days = span)
            thrs = (dn_end-dn+datetime.timedelta(days=1)).total_seconds()/3600.
            for hrs in linspace(0,thrs,thrs/3+1):
                t = (dn+datetime.timedelta(hours=hrs)).astimezone(pytz.utc)
                try:
                    kpi = pyglow.Point(t,0,0,250).kp
                except:
                    kpi = 0
                agraph.fill([t,t,t+datetime.timedelta(hours=3),t+datetime.timedelta(hours=3)],[0,kpi,kpi,0],color=kpc[int(kpi)],alpha=0.3,linewidth=1.0)

            agraph.set_ylim([0,9])
            agraph.xaxis.set_major_formatter(DateFormatter('%b %d'))
            agraph.set_yticks([0,3,6,9])
            agraph.set_title('%s %s %d' % (instrument, site_name, year))

        # Append the next file
        files.append('/rdata/airglow/fpi/results/%s_%s_%s.npz' % (instrument, site_name, dn.strftime('%Y%m%d')))

        ugraph = subplot2grid(grid_size, (mod(dn.day-1,span)+1, 0), rowspan=1, colspan=1)            
        vgraph = subplot2grid(grid_size, (mod(dn.day-1,span)+1, 1), rowspan=1, colspan=1)
        wgraph = subplot2grid(grid_size, (mod(dn.day-1,span)+1, 2), rowspan=1, colspan=1)
        bgraph = subplot2grid(grid_size, (mod(dn.day-1,span)+1, 3), rowspan=1, colspan=1)

        if os.path.isfile(files[-1]):

            try:
                # Calculate sunrise and sunset time
                npzfile = np.load(files[-1],allow_pickle=True)
                d = npzfile['FPI_Results'].ravel()[0]
                del npzfile.f
                npzfile.close()
        #        psst = ephem.Date(obs.previous_setting(ephem.Sun(), start=dn+datetime.timedelta(days=1))).datetime()
        #        nsrt = ephem.Date(obs.next_rising(ephem.Sun(), start=dn+datetime.timedelta(days=1))).datetime()
                psst = ephem.Date(obs.previous_setting(ephem.Sun(), start=d['sky_times'][0].astimezone(pytz.utc))).datetime()
                nsrt = ephem.Date(obs.next_rising(ephem.Sun(), start=d['sky_times'][0].astimezone(pytz.utc))).datetime()

                obs.date = psst
                moonset = ephem.Moon()
                moonset.compute(obs)

            #    print psst, 'Moon: ', moonset.alt, moonset.alt < obs.horizon
            #    print 'Previous moonrise: ', ephem.Date(obs.previous_rising(ephem.Moon(), start=psst)).datetime()
            #    print 'Next moonrise: ', ephem.Date(obs.next_rising(ephem.Moon(), start=psst)).datetime()

                if moonset.alt < obs.horizon:
                    # moon is not yet up
            #        print 'down at sunset'
                    rt = ephem.Date(obs.next_rising(ephem.Moon(), start=psst)).datetime()
                else:
                    # moon is up
            #        print 'up at sunset'
                    rt = ephem.Date(obs.previous_rising(ephem.Moon(), start=psst)).datetime()

                obs.date = nsrt
                moonrise = ephem.Moon()
                moonrise.compute(obs)
                if moonrise.alt > obs.horizon:
                    # moon is up at sunrise
            #        print 'up at sunrise'
                    st = ephem.Date(obs.next_setting(ephem.Moon(), start=nsrt)).datetime()
                else:
            #        print 'down at sunrise'
                    st = ephem.Date(obs.previous_setting(ephem.Moon(), start=nsrt)).datetime()

                if st < rt:
                    st = ephem.Date(obs.next_setting(ephem.Moon(), start=nsrt)).datetime()

                FPIDisplay.PlotDay(files[-1],directions=['East','West'],Doppler_Fig = fig, Doppler_Graph = ugraph, Temperature_Fig = fig, Temperature_Graph = bgraph, cull=True);
                ugraph.legend().set_visible(False)
                ugraph.set_ylabel('')
                if mod(dn.day-1,span) == 0:
                    ugraph.set_title('Zonal')
                else:
                    ugraph.set_title('')
                ugraph.set_xlabel('')
                ugraph.set_yticks(arange(-200,201,100))
                setp(ugraph.get_lines(),markersize=4)
                ugraph.fill_between([rt,st],[-200,-200],[200,200],alpha=0.1,linewidth=0,facecolor='k')
                ugraph.set_xlim([psst,nsrt])

                ind = isfinite(bin_U)
                center_dn = np.array([dn.astimezone(pytz.utc)+datetime.timedelta(hours = x) for x in center_time])
                center_dn = np.array([x.replace(tzinfo=sitetz) for x in center_dn])
                ugraph.fill_between(center_dn[ind],bin_U[ind]-bin_eU[ind],bin_U[ind]+bin_eU[ind],alpha=0.5,linewidth=0,facecolor='k')
                ugraph.set_xlim([psst, nsrt])

                FPIDisplay.PlotDay(files[-1],directions=['North','South'],Doppler_Fig = fig, Doppler_Graph = vgraph, Temperature_Fig = fig, Temperature_Graph = bgraph, cull=True);
                vgraph.legend().set_visible(False)
                vgraph.set_ylabel('')
                if mod(dn.day-1,span) == 0:
                    vgraph.set_title('Meridional')
                else:
                    vgraph.set_title('')
                vgraph.set_xlabel('')
                vgraph.set_yticks(arange(-200,201,100))
                setp(vgraph.get_lines(),markersize=4)
                vgraph.fill_between([rt,st],[-200,-200],[200,200],alpha=0.1,linewidth=0,facecolor='k')
                vgraph.set_xlim([psst,nsrt])

                ind = isfinite(bin_V)
                center_dn = np.array([dn.astimezone(pytz.utc)+datetime.timedelta(hours = x) for x in center_time])
                center_dn = np.array([x.replace(tzinfo=sitetz) for x in center_dn])
                vgraph.fill_between(center_dn[ind],bin_V[ind]-bin_eV[ind],bin_V[ind]+bin_eV[ind],alpha=0.5,linewidth=0,facecolor='k')
                vgraph.set_xlim([psst, nsrt])

                FPIDisplay.PlotDay(files[-1],directions=['Zenith'],Doppler_Fig = fig, Doppler_Graph = wgraph, Temperature_Fig = fig, Temperature_Graph = bgraph, cull=True);
                wgraph.legend().set_visible(False)
                wgraph.set_ylabel('')
                if mod(dn.day-1,span) == 0:
                    wgraph.set_title('Vertical')
                else:
                    wgraph.set_title('')
                wgraph.set_xlabel('')
                wgraph.set_yticks(arange(-200,201,100))
                setp(wgraph.get_lines(),markersize=4)
                wgraph.fill_between([rt,st],[-200,-200],[200,200],alpha=0.1,linewidth=0,facecolor='k')
                wgraph.set_xlim([psst, nsrt])


                bgraph.legend().set_visible(False)
                bgraph.set_ylabel('')
                if mod(dn.day-1,span) == 0:
                    bgraph.set_title('Temperature')
                else:
                    bgraph.set_title('')
                bgraph.set_xlabel('')
                bgraph.set_yticks(arange(500,1501,250))
                setp(bgraph.get_lines(),markersize=4)
                bgraph.fill_between([rt,st],[500,500],[1500,1500],alpha=0.1,linewidth=0,facecolor='k')
                bgraph.set_xlim([psst, nsrt])

                ind = isfinite(bin_T)
                center_dn = np.array([dn.astimezone(pytz.utc)+datetime.timedelta(hours = x) for x in center_time])
                center_dn = np.array([x.replace(tzinfo=sitetz) for x in center_dn])
                bgraph.fill_between(center_dn[ind],bin_T[ind]-bin_eT[ind],bin_T[ind]+bin_eT[ind],alpha=0.5,linewidth=0,facecolor='k')

                bgraph2 = bgraph.twinx()
                if center_dn[0].astimezone(pytz.utc).day == center_dn[-1].astimezone(pytz.utc).day:
                    bgraph2.set_ylabel(dn.astimezone(pytz.utc).strftime('%b %d %Y'))
                else:
                    bgraph2.set_ylabel('%s-%s' % (center_dn[0].astimezone(pytz.utc).strftime('%b %d'), center_dn[-1].astimezone(pytz.utc).strftime('%d %Y')))
                bgraph2.set_yticks([])

            except:
                error = 'Error'
                # No files, still create graphs, and label them if they are the top.
                if mod(dn.day-1,span) == 0:
                        ugraph.set_title('Zonal')
                        vgraph.set_title('Meridional')
                        wgraph.set_title('Vertical')
                        bgraph.set_title('Temperature')

                ugraph.yaxis.set_visible(False)
                vgraph.yaxis.set_visible(False)
                wgraph.yaxis.set_visible(False)
                bgraph.yaxis.set_visible(False)
                ugraph.xaxis.set_visible(False)
                vgraph.xaxis.set_visible(False)
                wgraph.xaxis.set_visible(False)
                bgraph.xaxis.set_visible(False)
        else:
            # No files, still create graphs, and label them if they are the top.
            if mod(dn.day-1,span) == 0:
                    ugraph.set_title('Zonal')
                    vgraph.set_title('Meridional')
                    wgraph.set_title('Vertical')
                    bgraph.set_title('Temperature')

            ugraph.yaxis.set_visible(False)
            vgraph.yaxis.set_visible(False)
            wgraph.yaxis.set_visible(False)
            bgraph.yaxis.set_visible(False)
            ugraph.xaxis.set_visible(False)
            vgraph.xaxis.set_visible(False)
            wgraph.xaxis.set_visible(False)
            bgraph.xaxis.set_visible(False)


        # Plot summary of when data are available on Kp plot
        cloudthreshold = -25.
        winderrorlimit = 50.
        temperrorlimit = 100.
        k_inst = 6

        for f in files:
            if os.path.isfile(f):
                npzfile = np.load(f,allow_pickle=True)
                d = npzfile['FPI_Results'].ravel()[0]
                del npzfile.f
                npzfile.close()

                ind = []
                try:
                    # Cloud Check
                    ind += list(np.where(d['Clouds']['mean'] > cloudthreshold)[0])
                except:
                    error = 'no cloudsensor'

                # Error Check
                ind += list(np.where(d['sigma_LOSwind'] > winderrorlimit)[0])
                ind += list(np.where(d['sigma_T'] > temperrorlimit)[0])
                # Direction Check
                ind += [x for x in range(len(d['direction'])) if 'Unknown' in d['direction'][x]]
                ind += [x for x in range(len(d['direction'])) if 'None' in d['direction'][x]]
                try:
                    # Plot laser start/stop times
                    agraph.plot([d['laser_times'][0].astimezone(pytz.utc), d['laser_times'][-1].astimezone(pytz.utc)],\
                                [k_inst, k_inst], 'k,-')
                except:
                    error = 'no laser'

                # Plot good X exposures as \
                good = list(set(range(len(d['direction'])))-set(ind))
                times = [x.astimezone(pytz.utc) for x in d['sky_times'][good]]
                agraph.plot(times,len(times)*[k_inst], 'b|')

        # If this is the last plot on the page, save it
        if mod(dn.day-1,span) == span-1:
            tight_layout()
            pdf_pages.savefig(fig)

    # Save the last page
    tight_layout()
    pdf_pages.savefig(fig)        
    pdf_pages.close()

def NetworkSummary(network, times, bin_time = np.arange(17,32,0.5),
                Tmin=500, Tmax=1500, Dmin=-200, Dmax=200, Imin=0, Imax=200, cloudy_temperature=-15.0,
                reference='laser'):
#
# Function to display a two dimensional summary of data.  The 
# y-axis is hours while the x-axis is the day.  The function
# stacks data in columns so daily/seasonal trends can be observed
#
# INPUTS:
#   network - the network to be analyzed (currently 'renoir' or 'nation')
#   times - a list of datetimes for plotting (x-axis)
# OPTIONAL INPUTS:
#   bin_time - array of times to bin data on to (default np.arange(17,32,0.5))
#   Tmin, Tmax - the min/max values for temperatures to be plotted (default 500,1500)
#   Dmin, Dmax - the min/max values for the Doppler shifts to be plotted (default -200, 200)
#   reference - the type of Doppler reference to use
# OUTPUTS:
#   (Temperature_Fig, Temperature_Ax) - references to figure and axis of temperature plot
#   (Zonal_Fig, Zonal_Ax) - references to figure and axis of zonal plot
#   (Meridional_Fig, Meridional_Ax) - references to figure and axis of meridional plot
#   (Intensity_Fig, Intesity_Ax) - references to the figure and axis for the airglow intensity plot
#
# HISTORY:
#   Written by Jonathan J. Makela on 2 Dec 2012
#   Added Intensity figures on 28 May 2013 (jjm)
#   Corrected to work with new file format on 22 Aug 2013
#   Extended to network capable on 4 Nov 2013

    # Create arrays to hold data
    N = np.size(times)

    all_T = np.zeros((len(bin_time)-1,N))*np.nan
    all_eT = np.zeros((len(bin_time)-1,N))*np.nan

    all_U = np.zeros((len(bin_time)-1,N))*np.nan
    all_eU = np.zeros((len(bin_time)-1,N))*np.nan

    all_V = np.zeros((len(bin_time)-1,N))*np.nan
    all_eV = np.zeros((len(bin_time)-1,N))*np.nan

    all_W = np.zeros((len(bin_time)-1,N))*np.nan
    all_eW = np.zeros((len(bin_time)-1,N))*np.nan

    all_I = np.zeros((len(bin_time)-1,N))*np.nan
    all_eI = np.zeros((len(bin_time)-1,N))*np.nan

    # Center of the bins (used in plotting)
    center_time = (bin_time[0:-1]+bin_time[1:])/2.

    # Work with the data from the times provided
    for myt in times:
        # Load the L2 data product
        t = dates.num2date(myt)
        L2data = FPIprocessLevel2.GetLevel2('renoir',datetime.datetime(t.year,t.month,t.day))

        doy = (t-dates.num2date(times[0])).days
        
        zonal_t = []
        zonal_u = []
        zonal_ue = []
        
        meridional_t = []
        meridional_v = []
        meridional_ve = []
        
        for look in L2data:
            # TODO: ALLOW USER TO CHOSE WHAT DIRECTIONS TO INCLUDE IN WHAT CALCULATION
            if('zenith' in str.lower(look.key)):
                if(len(look.w) > 0):
                    # Bin the temperature data
                    (bin_T, bin_eT) = FPI.bin_and_mean(look.t1,look.T,look.Te,bin_time)
                    all_T[:,doy-1] = bin_T
                    all_eT[:,doy-1] = bin_eT
                    
                    # Bin the intensity data
                    (bin_I, bin_eI) = FPI.bin_and_mean(look.t1,look.i,look.ie,bin_time)
            
                    # Save to array
                    all_I[:,doy-1] = bin_I
                    all_eI[:,doy-1] = bin_eI
            
                    (bin_dop, bin_edop) = FPI.bin_and_mean(look.t1,look.w,look.we,bin_time)
            
                    # Save to array
                    all_W[:,doy-1] = bin_dop
                    all_eW[:,doy-1] = bin_edop

            if 'east' in str.lower(look.key) or 'west' in str.lower(look.key) or 'cv' in str.lower(look.key):
                if(len(look.u) > 0):
                    zonal_t.extend(look.t1)
                    zonal_u.extend(look.u)
                    zonal_ue.extend(look.ue)
                    
            if 'north' in str.lower(look.key) or 'south' in str.lower(look.key) or 'cv' in str.lower(look.key):
                if(len(look.v) > 0):
                    meridional_t.extend(look.t1)
                    meridional_v.extend(look.v)
                    meridional_ve.extend(look.ve)

        # Bin all zonal measurements for this date
        if(len(zonal_t) > 0):
            zonal_t = np.array(zonal_t)
            zonal_u = np.array(zonal_u)
            zonal_ue = np.array(zonal_ue)
            (bin_dop, bin_edop) = FPI.bin_and_mean(zonal_t,zonal_u,zonal_ue,bin_time)
                
            # Save to array
            all_U[:,doy-1] = bin_dop
            all_eU[:,doy-1] = bin_edop
            
        # Bin all meridional measurements for this date
        if(len(meridional_t) > 0):
            meridional_t = np.array(meridional_t)
            meridional_v = np.array(meridional_v)
            meridional_ve = np.array(meridional_ve)
            (bin_dop, bin_edop) = FPI.bin_and_mean(meridional_t,meridional_v,meridional_ve,bin_time)
                
            # Save to array
            all_V[:,doy-1] = bin_dop
            all_eV[:,doy-1] = bin_edop
                    
    # Formatter for the local time axis
    def LT(x, pos):
        return ('%02d' % (np.mod(x,24)))

    ###########################################
    # Create temperature plot handle
    Temperature_Fig = plt.figure(figsize=(5,1.25))
    Temperature_Ax = Temperature_Fig.add_subplot(111)

    # Mask the data to be plotted
    masked_T = ma.masked_where((np.isnan(all_T)) | (all_eT > 100) | (all_T > Tmax), all_T)

    # Plot the data
    mappable_T = Temperature_Ax.pcolormesh(times,center_time,masked_T,vmin=Tmin,vmax=Tmax)

    # Fix y-axis labels
    Temperature_Ax.set_ylim(center_time[0],center_time[-1])
    Temperature_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Temperature_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Temperature_Ax.set_xlim(times[0],times[-1])
    Temperature_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Temperature_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Temperature_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Temperature_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Temperature_Fig.colorbar(mappable_T, ax=Temperature_Ax, ticks=[Tmin,(Tmin+Tmax)/2.0,Tmax])
    cbar.set_label('[K]')

    # Default title
    Temperature_Ax.set_title('Neutral Temperatures')

    ###########################################
    # Create zonal wind plot handle
    Zonal_Fig = plt.figure(figsize=(5,1.25))
    Zonal_Ax = Zonal_Fig.add_subplot(111)

    # Mask the data to be plotted
    masked_U = ma.masked_where((np.isnan(all_U)) | (all_eU > 25) | (abs(all_U) > Dmax*2), all_U)

    # Plot the data
    mappable_U = Zonal_Ax.pcolormesh(times,center_time,masked_U,vmin=Dmin,vmax=Dmax)

    # Fix y-axis labels
    Zonal_Ax.set_ylim(center_time[0],center_time[-1])
    Zonal_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Zonal_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Zonal_Ax.set_xlim(times[0],times[-1])
    Zonal_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Zonal_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Zonal_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Zonal_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Zonal_Fig.colorbar(mappable_U, ax=Zonal_Ax, ticks=[Dmin,(Dmin+Dmax)/2.0,Dmax])
    cbar.set_label('[m/s]')

    # Default title
    Zonal_Ax.set_title('Zonal Neutral Wind')

    ###########################################
    # Create meridional wind plot handle
    Meridional_Fig = plt.figure(figsize=(5,1.25))
    Meridional_Ax = Meridional_Fig.add_subplot(111)

    # Mask the data to be plotted
    masked_V = ma.masked_where((np.isnan(all_V)) | (all_eV > 25) | (abs(all_V) > Dmax*2), all_V)

    # Plot the data
    mappable_V = Meridional_Ax.pcolormesh(times,center_time,masked_V,vmin=Dmin,vmax=Dmax)

    # Fix y-axis labels
    Meridional_Ax.set_ylim(center_time[0],center_time[-1])
    Meridional_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Meridional_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Meridional_Ax.set_xlim(times[0],times[-1])
    Meridional_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Meridional_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Meridional_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Meridional_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Meridional_Fig.colorbar(mappable_V, ax=Meridional_Ax, ticks=[Dmin,(Dmin+Dmax)/2.0,Dmax])
    cbar.set_label('[m/s]')

    # Default title
    Meridional_Ax.set_title('Meridional Neutral Wind')

    ###########################################
    # Create vertical wind plot handle
    Vertical_Fig = plt.figure(figsize=(5,1.25))
    Vertical_Ax = Vertical_Fig.add_subplot(111)
    Dmin = -50.
    Dmax = 50.

    # Mask the data to be plotted
    masked_W = ma.masked_where((np.isnan(all_W)) | (all_eW > 25) | (all_W > Dmax*2), all_W)

    # Plot the data
    mappable_W = Vertical_Ax.pcolormesh(times,center_time,masked_W,vmin=Dmin,vmax=Dmax)

    # Fix y-axis labels
    Vertical_Ax.set_ylim(center_time[0],center_time[-1])
    Vertical_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Vertical_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Vertical_Ax.set_xlim(times[0],times[-1])
    Vertical_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Vertical_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Vertical_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Vertical_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Vertical_Fig.colorbar(mappable_W, ax=Vertical_Ax, ticks=[Dmin,(Dmin+Dmax)/2.0,Dmax])
    cbar.set_label('[m/s]')

    # Default title
    Vertical_Ax.set_title('Vertical Neutral Wind')

    ###########################################
    # Create intensity plot handle
    Intensity_Fig = plt.figure(figsize=(5,1.25))
    Intensity_Ax = Intensity_Fig.add_subplot(111)

    # Mask the data to be plotted
    masked_I = ma.masked_where((np.isnan(all_I)) | (all_eI > 100) | (all_I > Imax*2), all_I)

    # Plot the data
    mappable_I = Intensity_Ax.pcolormesh(times,center_time,masked_I,vmin=Imin,vmax=Imax)

    # Fix y-axis labels
    Intensity_Ax.set_ylim(center_time[0],center_time[-1])
    Intensity_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Intensity_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Intensity_Ax.set_xlim(times[0],times[-1])
    Intensity_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Intensity_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Intensity_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Intensity_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Intensity_Fig.colorbar(mappable_I, ax=Intensity_Ax, ticks=[Imin,(Imin+Imax)/2.0,Imax])
    cbar.set_label('[arb]')

    # Default title
    Intensity_Ax.set_title('630.0-nm Intensities')

    return (Temperature_Fig, Temperature_Ax), (Zonal_Fig, Zonal_Ax), (Meridional_Fig, Meridional_Ax), (Vertical_Fig, Vertical_Ax), (Intensity_Fig, Intensity_Ax)


def DataSummary(files, times, bin_time = np.arange(17,32,0.5),
                Tmin=500, Tmax=1500, Dmin=-200, Dmax=200, Imin=0, Imax=200, cloudy_temperature=-15.0,
                reference='Zenith'):
#
# Function to display a two dimensional summary of data.  The 
# y-axis is hours while the x-axis is the day.  The function
# stacks data in columns so daily/seasonal trends can be observed
#
# INPUTS:
#   files - a list of files to be included in the summary plot
#   times - a list of datetimes for plotting (x-axis)
# OPTIONAL INPUTS:
#   bin_time - array of times to bin data on to (default np.arange(17,32,0.5))
#   Tmin, Tmax - the min/max values for temperatures to be plotted (default 500,1500)
#   Dmin, Dmax - the min/max values for the Doppler shifts to be plotted (default -200, 200)
#   reference - the type of Doppler reference to use
# OUTPUTS:
#   (Temperature_Fig, Temperature_Ax) - references to figure and axis of temperature plot
#   (Zonal_Fig, Zonal_Ax) - references to figure and axis of zonal plot
#   (Meridional_Fig, Meridional_Ax) - references to figure and axis of meridional plot
#   (Intensity_Fig, Intesity_Ax) - references to the figure and axis for the airglow intensity plot
#
# HISTORY:
#   Written by Jonathan J. Makela on 2 Dec 2012
#   Added Intensity figures on 28 May 2013 (jjm)
#   Corrected to work with new file format on 22 Aug 2013

    # Sort the files given
    files.sort()

    # Create arrays to hold data
    N = np.size(times)

    all_T = np.zeros((len(bin_time)-1,N))*np.nan
    all_eT = np.zeros((len(bin_time)-1,N))*np.nan

    all_U = np.zeros((len(bin_time)-1,N))*np.nan
    all_eU = np.zeros((len(bin_time)-1,N))*np.nan

    all_V = np.zeros((len(bin_time)-1,N))*np.nan
    all_eV = np.zeros((len(bin_time)-1,N))*np.nan

    all_W = np.zeros((len(bin_time)-1,N))*np.nan
    all_eW = np.zeros((len(bin_time)-1,N))*np.nan

    all_I = np.zeros((len(bin_time)-1,N))*np.nan
    all_eI = np.zeros((len(bin_time)-1,N))*np.nan

    # Center of the bins (used in plotting)
    center_time = (bin_time[0:-1]+bin_time[1:])/2.

    # Work with the data from the files provided
    for f in files:
        # Load the file
        npzfile = np.load(f,allow_pickle=True)

        # Save data to FPI_Results and site dictionaries
        FPI_Results = npzfile['FPI_Results']
        FPI_Results = FPI_Results.reshape(-1)[0]
        site = npzfile['site']
        site = site.reshape(-1)[0]
        npzfile.close()

        # Calculate the day of year
        # TODO: GENERALIZE SO IT ISNT JUST PLOTTING OVER A YEAR
    #       doy = FPI_Results['sky_times'][0].timetuple().tm_yday

        doy = (FPI_Results['sky_times'][0]-dates.num2date(times[0])).days
        if (doy < 1) or (doy > N):
            continue
        
        # Find the zero offset of the Doppler shift
        (ref_Dop, e_ref_Dop) = FPI.DopplerReference(FPI_Results,reference=reference)

        # Calculate the vertical wind and interpolate it
        ind = FPI.all_indices('Zenith',FPI_Results['direction'])
        w = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]) # -1 because LOS is towards instrument
        sigma_w = FPI_Results['sigma_LOSwind'][ind]
        dt = []
        for x in FPI_Results['sky_times'][ind]:
            diff = (x - FPI_Results['sky_times'][0])
            dt.append(diff.seconds+diff.days*86400.)
        dt = np.array(dt)

        # Remove outliers
        ind = abs(w) < 200.
        
        if sum(ind) <= 1:
            # No good data, just use all ind
            ind = abs(w) > 0.

        if sum(ind) <= 1:
            continue

        # Interpolate
        w2 = interpolate.interp1d(dt[ind],w[ind],bounds_error=False,fill_value=0.0)
        sigma_w2 = interpolate.interp1d(dt[ind],sigma_w[ind],bounds_error=False,fill_value=0.0)
        dt = []
        
        for x in FPI_Results['sky_times']:
            diff = (x - FPI_Results['sky_times'][0])
            dt.append(diff.seconds+diff.days*86400.)
        w = w2(dt)
        sigma_w = sigma_w2(dt)

        # Fill in the arrays with binned data
        for x in np.unique(FPI_Results['direction']):
            # Use zenith for temperatures
            # TODO: GENERALIZE THIS TO ALLOW USER TO DETERMINE WHICH DIRECTION
            if x == 'Zenith':
                ind = FPI.all_indices(x, FPI_Results['direction'])

                # Grab the data
                st = FPI_Results['sky_times'][ind]
                T = FPI_Results['T'][ind]
                eT = FPI_Results['sigma_T'][ind]
                I = FPI_Results['skyI'][ind]
                eI = FPI_Results['sigma_skyI'][ind]
                dop = (FPI_Results['LOSwind'][ind]-ref_Dop[ind])
                # TODO: ERROR BASED ON BOTH LOS AND ZENITH REFERENCE
                edop = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)

                # Check if clouds are provided
                if 'Clouds' in FPI_Results.keys():
                    if FPI_Results['Clouds'] is not None:
                        clouds = FPI_Results['Clouds']['mean'][ind]
                        idx = clouds < cloudy_temperature
                        st = st[idx]
                        T = T[idx]
                        eT = eT[idx]
                        I = I[idx]
                        eI = eI[idx]
                        dop = dop[idx]
                        edop = edop[idx]

                # Find bad data points
                idx = (eT < 100) & (eT > 0)
                st = st[idx]
                T = T[idx]
                eT = eT[idx]
                I = I[idx]
                eI = eI[idx]
                dop = dop[idx]
                edop = edop[idx]

                if len(st) > 0:
                    # Bin the data
                    (bin_T, bin_eT) = FPI.bin_and_mean(st,T,eT,bin_time)

                    # Save to array
                    all_T[:,doy-1] = bin_T
                    all_eT[:,doy-1] = bin_eT

                    # Bin the intensity data
                    (bin_I, bin_eI) = FPI.bin_and_mean(st,I,eI,bin_time)

                    # Save to array
                    all_I[:,doy-1] = bin_I
                    all_eI[:,doy-1] = bin_eI

                    (bin_dop, bin_edop) = FPI.bin_and_mean(st,dop,edop,bin_time)

                    # Save to array
                    all_W[:,doy-1] = bin_dop
                    all_eW[:,doy-1] = bin_edop

            # Use east for zonal wind
            # TODO: GENERALIZE THIS TO ALLOW USER TO DETERMINE WHICH DIRECTION
            if x == 'East':
                ind = FPI.all_indices(x,FPI_Results['direction'])

                # Grab the data
                st = FPI_Results['sky_times'][ind]
                dop = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
                # TODO: ERROR BASED ON BOTH LOS AND ZENITH REFERENCE
                edop = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)

                # Check if clouds are provided
                if 'Clouds' in FPI_Results.keys():
                    if FPI_Results['Clouds'] is not None:
                        clouds = FPI_Results['Clouds']['mean'][ind]
                        idx = clouds < cloudy_temperature
                        st = st[idx]
                        dop = dop[idx]
                        edop = edop[idx]

                # Find bad data points
                idx = (edop < 50) & (edop > 0)
                st = st[idx]
                dop = dop[idx]
                edop = edop[idx]

                # Bin the data
                if len(st) > 0:
                    (bin_dop, bin_edop) = FPI.bin_and_mean(st,dop,edop,bin_time)

                    # Save to array
                    all_U[:,doy-1] = bin_dop
                    all_eU[:,doy-1] = bin_edop

            # Use north for meridional wind
            # TODO: GENERALIZE THIS TO ALLOW USER TO DETERMINE WHICH DIRECTION
            if x == 'North':
                ind = FPI.all_indices(x,FPI_Results['direction'])

                # Grab the data
                st = FPI_Results['sky_times'][ind]
                dop = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
                # TODO: ERROR BASED ON BOTH LOS AND ZENITH REFERENCE
                edop = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)

                # Check if clouds are provided
                if 'Clouds' in FPI_Results.keys():
                    if FPI_Results['Clouds'] is not None:
                        clouds = FPI_Results['Clouds']['mean'][ind]
                        idx = clouds < cloudy_temperature
                        st = st[idx]
                        dop = dop[idx]
                        edop = edop[idx]

                # Find bad data points
                idx = (edop < 50) & (edop > 0)
                st = st[idx]
                dop = dop[idx]
                edop = edop[idx]

                # Bin the data
                if len(st) > 0:
                    (bin_dop, bin_edop) = FPI.bin_and_mean(st,dop,edop,bin_time)

                    # Save to array
                    all_V[:,doy-1] = bin_dop
                    all_eV[:,doy-1] = bin_edop                # Bin the data

    # Formatter for the local time axis
    def LT(x, pos):
        return ('%02d' % (np.mod(x,24)))

    ###########################################
    # Create temperature plot handle
    Temperature_Fig = plt.figure(figsize=(5,1.25))
    Temperature_Ax = Temperature_Fig.add_subplot(111)

    # Mask the data to be plotted
    masked_T = ma.masked_where((np.isnan(all_T)) | (all_eT > 100) | (all_T > Tmax*2), all_T)

    # Plot the data
    mappable_T = Temperature_Ax.pcolormesh(times,center_time,masked_T,vmin=Tmin,vmax=Tmax)

    # Fix y-axis labels
    Temperature_Ax.set_ylim(center_time[0],center_time[-1])
    Temperature_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Temperature_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Temperature_Ax.set_xlim(times[0],times[-1])
    Temperature_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Temperature_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Temperature_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Temperature_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Temperature_Fig.colorbar(mappable_T, ax=Temperature_Ax, ticks=[Tmin,(Tmin+Tmax)/2.0,Tmax])
    cbar.set_label('[K]')

    # Default title
    Temperature_Ax.set_title('Neutral Temperatures')

    ###########################################
    # Create zonal wind plot handle
    Zonal_Fig = plt.figure(figsize=(5,1.25))
    Zonal_Ax = Zonal_Fig.add_subplot(111)

    # Mask the data to be plotted
    masked_U = ma.masked_where((np.isnan(all_U)) | (all_eU > 25) | (all_U > Dmax*2), all_U)

    # Plot the data
    mappable_U = Zonal_Ax.pcolormesh(times,center_time,masked_U,vmin=Dmin,vmax=Dmax)

    # Fix y-axis labels
    Zonal_Ax.set_ylim(center_time[0],center_time[-1])
    Zonal_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Zonal_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Zonal_Ax.set_xlim(times[0],times[-1])
    Zonal_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Zonal_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Zonal_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Zonal_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Zonal_Fig.colorbar(mappable_U, ax=Zonal_Ax, ticks=[Dmin,(Dmin+Dmax)/2.0,Dmax])
    cbar.set_label('[m/s]')

    # Default title
    Zonal_Ax.set_title('Zonal Neutral Wind')

    ###########################################
    # Create meridional wind plot handle
    Meridional_Fig = plt.figure(figsize=(5,1.25))
    Meridional_Ax = Meridional_Fig.add_subplot(111)

    # Mask the data to be plotted
    masked_V = ma.masked_where((np.isnan(all_V)) | (all_eV > 25) | (all_V > Dmax*2), all_V)

    # Plot the data
    mappable_V = Meridional_Ax.pcolormesh(times,center_time,masked_V,vmin=Dmin,vmax=Dmax)

    # Fix y-axis labels
    Meridional_Ax.set_ylim(center_time[0],center_time[-1])
    Meridional_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Meridional_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Meridional_Ax.set_xlim(times[0],times[-1])
    Meridional_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Meridional_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Meridional_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Meridional_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Meridional_Fig.colorbar(mappable_V, ax=Meridional_Ax, ticks=[Dmin,(Dmin+Dmax)/2.0,Dmax])
    cbar.set_label('[m/s]')

    # Default title
    Meridional_Ax.set_title('Meridional Neutral Wind')

    ###########################################
    # Create vertical wind plot handle
    Vertical_Fig = plt.figure(figsize=(5,1.25))
    Vertical_Ax = Vertical_Fig.add_subplot(111)
    Dmin = -50.
    Dmax = 50.

    # Mask the data to be plotted
    masked_W = ma.masked_where((np.isnan(all_W)) | (all_eW > 25) | (all_W > Dmax*2), all_W)

    # Plot the data
    mappable_W = Vertical_Ax.pcolormesh(times,center_time,masked_W,vmin=Dmin,vmax=Dmax)

    # Fix y-axis labels
    Vertical_Ax.set_ylim(center_time[0],center_time[-1])
    Vertical_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Vertical_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Vertical_Ax.set_xlim(times[0],times[-1])
    Vertical_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Vertical_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Vertical_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Vertical_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Vertical_Fig.colorbar(mappable_W, ax=Vertical_Ax, ticks=[Dmin,(Dmin+Dmax)/2.0,Dmax])
    cbar.set_label('[m/s]')

    # Default title
    Vertical_Ax.set_title('Vertical Neutral Wind')

    ###########################################
    # Create intensity plot handle
    Intensity_Fig = plt.figure(figsize=(5,1.25))
    Intensity_Ax = Intensity_Fig.add_subplot(111)

    # Mask the data to be plotted
    masked_I = ma.masked_where((np.isnan(all_I)) | (all_eI > 100) | (all_I > Imax*2), all_I)

    # Plot the data
    mappable_I = Intensity_Ax.pcolormesh(times,center_time,masked_I,vmin=Imin,vmax=Imax)

    # Fix y-axis labels
    Intensity_Ax.set_ylim(center_time[0],center_time[-1])
    Intensity_Ax.yaxis.set_major_formatter(FuncFormatter(LT))
    Intensity_Ax.set_ylabel('LT [hrs]')

    # Fix the x-axis labels
    Intensity_Ax.set_xlim(times[0],times[-1])
    Intensity_Ax.xaxis.set_major_locator(dates.MonthLocator(interval = 6))
    Intensity_Ax.xaxis.set_minor_locator(dates.MonthLocator())
    Intensity_Ax.xaxis.set_major_formatter(dates.DateFormatter('%b\n%Y'))
    Intensity_Ax.tick_params(axis='both', which='major', labelsize=8)

    # Draw colorbar
    cbar = Intensity_Fig.colorbar(mappable_I, ax=Intensity_Ax, ticks=[Imin,(Imin+Imax)/2.0,Imax])
    cbar.set_label('[arb]')

    # Default title
    Intensity_Ax.set_title('630.0-nm Intensities')

    return (Temperature_Fig, Temperature_Ax), (Zonal_Fig, Zonal_Ax), (Meridional_Fig, Meridional_Ax), (Vertical_Fig, Vertical_Ax), (Intensity_Fig, Intensity_Ax)

def PlotDay(f, full_clear=-30, full_cloud=-20,
            Tmin=300, Tmax=1500, Dmin=-200, Dmax=200,
            reference='laser',directions=None,
	    Temperature_Fig = None, Temperature_Graph = None, Doppler_Fig = None, Doppler_Graph = None,
        cull=False):
#
# Function to plot a single night's data for a single station
#
# INPUTS:
#
# OPTION OUTPUTS:
#
# OUTPUTS:
#
# HISTORY:
#   Written by Jonathan J. Makela on 3 Dec 2012
#   Removed the use of "ze_corr" on 17 Jul 2013 (bjh)
#   Defaulted reference to 'laser' on 22 Jul 2013 (bjh)
#   Added color/format information here, because plot-specific
#       information doesn't seem like it should belong to the instrument.
#   Removed zenith_times argument (bjh)

    # Read in the file
    npzfile = np.load(f,allow_pickle=True)
    FPI_Results = npzfile['FPI_Results']
    FPI_Results = FPI_Results.reshape(-1)[0]
    site = npzfile['site']
    site = site.reshape(-1)[0]
    del npzfile.f # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
    npzfile.close()
    
    # Assign markers and colors to Direction keys in a reasonable way.
    # All "similiar" direction keys should get the same marker, 
    # but different colors. (CV_ANN_UAO_1 and IN_ANN_UAO are "similiar")
    all_markers = ['s', 'd', '^', 'x', 'p', '*', '+', 'v', '>', '<']
    all_colors  = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    all_markers.reverse() # pop() starts from the end
    all_colors.reverse() 
    fmt = { 'East': {'Color': 'y', 'Marker': 'o'},
            'West': {'Color': 'm', 'Marker': 'o'},
            'North': {'Color': 'r', 'Marker': 'o'},
            'South': {'Color': 'g', 'Marker': 'o'},
            'Zenith': {'Color': 'k', 'Marker': 'o'},
            'Laser' : None
        }
    cvs = [d for d in site['Directions'].keys() if d not in fmt] # e.g., ['CV_ANN_UAO_1', 'IN_EKU_UAO', etc.]
    keys = list(set((['_'.join(c.split('_')[1:3]) for c in cvs]))) # e.g., ['ANN_UAO', 'EKU_UAO'] (no repeats)
    # make a unique marker for each key, and vary the color
    keymap = {}
    for key in keys:
        keymap[key] = (all_markers.pop(), list(all_colors)) # (unique marker, copy of color list)
    # assign
    for direc in cvs:
        # create the key
        key = '_'.join(direc.split('_')[1:3])
        mark, colors = keymap[key]
        # assign the color and marker
        fmt[direc] = {'Color': colors.pop(), 'Marker': mark} # vary the color, keep the marker

            
    
    
    if Temperature_Fig is None:
	# Add Temperature figure
   	Temperature_Fig = plt.figure()
    	Temperature_Graph = Temperature_Fig.add_subplot(111)

    if Doppler_Fig is None:
        # Add Doppler figure
        Doppler_Fig = plt.figure()
        Doppler_Graph = Doppler_Fig.add_subplot(111)
        
    (ref_Dop, e_ref_Dop) = FPI.DopplerReference(FPI_Results,reference=reference)

    # Calculate the vertical wind and interpolate it
    ind = FPI.all_indices('Zenith',FPI_Results['direction'])
    w = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]) # LOS is away from instrument
    sigma_w = FPI_Results['sigma_LOSwind'][ind]
    dt = []
    for x in FPI_Results['sky_times'][ind]:
        diff = (x - FPI_Results['sky_times'][0])
        dt.append(diff.seconds+diff.days*86400.)
    dt = np.array(dt)

    # Remove outliers
    ind = abs(w) < 200.

    if sum(ind) <= 1:
	# No good data, just use all ind
	    ind = np.array([True for i in range(len(w))]) # There has to be a clearer way to do this...

    if len(ind) == 0:
       raise Exception('%s: No Zenith look directions' % f)

    # Interpolate
    w2 = interpolate.interp1d(dt[ind],w[ind],bounds_error=False,fill_value=0.0)
    sigma_w2 = interpolate.interp1d(dt[ind],sigma_w[ind],bounds_error=False,fill_value=0.0)
    dt = []
    
    for x in FPI_Results['sky_times']:
        diff = (x - FPI_Results['sky_times'][0])
        dt.append(diff.seconds+diff.days*86400.)
    w = w2(dt)
    sigma_w = sigma_w2(dt)

    # Check if specific directions are requested
    if directions is None:
        # Create the plot with each direction plotted in a different colors
        l = [x for x in np.unique(FPI_Results['direction']) if x not in ['Unknown', 'Laser'] ]
    else:
        l = directions
        
    for x in l:
        ind = FPI.all_indices(x,FPI_Results['direction'])
        
        # Inital plot of bogus value to get the legend right
        Temperature_Graph.plot(-999,-999,color=fmt[x]['Color'], marker=fmt[x]['Marker'], label=x)

        if x == 'Zenith':
            Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind])
            Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2)
        else:
            Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
            Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)
        if x == 'South' or x == 'West':
            Doppler_Wind = -Doppler_Wind
            
        Doppler_Graph.plot(FPI_Results['sky_times'][ind],Doppler_Wind,color=fmt[x]['Color'],alpha=0.3,marker=None,label='_nolegend_')
        Doppler_Graph.plot(-999,-999,color=fmt[x]['Color'],marker=fmt[x]['Marker'],label=x)
        
	if ('Clouds' in FPI_Results.keys()) is False:
		FPI_Results['Clouds'] = None


        # Loop through each point (needed to be done this way to allow setting
        # the alpha value for each individual point)
        for (t,y,ey,z,ez,wq,tq) in zip(FPI_Results['sky_times'][ind], FPI_Results['T'][ind], FPI_Results['sigma_T'][ind],Doppler_Wind,Doppler_Error,FPI_Results['wind_quality_flag'][ind],FPI_Results['temp_quality_flag'][ind]):
        
            # Calculate the alpha value to be used for this point, based on quality flag
            w_alpha_val = 1.
            t_alpha_val = 1.
            if wq == 1:
                w_alpha_val = 0.5
            if tq == 1:
                t_alpha_val = 0.5
            if wq == 2:
                w_alpha_val = 0.2
            if tq == 2:
                t_alpha_val = 0.2
            # Plot the points and error bars on the graphs
            Temperature_Graph.errorbar(t,y,yerr=ey,alpha=t_alpha_val,color=fmt[x]['Color'],fmt=fmt[x]['Marker'],label='_nolegend_')
            Doppler_Graph.errorbar(t,z,yerr=ez,alpha=w_alpha_val,color=fmt[x]['Color'],fmt=fmt[x]['Marker'],label='_nolegend_')

    # Plotting font setup
    fontP = FontProperties() 
    fontP.set_size('10')
    fontsize = 10 # because legend() and xlabel() do it different ways
    
    # Format the Temperature plot               
    Temperature_Graph.set_ylim(Tmin,Tmax)
    Temperature_Graph.set_xlim(FPI_Results['sky_times'][0],FPI_Results['sky_times'][-1])
    Temperature_Graph.legend(ncol=4, prop=fontP)
    Temperature_Graph.xaxis.set_major_formatter(dates.DateFormatter('%H'))
    Temperature_Graph.set_ylabel('Neutral Temperature [K]', fontsize = fontsize)
    Temperature_Graph.set_xlabel('Universal Time', fontsize = fontsize)
    Temperature_Graph.grid(True)
    #bjh: How do I control font size? This works on non-date axes,
    # but throws a weird error here:
    #for tick in Doppler_Graph.xaxis.get_major_ticks():
    #    tick.label.set_fontsize(fontsize-2)

    # Format the Doppler plot
    Doppler_Graph.set_ylim(Dmin,Dmax)
    Doppler_Graph.set_xlim(FPI_Results['sky_times'][0],FPI_Results['sky_times'][-1])
    Doppler_Graph.legend(ncol=4, prop=fontP)
    Doppler_Graph.xaxis.set_major_formatter(dates.DateFormatter('%H'))
    Doppler_Graph.set_ylabel('Neutral Winds [m/s]', fontsize = fontsize)
    Doppler_Graph.set_xlabel('Universal Time', fontsize = fontsize)
    Doppler_Graph.plot([FPI_Results['sky_times'][0],FPI_Results['sky_times'][-1]],[0,0],'k--')
    Doppler_Graph.grid(True)
    #for tick in Doppler_Graph.xaxis.get_major_ticks():
    #    tick.label.set_fontsize(fontsize-2)

    # Mark the plots with the Doppler reference type
    Doppler_Fig.text(0.8,0.0,'Doppler ref: %s' % reference, horizontalalignment='right', verticalalignment='bottom', fontsize = fontsize-2)

    # If no cloud data, mark this on the plots
    if FPI_Results['Clouds'] is None:
        Temperature_Fig.text(0.0,0.0,'No cloud data', verticalalignment='bottom', fontsize = fontsize-2)
        Doppler_Fig.text(0.0,0.0,'No cloud data', verticalalignment='bottom', fontsize = fontsize-2)

    # Create title based on if start/stop time are on the same day or not
    if FPI_Results['sky_times'][0].day != FPI_Results['sky_times'][-1].day:
        Temperature_Graph.set_title(site['Abbreviation'] + ':' + (FPI_Results['sky_times'][0]).strftime('%d %b, %Y %H:%M LT') + ' - ' + (FPI_Results['sky_times'][-1]).strftime('%d %b, %Y %H:%M LT'), fontsize = fontsize)
        Doppler_Graph.set_title(site['Abbreviation'] + ':' + (FPI_Results['sky_times'][0]).strftime('%d %b, %Y %H:%M LT') + ' - ' + (FPI_Results['sky_times'][-1]).strftime('%d %b, %Y %H:%M LT'), fontsize = fontsize)
    else:
        Temperature_Graph.set_title(site['Abbreviation'] + ':' + (FPI_Results['sky_times'][0]).strftime('%d %b, %Y %H:%M') + ' - ' + (FPI_Results['sky_times'][-1]).strftime('%H:%M LT'), fontsize = fontsize)
        Doppler_Graph.set_title(site['Abbreviation'] + ':' + (FPI_Results['sky_times'][0]).strftime('%d %b, %Y %H:%M') + ' - ' + (FPI_Results['sky_times'][-1]).strftime('%H:%M LT'), fontsize = fontsize)

    return(Temperature_Fig, Temperature_Graph), (Doppler_Fig, Doppler_Graph)

def PlotDiagnosticDay(f, cloud_thresh = [-22.,-10.],\
			sky_quality_thresh =[-np.inf,-np.inf],\
			LASERPARAMINTERP = 'linear',\
			):
    '''
    Function to plot a single night's diagnotic file
    
    INPUTS:
       f                  - absolute path of the NPZ file.
       cloud_thresh       - [float,float], K. The two cloud sensor readings that indicate 
                             partially- and fully-cloudy. This affects the quality flag.
                          - None. Extracts the cloud thresholds from log file.
                             It must be located along with f within the same directory.
       sky_quality_thresh - [float,float]. The two intensity thresholds that indicate 
                             q1, q2 low sky intensity thresholds.
                          - None. Extract the thresholds from the instrument dictionary
                             saved within NPZ file.
       LASERPARAMINTERP   - 'linear' or 'spline'  interpolation used for instrument parameters
    
    OUTPUTS:
    
    HISTORY:
        5 Feb 2021, Compiled and tested by L. Navarro.
        
    '''
    import re
    # Read in the file
    npzfile = np.load(f,allow_pickle=True)
    FPI_Results = npzfile['FPI_Results']
    FPI_Results = FPI_Results.reshape(-1)[0]
    instrument = npzfile['instrument']
    instrument = instrument.reshape(-1)[0]
    site = npzfile['site']
    site = site.reshape(-1)[0]
    del npzfile.f # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
    npzfile.close()
    
    instrsitedate=os.path.basename(f)
    
#     instr_name ,_,datestr=instrsitedate.split("_")[:3]
#     nominal_dt=datetime.datetime.strptime(datestr[:8],'%Y%m%d')
#     site_name = fpiinfo.get_site_of(instr_name, nominal_dt)
#     site = fpiinfo.get_site_info(site_name, nominal_dt)
    
    #Reading all looking directions
    valid_az = np.array([site['Directions'][direc]['az'] for direc in site['Directions']])
    valid_ze = np.array([site['Directions'][direc]['ze'] for direc in site['Directions']])
    
    #Sky quality thresholds
    if sky_quality_thresh is None:
        sky_quality_thresh=instrument['skyI_quality_thresh']
    
    #Cloud cover thresholds
    if cloud_thresh is None:
        logname = os.path.dirname(f) +'/' + instrsitedate + '.log'
        h=open(logname,'r')
        lines=h.readlines()
        h.close()
        q1=-np.inf
        q2=-np.inf
        for line in lines:
            q1s=re.findall(r">(.\d*): W1T1",line)
            if len(q1s)>0:
                q1=np.float(q1s[0])
            q2s=re.findall(r">(.\d*): W2T1",line)
            if len(q2s)>0:
                q2=np.float(q2s[0])
            if np.isfinite(q1) and np.isfinite(q2):
                break
        cloud_thresh=[q1,q2]
    
    #Read output results
    sky_times=FPI_Results['sky_times']
    sky_redchi=FPI_Results['sky_chisqr']
    skyI=FPI_Results['skyI']
    all_az=FPI_Results['az']
    all_ze=FPI_Results['ze']
    direction=FPI_Results['direction']
    laser_times=FPI_Results['laser_times']
    laser_redchi=FPI_Results['laser_chisqr']
    laser_value=FPI_Results['laser_value']
    laser_stderr=FPI_Results['laser_stderr']
    reference=FPI_Results['reference']
    center=FPI_Results['center_pixel']
    
    #use laser or not
    uselaser='laser' in reference
    
    
    
    
    
    
    fig = plt.figure(dpi=300, figsize=(10,7.5)) # Figure for diagnostics to be drawn to
    
    fontP = FontProperties()
    fontP.set_size('small')
    
    ################ Plot center location ######################
    # only plot this one if uselaser and centerfinding succeeded for > 1 point
    if uselaser and center is not None and len(center)>0 and np.shape(center)!=(2,):
        lt0 = laser_times[0]
        tdiffvec = [laser_time - lt0 for laser_time in laser_times]
        dt = np.array([tdiff.seconds + tdiff.days*86400. for tdiff in tdiffvec])
        t = np.linspace(0, dt.max(), 500)
        datet = [lt0 + datetime.timedelta(seconds=ts) for ts in t]
        
        dt_laser = [tdiff.seconds for tdiff in tdiffvec]
        npoly = np.floor(len(dt_laser)/10) # use an adaptive degree polynomial to fit
        # Limit the maximum degree, because of numerical sensitivity for large orders.
        if npoly > 10:
            npoly = 10 
        pf_cx = np.polyfit(dt_laser,center[:,0],npoly)
        cx = np.poly1d(pf_cx)
        pf_cy = np.polyfit(dt_laser,center[:,1],npoly)
        cy = np.poly1d(pf_cy)
        
        
        
        xdev = center[:,0] - center[0,0]
        ydev = center[:,1] - center[0,1]
        xdevi = cx(t) - center[0,0]
        ydevi = cy(t) - center[0,1]

        ax = fig.add_subplot(421)
        ax.plot(datet, xdevi, 'b', label='x')
        ax.plot(laser_times, xdev, 'b.')
        ax.plot(datet, ydevi, 'r', label='y')
        ax.plot(laser_times, ydev, 'r.')
        ax.legend(loc='best', prop=fontP)
        ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        m = 0.5 # maximum deviation to show on plot
        if all(xdev < m) and all(xdev > -m) and all(ydev < m) and all(ydev > -m): # Set the x and y lims at +/- m 
            ax.set_ylim([-m, m])

        ax.set_ylabel('Center deviation\n[pixels]')
        #ax.set_xlabel('Universal Time, [hours]')
        ax.set_title(site['Abbreviation'] + ':' + \
            laser_times[0].strftime(' %d %b, %Y %H:%M LT') + ' - ' + laser_times[-1].strftime('%H:%M LT') )

        ax.grid(True)

    ####################### Laser Fit Chi^2 #######################
    
    if uselaser:
        ax = fig.add_subplot(422)
        ax.plot(laser_times, laser_redchi,'k.-')
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
        ax.set_ylabel('Laser Fit\nReduced Chi^2')
        ax.set_title(site['Abbreviation'] + ':' + \
                laser_times[0].strftime(' %d %b, %Y %H:%M LT') + ' - ' + laser_times[-1].strftime('%H:%M LT') )
        ax.grid(True)

    ####################### Spline fit for I #######################
    # Show the laser intensity halfway out in the spectrum
    if uselaser:
        r = 0.5 # r/rmax, where r is the radius of the radial bin at which to measure the intensity
        I = laser_value['I']
        a1 = laser_value['a1']
        a2 = laser_value['a2']
        Ir = I * ( 1 + a1*r + a2*r**2 )
        # assume errorbar at r=rmax/2 is approximately equal to that at r=0
        Ie = laser_stderr['I']
        ax = fig.add_subplot(423)
        ax.errorbar(laser_times,Ir,yerr=Ie,fmt='k.-')
        ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        ax.set_ylabel('Laser Intensity\n[counts]')
        ax.grid(True)

    ####################### Spline fit for t #######################
    # Show the spline fit for a certain parameter
    if uselaser:
        p=laser_value['t']
        w=laser_stderr['t']
        s=len(w)
        if LASERPARAMINTERP == 'spline':
            laser_spfit = interpolate.UnivariateSpline(np.array(dt),p, w=w, s=s)
        elif LASERPARAMINTERP == 'linear':
            laser_spfit = interpolate.interp1d(np.array(dt), p)
        
        dt = []
        for x in laser_times:
            diff = (x - lt0)
            dt.append(diff.seconds+diff.days*86400.)
        
        dt_vec = np.linspace(0, max(dt), 500)
        p_vec = laser_spfit(dt_vec)
        datet_vec = [laser_times[0] + datetime.timedelta(seconds=ts) for ts in dt_vec]

        ax = fig.add_subplot(425)
        ax.plot(datet_vec,p_vec,'k')
        ax.errorbar(laser_times,p,yerr=w,fmt='k.')
        ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        ax.set_ylabel('Etalon Gap, [m]')
        ax.grid(True)

    
    ##################### Look Direction Validation ####################### OUT

    # Plot the look directions
    all_ze = np.array(all_ze)
    all_az = np.array(all_az)
    valid_ze = np.array(valid_ze)
    valid_az = np.array(valid_az)

    # Flip az for negative ze angles
    idx = all_ze < 0
    all_ze[idx] = -all_ze[idx]
    all_az[idx] = all_az[idx] + 180
    idx = valid_ze < 0
    valid_ze[idx] = -valid_ze[idx]
    valid_az[idx] = valid_az[idx] + 180

    az_rad = all_az * np.pi/180.0
    valid_az_rad = valid_az * np.pi/180.0

    ax = fig.add_subplot(428, projection='polar')
    ax.plot(valid_az_rad, valid_ze, 'kx', label = 'valid')
    valid_idx = [d is not 'Unknown' for d in direction]
    invalid_idx = [not i for i in valid_idx]
    ax.plot(az_rad[np.where(valid_idx)], all_ze[np.where(valid_idx)], 'k.', label = 'actual')
    ax.plot(az_rad[np.where(invalid_idx)], all_ze[np.where(invalid_idx)], 'r.', label = 'unrecognized')
    # Now make it look like a cardinal plot, not a math plot
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi/2)
    ax.set_rmax(90.)
    ax.set_rgrids([30,60])
    ax.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.,numpoints=1,prop=fontP)
    # Indicate any non-displayed points (zenith > 90)
    nnshown = sum(abs(all_ze) > 90.)
    if nnshown > 0:
        ax.set_xlabel('%03i points not shown (|ze| > 90)' % nnshown)
    
    ##################### Sky Fit Chi^2 ####################### 
    
    ax = fig.add_subplot(424)
    ax.plot(sky_times, sky_redchi,'k.-')
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
    ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
    ax.set_ylabel('Sky Fit\nReduced Chi^2')
    ax.grid(True)

    ####################### Plot of skyI ####################### 
    ax = fig.add_subplot(426)
    for direc in list(set(direction)):
        # account for different exposure times
        I = np.array([si for (si,d) in zip(skyI, direction) if d == direc])
        t = np.array([si for (si,d) in zip(sky_times, direction) if d == direc])
        ax.semilogy(t, I, '.-', label=direc)
    tp0 = sky_times[0] - datetime.timedelta(hours=0.5)
    tp1 = sky_times[-1] + datetime.timedelta(hours=0.5)
    ax.semilogy([tp0, tp1],[sky_quality_thresh[0], sky_quality_thresh[0]],'k--',lw=0.5,label='qual thresh (q=1)')
    ax.semilogy([tp0, tp1],[sky_quality_thresh[1], sky_quality_thresh[1]],'k--',lw=0.5,label='qual thresh (q=2)')
    
    ax.set_xlim([tp0, tp1])
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
    ax.set_ylabel('Line Intensity\n[arbitrary]')
    ax.set_xlabel('Universal Time')
    ax.legend(loc='best', prop={'size':6}, numpoints=1, ncol=5, framealpha=0.5)
    ax.grid(True)
    
    ####################### Plot of cloud cover #######################
    ax = fig.add_subplot(427)
    
    if FPI_Results['Clouds'] is not None:
        ct = FPI_Results['sky_times']
        cloud = FPI_Results['Clouds']['mean']
        ax.plot(ct, cloud, 'k.-')
        ax.plot([ct[0],ct[-1]], [cloud_thresh[0],cloud_thresh[0]], 'k--')
        ax.plot([ct[0],ct[-1]], [cloud_thresh[1],cloud_thresh[1]], 'k--')
        ax.set_xlim([ct[0] - datetime.timedelta(hours=0.5), ct[-1] + datetime.timedelta(hours=0.5)])
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        ax.set_ylim([-50,0])
        ax.set_ylabel('Cloud indicator\n[degrees C]')
        ax.grid(True)
        ax.set_xlabel('Universal Time, [hours]')
    
    fig.tight_layout()
    
    return fig




def CompareData(files, full_clear=-30, full_cloud=-20,
            Tmin=500, Tmax=1500, Dmin=-200, Dmax=200,
            reference='Laser',directions=None,displayhours=False,
	    Temperature_Fig = None, Temperature_Graph = None, Doppler_Fig = None, Doppler_Graph = None):
#
# Function to plot a single night's data for a single station
#
# INPUTS:
#
# OPTION OUTPUTS:
#
# OUTPUTS:
#
# HISTORY:
#   Written by Jonathan J. Makela on 3 Dec 2012
#   Removed the use of "ze_corr" on 17 Jul 2013 (bjh)
#   Defaulted reference to 'Laser' on 22 Jul 2013 (bjh)
#   Added color/format information here, because plot-specific
#       information doesn't seem like it should belong to the instrument.

    # Assign markers and colors to Direction keys in a reasonable way.
    # All "similiar" direction keys should get the same marker, 
    # but different colors. (CV_ANN_UAO_1 and IN_ANN_UAO are "similiar")
    all_markers = ['s', 'd', '^', 'x', 'p', '*', '+', 'v', '>', '<']
    all_colors  = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    all_markers.reverse() # pop() starts from the end
    all_colors.reverse() 

    # start and stop time of entire dataset
    t0 = None
    t1 = None
    
    for f in files:
        # Read in the file
        npzfile = np.load(f,allow_pickle=True)
        FPI_Results = npzfile['FPI_Results']
        FPI_Results = FPI_Results.reshape(-1)[0]
        site = npzfile['site']
        site = site.reshape(-1)[0]
        npzfile.close()

        # Get the instrument and site name from the filename
        inst_name, site_name = os.path.basename(f).split('_')[0:2]

        if Temperature_Fig is None:
            # Add Temperature figure
            Temperature_Fig = plt.figure()
            Temperature_Graph = Temperature_Fig.add_subplot(111)

        if Doppler_Fig is None:
            # Add Doppler figure
            Doppler_Fig = plt.figure()
            Doppler_Graph = Doppler_Fig.add_subplot(111)

	if len(all_colors) == 0:
            all_markers = ['s', 'd', '^', 'x', 'p', '*', '+', 'v', '>', '<']
            all_colors  = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
            all_markers.reverse() # pop() starts from the end
	    all_colors.reverse()

        fmt = {'Color': all_colors.pop(), 'Marker': all_markers.pop()}
                
        (ref_Dop, e_ref_Dop) = FPI.DopplerReference(FPI_Results,reference=reference)

        # Calculate the vertical wind and interpolate it
        ind = FPI.all_indices('Zenith',FPI_Results['direction'])
        ##w = -1. * (FPI_Results['LOSwind'][ind]-ref_Dop[ind]) # -1 because LOS is towards instrument
	w = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]) # LOS is away from instrument
        sigma_w = FPI_Results['sigma_LOSwind'][ind]
        dt = []
        for x in FPI_Results['sky_times'][ind]:
            diff = (x - FPI_Results['sky_times'][0])
            dt.append(diff.seconds+diff.days*86400.)
        dt = np.array(dt)

        # Remove outliers
        ind = abs(w) < 200.

        if sum(ind) <= 1:
            # No good data, just use all ind
                ind = abs(w) > 200.

        if len(ind) == 0:
           raise Exception('%s: No Zenith look directions' % f)

        # Interpolate
        w2 = interpolate.interp1d(dt[ind],w[ind],bounds_error=False,fill_value=0.0)
        sigma_w2 = interpolate.interp1d(dt[ind],sigma_w[ind],bounds_error=False,fill_value=0.0)
        dt = []
        
        for x in FPI_Results['sky_times']:
            diff = (x - FPI_Results['sky_times'][0])
            dt.append(diff.seconds+diff.days*86400.)
        w = w2(dt)
        sigma_w = sigma_w2(dt)

        # Check if specific directions are requested
        if directions is None:
            # Create the plot with each direction plotted in a different colors
            l = [x for x in np.unique(FPI_Results['direction']) if x not in ['Unknown', 'Laser'] ]
        else:
            l = directions
            
        for x in l:
            ind = FPI.all_indices(x,FPI_Results['direction'])

            temp0 = np.min(FPI_Results['sky_times'][ind])
            temp1 = np.max(FPI_Results['sky_times'][ind])
            if t0 is None:
                t0 = temp0
                t1 = temp1
            else:
                if temp0 < t0:
                    t0 = temp0
                if temp1 > t1:
                    t1 = temp1
                    
            # Inital plot of bogus value to get the legend right
            Temperature_Graph.plot(-999,-999,color=fmt['Color'], marker=fmt['Marker'], label='%s, %s' % (inst_name,x))

	    if x == 'Zenith':
		Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind])
                Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2)
            else:
                Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
                Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)
            if x == 'South' or x == 'West':
                Doppler_Wind = -Doppler_Wind

#            if x == 'Zenith':
#                Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind])
#                Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2)
#            else:
#                Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
#                Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)
#            if x == 'North' or x == 'East':
#                Doppler_Wind = -Doppler_Wind
                
#            Doppler_Graph.plot(FPI_Results['sky_times'][ind],Doppler_Wind,color=fmt['Color'],alpha=0.3,marker=None,label='_nolegend_')
            Doppler_Graph.plot(-999,-999,color=fmt['Color'],marker=fmt['Marker'],label='%s, %s' % (inst_name,x))
            
            if ('Clouds' in FPI_Results.keys()) is False:
                    FPI_Results['Clouds'] = None

            if FPI_Results['Clouds'] is None:
                # No cloud data so create plot without alpha
                for (t,y,ey,z,ez) in zip(FPI_Results['sky_times'][ind], FPI_Results['T'][ind], FPI_Results['sigma_T'][ind],Doppler_Wind,Doppler_Error):

		    if (ez < 10.) and (abs(z) < 200.):
			if displayhours:
			    t = (t.hour * 3600 + t.minute * 60 + t.second)/3600.
			    if t > 18:
				t = t-24.

                        # Plot the points and error bars on the graphs
                        Temperature_Graph.errorbar(t,y,yerr=ey,color=fmt['Color'],fmt=fmt['Marker'],label='_nolegend_')
                        Doppler_Graph.errorbar(t,z,yerr=ez,color=fmt['Color'],fmt=fmt['Marker'],label='_nolegend_')
            else:
                # Loop through each point (needed to be done this way to allow setting
                # the alpha value for each individual point)
                for (t,y,ey,z,ez,sky_temp) in zip(FPI_Results['sky_times'][ind], FPI_Results['T'][ind], FPI_Results['sigma_T'][ind],Doppler_Wind,Doppler_Error,FPI_Results['Clouds']['mean'][ind]):
                
                    # Calculate the alpha value to be used for this point.  Scaled linearly
                    # between [.1,1] corresponding to a sky temp range of [full_cloud, full_clear]
                    alpha_val = 1-max(0,min(((sky_temp-full_clear)/(full_cloud-full_clear)),.9))
                
                    # Plot the points and error bars on the graphs
                    Temperature_Graph.errorbar(t,y,yerr=ey,alpha=alpha_val,color=fmt['Color'],fmt=fmt['Marker'],label='_nolegend_')
                    Doppler_Graph.errorbar(t,z,yerr=ez,alpha=alpha_val,color=fmt['Color'],fmt=fmt['Marker'],label='_nolegend_')

    # Plotting font setup
    fontP = FontProperties() 
#    fontP.set_size('10')
#    fontsize = 18 # because legend() and xlabel() do it different ways
    
    # Format the Temperature plot               
    Temperature_Graph.set_ylim(Tmin,Tmax)
    if displayhours:
        t0 = np.floor((t0.hour * 3600 + t0.minute * 60 + t0.second)/3600.) - 24.
        t1 = np.ceil((t1.hour * 3600 + t1.minute * 60 + t1.second)/3600.)

    Temperature_Graph.set_xlim(t0,t1)

    #Temperature_Graph.set_xlim(FPI_Results['sky_times'][0],FPI_Results['sky_times'][-1])
    Temperature_Graph.legend(ncol=4, prop=fontP)
    if not(displayhours):
        Temperature_Graph.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        Temperature_Graph.set_ylabel('Neutral Temperature [K]')#, fontsize = fontsize)
    Temperature_Graph.set_xlabel('Universal Time')#, fontsize = fontsize)
    Temperature_Graph.grid(True)

    # Format the Doppler plot
    Doppler_Graph.set_ylim(Dmin,Dmax)
    Doppler_Graph.set_xlim(t0,t1)
    #Doppler_Graph.set_xlim(FPI_Results['sky_times'][0],FPI_Results['sky_times'][-1])
    Doppler_Graph.legend(ncol=4, prop=fontP)
    Doppler_Graph.xaxis.set_major_formatter(dates.DateFormatter('%H'))
    Doppler_Graph.set_ylabel('Neutral Winds [m/s]')#, fontsize = fontsize)
    if not(displayhours):
        Doppler_Graph.set_xlabel('Universal Time')#, fontsize = fontsize)
        Doppler_Graph.plot([FPI_Results['sky_times'][0],FPI_Results['sky_times'][-1]],[0,0],'k--')
    Doppler_Graph.grid(True)

    # Mark the plots with the Doppler reference type
    Doppler_Fig.text(1.0,0.0,'Doppler ref: %s' % reference, horizontalalignment='right', verticalalignment='bottom')

    # If no cloud data, mark this on the plots
    if FPI_Results['Clouds'] is None:
        Temperature_Fig.text(0.0,0.0,'No cloud data', verticalalignment='bottom')#, fontsize = fontsize)
        Doppler_Fig.text(0.0,0.0,'No cloud data', verticalalignment='bottom')#, fontsize = fontsize)

    # Create title based on if start/stop time are on the same day or not
    if FPI_Results['sky_times'][0].day != FPI_Results['sky_times'][-1].day:
        Temperature_Graph.set_title(site['Abbreviation'] + ':' + (FPI_Results['sky_times'][0]).strftime('%d %b, %Y %H:%M LT') + ' - ' + (FPI_Results['sky_times'][-1]).strftime('%d %b, %Y %H:%M LT'))#, fontsize = fontsize)
        Doppler_Graph.set_title(site['Abbreviation'] + ':' + (FPI_Results['sky_times'][0]).strftime('%d %b, %Y %H:%M LT') + ' - ' + (FPI_Results['sky_times'][-1]).strftime('%d %b, %Y %H:%M LT'))#, fontsize = fontsize)
    else:
        Temperature_Graph.set_title(site['Abbreviation'] + ':' + (FPI_Results['sky_times'][0]).strftime('%d %b, %Y %H:%M') + ' - ' + (FPI_Results['sky_times'][-1]).strftime('%H:%M LT'))#, fontsize = fontsize)
        Doppler_Graph.set_title(site['Abbreviation'] + ':' + (FPI_Results['sky_times'][0]).strftime('%d %b, %Y %H:%M') + ' - ' + (FPI_Results['sky_times'][-1]).strftime('%H:%M LT'))#, fontsize = fontsize)

    return(Temperature_Fig, Temperature_Graph), (Doppler_Fig, Doppler_Graph)

def Map(sites,m,times,MapTemperatures=False, MapDopplers=False, avg_time=15, title_stub=''):
# Function to produce a map of values from several FPIs
#
# INPUTS:
#
# OPTION INPUTS:
#


# OUTPUTS:
#
# HISTORY:
#   Written on 6 December 2012 by Jonathan J. Makela

    myTemperature = None
    Tsize = None
    myDoppler = None
    Dalpha = None

    for t0 in times:
        # draw coastlines and fill the continents (alpha value is set to partial transparancy because for some
        # reason the fill is over top the quivers used to denote the wind vectors
        m.drawcoastlines()
        #m.fillcontinents(alpha=.5)
        m.drawstates()
        
        # Define the q vector as None (needed below in case no values are found)
        q = None
        
        for site in sites:
            # Find out the cloud status at this site at this time
            if np.size(sites[site]['bw_date']) > 0:
                closest = sorted(sites[site]['bw_date'], key=lambda d:abs(t0-d))[0]
                sky_temp = sites[site]['bw_sky'][(sites[site]['bw_date'] == closest).nonzero()][0]
            else:
                sky_temp = -999
                        
            # convert to map projection coords.
            # Note that lon,lat can be scalars, lists or numpy arrays.
            xpt,ypt = m(sites[site]['lon'],sites[site]['lat'])
            # convert back to lat/lon
            lonpt, latpt = m(xpt,ypt,inverse=True)
            if sky_temp < -25.0:
                m.plot(xpt,ypt,'bo')  # plot a blue dot there
            else:
                m.plot(xpt,ypt,'bx')
                
            # put some text next to the dot, offset a little bit
            # (the offset is in map projection coordinates)
            plt.text(xpt+50000,ypt+50000,'%s' % (site))

            # Load a npz file with processed data and read required variables
            npzfile = np.load(sites[site]['fname'],allow_pickle=True)
            FPI_Results = npzfile['FPI_Results']
            FPI_Results = FPI_Results.reshape(-1)[0]
            
            direction = []
            cmin = 500
            cmax = 1200
            dmin = -200
            dmax = 200
            
            xt,yt = m(0,0)
            m.scatter(xt,yt,c=500,s=250,vmin=cmin,vmax=cmax, edgecolor='none')

            # Find the zero offset
            ind = FPI.all_indices('Zenith',FPI_Results['direction'])

            # Grab times, convert to decimal hour
            zt_h = FPI.dt2h(FPI_Results['sky_times'][ind])
            all_h = FPI.dt2h(FPI_Results['sky_times'])
    
            # Grab the zenith Doppler values
            z = FPI_Results['LOSwind'][ind]

            # Perform interpolation
            if len(z) > 0:
                ref_Dop = np.interp(all_h, zt_h, z)
            else:
                ref_Dop = np.nan*all_h

            for x in np.unique(FPI_Results['direction']):
                if x != 'Unknown':
                    ind = FPI.all_indices(x,FPI_Results['direction'])
                    
                    # Find times within 15 minutes of a requested time
                    lt = FPI_Results['sky_times'][ind]

                    # Grab temperatures
                    if MapTemperatures:
                        Temperature = FPI_Results['T'][ind]
                        e_Temperature = FPI_Results['sigma_T'][ind]

                    # Grab Dopplers
                    if MapDopplers:
                        Doppler = (FPI_Results['LOSwind'][ind]-ref_Dop[ind])/FPI_Results['ze_corr'][ind]
                        e_Doppler = FPI_Results['sigma_LOSwind'][ind]

                    # Calculate time differece from the requested time
                    dt = t0 - lt

                    # Take the weighted average of the Temperature and Wind
                    ind = (dt <= datetime.timedelta(minutes = avg_time)) & (dt >= datetime.timedelta(minutes = 0))

                    # Calculate the average values within the window
                    if ind.any():
                        
                        if sum(ind == True) > 1:
                            if MapTemperatures:
                                print Temperature[ind]
                                print abs(1/e_Temperature[ind])
                                myTemperature, eMyTemperature = FPI.weighted_avg_and_std(Temperature[ind], abs(1/e_Temperature[ind]))
                        
                            if MapDopplers:
                                myDoppler, eMyDoppler = FPI.weighted_avg_and_std(Doppler[ind], abs(1/e_Doppler[ind]))
                        else:
                            if MapTemperatures:
                                myTemperature = Temperature[ind]
                                eMyTemperature = e_Temperature[ind]

                            if MapDopplers:
                                myDoppler = Doppler[ind]
                                eMyDoppler = e_Doppler[ind]

                        if MapTemperatures:
                            if (myTemperature < cmax) & (myTemperature > cmin) & (eMyTemperature < 100) & (sky_temp < -25.0):
                                sites[site]['LookDirection'][x][2] = t0
                                sites[site]['LookDirection'][x][3] = myTemperature
                                Tsize = 250.0
                            else:
                                # Bad point, use last value but modify based on how old the
                                # measurement is
                                if sites[site]['LookDirection'][x][2] is not None:
                                    myTemperature = sites[site]['LookDirection'][x][3]
                                    Tsize = 250.*max([1 - (t0 - sites[site]['LookDirection'][x][2]).seconds/3600., 0])
                                else:
                                    myTemperature = 0
                                    Tsize = 0

                        if MapDopplers:
                            if (myDoppler < dmax) & (myDoppler > dmin) & (eMyDoppler < 100) & (sky_temp < -25.0):
                                sites[site]['LookDirection'][x][4] = t0
                                sites[site]['LookDirection'][x][5] = myDoppler
                                Dalpha = 1.0
                            else:
                                # Bad point, use last value but modify based on how old the
                                # measurement is
                                if sites[site]['LookDirection'][x][4] is not None:
                                    myDoppler = sites[site]['LookDirection'][x][5]
                                    Dalpha = 1.0*max([1 - (t0 - sites[site]['LookDirection'][x][2]).seconds/3600., 0])
                                else:
                                    myDoppler = -999
                                    Dalpha = 0.0
                    else:
                        # No data for this site/direction, use last value but modify based on how old the
                        # measurement is
                        if MapTemperatures:
                            if sites[site]['LookDirection'][x][2] is not None:
                                myTemperature = sites[site]['LookDirection'][x][3]
                                Tsize = 250.*max([1 - (t0 - sites[site]['LookDirection'][x][2]).seconds/3600., 0])
                            else:
                                myTemperature = 0
                                Tsize = 0

                        if MapDopplers:
                            if sites[site]['LookDirection'][x][4] is not None:
                                myDoppler = sites[site]['LookDirection'][x][5]
                                Dalpha = 1.0*max([1 - (t0 - sites[site]['LookDirection'][x][2]).seconds/3600., 0])
                            else:
                                myDoppler = -999
                                Dalpha = 0.0
                    
                    # The location at which to plot the current data                    
                    lon = sites[site]['LookDirection'][x][1]
                    lat = sites[site]['LookDirection'][x][0]

                    # Plot this value on the map
                    xpt,ypt = m(lon,lat)

                    if MapTemperatures:
                        m.scatter(xpt,ypt,c=myTemperature,s=Tsize,vmin=cmin,vmax=cmax, edgecolor='none')

                    if MapDopplers:
                        if (x == 'East' or x == 'West') and (myDoppler != -999):
                            u = myDoppler
                            v = 0
                    
                            # Draw quivers
                            q = m.quiver(xpt,ypt,u,v,angles='uv',scale=1000, color='black',alpha=Dalpha)
                        elif (x == 'North' or x == 'South') and (myDoppler != -999):
                            u = 0
                            v = myDoppler
                    
                            # Draw quivers
                            q = m.quiver(xpt,ypt,u,v,angles='uv',scale=1000, color='black',alpha=Dalpha)

        cb = m.colorbar()
        cb.set_alpha(1)
        cb.set_label('[K]')
            
        if q is not None:
            x,y = m(-77.0,31.0)
            p = plt.quiverkey(q,x,y,50,"50 m/s",coordinates='data',color='r')

        # draw parallels and meridians.
        # label parallels on left and bottom
        # meridians on bottom and left
        # labels = [left,right,top,bottom]
        parallels = np.arange(0.,81,5.)
        m.drawparallels(parallels,labels=[True,False,False,True])
        
        meridians = np.arange(-100.,-60.,5.)
        m.drawmeridians(meridians,labels=[True,False,False,True])

        plt.title('%s %s UT' % (title_stub, t0.strftime('%d %b %Y %H:%M')))
        plt.savefig('NATION_%s.png' % t0.strftime('%Y%m%d_%H%M%S'))
        plt.clf()

def DisplayRaw(f, cmin=None, cmax=None):
#
# Function to display information about a raw FPI image and plot
# the ring pattern.
#
# INPUTS:
#   f - full path name to image to display
#
# OPTIONAL INPUTS:
#   cmin, cmax - color limits (min, max) to scale the resultant plot by
#
# OUTPUT:
#   p - a reference to the plot created
#
# HISTORY:
#   Written by Jonathan J. Makela on 19 Dec 2012

    d = FPI.ReadIMG(f)

    print 'Filename: %s' % f
    print 'Exposure time: %.1f s' % d.info['ExposureTime']
    print 'CCD Temperature: %.1f C' % d.info['CCDTemperature']
    print 'SkyScanner Direction (az, ze): (%.1f, %.1f)' % (d.info['azAngle'], d.info['zeAngle'])
    
    p = plt.matshow(np.reshape(d.getdata(),d.size))
    if cmin is not None:
        plt.clim([cmin,cmax])
        plt.title(d.info['LocalTime'].strftime('%Y-%m-%d %H:%M:%S LT'))
        plt.colorbar(ax=p.axes,orientation='vertical',shrink=0.8)

        p.axes.get_xaxis().set_visible(False)
        p.axes.get_yaxis().set_visible(False)

    return p
