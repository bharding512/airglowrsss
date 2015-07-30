
# Functions to process FPI data on remote2

import matplotlib
matplotlib.use('AGG')

import FPI
import glob
from lmfit import Parameters
#import logging
from optparse import OptionParser
import MySQLdb as mdb
import datetime
import numpy as np
import os
import BoltwoodSensor
import X300Sensor
import FPIDisplay
import pytz
import multiprocessing
import subprocess
import re
import fpiinfo
import matplotlib.pyplot as plt
import traceback
import FPIResults
import FPIwindfield
from matplotlib import dates
import shutil

# For quick testing of these modules, force a reload
reload(FPIDisplay)
reload(FPI)
reload(fpiinfo)

   
    
def quality_hack(instr_name, year, doy, FPI_Results, logfile):
    '''
    This is a hack function, giving us a place to document and implement
    any manual corrections that are needed for specific dates.
    Ideally, this function should do nothing, and our processing should
    elegantly handle all fault cases. Unfortunately, this is not the case.
    '''
    # ('minime08',2012,319) Doppler reference appears to fail. Only use the
    # first half of the night for Doppler referencing.
    if instr_name == 'minime08' and year == 2012 and doy == 319:
        dref,drefe = FPI.DopplerReference(FPI_Results, reference='laser',AVERAGING_TIME=[17.,0.])
        FPI_Results['LOSwind'] = FPI_Results['LOSwind'] - dref
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                      'FPIprocess.quality_hack(): Manual intervention to tweak Doppler reference.\n')
                      
    # ('minime08', 2013, 274): try to fix laser drift
    if instr_name == 'minime08' and year == 2013 and doy == 274:
        t = FPI_Results['sky_times']
        tsec = np.array([(ti - t[0]).total_seconds() for ti in t])
        drift_corr = -100.0 * (tsec[-1] - tsec)/(tsec[-1] - tsec[0])
        FPI_Results['LOSwind'] = FPI_Results['LOSwind'] + drift_corr
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                      'FPIprocess.quality_hack(): Manual intervention to fix laser drift.\n')
    
    return FPI_Results
    
    

    
    
def process_instr(instr_name ,year, doy, reference='laser', use_npz = False, zenith_times=[17.,7.], wind_err_thresh=100., temp_err_thresh=100., cloud_thresh = [-22.,-10.]):
    '''
    Process all the data from the instrument with name instr_name on
    the day specified by (year, doy). This function looks in the appropriate
    directory on remote2, analyzes the laser and sky images to obtain line of
    sight winds, temperatures, etc., and saves them in a .npz file in the
    appropriate place. It also generates summary plots and sends them to the
    website. Finally, it generates summary ASCII files for Madrigal.
    INPUTS:
        instr_name - 'minime01', 'minime02', etc.
        year - int
        doy - int (1-indexed, right?)
    OPTIONAL INPUTS:
        reference - str, 'laser' or 'zenith'. Passed on to FPI.ParameterFit(...)
        use_npz - bool, if True, the call to FPI.ParameterFit(...) will be skipped,
                  and the previously-analyzed data in the saved npz file will be used
                  instead. Cloud and X300 data will be reloaded and added to the npz file.
                  New wind and temperature plots will be sent to the website,
                  along with the updated log file.
        zenith_times - 2x1 array, The local times outside of which the zenith wind should
                       be ignored for use as a reference.
        wind_err_thresh - float, m/s. Samples with a fit error above this should get a quality
                            flag of 2.
        temp_err_thresh - float, K. Samples with a fit error above this should get a quality
                            flag of 2.      
        cloud_thresh - [float,float], K. The two cloud sensor readings that indicate 
                        partially- and fully-cloudy. This affects the quality flag.              
    OUTPUTS:
        warnings - str - If this script believes a manual check of the data
                   is a good idea, a message will be returned in this string.
                   If not, the empty string will be returned. For now, the
                   whole log will be the message. 
        
    '''
    

    # Define constants that do not depend on the site
    direc_tol = 10.0 # tolerance in degrees to recognize a look direction with
    fpi_dir =               '/rdata/airglow/fpi/'
    bw_dir =                '/rdata/airglow/templogs/cloudsensor/'
    x300_dir =              '/rdata/airglow/templogs/x300/'
    results_stub =          '/rdata/airglow/fpi/results/'
    temp_plots_stub =       '/rdata/airglow/fpi/results/temporary_plots/' # where to save pngs on remote2
    scp_user =              'data@airglow.ece.illinois.edu' # User to scp to
    db_image_stub =         'SummaryImages/' # relative path from web server directory on airglow
    db_log_stub =           'SummaryLogs/'
    web_images_stub =       '/data/SummaryImages/' # absolute location on airglow of summary images
    web_logs_stub =         '/data/SummaryLogs/' # absolute location on airglow of summary logs
    madrigal_stub =         '/rdata/airglow/database/' # where reduced ascii txt and png files are saved for Madrigal
    share_stub    =         '/rdata/airglow/share/' # where to save a copy of the .npz file for sharing with collaborators
    
    notify_the_humans = False # default

    nominal_dt = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1)

    # Import the site information
    site_name = fpiinfo.get_site_of(instr_name, nominal_dt)
    site = fpiinfo.get_site_info(site_name)
    # Import the instrument information
    instrument = fpiinfo.get_instr_info(instr_name, nominal_dt)

    # Create "minime05_uao_20130729" string
    datestr = nominal_dt.strftime('%Y%m%d')
    instrsitedate = instr_name + '_' + site_name + '_' + datestr



    # Construct the directories to the relevant data products
    data_stub = fpi_dir + instr_name + '/' + site_name + '/'
    bw_dir_stub =       bw_dir + site_name + '/'
    bw_name_stub =      'Cloud_' + site_name + '_'
    x300_dir_stub =     x300_dir + site_name + '/'
    x300_name_stub =    'TempL_' + site_name + '_'
    

    # Determine the times between which files will be accepted.
    # Use local-time noon on doy to local-time noon on doy+1
    local = pytz.timezone(site['Timezone'])
    start_dt = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1, hours = 12)
    stop_dt  = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1, hours = 36)
    start_dt = local.localize(start_dt)
    stop_dt  = local.localize(stop_dt)

    
    # Create a list of relevant filenames by searching in adjacent
    # days' directories also, since we might have problems with UT/LT,
    # time zones, etc.
    laser_fns = []
    sky_fns = []
    for day_offset in [-1, 0, 1]: # search folders around the day of interest
        dir_dt = nominal_dt + datetime.timedelta(days = day_offset)
        # Create the YYYYMMDD date format
        yearstr = dir_dt.strftime('%Y')
        datestr = dir_dt.strftime('%Y%m%d')
        # Create the directory name for the data to be processed
        data_dir = data_stub + yearstr + '/' + datestr + '/'
        # Look in the data directory, and grab the files if they were
        # taken between the start and stop times.
        # First, lasers:
        fns_all = glob.glob(data_dir + '*L*.img')
        for fn in fns_all:
            d = FPI.ReadIMG(fn)
            dtime = local.localize(d.info['LocalTime'])
            if dtime > start_dt and dtime < stop_dt:
                laser_fns.append(fn)
        # Second, sky:
        fns_all = glob.glob(data_dir + '*X*.img')
        for fn in fns_all:
            d = FPI.ReadIMG(fn)
            dtime = local.localize(d.info['LocalTime'])
            if dtime > start_dt and dtime < stop_dt:
                sky_fns.append(fn) 
    laser_fns.sort()
    sky_fns.sort()


    if not laser_fns and not sky_fns:
        raise Exception('No %s data found between %s and %s. \n' % (instr_name, str(start_dt), str(stop_dt)))

    datestr = nominal_dt.strftime('%Y%m%d')
    
    # Open a logfile
    logname = results_stub + instrsitedate + '.log'
    if use_npz: # append to the previous log file
        logfile = open(logname,'a') # overwrite previous log
        logfile.write('\n' + datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Rerunning processing to obtain new cloud and temperature data, and updating plots.\n')
    else: # create a new log file
        logfile = open(logname,'w') # overwrite previous log
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Logfile Created\n')
    

    if not laser_fns and sky_fns and reference=='laser': # This is not a big deal if reference=='zenith'
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No %s laser data found between %s and %s. Sky data found. <BADLASER> \n' % (instr_name, str(start_dt), str(stop_dt)))
        logfile.close()
        raise Exception('No %s laser data found between %s and %s. Sky data found.\n' % (instr_name, str(start_dt), str(stop_dt)))
    

    npzname = results_stub + instrsitedate + '.npz' # the name of the npz file to save to
    Diagnostic_Fig = plt.figure(dpi=300, figsize=(8,6)) # Figure for diagnostics to be drawn to   
    if use_npz: # Load previously-analyzed results
        npzfile = np.load(npzname)
        FPI_Results = npzfile['FPI_Results']
        FPI_Results = FPI_Results.reshape(-1)[0]
        del npzfile.f # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
        npzfile.close()
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                'Laser and sky image skipped. Loaded npz file: %s\n' % npzname.split('/')[-1])
    else:
        # Try to analyze the data
        try:
            # Run the analysis
            (FPI_Results, notify_the_humans) = FPI.ParameterFit(instrument, site, laser_fns, sky_fns, \
                    direc_tol=direc_tol, N=instrument['N'], N0=instrument['N0'], N1=instrument['N1'], \
                                logfile=logfile, diagnostic_fig=Diagnostic_Fig, reference=reference, \
                                zenith_times=zenith_times)
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                'Laser and sky image analysis complete.\n')
        except: 
            # FPI.ParameterFit crashed. For now, write the log and re-raise the Exception.
            tracebackstr = traceback.format_exc()
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error analyzing %s. Traceback listed below.\n-----------------------------------\n%s\n-----------------------------------\n' % (datestr,tracebackstr))
            notify_the_humans = True # This is redundant, since raising will ensure humans are notified.
    #        Diagnostic_Fig = plt.figure()
    #        Diagnostic_Fig.text(0.5, 0.5, 'Error - see log file', fontsize = 20, ha='center')
            logfile.close()
            raise

        # Do any necessary manual corrections, as per the quality_hack function
        FPI_Results = quality_hack(instr_name, year, doy, FPI_Results, logfile)


    # Grab the SVN number
    p = subprocess.Popen('svnversion /usr/local/share/airglowrsss/Python/modules/',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdout,stderr) = p.communicate()
    sv = re.split(':|\n', stdout)[0]

    # Save the SVN version number
    FPI_Results['SVNRevision'] = sv
    
    # Try to load the Boltwood and X300 data. This is inside a try block so that
    # the npz will still save if there is a problem.
    FPI_Results['Clouds'] = None # defaults
    FPI_Results['Dome']   = None
    FPI_Results['Inside'] = None
    try: 
        # call modules from BoltWood for the two days required and combine the data
        if (FPI_Results['sky_times'][0].strftime("%Y%m%d") == FPI_Results['sky_times'][-1].strftime("%Y%m%d")):
            bw_date, bw_sky, bw_amb = BoltwoodSensor.ReadTempLog(
                    bw_dir_stub + bw_name_stub + FPI_Results['sky_times'][0].strftime("%Y%m%d") + ".txt",
                    tz=site['Timezone']
                    )
            x_date, x_dome, x_in = X300Sensor.ReadTempLog(
                    x300_dir_stub + x300_name_stub + FPI_Results['sky_times'][0].strftime("%Y%m%d") + ".txt",
                    tz=site['Timezone']
                    )
        else:
            bw_date1, bw_sky1, bw_amb1 = BoltwoodSensor.ReadTempLog(
                    bw_dir_stub + bw_name_stub + FPI_Results['sky_times'][0].strftime("%Y%m%d") + ".txt",
                    tz=site['Timezone'])
            bw_date2, bw_sky2, bw_amb2 = BoltwoodSensor.ReadTempLog(
                    bw_dir_stub + bw_name_stub + FPI_Results['sky_times'][-1].strftime("%Y%m%d") + ".txt", 
                    tz=site['Timezone']
                    )
            bw_date = np.hstack((bw_date1,bw_date2))
            bw_sky = np.hstack((bw_sky1, bw_sky2))
            bw_amb = np.hstack((bw_amb1, bw_amb2))
                
            x_date1, x_dome1, x_in1 = X300Sensor.ReadTempLog(
                    x300_dir_stub + x300_name_stub + FPI_Results['sky_times'][0].strftime("%Y%m%d") + ".txt", 
                    tz=site['Timezone']
                    )
            x_date2, x_dome2, x_in2 = X300Sensor.ReadTempLog(
                    x300_dir_stub + x300_name_stub + FPI_Results['sky_times'][-1].strftime("%Y%m%d") + ".txt", 
                    tz=site['Timezone']
                    )
            x_date = np.hstack((x_date1,x_date2))
            x_dome = np.hstack((x_dome1,x_dome2))
            x_in = np.hstack((x_in1,x_in2))
            
        # Create a data structure containing the sensor temperatures
        Clouds = np.nan*np.zeros((np.size(FPI_Results['sky_times']),3))
        Dome = np.nan*np.zeros((np.size(FPI_Results['sky_times']),3))
        Inside = np.nan*np.zeros((np.size(FPI_Results['sky_times']),3))
        count = 0
        for (t, dt) in zip(FPI_Results['sky_times'], FPI_Results['sky_intT']):
            # Start and stop time of image data
            t0 = np.array(t                                   )
            t1 = np.array(t + datetime.timedelta(seconds = dt))
            # Find the corresponding times in the Boltwood data
            i = ((bw_date >= t0) & (bw_date <= t1))
            if i.any():
                Clouds[count,:] = [bw_sky[i].mean(), bw_sky[i].max(), bw_sky[i].min()]
            else:
                Clouds[count,:] = [np.nan,np.nan,np.nan]
                  
            # Find the corresponding times in the X300 data
            i = ((x_date >= t0) & (x_date <= t1))
            if i.any():
                Dome[count,:] = [x_dome[i].mean(), x_dome[i].max(), x_dome[i].min()]
                Inside[count,:] = [x_in[i].mean(), x_in[i].max(), x_in[i].max()]
            else:
                Dome[count,:] = [np.nan,np.nan,np.nan]
                Inside[count,:] = [np.nan,np.nan,np.nan]
                
            count = count+1

        if(np.size(bw_date) == 0):
            FPI_Results['Clouds'] = None
        else:
            FPI_Results['Clouds'] = {'mean': Clouds[:,0], 'max': Clouds[:,1], 'min': Clouds[:,2]}
        
        if(np.size(x_date) == 0):
            FPI_Results['Dome'] = None
            FPI_Results['Inside'] = None
        else:
            FPI_Results['Dome'] = {'mean': Dome[:,0], 'max': Dome[:,1], 'min': Dome[:,2]}
            FPI_Results['Inside'] = {'mean': Inside[:,0], 'max': Inside[:,1], 'min': Inside[:,2]}

        # Write things to the log
        if FPI_Results['Clouds'] is None:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No Boltwood cloud sensor data found.\n')
            c = np.nan*np.zeros(len(FPI_Results['LOSwind']))
            FPI_Results['Clouds'] = {'mean': c, 'max': c, 'min': c}
            notify_the_humans = True
        else:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Found and loaded Boltwood cloud sensor data.\n')
        if FPI_Results['Dome'] is None:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No X300 temperature sensor data found.\n')
            #notify_the_humans = True
        else:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Found and loaded X300 temperature sensor data.\n')
    except: # There was an error in the Boltwood or X300 code. Write the log but continue.
        c = np.nan*np.zeros(len(FPI_Results['LOSwind']))
        FPI_Results['Clouds'] = {'mean': c, 'max': c, 'min': c}
        tracebackstr = traceback.format_exc()
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error obtaining Boltwood or X300 data for %s. Traceback listed below. Analysis will continue without these data.\n-----------------------------------\n%s\n-----------------------------------\n' % (datestr,tracebackstr))
        notify_the_humans = True
        
        
    # Helper function to determine if laser wavelength drift is occuring
    def laser_is_drifting():
        '''
        Return True if the laser is suspected to be drifting.
        Return False if zenith reference was used.
        
        The laser is "suspected to be drifting" if the following 2 criteria are satisfied:
        1) The vertical wind at the beginning of the night and the end of the night are different by > 30 m/s.
        2) The laser intensity varies by more than 20% from the median across the night.
        '''
        fpir = FPI_Results
        if fpir['reference']=='zenith':
            return False
        
        direction = fpir['direction']
        LOSwind = fpir['LOSwind']
            
        direc = 'Zenith'
        w = np.array([si for (si,d) in zip(LOSwind, direction) if d == direc])

        lasI = fpir['laser_value']['I']
        lasIe = fpir['laser_stderr']['I']

        # Check if laser varies by more than 20 %
        las_flag = sum(abs(lasI - np.median(lasI))/np.median(lasI) > 0.2) > 2
        # Check if vertical wind drifts by more than 30 m/s
        wstart = np.median(w[:5])
        wend   = np.median(w[-5:])
        w_flag = abs(wend-wstart) > 30.
        
        return w_flag and las_flag
    
    
    ######### Quality Flags #########
    # Compute quality flags and write any tripped flags to log
    [t_flag_manual, w_flag_manual] = fpiinfo.get_bad_data_flags(instr_name, nominal_dt)
    laser_drift = laser_is_drifting()
    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + ' Quality Flag Calculation:\n')
    logfile.write('\t\tManual Wind Flag: %i\n' % w_flag_manual)
    logfile.write('\t\tManual Temperature Flag: %i\n' % t_flag_manual)
    logfile.write('\t\tLaser Drift: %s\n' % laser_drift)
    logfile.write('\t\tSample-specific flags:\n')
    wind_quality_flag = np.zeros(len(FPI_Results['sky_times']))
    temp_quality_flag = np.zeros(len(FPI_Results['sky_times']))
    for ii in range(len(FPI_Results['sky_times'])): # Loop over samples
        logfile.write('\t\t%s ' % FPI_Results['sky_fns'][ii].split('/')[-1])
        t_flag = 0 # default
        w_flag = 0 # default
        if (FPI_Results['skyI'][ii] < instrument['skyI_quality_thresh']):
            # The sky brightness is low enough that OH is probably an issue.
            t_flag = 1
            w_flag = 1
            logfile.write('[skyI low W1T1] ')
        if (reference == 'zenith'): # Zenith Processing
            t_flag = 1
            w_flag = 1
            logfile.write('[zen ref W1T1] ')
        if laser_drift:
            w_flag = 1
            logfile.write('[laser drift W1] ')
        c = FPI_Results['Clouds']['mean'][ii]
        if np.isnan(c): # There isn't a cloud sensor, or it isn't working
            t_flag = 1
            w_flag = 1
            logfile.write('[no cloud data W1T1] ')
        if c > cloud_thresh[0]: # It's at least a little bit cloudy. Caution.
            t_flag = 1
            w_flag = 1
            logfile.write('[cloud>%.0f: W1T1] '%cloud_thresh[0])
        if c > cloud_thresh[1]: # It's definitely cloudy
            t_flag = 1
            w_flag = 2
            logfile.write('[cloud>%.0f: W2T1] '%cloud_thresh[1])
        w_fit_err = FPI_Results['sigma_fit_LOSwind'][ii]
        if (w_fit_err > wind_err_thresh): # Sample is not trustworthy at all
            w_flag = 2
            logfile.write('[sigma_wind large: W2] ')
        t_fit_err = FPI_Results['sigma_T'][ii]
        if (t_fit_err > temp_err_thresh): # Sample is not trustworthy at all
            w_flag = 2
            logfile.write('[sigma_T large: T2] ')
        # A manual override can increase the calculated flag, but not decrease
        t_flag = max(t_flag,t_flag_manual)
        w_flag = max(w_flag,w_flag_manual)
        wind_quality_flag[ii] = w_flag
        temp_quality_flag[ii] = t_flag
        logfile.write('[FINAL: W%iT%i]\n' % (w_flag,t_flag))
    FPI_Results['wind_quality_flag'] = wind_quality_flag
    FPI_Results['temp_quality_flag'] = temp_quality_flag
            


    # Save the results
    np.savez(npzname, FPI_Results=FPI_Results, site=site, instrument=instrument)
    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Results saved to %s\n' % npzname)
    if site['share']: # save a copy of the npz file in a separate folder
        npznameshare = share_stub + site_name + '/' + instrsitedate + '.npz'
        np.savez(npznameshare, FPI_Results=FPI_Results, site=site, instrument=instrument)
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Results also saved to %s\n' % npznameshare)
    if instrument['send_to_madrigal']: # save the summary ASCII file to send to the Madrigal database
        asciiname = madrigal_stub + instrsitedate + '.txt'
        FPIResults.CreateL1ASCII(npzname, asciiname)
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'ASCII results saved to %s\n' % asciiname)

    # Load the results
    npzfile = np.load(npzname)
    FPI_Results = npzfile['FPI_Results']
    FPI_Results = FPI_Results.reshape(-1)[0]
    site = npzfile['site']
    site = site.reshape(-1)[0]
    del npzfile.f # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
    npzfile.close()
    
    # Try to make plots
    try: 
        # Plot some quick-look single-station data (LOSwinds and Temps)
        (Temperature_Fig, Temperature_Graph), (Doppler_Fig, Doppler_Graph) = FPIDisplay.PlotDay(npzname, reference = reference, Zenith_Times=zenith_times)
        
        # Add Level 1 diagnostics to diagnostics fig
        ax = Diagnostic_Fig.add_subplot(427) # TODO: generalize diagnostic figure generation?
        if FPI_Results['Clouds'] is not None:
            ct = FPI_Results['sky_times']
            cloud = FPI_Results['Clouds']['mean']
            ax.plot(ct, cloud, 'k.-')
            ax.plot([ct[0],ct[-1]], [cloud_thresh[0],cloud_thresh[0]], 'k--')
            ax.plot([ct[0],ct[-1]], [cloud_thresh[1],cloud_thresh[1]], 'k--')
            ax.set_xlim([ct[0] - datetime.timedelta(hours=0.5), ct[-1] + datetime.timedelta(hours=0.5)])
            ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
            ax.set_ylim([-50,0])
            ax.set_ylabel('Cloud indicator [degrees C]')
            ax.grid(True)
            ax.set_xlabel('Universal Time, [hours]')
        #obj1.plot_diagnostics(ax)
        Diagnostic_Fig.tight_layout()
        # Concatenate level0 and level1 logs
        #logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Level 1 processing log:\n' + \
        #              '--------------------------------------------------\n' + obj1.log + \
        #              '\nEnd Level 1 processing log\n' + \
        #              '--------------------------------------------------\n')
    except: 
        # Summary plots crashed. We still want to send the diagnostics and log to the
        # website, so continue on with dummy plots.
        tracebackstr = traceback.format_exc()
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error making plots for %s. Traceback listed below. Dummy plots sent to website.\n-----------------------------------\n%s\n-----------------------------------\n' % (datestr,tracebackstr))
        notify_the_humans = True
        Temperature_Fig = plt.figure()
        Temperature_Fig.text(0.5, 0.5, 'Plotting error - see log file', fontsize = 20, ha='center')
        Doppler_Fig = plt.figure()
        Doppler_Fig.text(0.5, 0.5, 'Plotting error - see log file', fontsize = 20, ha='center')
        
        
    # Create (or update) the level 3 wind field map, if use_npz=False (so that use_npz=True is fast)
    network_name = site['Network']
    gif_fn = None # filename of the gif to send to the website
    if network_name in  ['nation','peru'] and not use_npz:
        try:
            wf = FPIwindfield.WindField([network_name], year, doy, timestep_fit = 1800)
            if len(wf.instr_names) >= 3: # If there are 3 sites, try the inversion
                losfig = wf.plot_los_winds()
                plt.close(losfig)
                wf.run_inversion()
                gif_fn = wf.make_quiver_gif(show_vert_wind=True)
                shutil.copy(gif_fn, '/home/bhardin2/public_html/windfield_movies/')
                quicklook_fig, quicklook_fn = wf.make_quicklook(show_vert_wind=True)
                shutil.copy(quicklook_fn, '/home/bhardin2/public_html/windfield_quicklooks/')
                plt.close(quicklook_fig)
        except:
            tracebackstr = traceback.format_exc()
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error creating windfield quicklook plot. Traceback listed below.\n-----------------------------------\n%s\n-----------------------------------\n' % (tracebackstr,))
            notify_the_humans = True
    
    
    # Update the database with the summary images and log file
    summary_figs = [Diagnostic_Fig,
                    Temperature_Fig,
                    Doppler_Fig,]
    summary_fns = [instrsitedate + '_diagnostics.png',
                   instrsitedate + '_temperature.png', # e.g., 'minime05_uao_20130107_temperature.png'
                   instrsitedate + '_winds.png',]
    db_ids = [instrument['sql_diagnostics_id'],
              instrument['sql_temperatures_id'],
              instrument['sql_winds_id'], ]
    if use_npz: # don't update diagnostics fig
        summary_figs.pop(0)
        summary_fns.pop(0)
        db_ids.pop(0)

    site_id = site['sql_id']
    utc = pytz.utc # Define timezones
    # Start and stop time of observations
    d = FPI.ReadIMG(sky_fns[0])
    dtime = local.localize(d.info['LocalTime'])
    startut = dtime.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')
    d = FPI.ReadIMG(sky_fns[-1])
    dtime = local.localize(d.info['LocalTime'])
    stoput = dtime.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')  
    # Open the database (see http://zetcode.com/databases/mysqlpythontutorial/ for help)
    # Read the user and password from a file.
    con = mdb.connect(host='airglow.ece.illinois.edu', db='webdatabase', read_default_file="~/.my.cnf")
    cur = con.cursor()  
    # Create the summary images
    for fig, fn, db_id in zip(summary_figs, summary_fns, db_ids):
        fig.savefig(temp_plots_stub + fn) # save it in the remote2 data dir
        flag = subprocess.call(['scp', temp_plots_stub + fn, scp_user + ':' + web_images_stub + fn]) # send to airglow
        if flag != 0: # Sending png to airglow failed
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error sending %s to airglow server for displaying on website.\n' % fn)
            notify_the_humans = True            
        # If we are supposed to put this png in Madrigal, put it in the Madrigal database directory
        send_to_madrigal = instrument['send_to_madrigal'] and 'diagnostic' not in fn # don't send diag. fig
        if send_to_madrigal:
            fig.savefig(madrigal_stub + fn) # save the figure to the madrigal directory
        os.remove(temp_plots_stub + fn) # remove the png
        # update the database
        # Send winds png. First find out if the entry is in there (i.e., we are just updating the png)
        sql_cmd = 'SELECT id FROM DataSet WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (site_id, db_id, startut)
        cur.execute(sql_cmd)
        rows = cur.fetchall()
        log_fn = db_log_stub + instrsitedate + '_log.log'
        if len(rows) == 0: # Create the entry
            sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime, SummaryImage, LogFile) VALUES(%d, %d, \"%s\", \"%s\", \"%s\", \"%s\")' % (site_id, db_id, startut, stoput, db_image_stub + fn, log_fn)
            cur.execute(sql_cmd)
        else: # Entry exists. Update it.
            sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\",LogFile=\"%s\" WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (db_image_stub + fn, log_fn, site_id, db_id, startut)
            cur.execute(sql_cmd)
        # close figure to save memory
        plt.close(fig)
            
    # Send level 3 windfield gif to website, if we made one
    if gif_fn is not None:
        # Since the start/stop time of a network is ill-defined, just make up a start time: 21 UT on the 1st night, and 12 UT on the second night. This will have to be updated for networks other than NATION.
        d = FPI.ReadIMG(sky_fns[0])
        # Define midnight on the first night
        dtime = d.info['LocalTime'] - datetime.timedelta(hours=12)
        midnight = dtime - datetime.timedelta(hours=dtime.hour, minutes=dtime.minute, seconds=dtime.second)
        start = midnight + datetime.timedelta(hours=21)
        stop  = midnight + datetime.timedelta(hours=36)
        startut = start.strftime('%Y-%m-%d %H:%M:%S')
        stoput = stop.strftime('%Y-%m-%d %H:%M:%S')
        network_info = fpiinfo.get_network_info(network_name)
        network_id   = network_info['sql_id']
        gif_id       = network_info['quicklook_gif_id']
        flag = subprocess.call(['scp', gif_fn, scp_user + ':' + web_images_stub + gif_fn.split('/')[-1]]) # send to airglow
        if flag != 0: # Sending png to airglow failed
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error sending "%s" to airglow server for displaying on website.\n' % gif_fn)
            notify_the_humans = True
        # Register the gif. First find out if the entry is in there (i.e., we are just updating it)
        sql_cmd = 'SELECT id FROM DataSet WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (network_id, gif_id, startut)
        cur.execute(sql_cmd)
        rows = cur.fetchall()
        if len(rows) == 0: # Create the entry
            sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime, SummaryImage) VALUES(%d, %d, \"%s\", \"%s\", \"%s\")' % (network_id, gif_id, startut, stoput, db_image_stub + gif_fn.split('/')[-1])
            cur.execute(sql_cmd)
        else: # Entry exists. Update it.
            sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\" WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (db_image_stub + gif_fn.split('/')[-1], network_id, gif_id, startut)
            cur.execute(sql_cmd)
            

    # Send the log file
    logfile.close()
    subprocess.call(['scp', logname, scp_user + ':' + web_logs_stub + instrsitedate + '_log.log'])
    
    # Close the connection to airglow database
    con.close()

    # Return log as string if an email to humans is suggested.
    s = ''
    if notify_the_humans:
        with open (logname, "r") as myfile:
            s = myfile.read() # return the whole thing as a string, including new-lines
    return s



def process_site(site_name, year, doy, reference='laser'):
    '''
    Process all the data from all the instruments at the site with name
    site_name on the date specified with (year, doy). This function just
    calls process_instr(...) once per instrument, and ignores the output
    of process_instr(...).
    INPUTS:
        site_name - 'uao', 'eku', etc.
        year - int
        doy - int (1-indexed, right?)
    OPTIONAL INPUTS:
        reference - str, 'laser' or 'zenith'. Passed on to FPI.ParameterFit(...)
    HISTORY:
        19 Jul 2013 bjh - new version for remote2
    '''
    process_dn = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1)
    instr_names = fpiinfo.get_instr_at(site_name, process_dn)
    if len(instr_names) == 0:
        raise Exception('No instruments found at site "%s" on %s' % (site_name, process_dn))
    for instr_name in instr_names:
        print 'Starting (%s, %s, %s, %s, reference=%s)' % (instr_name, site_name, year, doy, reference)
        process_instr(instr_name, year, doy, reference)
    


def multiprocess_instr(arg_list, num_processes=16):
    '''
    INPUT: arg_list - a list of [instr, year, doy] that will be cycled through.
                      Optional extra inputs: [ ... , reference, use_npz] 
                      (see process_instr function for details)
           num_process - the maximum number of processes to work 
                         (recommend <= 16 for remote2).
    '''
    # Create a list to store the running threads
    threads = []

    # Loop to perform the analysis
    while threads or arg_list:
        # Check if we should start a new process
        if len(threads) < num_processes and arg_list:
            # Grab the argument list to operate on
            args = arg_list.pop()
            print 'Starting %s' % str(args)
            # Set up the thread
            p = multiprocessing.Process(target=process_instr, args=args)
            p.start()
            #print args
            threads.append(p)
        else:
            # No more free processes.  Check for one to end
            for thread in threads:
                if not thread.is_alive():
                    threads.remove(thread)

def process_zenith_ref(num_processes=16):
    '''
    Search through all the log files for the "<BADLASER>" string. Re-run the nights
    with this tag using reference='zenith'. Also run any nights hard-coded below.
    OPTIONAL INPUT:
       num_process - the maximum number of processes to work 
                         (recommend <= 16 for remote2).
    HISTORY:
       29 July 2013 bjh - created
    '''

    #################### HARD CODE BAD LASER DAYS HERE #######################
    # If there are bad laser days that are not automatically caught by this
    # function, hard code them here, in the form:
    # [instr_name, year, doy, 'zenith']. All sequences written here will be run.
    # One way to obtain this list is with the show_error_logs() function below.
    arg_list = []
    # For example:
    # arg_list.append(['minime02', 2013, 8, 'zenith'])
    # arg_list.extend([['minime05', 2013, i, 'zenith'] for i in range(203,209)])
    #
    # Add the days when the Morocco FPI didn't have a laser.
    arg_list.extend([['minime03', 2013, i, 'zenith'] for i in range(311,351)])
    # The following were obtained from FPIprocess.show_error_logs().
    # Most of them are just bad data, but we'll try anyway for completeness.
    arg_list.extend([
        ['minime02', 2011, 301, 'zenith'],
        ['minime02', 2011, 302, 'zenith'],
        ['minime02', 2011, 303, 'zenith'],
        ['minime02', 2011, 304, 'zenith'],
        ['minime02', 2011, 305, 'zenith'],
        ['minime02', 2011, 307, 'zenith'],
        ['minime02', 2011, 310, 'zenith'],
        ['minime02', 2011, 311, 'zenith'],
        ['minime02', 2011, 312, 'zenith'],
        ['minime02', 2011, 313, 'zenith'],
        ['minime02', 2011, 314, 'zenith'],
        ['minime02', 2011, 315, 'zenith'],
        ['minime02', 2011, 316, 'zenith'],
        ['minime02', 2011, 322, 'zenith'],
        ['minime02', 2011, 325, 'zenith'],
        ['minime02', 2011, 339, 'zenith'],
        ['minime02', 2011, 340, 'zenith'],
        ['minime02', 2011, 342, 'zenith'],
        ['minime02', 2011, 349, 'zenith'],
        ['minime02', 2012,  81, 'zenith'],
        ['minime02', 2012,  90, 'zenith'],
        ['minime02', 2012,  95, 'zenith'],
        ['minime02', 2012,  97, 'zenith'],
        ['minime02', 2012, 110, 'zenith'],
        ['minime02', 2012, 111, 'zenith'],
        ['minime02', 2012, 112, 'zenith'],
        ['minime02', 2012, 113, 'zenith'],
        ['minime02', 2012, 114, 'zenith'],
        ['minime02', 2012, 115, 'zenith'],
        ['minime02', 2012, 116, 'zenith'],
        ['minime02', 2012, 117, 'zenith'],
        ['minime02', 2012, 118, 'zenith'],
        ['minime02', 2012, 119, 'zenith'],
        ['minime02', 2012, 120, 'zenith'],
        ['minime02', 2012, 121, 'zenith'],
        ['minime02', 2012, 122, 'zenith'],
        ['minime02', 2012, 123, 'zenith'],
        ['minime02', 2012, 124, 'zenith'],
        ['minime02', 2012, 125, 'zenith'],
        ['minime02', 2012, 126, 'zenith'],
        ['minime02', 2012, 127, 'zenith'],
        ['minime02', 2012, 128, 'zenith'],
        ['minime02', 2012, 129, 'zenith'],
        ['minime02', 2012, 130, 'zenith'],
        ['minime02', 2012, 131, 'zenith'],
        ['minime02', 2012, 132, 'zenith'],
        ['minime02', 2012, 133, 'zenith'],
        ['minime02', 2012, 134, 'zenith'],
        ['minime02', 2012, 135, 'zenith'],
        ['minime02', 2012, 139, 'zenith'],
        ['minime02', 2012, 140, 'zenith'],
        ['minime02', 2012, 141, 'zenith'],
        ['minime02', 2012, 142, 'zenith'],
        ['minime02', 2012, 143, 'zenith'],
        ['minime02', 2012, 145, 'zenith'],
        ['minime02', 2012, 146, 'zenith'],
        ['minime02', 2012, 147, 'zenith'],
        ['minime02', 2012, 148, 'zenith'],
        ['minime02', 2012, 149, 'zenith'],
        ['minime02', 2012, 150, 'zenith'],
        ['minime02', 2012, 151, 'zenith'],
        ['minime02', 2012, 173, 'zenith'],
        ['minime02', 2012, 212, 'zenith'],
        ['minime02', 2012, 213, 'zenith'],
        ['minime02', 2012, 214, 'zenith'],
        ['minime02', 2012, 215, 'zenith'],
        ['minime02', 2012, 216, 'zenith'],
        ['minime02', 2012, 292, 'zenith'],
        ['minime02', 2012, 304, 'zenith'],
        ['minime02', 2012, 366, 'zenith'],
        ['minime02', 2013,   1, 'zenith'],
        ['minime02', 2013,   2, 'zenith'],
        ['minime02', 2013,   3, 'zenith'],
        ['minime02', 2013,   4, 'zenith'],
        ['minime02', 2013,   5, 'zenith'],
        ['minime02', 2013,   6, 'zenith'],
        ['minime02', 2013,   7, 'zenith'],
        ['minime02', 2013,   8, 'zenith'],
        ['minime02', 2013,   9, 'zenith'],
        ['minime02', 2013,  93, 'zenith'],
        ['minime02', 2013, 215, 'zenith'],
        ['minime02', 2007, 297, 'zenith'],
        ['minime06', 2011, 172, 'zenith'],
        ['minime06', 2011, 281, 'zenith'],
        ['minime06', 2012, 154, 'zenith'],
        ['minime06', 2012, 362, 'zenith'],
        ['minime06', 2013, 101, 'zenith'],
        ['minime07', 2012, 234, 'zenith'],
        ['minime07', 2012, 286, 'zenith'],
        ['minime07', 2012, 287, 'zenith'],
        ['minime07', 2012, 300, 'zenith'],
        ['minime07', 2012, 301, 'zenith'],
        ['minime07', 2012, 302, 'zenith'],
        ['minime07', 2012, 303, 'zenith'],
        ['minime07', 2012, 308, 'zenith'],
        ['minime07', 2012, 322, 'zenith'],
        ['minime07', 2012, 325, 'zenith'],
        ['minime07', 2012, 334, 'zenith'],
        ['minime07', 2013, 107, 'zenith'],
        ['minime07', 2013, 108, 'zenith'],
        ['minime07', 2013, 109, 'zenith'],
        ['minime07', 2013, 111, 'zenith'],
        ['minime07', 2013, 212, 'zenith'],
        ]
        )

    ##########################################################################

    
    logfile_stub = '/rdata/airglow/fpi/results/' # location of log files

    # Search for logfiles with the string <BADLASER>
    key = '<BADLASER>'
    logfiles = glob.glob(logfile_stub + '*.log')
    logfiles_to_rerun = []
    for log in logfiles:
         with open(log) as f:
            if key in f.read():
                logfiles_to_rerun.append(log)
    logfiles_to_rerun.sort()

    
    # From these filenames, find the instruments, years, and doys to re-run
    for log in logfiles_to_rerun:
        # Parse the filename to get instr_name, year, and doy
        (instr_name, site_name, datestr) = log.split('/')[-1].split('.')[0].split('_')
        year  = int(datestr[0:4])
        month = int(datestr[4:6])
        day   = int(datestr[6:8])
        dt = datetime.datetime(year, month, day)
        doy = ((dt - datetime.datetime(year, 1, 1)) + datetime.timedelta(days=1)).days
        
        arg_list.append([instr_name, year, doy, 'zenith'])

   
    # Re-run the analysis for each (instr_name, year, doy) using reference='zenith'
    print 'Found %i days to re-process.' % len(arg_list)
    multiprocess_instr(arg_list, num_processes)

def process_zenith_ref_no_multi():
    '''
    First, run all the nights hardcoded below using reference='zenith'.
    Then, search through all the log files for the "<BADLASER>" string. Re-run the nights
    with this tag using reference='zenith'.
    HISTORY:
       29 July 2013 bjh - created
       11 Feb  2014 bjh - modified for single-core processing
    '''

    #################### HARD CODE BAD LASER DAYS HERE #######################
    # If there are bad laser days that are not automatically caught by this
    # function, hard code them here, in the form:
    # [instr_name, year, doy, 'zenith']. All sequences written here will be run.
    # One way to obtain this list is with the show_error_logs() function below.
    arg_list = []
    # For example:
    # arg_list.append(['minime02', 2013, 8, 'zenith'])
    # arg_list.extend([['minime05', 2013, i, 'zenith'] for i in range(203,209)])
    #
    # Add the days when the Morocco FPI didn't have a laser.
    arg_list.extend([['minime03', 2013, i, 'zenith'] for i in range(311,351)])
    # The following were obtained from FPIprocess.show_error_logs().
    # Most of them are just bad data, but we'll try anyway for completeness.
    '''
    arg_list.extend([
        ['minime02', 2011, 301, 'zenith'],
        ['minime02', 2011, 302, 'zenith'],
        ['minime02', 2011, 303, 'zenith'],
        ['minime02', 2011, 304, 'zenith'],
        ['minime02', 2011, 305, 'zenith'],
        ['minime02', 2011, 307, 'zenith'],
        ['minime02', 2011, 310, 'zenith'],
        ['minime02', 2011, 311, 'zenith'],
        ['minime02', 2011, 312, 'zenith'],
        ['minime02', 2011, 313, 'zenith'],
        ['minime02', 2011, 314, 'zenith'],
        ['minime02', 2011, 315, 'zenith'],
        ['minime02', 2011, 316, 'zenith'],
        ['minime02', 2011, 322, 'zenith'],
        ['minime02', 2011, 325, 'zenith'],
        ['minime02', 2011, 339, 'zenith'],
        ['minime02', 2011, 340, 'zenith'],
        ['minime02', 2011, 342, 'zenith'],
        ['minime02', 2011, 349, 'zenith'],
        ['minime02', 2012,  81, 'zenith'],
        ['minime02', 2012,  90, 'zenith'],
        ['minime02', 2012,  95, 'zenith'],
        ['minime02', 2012,  97, 'zenith'],
        ['minime02', 2012, 110, 'zenith'],
        ['minime02', 2012, 111, 'zenith'],
        ['minime02', 2012, 112, 'zenith'],
        ['minime02', 2012, 113, 'zenith'],
        ['minime02', 2012, 114, 'zenith'],
        ['minime02', 2012, 115, 'zenith'],
        ['minime02', 2012, 116, 'zenith'],
        ['minime02', 2012, 117, 'zenith'],
        ['minime02', 2012, 118, 'zenith'],
        ['minime02', 2012, 119, 'zenith'],
        ['minime02', 2012, 120, 'zenith'],
        ['minime02', 2012, 121, 'zenith'],
        ['minime02', 2012, 122, 'zenith'],
        ['minime02', 2012, 123, 'zenith'],
        ['minime02', 2012, 124, 'zenith'],
        ['minime02', 2012, 125, 'zenith'],
        ['minime02', 2012, 126, 'zenith'],
        ['minime02', 2012, 127, 'zenith'],
        ['minime02', 2012, 128, 'zenith'],
        ['minime02', 2012, 129, 'zenith'],
        ['minime02', 2012, 130, 'zenith'],
        ['minime02', 2012, 131, 'zenith'],
        ['minime02', 2012, 132, 'zenith'],
        ['minime02', 2012, 133, 'zenith'],
        ['minime02', 2012, 134, 'zenith'],
        ['minime02', 2012, 135, 'zenith'],
        ['minime02', 2012, 139, 'zenith'],
        ['minime02', 2012, 140, 'zenith'],
        ['minime02', 2012, 141, 'zenith'],
        ['minime02', 2012, 142, 'zenith'],
        ['minime02', 2012, 143, 'zenith'],
        ['minime02', 2012, 145, 'zenith'],
        ['minime02', 2012, 146, 'zenith'],
        ['minime02', 2012, 147, 'zenith'],
        ['minime02', 2012, 148, 'zenith'],
        ['minime02', 2012, 149, 'zenith'],
        ['minime02', 2012, 150, 'zenith'],
        ['minime02', 2012, 151, 'zenith'],
        ['minime02', 2012, 173, 'zenith'],
        ['minime02', 2012, 212, 'zenith'],
        ['minime02', 2012, 213, 'zenith'],
        ['minime02', 2012, 214, 'zenith'],
        ['minime02', 2012, 215, 'zenith'],
        ['minime02', 2012, 216, 'zenith'],
        ['minime02', 2012, 292, 'zenith'],
        ['minime02', 2012, 304, 'zenith'],
        ['minime02', 2012, 366, 'zenith'],
        ['minime02', 2013,   1, 'zenith'],
        ['minime02', 2013,   2, 'zenith'],
        ['minime02', 2013,   3, 'zenith'],
        ['minime02', 2013,   4, 'zenith'],
        ['minime02', 2013,   5, 'zenith'],
        ['minime02', 2013,   6, 'zenith'],
        ['minime02', 2013,   7, 'zenith'],
        ['minime02', 2013,   8, 'zenith'],
        ['minime02', 2013,   9, 'zenith'],
        ['minime02', 2013,  93, 'zenith'],
        ['minime02', 2013, 215, 'zenith'],
        ['minime02', 2007, 297, 'zenith'],
        ['minime06', 2011, 172, 'zenith'],
        ['minime06', 2011, 281, 'zenith'],
        ['minime06', 2012, 154, 'zenith'],
        ['minime06', 2012, 362, 'zenith'],
        ['minime06', 2013, 101, 'zenith'],
        ['minime07', 2012, 234, 'zenith'],
        ['minime07', 2012, 286, 'zenith'],
        ['minime07', 2012, 287, 'zenith'],
        ['minime07', 2012, 300, 'zenith'],
        ['minime07', 2012, 301, 'zenith'],
        ['minime07', 2012, 302, 'zenith'],
        ['minime07', 2012, 303, 'zenith'],
        ['minime07', 2012, 308, 'zenith'],
        ['minime07', 2012, 322, 'zenith'],
        ['minime07', 2012, 325, 'zenith'],
        ['minime07', 2012, 334, 'zenith'],
        ['minime07', 2013, 107, 'zenith'],
        ['minime07', 2013, 108, 'zenith'],
        ['minime07', 2013, 109, 'zenith'],
        ['minime07', 2013, 111, 'zenith'],
        ['minime07', 2013, 212, 'zenith'],
        ]
        )
    '''
    ##########################################################################

    
    logfile_stub = '/rdata/airglow/fpi/results/' # location of log files

    # Search for logfiles with the string <BADLASER>
    key = '<BADLASER>'
    logfiles = glob.glob(logfile_stub + '*.log')
    logfiles_to_rerun = []
    count = -1
    for log in logfiles:
        with open(log) as f:
            count += 1
            if np.mod(count,100)==0:
                print '%i/%i' % (count,len(logfiles))
            if key in f.read():
                logfiles_to_rerun.append(log)
    logfiles_to_rerun.sort()

    
    # From these filenames, find the instruments, years, and doys to re-run
    for log in logfiles_to_rerun:
        # Parse the filename to get instr_name, year, and doy
        (instr_name, site_name, datestr) = log.split('/')[-1].split('.')[0].split('_')
        year  = int(datestr[0:4])
        month = int(datestr[4:6])
        day   = int(datestr[6:8])
        dt = datetime.datetime(year, month, day)
        doy = ((dt - datetime.datetime(year, 1, 1)) + datetime.timedelta(days=1)).days
        
        arg_list.append([instr_name, year, doy, 'zenith'])

   
    # Re-run the analysis for each (instr_name, year, doy) using reference='zenith'
    print 'Found %i days to re-process.' % len(arg_list)

    for arg in arg_list:
        print 'Starting %s_%i_%i' % (arg[0], arg[1], arg[2])
        try:
            process_instr(arg[0], arg[1], arg[2], arg[3])
        except Exception as e:
            print 'Failed: %s' % e


def show_error_logs():
    '''
    Print the last few lines of logs from all days where a log file was 
    created but no npz file was generated.
    Don't show logs if they contain the strings specified below.
    The idea is to catch bugs in the processing code to make sure
    that we are not overlooking good data.
    Also, at the end, print a list of args that can be easily copied to
    the process_zenith_ref function.
    '''

    keys_to_ignore = ['<BADLASER>', 'data found between']
    N = 20 # how many lines to show from end of log file
    
    logfile_stub = '/rdata/airglow/fpi/results/' # location of log files

    logfiles = glob.glob(logfile_stub + '*.log')
    npzfiles = glob.glob(logfile_stub + '*.npz')
    logfiles.sort()
    
    argstr = '[\n'
    for log in logfiles:
        # Parse filename to get year and doy
        (instr_name, site_name, datestr) = log.split('/')[-1].split('.')[0].split('_')
        year  = int(datestr[0:4])
        month = int(datestr[4:6])
        day   = int(datestr[6:8])
        dt = datetime.datetime(year, month, day)
        doy = ((dt - datetime.datetime(year, 1, 1)) + datetime.timedelta(days=1)).days
        # See if the npz file exists
        # Get npz filename
        npzfn = log[:-4] + '.npz'
        if npzfn not in npzfiles:
            with open(log) as f:
                loglines = f.readlines()
                logstr = ''.join(loglines)
                show = True
                for key in keys_to_ignore:
                    if key in logstr:
                        show = False
                        break
                if show:
                    print('===============================================')
                    print log + '\n'
                    M = np.min([N, len(loglines)])
                    for line in loglines[-M:]:
                        print line[:-1] # Don't print newline
                    print('===============================================')
                    argstr += '[\'%s\', %04s, %03s, \'zenith\'],\n' % (instr_name, year, doy)
    argstr += ']'
    print argstr
    return argstr

def find_unprocessed_days( years = range(2005,2020) , doys = range(1,367) ):
    '''
    Print and return all days that do not have corresponding log files.
    See also: show_error_logs().
    TODO: make this a bit more informative by actually crawling through 
    our database instead of through all instruments, sites, and days. Right
    now, data outside the proper directory structure will be ignored, as will
    sites and instruments not in fpiinfo (or spelled wrong)
    '''
    datastub = '/rdata/airglow/fpi'
    unproc = []
    for site_name in fpiinfo.get_all_sites_info().keys():
        for instr_name in fpiinfo.get_all_instr_names():
            for year in years:
                for doy in doys:
                    # Look in the directory for sky images
                    dt = datetime.datetime(year,1,1) + datetime.timedelta(days=doy-1)
                    datestr = dt.strftime('%Y%m%d')
                    path = '%s/%s/%s/%i/%s' % (datastub, instr_name, site_name, year, datestr)
                    files = glob.glob('%s/*X*.img' % path)
                    if files:
                        # Read the first file and find the appropriate "processing" doy
                        d = FPI.ReadIMG(files[0])
                        t = d.info['LocalTime'] - datetime.timedelta(hours=12) # find doy at beginning of night
                        procyear = t.year
                        procdoy = (t - datetime.datetime(t.year,1,1)).days + 1
                        procdatestr = t.strftime('%Y%m%d')
                        # Search for the corresponding log file
                        logfile = '%s/%s/%s_%s_%s.log' % (datastub, 'results', instr_name, site_name, procdatestr)
                        exists = os.path.isfile(logfile)
                        if not exists:
                            unproc.append((instr_name, procyear, procdoy, path, files[0]))
                        #print site_name,year,datestr,exists
    print unproc
    return unproc
                        
                        


def load_level0(instr_name, year, doy):
    '''
    Return the contents of the npz file specified by the arguments.
    Access contents with, e.g.,
        output['FPI_Results']
        output['instrument']
        output['site']
    '''
    # Load FPI_Results
    fpi_results_dir =  '/rdata/airglow/fpi/results/'
    process_dn = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1)
    site_name = fpiinfo.get_site_of(instr_name, process_dn)
    datestr = process_dn.strftime('%Y%m%d')
    instrsitedate = instr_name + '_' + site_name + '_' + datestr
    npzfn = fpi_results_dir + instrsitedate + '.npz'
    npzfile = np.load(npzfn)
    
    # Rearrange so it actually looks like a dictionary
    npzdict = {}
    for key in npzfile:
        npzdict[key] = npzfile[key].item()
   
    del npzfile.f 
    npzfile.close()
    return npzdict

















