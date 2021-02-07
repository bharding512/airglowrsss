
# Functions to process FPI data on remote2

import matplotlib
matplotlib.use('AGG')

import FPI
import glob
from lmfit import Parameters
from optparse import OptionParser
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
        dref,drefe = FPI.DopplerReference(FPI_Results, reference='laser')
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

    # Correct for Morocco timezone problems during Ramadan.
    # Daylight savings time paused during Ramadan, but the computer didn't know this.
    # After 2015, we switched the computer to UTC, so this won't happen again.
    if instr_name == 'minime03':
        t = FPI_Results['sky_times']
        tutc = np.array([ti.astimezone(pytz.utc).replace(tzinfo=None) for ti in t])

        # 2015
        # sky images
        idx = (tutc > datetime.datetime(2015,6,13,1)) & (tutc < datetime.datetime(2015,7,18,2))
        t[idx] = t[idx] - datetime.timedelta(hours=1)
        dst_secs = np.array([ti.dst().total_seconds() for ti in t])
        idx2 = (tutc >= datetime.datetime(2015,7,18,2)) & (tutc < datetime.datetime(2015,7,18,3)) & (dst_secs == 0.)
        t[idx2] = t[idx2] - datetime.timedelta(hours=1)
        if sum(idx) > 0:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                      'FPIprocess.quality_hack(): Manual intervention to fix timing offset during ' + \
                      '2015 Ramadan DST pause.\n')
        # laser images
        t = FPI_Results['laser_times']
        tutc = np.array([ti.astimezone(pytz.utc).replace(tzinfo=None) for ti in t])
        idx = (tutc > datetime.datetime(2015,6,13,1)) & (tutc < datetime.datetime(2015,7,18,2))
        t[idx] = t[idx] - datetime.timedelta(hours=1)
        dst_secs = np.array([ti.dst().total_seconds() for ti in t])
        idx2 = (tutc >= datetime.datetime(2015,7,18,2)) & (tutc < datetime.datetime(2015,7,18,3)) & (dst_secs == 0.)
        t[idx2] = t[idx2] - datetime.timedelta(hours=1)


        # 2014
        # sky images
        idx = (tutc > datetime.datetime(2014,6,28,1)) & (tutc < datetime.datetime(2014,8,2,2))
        t[idx] = t[idx] - datetime.timedelta(hours=1)
        dst_secs = np.array([ti.dst().total_seconds() for ti in t])
        idx2 = (tutc >= datetime.datetime(2014,8,2,2)) & (tutc < datetime.datetime(2014,8,2,3)) & (dst_secs == 0.)
        t[idx2] = t[idx2] - datetime.timedelta(hours=1)
        if sum(idx) > 0:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                      'FPIprocess.quality_hack(): Manual intervention to fix timing offset during ' + \
                      '2014 Ramadan DST pause.\n')
        # laser images
        t = FPI_Results['laser_times']
        tutc = np.array([ti.astimezone(pytz.utc).replace(tzinfo=None) for ti in t])
        idx = (tutc > datetime.datetime(2014,6,28,1)) & (tutc < datetime.datetime(2014,8,2,2))
        t[idx] = t[idx] - datetime.timedelta(hours=1)
        dst_secs = np.array([ti.dst().total_seconds() for ti in t])
        idx2 = (tutc >= datetime.datetime(2014,8,2,2)) & (tutc < datetime.datetime(2014,8,2,3)) & (dst_secs == 0.)
        t[idx2] = t[idx2] - datetime.timedelta(hours=1)

    return FPI_Results





def createL1ASCII(NPZ,OUT):
    '''
    Summary:
        Function to read in processed FPI npz file and output a text file of level 1 results.
        Level 1 results are line-of-sight winds and temperatures for each look direction.

    Inputs:
        NPZ = full file name and path of npz file
        OUT = full file name and path of ASCII results output

    History:
        2/21/13 -- Written by Daniel J. Fisher (dfisher2@illinois.edu)
        2/16/15 -- Modified for complete npz file
    '''

    # Load in FPI processed Data
    data = np.load(NPZ,allow_pickle=True)
    FPI_Results = data['FPI_Results'].ravel()[0]
    reference = data['FPI_Results'].ravel()[0]['reference']
    #instr = data['instrument'].ravel()[0]['Abbreviation']
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
        sky_temp = np.ones(len(timeywimey))*-999.
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
        note.write('UTCTime____________  ____Az  ___Ze  ______T  _T_Sig  ___Wind  _W_Sig  FitSig  CalSig  _____I  _I_Sig  ___Bkgd B_Sig  Sec  Chisqr  CldInd TF WF  __Ref  ______Wl  Vers\n')

        for a_tw, a_az, a_ze, a_t, a_e_t, a_w, a_e_w, a_ef, a_ec, a_i, a_e_i, a_b, a_e_b, a_it, a_cs, a_dt,a_tf,a_wf in zip(timeywimey, az, ze, temps, e_temps, winds, e_winds, e_fit, e_cal, i, e_i, b, e_b, inttime, chisqr, sky_temp,temp_flag,wind_flag):
            dn = a_tw.astimezone(pytz.utc)
            utctime = dn.strftime("%Y-%m-%d %H:%M:%S")
            line = "{:19s}  {:6.1f}  {:5.1f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:6.2f}  {:6.2f}  {:6.4f}  {:6.4f}  -999.00  -999  {:3.0f}  {:6.2f}  {:6.1f}  {:1.0f}  {:1.0f}  {:6s}  {:7.1e}  {:5s}\n".format(utctime, a_az, a_ze, a_t, a_e_t, a_w, a_e_w, a_ef, a_ec, a_i, a_e_i, a_it, a_cs, a_dt, a_tf, a_wf, reference, wavelength, version)
            #line = "{:19s}  {:6.1f}  {:5.1f}  {:7.2f}  {:6.2f}  {:7.2f}  {:6.2f}  {:6.1f}  {:4.2f}  {:6.1f}  {:4.2f}  {:3.0f}  {:6.2f}  {:6.1f}  {:1d}  {:1d}  {:6s}  {:7.1e}  {:5s}\n".format(utctime, a_az, a_ze, a_t, a_e_t, a_w, aax.xaxis_date()_e_w, a_i, a_e_i, a_b, a_e_b, a_it, a_cs, a_dt, t_flag, w_flag, reference, wavelength, version)
            note.write(line)

    note.closed




def get_all_laser_images(direc):
    '''
    Return all laser images in the specified directory, as a list of strings.
    Laser images are those of the following forms:

    L20091103001.img
    UAO_L_20091103_210000_001.img

    Return empty list if none are found.
    '''

    fns_1 = glob.glob(direc + '/L[0-9]*.img')
    fns_2 = glob.glob(direc + '/*_L_*.img')

    return fns_1 + fns_2



def get_all_sky_images(direc):
    '''
    Return all sky (i.e., airglow) images in the specified directory, as a list of strings.
    Sky images are those of the following forms:

    X20091103001.img
    UAO_X_20091103_210000_001.img

    Return empty list if none are found.
    '''

    fns_1 = glob.glob(direc + '/X[0-9]*.img')
    fns_2 = glob.glob(direc + '/*_X_*.img')

    return fns_1 + fns_2




def process_instr(instr_name ,year, doy, reference='laser', use_npz = False,
                  wind_err_thresh=100., temp_err_thresh=100., cloud_thresh = [-22.,-10.],
                  send_to_website=False, enable_share=False, send_to_madrigal=False,
                  enable_windfield_estimate=False,
                  fpi_dir='/rdata/airglow/fpi/', bw_dir='/rdata/airglow/templogs/cloudsensor/',
                  x300_dir='/rdata/airglow/templogs/x300/', results_stub='/rdata/airglow/fpi/results/',
):
    '''
    Process all the data from the instrument with name instr_name on
    the day specified by (year, doy). This function looks in the appropriate
    directory, analyzes the laser and sky images to obtain line of
    sight winds, temperatures, etc., and saves them in a .npz file in the
    appropriate place. It also generates summary plots and sends them to the
    website. Finally, it generates summary ASCII files for Madrigal.
    INPUTS:
        instr_name - 'minime01', 'minime02', etc.
        year - int
        doy - int (i.e., Jan 1 is doy=1)
    OPTIONAL INPUTS:
        reference - str, 'laser' or 'zenith'. Passed on to FPI.ParameterFit(...)
        use_npz - bool, if True, the call to FPI.ParameterFit(...) will be skipped,
                  and the previously-analyzed data in the saved npz file will be used
                  instead. Cloud and X300 data will be reloaded and added to the npz file.
                  New wind and temperature plots will be sent to the website,
                  along with the updated log file.
        wind_err_thresh - float, m/s. Samples with a fit error above this should get a quality
                            flag of 2.
        temp_err_thresh - float, K. Samples with a fit error above this should get a quality
                            flag of 2.
        cloud_thresh - [float,float], K. The two cloud sensor readings that indicate
                        partially- and fully-cloudy. This affects the quality flag.
        send_to_website - bool. If True, connect to the database at Illinois and update
                          the website. This will only work if run at Illinois.
        enable_share - bool. If True, and if site['share']==True, then save a copy of the npz
                       file elsewhere on remote2 for transfer to other institutions.
        send_to_madrigal - bool. If True, and if instrument['send_to_madrigal']==True,
                       then create an ASCII file of the results and save it to the appropriate
                       place at Illinois for Madrigal to read. This only works at Illinois.
        enable_windfield_estimate - bool. If True, and if this site is part of the 'nation'
                       or 'peru' networks, then the wind field estimation analysis will be
                       run. Right now, this only works at Illinois.
        fpi_dir - str. the base directory where the data are stored. It is assumed that
                  data are organized in a folder like: <fpi_dir>/minime01/car/2015/20150810/
        bw_dir - str. the directory where the Boltwood Cloud Sensor files are stored. Set
                      this to empty string ('') if no sensor exists
        x300_dir - str. the directory where the X300 temperature sensor files are stored.
                       Set this to empty string ('') if no sensor exists
        results_stub - str. the directory where the results will be saved. A .npz file
                       and a .log file will be created.

    OUTPUTS:
        warnings - str - If this script believes a manual check of the data
                   is a good idea, a message will be returned in this string.
                   If not, the empty string will be returned. For now, the message
                   is the entire log.

    '''


    # Define constants that do not depend on the site
    direc_tol = 10.0 # tolerance in degrees to recognize a look direction with
    ccd_temp_thresh = -60. # sky exposures with a CCD temp above this will get a quality flag.

    # Information about sending data to the website. Only used if send_to_website==True.
    temp_plots_stub= '/rdata/airglow/fpi/results/temporary_plots/' #where to save png files
    scp_user       = 'airglowgroup@webhost.engr.illinois.edu'
    db_image_stub  =        'SummaryImages/' # relative path from web server directory on airglow
    db_log_stub =           'SummaryLogs/'
    web_images_stub =       '/home/airglowgroup/data/SummaryImages/' # absolute location on airglow of summary images
    web_logs_stub =         '/home/airglowgroup/data/SummaryLogs/' # absolute location on airglow of summary logs
    # Information about sending to Madrigal. Only used if send_to_madrigal==True
    madrigal_stub =         '/rdata/airglow/database/' # where reduced ascii txt and png files are saved for Madrigal
    # Information about sending to partner institutions. Only used if enable_share==True
    share_stub    =         '/rdata/airglow/share/' # where to save a copy of the .npz file for sharing with collaborators

    notify_the_humans = False # default

    nominal_dt = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1)

    # Import the site information
    site_name = fpiinfo.get_site_of(instr_name, nominal_dt)
    site = fpiinfo.get_site_info(site_name, nominal_dt)
    # Import the instrument information
    instrument = fpiinfo.get_instr_info(instr_name, nominal_dt)

    # Create "minime05_uao_20130729" string
    datestr = nominal_dt.strftime('%Y%m%d')
    instrsitedate = instr_name + '_' + site_name + '_' + datestr

    # Construct the directories to the relevant data products
    data_stub      = fpi_dir + instr_name + '/' + site_name + '/'
    bw_dir_stub    = bw_dir + site_name + '/'
    bw_name_stub   = 'Cloud_' + site_name + '_'
    x300_dir_stub  = x300_dir + site_name + '/'
    x300_name_stub = 'TempL_' + site_name + '_'

    # Determine the times between which files will be accepted.
    # Define start and stop times as solar noon on doy, and on doy+1
    site_lon = np.mod(site['Location'][1]+180,360)-180
    start_dt = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1, hours = 12-24*site_lon/360.)
    stop_dt = start_dt + datetime.timedelta(hours=24)
    start_dt = pytz.utc.localize(start_dt)
    stop_dt = pytz.utc.localize(stop_dt)

    # Create a list of relevant filenames by searching in adjacent
    # days' directories also, since we might have problems with UT/LT,
    # time zones, etc.
    laser_fns = []
    sky_fns = []
    local = pytz.timezone(site['Timezone'])

    if not use_npz:
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
            fns_all = get_all_laser_images(data_dir)
            for fn in fns_all:
                d = FPI.ReadIMG(fn)
                dtime = local.localize(d.info['LocalTime'])
                if dtime > start_dt and dtime < stop_dt:
                    laser_fns.append(fn)
            # Second, sky:
            fns_all = get_all_sky_images(data_dir)
            for fn in fns_all:
                d = FPI.ReadIMG(fn)
                dtime = local.localize(d.info['LocalTime'])
                if dtime > start_dt and dtime < stop_dt:
                    sky_fns.append(fn)
        laser_fns.sort()
        sky_fns.sort()

    # Uncomment these for rapid testing of new code
    #laser_fns = laser_fns[::6]
    #sky_fns = [sky_fns[20],sky_fns[78],sky_fns[119]]

    if not laser_fns and not sky_fns and not use_npz:
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
    Diagnostic_Fig = plt.figure(dpi=300, figsize=(10,7.5)) # Figure for diagnostics to be drawn to
    if use_npz: # Load previously-analyzed results
        npzfile = np.load(npzname,allow_pickle=True)
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
                                logfile=logfile, diagnostic_fig=Diagnostic_Fig, reference=reference)
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


    # Grab the SVN revision number, so we know what code was used to process this day.
    svndir = '/'.join(FPI.__file__.split('/')[:-1]) # get the directory of FPI.py
    p = subprocess.Popen('svnversion %s'%svndir,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdout,stderr) = p.communicate()
    sv = re.split(':|\n', stdout)[0]

    # Save the SVN version number
    FPI_Results['SVNRevision'] = sv

    # Try to load the Boltwood and X300 data. This is inside a try block so that
    # the npz will still save if there is a problem.
    print( np.any(np.isfinite(FPI_Results['Clouds']['mean'])) )
    if not use_npz:# Do not overwrite
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
        print(np.any(np.isfinite(Clouds['mean'])))
        if(np.size(bw_date) == 0):
            if not use_npz:#Do not overwrite
                FPI_Results['Clouds'] = None
        else:
            FPI_Results['Clouds'] = {'mean': Clouds[:,0], 'max': Clouds[:,1], 'min': Clouds[:,2]}

        if(np.size(x_date) == 0):
            if not use_npz:#Do not overwrite
                FPI_Results['Dome'] = None
                FPI_Results['Inside'] = None
        else:
            FPI_Results['Dome'] = {'mean': Dome[:,0], 'max': Dome[:,1], 'min': Dome[:,2]}
            FPI_Results['Inside'] = {'mean': Inside[:,0], 'max': Inside[:,1], 'min': Inside[:,2]}

        # Write things to the log
        if FPI_Results['Clouds'] is None: # if it wasn't found
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No Boltwood cloud sensor data found.\n')
            c = np.nan*np.zeros(len(FPI_Results['LOSwind']))
            if not use_npz:#Do not overwrite
                FPI_Results['Clouds'] = {'mean': c, 'max': c, 'min': c}
            if bw_dir: # it should have been found
                notify_the_humans = True
        else:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Found and loaded Boltwood cloud sensor data.\n')
        if FPI_Results['Dome'] is None:# if it wasn't found but should have been
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No X300 temperature sensor data found.\n')
            #notify_the_humans = True # No need to notify about X300's.
        else:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Found and loaded X300 temperature sensor data.\n')
    except: # There was an error in the Boltwood or X300 code. Write the log but continue.
        c = np.nan*np.zeros(len(FPI_Results['LOSwind']))
        if not use_npz:#Do not overwrite
            FPI_Results['Clouds'] = {'mean': c, 'max': c, 'min': c}
        tracebackstr = traceback.format_exc()
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Unknown error obtaining Boltwood or X300 data for %s. Traceback listed below. Analysis will continue without these data.\n-----------------------------------\n%s\n-----------------------------------\n' % (datestr,tracebackstr))
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
        # The order here is important. It goes from least concerning to most
        # concerning.
        if (FPI_Results['skyI'][ii] < instrument['skyI_quality_thresh'][0]):
            # The sky brightness is low enough that OH may be an issue.
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
        if FPI_Results['sky_ccd_temperature'][ii] > ccd_temp_thresh:
            t_flag = 1
            w_flag = 1
            logfile.write('[ccdtemp>%.0f: W1T1] '%ccd_temp_thresh)
        if c > cloud_thresh[1]: # It's definitely cloudy
            t_flag = 1
            w_flag = 2
            logfile.write('[cloud>%.0f: W2T1] '%cloud_thresh[1])
        if (FPI_Results['skyI'][ii] < instrument['skyI_quality_thresh'][1]):
            # The sky brightness is low enough that OH is almost certainly issue.
            t_flag = 2
            w_flag = 2
            logfile.write('[skyI low W2T2] ')
        w_fit_err = FPI_Results['sigma_fit_LOSwind'][ii]
        if (w_fit_err > wind_err_thresh): # Sample is not trustworthy at all
            w_flag = 2
            t_flag = 2
            logfile.write('[sigma_wind large: W2T2] ')
        t_fit_err = FPI_Results['sigma_T'][ii]
        if (t_fit_err > temp_err_thresh): # Sample is not trustworthy at all
            t_flag = 2
            w_flag = 2
            logfile.write('[sigma_T large: W2T2] ')
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
    if enable_share and site['share']: # save a copy of the npz file in a separate folder
        npznameshare = share_stub + site_name + '/' + instrsitedate + '.npz'
        np.savez(npznameshare, FPI_Results=FPI_Results, site=site, instrument=instrument)
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Results also saved to %s\n' % npznameshare)
    if send_to_madrigal and instrument['send_to_madrigal']: # save the summary ASCII file to send to the Madrigal database
        asciiname = madrigal_stub + instrsitedate + '.txt'
        createL1ASCII(npzname, asciiname)
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'ASCII results saved to %s\n' % asciiname)

    # Load the results
    npzfile = np.load(npzname,allow_pickle=True)
    FPI_Results = npzfile['FPI_Results']
    FPI_Results = FPI_Results.reshape(-1)[0]
    site = npzfile['site']
    site = site.reshape(-1)[0]
    del npzfile.f # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
    npzfile.close()

    # Try to make plots
    try:
        # Plot some quick-look single-station data (LOSwinds and Temps)
        if use_npz:
            Diagnostic_Fig=FPIDisplay.PlotDiagnosticDay(npzname, cloud_thresh, instrument['skyI_quality_thresh'], LASERPARAMINTERP='linear')
        else:
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
                ax.set_ylabel('Cloud indicator\n[degrees C]')
                ax.grid(True)
                ax.set_xlabel('Universal Time, [hours]')
            Diagnostic_Fig.tight_layout()

        (Temperature_Fig, Temperature_Graph), (Doppler_Fig, Doppler_Graph) = FPIDisplay.PlotDay(npzname, reference = reference)
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

    gif_fn = None # filename of the gif to send to the website
    if enable_windfield_estimate:
        # Create (or update) the level 3 wind field map, if use_npz=False (so that use_npz=True is fast)
        network_name = site['Network']

        if network_name in  ['nation','peru'] and not use_npz:
            try:
                import FPIwindfield
                wf = FPIwindfield.WindField([network_name], year, doy, timestep_fit = 1800)
                if len(wf.instr_names) >= 3: # If there are 3 sites, try the inversion
                    losfig = wf.plot_los_winds()
                    plt.close(losfig)
                    wf.run_inversion()
                    if wf.success:
                        gif_fn = wf.make_quiver_gif(show_vert_wind=True)
                        shutil.copy(gif_fn, '/home/bhardin2/public_html/windfield_movies/')
                        quicklook_fig, quicklook_fn = wf.make_quicklook(show_vert_wind=True)
                        shutil.copy(quicklook_fn, '/home/bhardin2/public_html/windfield_quicklooks/')
                        plt.close(quicklook_fig)
                        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Wind field estimation successful. Created gif and quicklook plot.\n')
                    else:
                        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Wind field estimation not successful. No plots created.\n')
            except:
                tracebackstr = traceback.format_exc()
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error creating windfield quicklook plot. Traceback listed below.\n-----------------------------------\n%s\n-----------------------------------\n' % (tracebackstr,))
                notify_the_humans = True

    # Save files in temporary folder.
    summary_figs = [Diagnostic_Fig,
                    Temperature_Fig,
                    Doppler_Fig,]
    summary_fns = [instrsitedate + '_diagnostics.png',
                   instrsitedate + '_temperature.png', # e.g., 'minime05_uao_20130107_temperature.png'
                   instrsitedate + '_winds.png',]
    db_ids = [instrument['sql_diagnostics_id'],
              instrument['sql_temperatures_id'],
              instrument['sql_winds_id'], ]
    #if use_npz: # don't update diagnostics fig
    #    summary_figs.pop(0)
    #    summary_fns.pop(0)

    for fig, fn in zip(summary_figs, summary_fns):
        fig.savefig(temp_plots_stub + fn, bbox_inches='tight') # save it in the remote2 data dir
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Created %s.\n' % fn)
        plt.close(fig)

    if send_to_website:
        import MySQLdb as mdb
	from sshtunnel import SSHTunnelForwarder

	def query(sql_cmd):
        	with SSHTunnelForwarder(
                	('webhost.engr.illinois.edu', 22),
                	ssh_username='airglowgroup',
                	ssh_private_key='/home/airglow/.ssh/id_rsa',
                	remote_bind_address=('127.0.0.1', 3306)
        	) as server:
                	con = mdb.connect(host='127.0.0.1', db='airglowgroup_webdatabase', port=server.local_bind_port, read_default_file="/home/airglow/.my.cnf")
                	cur = con.cursor()
                	cur.execute(sql_cmd)
	                rows = cur.fetchall()
        	        return rows

        site_id = site['sql_id']
        utc = pytz.utc # Define timezones
        # Start and stop time of observations.
        #d = FPI.ReadIMG(sky_fns[0])
        #dtime = local.localize(d.info['LocalTime'])
        #startut = dtime.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')
        #d = FPI.ReadIMG(sky_fns[-1])
        #dtime = local.localize(d.info['LocalTime'])
        #stoput = dtime.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')
        # Open the database (see http://zetcode.com/databases/mysqlpythontutorial/ for help)
        # Read the user and password from a file.
#        con = mdb.connect(host='webhost.engr.illinois.edu', db='airglowgroup_webdatabase', read_default_file="/home/airglow/.my.cnf")
#        cur = con.cursor()
        # Send summary images to server
        for fn, db_id in zip(summary_fns, db_ids):
            flag = subprocess.call(['scp', temp_plots_stub + fn, scp_user + ':' + web_images_stub + fn]) # send to airglow
            if flag != 0: # Sending png to airglow failed
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error sending %s to airglow server for displaying on website.\n' % fn)
                notify_the_humans = True
            # If we are supposed to put this png in Madrigal, put it in the Madrigal database directory
            if send_to_madrigal and instrument['send_to_madrigal']:
                shutil.copy(temp_plots_stub + fn ,madrigal_stub) # save the figure to the madrigal directory
            # update the database
            # Send png. First find out if the entry is in there (i.e., we are just updating the png)
            sql_cmd = 'SELECT id FROM DataSet WHERE SummaryImage = \"%s\"' % (db_image_stub + fn)
            #sql_cmd = 'SELECT id FROM DataSet WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (site_id, db_id, startut)
	    rows = query(sql_cmd)
#            cur.execute(sql_cmd)
#            rows = cur.fetchall()
            log_fn = db_log_stub + instrsitedate + '_log.log'
            if len(rows) == 0: # Create the entry
                # Start and stop time of observations
                if use_npz:
                    d0,dn=FPI_Results['sky_times'][0],FPI_Results['sky_times'][-1]
                else:
                    d=FPI.ReadIMG(sky_fns[0])
                    d0 = local.localize(d.info['LocalTime'])
                    d = FPI.ReadIMG(sky_fns[-1])
                    dn = local.localize(d.info['LocalTime'])
                
                startut = d0.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')
                stoput = dn.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')
                
                sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime, SummaryImage, LogFile) VALUES(%d, %d, \"%s\", \"%s\", \"%s\", \"%s\")' % (site_id, db_id, startut, stoput, db_image_stub + fn, log_fn)
		logfile.write(sql_cmd + '\n')
		query(sql_cmd)
#                cur.execute(sql_cmd)
            else: # Entry exists. Update it.
                sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\",LogFile=\"%s\" WHERE id = %d' % (db_image_stub + fn, log_fn, rows[0][0])
               # sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\",LogFile=\"%s\" WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (db_image_stub + fn, log_fn, site_id, db_id, startut)
#                cur.execute(sql_cmd)
		logfile.write(sql_cmd + '\n')
		query(sql_cmd)

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
#            cur.execute(sql_cmd)
#            rows = cur.fetchall()
	    row = query(sql_cmd)
            if len(rows) == 0: # Create the entry
	    #if True:
                sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime, SummaryImage) VALUES(%d, %d, \"%s\", \"%s\", \"%s\")' % (network_id, gif_id, startut, stoput, db_image_stub + gif_fn.split('/')[-1])
		logfile.write(sql_cmd + '\n')
		query(sql_cmd)
#                cur.execute(sql_cmd)
            else: # Entry exists. Update it.
                sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\" WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (db_image_stub + gif_fn.split('/')[-1], network_id, gif_id, startut)
		query(sql_cmd)
#                cur.execute(sql_cmd)



    logfile.close()
    if send_to_website:
        # Send the log file
        subprocess.call(['scp', logname, scp_user + ':' + web_logs_stub + instrsitedate + '_log.log'])
        # Close the connection to airglow database
#        con.close()

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
        doy - int (1-indexed)
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





def process_directory(data_direc, results_direc, instrument, site, reference='laser',
                      wind_err_thresh=100., temp_err_thresh=100., make_plots=True, Tmin=300., Tmax=1500.,
                      bw_cloud_file_1 = '', bw_cloud_file_2 = '', cloud_thresh = [-22.,-10.], horizon_cutoff=-6):
    '''
    Process all the data from the specified folder. This function analyzes
    the laser and sky images to obtain line of sight winds,
    temperatures, etc., and saves them in a .npz file in the
    directory results_direc. It also generates summary plots and saves them in
    the directory results_direc. It does *not* look for X300 or Boltwood Cloud
    Sensor data. If this is running on remote2, you should use process_instr
    instead.
    INPUTS:
        data_direc - directory containing the FPI data to be processed
        results_direc - directory to save the npz file and summary images
        instrument - dictionary containing necessary instrument parameters
        site - dictionary containing necessary site parameters
    OPTIONAL INPUTS:
        reference - str, 'laser' or 'zenith'. Passed on to FPI.ParameterFit(...)
        wind_err_thresh - float, m/s. Samples with a fit error above this should get a quality
                            flag of 2.
        temp_err_thresh - float, K. Samples with a fit error above this should get a quality
                            flag of 2.
        make_plots      - whether to make and save summary graphics
        Tmin, Tmax      - float, K. Sets the scale of the y axis on the temperature plot
        bw_cloud_file_1 - location of the Boltwood cloud sensor file relevant to the pre-midnight data
                          (if '', don't try to ingest cloud data at all)
        bw_cloud_file_2 - location of the Boltwood cloud sensor file relevant to the post-midnight data
                          (if '', ignore the second file)
        cloud_thresh    - [float,float], K. The two cloud sensor readings that indicate
                          partially- and fully-cloudy. This affects the quality flag.
        horizon_cutoff  - int, deg. For solar elevation angles above this, the sun is considered "up" and
                          data are discarded.

    OUTPUTS:
        warnings - str - If this script believes a manual check of the data
                   is a good idea, a message will be returned in this string.
                   If not, the empty string will be returned. For now, the message
                   is the entire log.

    '''


    # Define constants that do not depend on the site
    direc_tol = 10.0 # tolerance in degrees to recognize a look direction with
    ccd_temp_thresh = -60. # sky exposures with a CCD temp above this will get a quality flag.

    notify_the_humans = False # default

    # Find laser and sky files
    laser_fns = get_all_laser_images(data_direc)
    sky_fns   = get_all_sky_images(data_direc)
    local = pytz.timezone(site['Timezone'])
    laser_fns.sort()
    sky_fns.sort()


    if not laser_fns and not sky_fns:
        raise Exception('No data found.\n')

    # Open the first image to get the date
    d = FPI.ReadIMG(sky_fns[0])
    t0 = local.localize(d.info['LocalTime'])
    datestr = t0.strftime('%Y_%m_%d_%Hh')
    instrdatestr = instrument['name'] + '_' + datestr

    # Open a logfile
    logname = results_direc + instrdatestr + '.log'
    logfile = open(logname,'w') # overwrite previous log
    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Logfile Created\n')


    if not laser_fns and sky_fns and reference=='laser': # This is not a big deal if reference=='zenith'
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No %s laser data found between %s and %s. Sky data found. <BADLASER> \n' % (instr_name, str(start_dt), str(stop_dt)))
        logfile.close()
        raise Exception('No %s laser data found between %s and %s. Sky data found.\n' % (instr_name, str(start_dt), str(stop_dt)))


    npzname = results_direc + instrdatestr + '.npz' # the name of the npz file to save to
    Diagnostic_Fig = plt.figure(dpi=300, figsize=(10,7.5)) # Figure for diagnostics to be drawn to
    try:
        # Run the analysis
        (FPI_Results, notify_the_humans) = FPI.ParameterFit(instrument, site, laser_fns, sky_fns, \
                direc_tol=direc_tol, N=instrument['N'], N0=instrument['N0'], N1=instrument['N1'], \
                            logfile=logfile, diagnostic_fig=Diagnostic_Fig, reference=reference, horizon_cutoff=horizon_cutoff)
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
            'Laser and sky image analysis complete.\n')
    except:
        # FPI.ParameterFit crashed. For now, write the log and re-raise the Exception.
        tracebackstr = traceback.format_exc()
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error analyzing %s. Traceback listed below.\n-----------------------------------\n%s\n-----------------------------------\n' % (instrdatestr,tracebackstr))
        notify_the_humans = True # This is redundant, since raising will ensure humans are notified.
    #        Diagnostic_Fig = plt.figure()
    #        Diagnostic_Fig.text(0.5, 0.5, 'Error - see log file', fontsize = 20, ha='center')
        logfile.close()
        raise


    # Grab the SVN revision number, so we know what code was used to process this day.
    svndir = '/'.join(FPI.__file__.split('/')[:-1]) # get the directory of FPI.py
    p = subprocess.Popen('svnversion %s'%svndir,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdout,stderr) = p.communicate()
    sv = re.split(':|\n', stdout)[0]

    # Save the SVN version number
    FPI_Results['SVNRevision'] = sv


    # Try to load the Boltwood data. This is inside a try block so that
    # the npz will still save if there is a problem.
    FPI_Results['Clouds'] = None # defaults
    FPI_Results['Dome']   = None # assume no X300
    FPI_Results['Inside'] = None # assume no X300
    try:
        # call modules from BoltWood for the two days required and combine the data
        if bw_cloud_file_1 != '': # Look for boltwood data
            if bw_cloud_file_2 == '': # only use one file
                bw_date, bw_sky, bw_amb = BoltwoodSensor.ReadTempLog(
                        bw_cloud_file_1,
                        tz=site['Timezone']
                        )

            else: # use two files
                bw_date1, bw_sky1, bw_amb1 = BoltwoodSensor.ReadTempLog(
                        bw_cloud_file_1,
                        tz=site['Timezone'])
                bw_date2, bw_sky2, bw_amb2 = BoltwoodSensor.ReadTempLog(
                        bw_cloud_file_2,
                        tz=site['Timezone']
                        )
                bw_date = np.hstack((bw_date1,bw_date2))
                bw_sky = np.hstack((bw_sky1, bw_sky2))
                bw_amb = np.hstack((bw_amb1, bw_amb2))


            # Create a data structure containing the sensor temperatures
            Clouds = np.nan*np.zeros((np.size(FPI_Results['sky_times']),3))
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


                count = count+1

            if(np.size(bw_date) == 0):
                FPI_Results['Clouds'] = None
            else:
                FPI_Results['Clouds'] = {'mean': Clouds[:,0], 'max': Clouds[:,1], 'min': Clouds[:,2]}


        # Write things to the log
        if FPI_Results['Clouds'] is None: # if it wasn't found
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No Boltwood cloud sensor data found.\n')
            c = np.nan*np.zeros(len(FPI_Results['LOSwind']))
            FPI_Results['Clouds'] = {'mean': c, 'max': c, 'min': c}

        else:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Found and loaded Boltwood cloud sensor data.\n')

    except: # There was an error in the Boltwood code. Write the log but continue.
        c = np.nan*np.zeros(len(FPI_Results['LOSwind']))
        FPI_Results['Clouds'] = {'mean': c, 'max': c, 'min': c}
        tracebackstr = traceback.format_exc()
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Unknown error obtaining Boltwood data for %s. Traceback listed below. Analysis will continue without these data.\n-----------------------------------\n%s\n-----------------------------------\n' % (datestr,tracebackstr))
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
    laser_drift = laser_is_drifting()
    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + ' Quality Flag Calculation:\n')
    logfile.write('\t\tLaser Drift: %s\n' % laser_drift)
    logfile.write('\t\tSample-specific flags:\n')
    wind_quality_flag = np.zeros(len(FPI_Results['sky_times']))
    temp_quality_flag = np.zeros(len(FPI_Results['sky_times']))
    for ii in range(len(FPI_Results['sky_times'])): # Loop over samples
        logfile.write('\t\t%s ' % FPI_Results['sky_fns'][ii].split('/')[-1])
        t_flag = 0 # default
        w_flag = 0 # default
        if (FPI_Results['skyI'][ii] < instrument['skyI_quality_thresh'][0]):
            # The sky brightness is low enough that OH may be an issue.
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
        if FPI_Results['sky_ccd_temperature'][ii] > ccd_temp_thresh:
            t_flag = 1
            w_flag = 1
            logfile.write('[ccdtemp>%.0f: W1T1] '%ccd_temp_thresh)
        if c > cloud_thresh[1]: # It's definitely cloudy
            t_flag = 1
            w_flag = 2
            logfile.write('[cloud>%.0f: W2T1] '%cloud_thresh[1])
        if (FPI_Results['skyI'][ii] < instrument['skyI_quality_thresh'][1]):
            # The sky brightness is low enough that OH is almost certainly issue.
            t_flag = 2
            w_flag = 2
            logfile.write('[skyI low W2T2] ')
        w_fit_err = FPI_Results['sigma_fit_LOSwind'][ii]
        if (w_fit_err > wind_err_thresh): # Sample is not trustworthy at all
            w_flag = 2
            t_flag = 2
            logfile.write('[sigma_wind large: W2T2] ')
        t_fit_err = FPI_Results['sigma_T'][ii]
        if (t_fit_err > temp_err_thresh): # Sample is not trustworthy at all
            t_flag = 2
            w_flag = 2
            logfile.write('[sigma_T large: W2T2] ')
        # A manual override can increase the calculated flag, but not decrease
        #t_flag = max(t_flag,t_flag_manual)
        #w_flag = max(w_flag,w_flag_manual)
        wind_quality_flag[ii] = w_flag
        temp_quality_flag[ii] = t_flag
        logfile.write('[FINAL: W%iT%i]\n' % (w_flag,t_flag))
    FPI_Results['wind_quality_flag'] = wind_quality_flag
    FPI_Results['temp_quality_flag'] = temp_quality_flag


    # Save the results
    np.savez(npzname, FPI_Results=FPI_Results, site=site, instrument=instrument)
    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Results saved to %s\n' % npzname)

    # Load the results
    npzfile = np.load(npzname,allow_pickle=True)
    FPI_Results = npzfile['FPI_Results']
    FPI_Results = FPI_Results.reshape(-1)[0]
    site = npzfile['site']
    site = site.reshape(-1)[0]
    del npzfile.f # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
    npzfile.close()

    # Try to make plots
    if make_plots:
        try:
            # Plot some quick-look single-station data (LOSwinds and Temps)
            (Temperature_Fig, Temperature_Graph), (Doppler_Fig, Doppler_Graph) = FPIDisplay.PlotDay(npzname, reference = reference, Tmin=Tmin, Tmax=Tmax)

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
                ylim_kwargs={'ymax':0}
                if np.nanmedian(cloud)>-50:
                    ylim_kwargs['ymin']=-50
                ax.set_ylim(**ylim_kwargs)
                ax.set_ylabel('Cloud indicator\n[degrees C]')
                ax.grid(True)
                ax.set_xlabel('Universal Time, [hours]')
            Diagnostic_Fig.tight_layout()

            Temperature_Fig.savefig(results_direc + instrdatestr + '_temp.png')
            Doppler_Fig.savefig(    results_direc + instrdatestr + '_wind.png')
            Diagnostic_Fig.savefig( results_direc + instrdatestr + '_diagnostic.png')
            plt.close(Doppler_Fig)
            plt.close(Temperature_Fig)
            plt.close(Diagnostic_Fig)

        except:
            # Summary plots crashed. We still want to send the diagnostics and log to the
            # website, so continue on with dummy plots.
            tracebackstr = traceback.format_exc()
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Error making plots for %s. Traceback listed below. Dummy plots sent to website.\n-----------------------------------\n%s\n-----------------------------------\n' % (instrdatestr,tracebackstr))
            notify_the_humans = True
            Temperature_Fig = plt.figure()
            Temperature_Fig.text(0.5, 0.5, 'Plotting error - see log file', fontsize = 20, ha='center')
            Doppler_Fig = plt.figure()
            Doppler_Fig.text(0.5, 0.5, 'Plotting error - see log file', fontsize = 20, ha='center')

    logfile.close()

    # Return log as string if an email to humans is suggested.
    s = ''
    if notify_the_humans:
        with open (logname, "r") as myfile:
            s = myfile.read() # return the whole thing as a string, including new-lines
    return s






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
            kwargs={}
            i_argdict=[i for i,arg in enumerate(args) if isinstance(arg,dict)]
            if len(i_argdict)>0:
                kwargs=args.pop(i_argdict[0])
            # Set up the thread
            p = multiprocessing.Process(target=process_instr, args=args, kwargs=kwargs)
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
    # If instrument not at a site at this time, return IOError, file not found, for consistency
    # with the error that occurs if the instrument didn't take data on that night
    if site_name is None:
        raise IOError('File Not Found for %s_%s_%03i' % (instr_name, year, doy))
    datestr = process_dn.strftime('%Y%m%d')
    instrsitedate = instr_name + '_' + site_name + '_' + datestr
    npzfn = fpi_results_dir + instrsitedate + '.npz'
    npzfile = np.load(npzfn,allow_pickle=True)

    # Rearrange so it actually looks like a dictionary
    npzdict = {}
    for key in npzfile:
        npzdict[key] = npzfile[key].item()

    del npzfile.f
    npzfile.close()
    return npzdict

















