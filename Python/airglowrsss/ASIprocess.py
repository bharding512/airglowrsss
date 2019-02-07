import matplotlib
matplotlib.use('AGG')

import time
from mpl_toolkits.basemap import Basemap
import ASIDisplay
import ASI
import glob as glob
import datetime as datetime
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from PIL import Image
import MySQLdb as mdb
import os
import pytz
import traceback
import multiprocessing
import asiinfo
import TifImagePlugin
#import pdb

def process_instr(inst,year,doy,do_DB=True):
    # Function to process a single day of ASI (all-sky imaging) data.  A
    # movie and a keogram will be generated.  If requested, the information
    # will be added into the airglow database
    #
    # INPUTS:
    #   inst - instrument code to process
    #   year, doy - year and day of year of date to be processed
    #
    # OPTIONAL INPUTS:
    #   do_DB - put information into the airglow SQL database (default=True)
    #
    # OUTPUTS:
    #   warnings - warning message if processing crashed (a traceback)
    #              If no issues, an empty string will be returned. 
    #              TODO: Mimic FPIprocess warnings
    #
    # HISTORY:
    #   Written by Jonathan J. Makela on 15 July 2013
    #   Reworked to use instruments as the input, rather than site on 22 July 2013

    # Create the process_dn date
    process_dn = datetime.datetime(year,1,1) + datetime.timedelta(days = doy -1)
    warnings = ''

    # Create needed string variables
    datestr = process_dn.strftime('%Y%m%d')
    
    # User to scp to
    scp_user = 'airglowgroup@webhost.engr.illinois.edu'
    web_stub = 'SummaryImages/'
    file_stub = '/home/airglowgroup/data/SummaryImages/'

    inst = inst.lower()

    # Find out where the instrument was
    site = asiinfo.get_site_of(inst,process_dn)

    # Get the info for the instrument on this date
    inst_info = asiinfo.get_instr_info(inst,process_dn)
    site_info = asiinfo.get_site_info(site)
    filters = inst_info['filters']
    filter_names = inst_info['filter_names']
    unwarp_ht = inst_info['unwarp_ht']
    inst_id = inst_info['sql_inst_id']
    horz = inst_info['horizon']
    t_lat = inst_info['t_lat']
    t_lon = inst_info['t_lon']
    nfile = inst_info['cal_file']
    kernel_size = inst_info['kernel_size']
    site_id = site_info['sql_id']
    site_name = site_info['Name']
    show_countries = site_info['borders']
    dark_flag = inst_info['ignore_dark']
    
    # Load in the data from the requested calibration file
    npzfile = np.load(nfile)
    el = npzfile['el']
    az = npzfile['az']
    rx_lat = npzfile['rx_lat']
    rx_lon = npzfile['rx_lon']
    if 'nf' in inst: #get CNFI's mask that used to be on el/az to apply to data
        maskDat = npzfile['mask']
    else:
        maskDat = None
    npzfile.close()

    # Determine the times between which files will be accepted.
    # Use local-time noon on doy to local-time noon on doy+1
    local = pytz.timezone(site_info['Timezone'])
    start_dt = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1, hours = 12)
    stop_dt  = datetime.datetime(year,1,1) + datetime.timedelta(days = doy-1, hours = 36)
    start_dt = local.localize(start_dt)
    stop_dt  = local.localize(stop_dt)
    
    # Location of files
    data_stub = '/rdata/airglow/imaging/' + inst + '/' + site + '/'
    warn_flag = True

    # Loop through each requested filter
    for (fils,fil_name, un_ht, in_id, target_lat, target_lon) in zip(filters, filter_names, unwarp_ht, inst_id, t_lat, t_lon):
    
        # Variable to hold the files that were created
        files_created = []
       
        try:
            con = None

            # Create a list of relevant filenames by searching in adjacent
            # days' directories also, since we might have problems with UT/LT,
            # time zones, etc.
            files = []
            temps = []
            darks = []
            for day_offset in [-1, 0, 1]: # search folders around the day of interest
                dir_dt = process_dn + datetime.timedelta(days = day_offset)
                # Create the YYYYMMDD date format
                yearstr = dir_dt.strftime('%Y')
                doystr = dir_dt.strftime('%j')
                # Create the directory name for the data to be processed
                data_dir = data_stub + yearstr + '/' + doystr + '/'
                # Look in the data directory, and grab the files if they were
                # taken between the start and stop times.

                #fns = glob.glob(data_dir + site.upper() + '_' + fils + '*.tif')
                #dks = glob.glob(data_dir + site.upper() + '_[dD][aA][rR][kK]_*.tif')
                fns = glob.glob(data_dir + '*_' + fils + '*.tif')
                dks = glob.glob(data_dir + '*_[dD][aA][rR][kK]_*.tif')
                if len(fns) == 0:
                    # Do old data format if nothing is returned.
                    fns = glob.glob(data_dir + fils + '*.tif')
                    dks = glob.glob(data_dir + '[dD][aA][rR][kK]_*.tif')
                for fn in fns:
                    d = Image.open(fn)
                    dtime = local.localize(d.info['LocalTime'])
                    if dtime > start_dt and dtime < stop_dt:
                        files.append(fn)
                        temps.append(d.info['CCDTemperature'])
                for dk in dks:
                    d = Image.open(dk)
                    dtime = local.localize(d.info['LocalTime'])
                    if dtime > start_dt and dtime < stop_dt:
                        darks.append(dk) 
            files.sort()
            darks.sort()
            if len(darks) is 0 or dark_flag:
                darks = None

        # Go through images for each filter
            if len(files) is not 0 :
                files.sort()
                warn_flag = False

                # Warning if CCD is too hot
                cool_point = asiinfo.get_instr_info(inst,process_dn)['ccd_temp_set']
                over_limit = sum(np.array(temps) >= cool_point+3)
                if over_limit > 0:
                    warnings = warnings + 'CCD OVERHEATING: %s- %03i of %03i images over set limit\n'%(fils,over_limit,len(files))

                # Output names
                movie_name = inst + '_' + site + '_' + datestr + '_' + fils + 'movie.avi'
                keo_name = inst + '_' + site + '_' + datestr + '_' + fils + 'keogram.png'

                # Convert from the el/az to lat/lon for the requested unwarp height
                lat,lon = ASI.ConvertAzEl2LatLon(az,el,un_ht,rx_lat,rx_lon,horizon=horz)
                min_lat = np.min(lat[np.isfinite(lat)])
                max_lat = np.max(lat[np.isfinite(lat)])
                min_lon = np.min(lon[np.isfinite(lon)])
                max_lon = np.max(lon[np.isfinite(lon)])
                
                # Create the map
                m = Basemap(llcrnrlon=min_lon-0.5,llcrnrlat=min_lat-0.5,urcrnrlon=max_lon+0.5,urcrnrlat=max_lat+0.5,
                            projection='merc', area_thresh=1000,
                            resolution='i')

                # create mask to mask image data
                if maskDat is not None:
                    mask = np.isnan(lat) | maskDat
                else:
                    mask = np.isnan(lat)
                
                # Get lat and lon without nans so we can plot
                lat1,lon1 = ASI.ConvertAzEl2LatLon(az,el,un_ht,rx_lat,rx_lon,horizon=-10.0) #-10 so there are no nans
                
                # Create the movie
                ASI.MapMovie(files,m,lat1,lon1,mask,movie_name=movie_name,darks=darks,sitename=site_name,kernel_size=kernel_size, displayCountries=show_countries,filt=fil_name)

                # SCP the movie over to airglow (call was Popen)
                subprocess.call(['scp',movie_name,scp_user + ':' + file_stub])

                files_created.append(movie_name)

                # Create the keogram
                #pdb.set_trace()
                f = ASIDisplay.Keogram(files,lat,lon,target_lat,target_lon,sitename=site_name,filt=fil_name,darks=darks)
                f.savefig(keo_name)

                subprocess.call(['scp', keo_name, scp_user + ':' + file_stub])
        
                files_created.append(keo_name)

                # Add info to the database, if requested
                if do_DB:
                    # Start and stop time of observations
                    d = Image.open(files[0])
                    startut = d.info['UniversalTime']
                    d = Image.open(files[-1])
                    stoput = d.info['UniversalTime']

                    # Open and populate the database (see http://zetcode.com/databases/mysqlpythontutorial/ for help)
                    con = mdb.connect(host='webhost.engr.illinois.edu',db='airglowgroup_webdatabase',read_default_file='/home/airglow/.my.cnf')
                    cur = con.cursor()

                    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
                    sql_cmd = 'SELECT id FROM DataSet WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (site_id, in_id, startut)
                    cur.execute(sql_cmd)
                    rows = cur.fetchall()
                    if len(rows) == 0:
                        sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime,SummaryImage, SummaryMovie) VALUES(%d, %d, \"%s\", \"%s\", \"%s\", \"%s\")' % (site_id, in_id, startut, stoput, web_stub + keo_name, web_stub + movie_name)
                        cur.execute(sql_cmd)
                    else:
                        # Entry exists.  Update it
                        sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\", SummaryMovie=\"%s\" WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (web_stub + keo_name, web_stub + movie_name, site_id, in_id, startut)
                        cur.execute(sql_cmd)
            else:
                print 'No files used for '+fils+' on this day.'
            
        except:
            print 'Something bad happenend...'
            warnings = warnings + traceback.format_exc() + '\n'
            
        finally:
            if con:
                con.close()

            # Delete the files that were created
            for f in files_created:
                os.remove(f)
    
    if warn_flag:
        warnings = warnings + 'No files in folders!\n'
                
    return(warnings)
    

def ASI_multiprocess(arg_list, num_processes):
    # INPUT: arg_list - a list of [site, year, doy] that will be cycled through.
    #        num_process - the maximum number of processes to work.

    # Create a list to store the running threads
    threads = []

    # Loop to perform the analysis
    while threads or arg_list:
        # Check if we should start a new process
        if len(threads) < num_processes:
            # Grab the argument list to operate on
            args = arg_list.pop()

            # Set up the thread
            p = multiprocessing.Process(target=process_instr, args=args)
            p.start()
            print args
            threads.append(p)
        else:
            # No more free processes.  Check for one to end
            for thread in threads:
                if not thread.is_alive():
                    threads.remove(thread)
