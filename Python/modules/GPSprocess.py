import matplotlib
matplotlib.use('AGG')

#from mpl_toolkits.basemap import Basemap
import ScintMonDisplay
import glob as glob
import datetime as datetime
import pytz
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import Image
import MySQLdb as mdb
import os
import re
import tarfile
import multiprocessing
import Emailer
import gpsinfo
import CASES as cases

def process_instr(inst,year,doy,do_DB=True):
    # Function to process a single day of GPS data.  A
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
    # HISTORY:
    #   Written by Daniel J. Fisher on 05 Feb 2014

    # Create the process_dn date
    process_dn = datetime.datetime(year,1,1) + datetime.timedelta(days = doy -1)

    instr = inst[0:-2]

    # Find out where the instrument was
    site = gpsinfo.get_site_of(inst,process_dn)

    # Get the info for the instrument on this date
    inst_info = gpsinfo.get_instr_info(inst,process_dn)
    site_info = gpsinfo.get_site_info(site)
    inst_id = inst_info['sql_inst_id']
    horz = inst_info['horizon']
    site_id = site_info['sql_id']
    site_name = site_info['Name']
    timezone = site_info['Timezone']

    # Create needed string variables
    year = process_dn.strftime('%Y')
    month = process_dn.strftime('%m')
    day = process_dn.strftime('%d')
    date = process_dn.strftime('%y%m%d')

    # Create the data directory name for the data to be processed
    data_dir = '/rdata/airglow/gps/' + inst + '/' + site + '/' + year + '/'
    os.chdir(data_dir)

    # Do Scintmon S4 Processing
    if instr == 'scintmo':
        # Move data to data folder
        raw_dir  = data_dir + 'raw_data/'
        os.system('mv ' +raw_dir+date + '*.[f,o,n]* ' + data_dir)
        files = glob.glob(date + '*.fsl')
        files.sort()
        # Go through all s4 data per day
        for f in files:
            s4filename = inst+f[-5]+'_'+site+'_'+year+month+day+'_s4.png'
            os.system('chmod 777 ' + data_dir+f[0:-4]+'*')
            print 'Getting lsum4...'
            os.system('/usr/local/bin/lsum4 -n ' + f[0:-4])
            print 'Move raw data...'
            os.system('chmod 750 '+data_dir+f[0:-4]+'*')

            sumfile = glob.glob(f[0:-4]+'.sum')
            if(sumfile):
                # change ownership
                os.system('chmod 750 '+sumfile[0]) 
                os.system('chown airglow.gps '+sumfile[0])
                # get output image and times
                startut, stoput, s4fig = ScintMonDisplay.PlotDay(sumfile[0],horz)
                s4fig.savefig(data_dir+s4filename)
                # Add info to the database, if requested
                if do_DB:
                    print 'Do databasing...'
                    database([s4filename],data_dir,startut,stoput,[inst_id],site_id)
                os.system('rm -f ' + data_dir+s4filename)
            else:
                subject = "!!! Processing Error (\'" + inst +'\','+year+','+str(doy)+') @ ' + site 
                print subject
                Emailer.emailerror(subject, 'lsum4 problem - move rawdata back and reprocess?')


    # Do Scinda S4/TEC Processing
    elif instr == 'scinda':
        try:
            print 'Running wintep-p...'
            os.system('/usr/local/bin/wintec-p-daily.pl ' + date+'*.scn')
            print 'Move raw data...'
            # Standardize filenames
            s4filename = inst+'_'+site+'_'+year+month+day+'_s4.png'
            tecfilename = inst+'_'+site+'_'+year+month+day+'_tec.png'
            os.rename(data_dir+'plots/S'+date+'.PNG',data_dir+s4filename)
            os.rename(data_dir+'plots/T'+date+'.PNG',data_dir+tecfilename)
            # Get times
            files = glob.glob(data_dir+date+'*.scn')
            files.sort()
            local = pytz.timezone(timezone)
            # find start time
            data = open(files[0],'r')
            r = re.compile('[ \t\n\r:]+')
            line = r.split(data.readline())
            dtime = datetime.datetime(2000+int(line[1]),int(line[2]),int(line[3]))+datetime.timedelta(seconds=line[4])
            dtime = local.localize(dtime)
            startut = dtime.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')
            data.close()
            # find end time
            r = re.compile('[ \t\n\r:]+')
            for zelda in reversed(open(files[-1]).readlines()):
		line = r.split(zelda)
                if line[0] == 'T':
                    dtime = datetime.datetime(2000+int(line[1]),int(line[2]),int(line[3]))+datetime.timedelta(seconds=line[4])
                    dtime = local.localize(dtime)
                    break
            # Add info to the database, if requested
            if do_DB:
                print 'Do databasing...'
                database([s4filename,tecfilename],data_dir,startut,stoput,inst_id,site_id)
            os.system('rm -f ' + data_dir+s4filename)
            os.system('rm -f ' + data_dir+tecfilename)
        except:
            subject = "!!! Processing Error (\'" + inst +'\','+year+','+str(doy)+') @ ' + site 
            print subject
            Emailer.emailerror(subject, 'wintec problem')


    # Do Cases Processing
    elif instr == 'cases':

        # NOTE: Python code (originally from /programs/plot_CASES_day.py
        try:
            # Paths to use
            name = inst + '_' + site
            datapath = '/rdata/airglow/gps/'+inst+'/'+site+'/streaming/'
            pngpath = '/rdata/airglow/gps/results/'
            log_fname = '/rdata/airglow/gps/results/'+inst+'_'+site+'.log'
            #log_cmd = '/usr/bin/perl /usr/local/share/Python/programs/load_'+inst+'_'+site+'_logfile.pl'

            # Create the filename of the data to be parsed
            fname = '{:s}dataout_{:s}_{:03d}.bin'.format(datapath,year,doy) 

            # Run binflate
            os.popen('/rdata/airglow/gps/cases01/hka/streaming/binflate -i ' + fname)

            # Load txinfo.log file
            txfname = 'txinfo.log'
            txprn,txdn, el, az, txsystem = cases.load_txinfo(txfname)

            # Load scint.log file
            s4fname = 'scint.log'
            s4prn, s4dn, s4, s4system = cases.load_scint(s4fname)
            startut = s4dn[0].strftime('%Y-%m-%d %H:%M:%S')
            stoput = s4dn[-1].strftime('%Y-%m-%d %H:%M:%S')

            # Create plots
            dn = datetime.date(int(year),1,1)+datetime.timedelta(days=doy-1)
            s4fname = '{:s}_{:s}H.png'.format(name, dn.strftime('%y%m%d'))
            cases.plot_s4summary(txprn,txdn,el,az,txsystem,s4prn,s4dn,s4,s4system,pngpath+s4fname)

            # Write the logfile
            fid = open(log_fname,'w')

            fid.writelines('Site,Instrument,StartUTTime,StopUTTime,SummaryImage,MovieFile\n')
            fid.writelines('{:d},{:d},{:s},{:s},{:s}_{:s}H.png'.format(site_id,inst_id,
                           startut, stoput, name, dn.strftime('%y%m%d')))
            fid.close()

            # Load the log file into the database
            #os.popen(log_cmd)
            if do_DB:
                database([s4fname],pngpath,startut,stoput,[inst_id],site_id)

            # Move the data
            os.popen('mv {:s}dataout*_{:s}_{:03d}.bin .'.format(datapath,year,doy
                ) )
            tar = tarfile.open('{:s}_{:s}.tgz'.format(name, dn.strftime('%y%m%d')), 'w:gz')

            tar.add('dataout_{:s}_{:03d}.bin'.format(year,doy))

            if os.path.exists('dataoutiq_{:s}_{:03d}.bin'.format(year,doy)):
                tar.add('dataoutiq_{:s}_{:03d}.bin'.format(year,doy))

            tar.close()

            # Clean up files
            #os.popen('tar czvf {:s}.tgz dataout*_{:s}_{:03d}.bin'.format(dn.strftime('%y%m%d'), year, DOY))
            os.popen('rm dataout*_{:s}_{:03d}.bin'.format(year,doy))
            #os.popen('mv {:s}.tgz {:s}{:s}'.format(dn.strftime('%y%m%d'),data_dir,year))
            os.popen('rm channel.log')
            os.popen('rm iono.log')
            os.popen('rm navsol.log')
            os.popen('rm scint.log')
            os.popen('rm txinfo.log')

        except:
            subject = "!!! Processing Error (\'" + inst +'\','+year+','+str(doy)+') @ ' + site
            print subject
            Emailer.emailerror(subject, 'Cases problem')


    # Send Error if Neither
    else:
        subject = "!!! GPS Problem " + date
        print subject
        Emailer.emailerror(subject, 'Something weird has happened: '+inst)

    os.chdir('/rdata/airglow/rx/')


def database(filenames,data_dir,startut,stoput,inst_id,site_id,send_to_madrigal=False):
    # Function to add files to the airglow database properly.
    # TODO: future madrigal support
    #
    # INPUTS:
    #   filenames - List of filenames to add to database
    #   data_dir - path to filenames
    #   startut - Starttime of gps data
    #   stoput - Endtime of gps data
    #   inst_id - sql id of the instruments used
    #   site_id - sql id of the location
    #
    # OPTIONAL:
    #   send_to_madrigal - True or Flase to send data to madrigal folder
    #
    # HISTORY:
    #   Written by Daniel J. Fisher on 27 Feb 2014

    # User to scp to
    scp_user = 'airglowgroup@webhost.engr.illinois.edu'
    web_stub = 'SummaryImages/'
    file_stub = '/home/airglowgroup/data/SummaryImages/'
    madgrigal_stub = '/rdata/airglow/database/'

    # Open the database (see http://zetcode.com/databases/mysqlpythontutorial/ for help)
    # Read the user and password from a file.
    con = mdb.connect(host='webhost.engr.illinois.edu', db='airglowgroup_webdatabase', read_default_file="/home/airglow/.my.cnf")
    cur = con.cursor()
    for fn,db_id in zip(filenames, inst_id):
        # Create the summary images
        subprocess.call(['scp', data_dir+fn, scp_user + ':' + file_stub + fn]) # send to airglow
        '''
        # If we are supposed to put this png in Madrigal, put it in the Madrigal database directory
        if send_to_madrigal:
            os.rename(data_dir+fn, madrigal_stub + fn) # send to madrigal
        '''
        # update the database
        # First find out if the entry is in there (i.e., we are just updating the png)
        sql_cmd = 'SELECT id FROM DataSet WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (site_id, db_id, startut)
        cur.execute(sql_cmd)
        rows = cur.fetchall()
        if len(rows) == 0: # Create the entry
            sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime, SummaryImage) VALUES(%d, %d, \"%s\", \"%s\", \"%s\")' % (site_id, db_id, startut, stoput, web_stub+fn)
            cur.execute(sql_cmd)
        else: # Entry exists. Update it.
            sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\" WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (web_stub+fn, site_id, db_id, startut)
            cur.execute(sql_cmd)
    con.close()


#def plot_cases_day(year,doy):


def GPS_multiprocess(arg_list, num_processes):
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
