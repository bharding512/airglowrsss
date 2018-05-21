#!/usr/bin/python
'''
Run program to unzip and sort data onto airglow's remote2 server

History: 25 Aug 2011 - initial script written (PERL)
         16 Jul 2013 - rewritten to PYTHON
         12 Feb 2014 - Updated to v3.0 - txtcheck

Written by Daniel J. Fisher (dfisher2@illinois.edu)
'''

# Import required modules
import os
import sys
import traceback
import tarfile
from glob import glob
import numpy as np
import datetime as dt
import pytz
from optparse import OptionParser
from collections import defaultdict
import Emailer
from Zipper import activeinstruments
import matplotlib
matplotlib.use('AGG')
import fpiinfo
import asiinfo
import gpsinfo
# Data Processing modules
import FPI
from PIL import Image
import FPIprocess
import ASIprocess
import GPSprocess


def instrumentcode():
    '''
    Summary:
        code = instrumentcode()
        Is the dictionary of all full instruments names given the 3 letter abbreviation.
            Contains pointers for cloud & templ data.
    
    Outputs:
        code =  dictionary of all full names from abbr

    History:
        2/25/14 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # Instrument Dictionary  TODO: make all 3 letters or use directory structure
    code = defaultdict(str)
    code['fpi'] = 'minime'
    code['asi'] = 'casi'
    code['nfi'] = 'cnfi'
    code['pic'] = 'picasso'
    code['sky'] = 'skymon'
    code['swe'] = 'swenson'
    code['cas'] = 'cases'
    code['tec'] = 'scinda'
    code['scn'] = 'scintmon'
    code['bwc'] = 'cloudsensor'
    code['x3t'] = 'x300'
    # Pointers for them too
    code['cloud'] = 'bwc'
    code['templ'] = 'x3t'

    return(code)


def sortinghat(dir_data,f):
    '''
    Summary
    -------
        result = sortinghat(dir_data,f,code):
        Concatinates/Unzips and moves data to folder. Returns errors and filenames.

    Inputs
    ------
        dir_data = directory that data will be placed in
        f = glob'd info file that contains name, parts, and size

    Outputs
    -------
        result = names - successful sort, else - error

    History
    -------
        7/18/13 -- Written by DJF (dfisher2@illionis.edu)
    '''
    
    result = 0
    code = instrumentcode()
    # Read info file
    info = open(f, 'r')
    zelda = info.readline().rstrip().split('.tar.gz',1)[0]+'.tar.gz'
    parts = int(info.readline())
    tsize = int(info.readline())
    time = dt.datetime.strptime(info.readline()[:19],'%Y-%m-%d %H:%M:%S')
    info.close()

    # Parse Name
    site = f[6:9]
    instr= f[0:3]
    inum = f[3:5]
    year = f[12:14]
    mon  = f[14:16]
    day  = f[16:18]
    emails = activeinstruments()[site][instr][inum]['email']

    # Check all parts present
    print 'Parts:',len(glob(zelda + '*')),'/',parts
    ## Case 1 - no data created last night
    if(tsize ==0 and parts == 0):
        ## Emails Warning that system is down!
        print '!!! No Data Collected'
        if instr in ['asi','nfi','pic','sky','swe']:
            msg = "%s%s down at %s!\nIs it a full moon?\nInternet & Sortinghat are working, is it an instrument/PC issue." %(code[instr],inum,site)
        else:
            msg = "%s%s down at %s!\nInternet & Sortinghat are working, this is an instrument/PC issue." %(code[instr],inum,site)
        subject = "!!! No data collected on %02s-%02s-%02s" %(mon,day,year)
        Emailer.emailerror(emails,subject,msg)
        # Move info file to tracking folder
        os.system('mv ' + f + ' ./tracking')
    ## Case 2 - all parts sent over in rx
    elif len(glob(zelda + '*')) == parts:
        # Check that folder to data exists
        try:
            os.makedirs(dir_data)
            os.system('chmod u+rwx,go+rX,go-w ' +dir_data)
            print "!!! Folder Created - Verify..."
        except OSError:
            print '!!! Folders Exist... moving on'
        # Concatinate the split files
        oscar = glob(zelda+'*')
        oscar.sort()
        print "\n".join(str(x) for x in oscar)
        os.system('cat ' + zelda + '* > temp.tar.gz')
        os.system('chmod 770 temp.tar.gz')
        # Check that size is correct
        statinfo = os.stat("temp.tar.gz")
        print 'Sizes:',statinfo.st_size,'/',tsize
        if statinfo.st_size == tsize:
            # Untar the gunzip files 
            tar = tarfile.open("temp.tar.gz","r:gz")
            try:
                result = tar.getnames()
                tar.extractall(dir_data)
                os.system('mv ' + f + ' ./tracking')
            except:
                age = (dt.datetime.utcnow()-time).total_seconds()/3600.0
                subject = "!!! Extract Error on %02s-%02s-%02s" %(mon,day,year)
                print subject
                msg = "%s%s issue at %s!\nThis file will not untar.\nBad Zip? Try -p %i" %(code[instr],inum,site,age/24)
                Emailer.emailerror(emails,subject,msg)
                result = []
            tar.close()
        else:
            print '!!! Waiting for complete parts...'
    ## Case 3 - all parts not yet sent
    else:
        print '!!! Waiting for parts...'
    os.system('rm -f temp.tar.gz')
    return(result)
    
    
def makeinfo(f):
    '''
    Summary
    -------
        makeinfo(f):
        Takes f (.tar.gz) and creates an appropriate txt info file.

    Inputs
    ------
        f = name of an unzipped tarball file (not standard Zip->Send->Sort)

    History
    -------
        2/13/14 -- Written by DJF (dfisher2@illinois.edu)
    '''
    
    now = dt.datetime.utcnow()
    statinfo = os.stat(f)
    f = f.lower()
    if len(f) == 25:
        check = open(f[0:-7]+'.txt', 'w')
    else:
        check = open(f[0:-4]+'.txt', 'w')
    check.write(f + '\n' + '1\n' + str(statinfo.st_size) +'\n' + str(now) + '\n999')
    check.close()
    return
    

def sorter(san,pgm):
    '''
    Summary
    -------
        sorter(san,core):
        program that runs through all rx data for one site to sort it.

    Inputs
    ------
        san = site abbreviated name (e.g. uao)
        pgm = program running (for pid checking)

    History
    -------
        2/13/14 -- Written by DJF (dfisher2@illinois.edu)
        6/11/15 -- Modified for site-by-site multicore
    '''
    
    # The location of programs (needed since running in crontab)
    dir_local = '/rdata/airglow/'
    #dir_script = '/usr/local/share/airglowrsss/Python/Programs/'
    #python = '/usr/local/python/'
    dir_share = '/rdata/airglow/share/'
    
    # Close Program if already running (just in case...)
    pid = str(os.getpid())
    pidfile = "/tmp/Sorter_%s.pid"%pgm
    if os.path.isfile("/tmp/Sorterall.pid"):
        print "Sorterall exists, must finish before calling new batch"
        sys.exit()
    elif os.path.isfile(pidfile):
        print "%s already exists, exiting" % pidfile
        sys.exit()
    else:
        file(pidfile, 'w').write(pid)
        
    # Start Write to Log File
    print "\n!!!!!!!!!!!!!!!!!!!!"
    print '!!! BEGIN TIMESTAMP:',dt.datetime.now()

    # Load instrument dictionary
    code = instrumentcode()
    # Set order so bwc & x3t process first
    ids = ['Cloud','TempL']
    ids.extend(code.keys())
   
    
    # TRY YOUR HARDEST
    try:
        # Get Data in RX folder
        os.chdir(dir_local+'rx/')
        #os.system('chmod 774 *') No longer have permissions, tx sends as 774.
        # Get info files for non-standard (Zip->Send->Sort) data
        rxfiles = ["fpi04_kaf","cas01_hka"]
        for x in rxfiles:
            for files in glob(x + '*'):
                makeinfo(files)
        # Go through txt files to Sort data
        for i in ids:
            # do everything if all, else limit to a site.
            if pgm == 'all':
                search = '*.txt'
            else:
                search = '*'+san+'*.txt'
            for f in glob(i + search):
                # Get information for assembling & sos this unAmericarting file
                name = f[0:18]              # name         = IIIII_SSS_YYYYMMDD
                instrument = f[0:5].lower() # instrument   = IIIII
                instr = f[0:3].lower()      # instrument   = III__
                inum = f[3:5]               # instrument # = ___II
                site = f[6:9].lower()       # site         = SSS
                # FOR OLDER FILES THAT DID DOY, ELSE STANDARD DAY
                if f[17] in ['.']:
                    date = f[10:17]             # date         = YYYYDDD
                    dates = f[12:17]            # dates        = YYDDD
                    year = int(f[10:14])        # year         = YYYY
                    doy = int(f[14:17])         # doy          = DDD
                    dn = dt.datetime(year,1,1)+dt.timedelta(days = doy-1)
                    month = dn.timetuple().tm_mon
                    day = dn.timetuple().tm_mday
                else:
                    date = f[10:18]             # date         = YYYYMMDD
                    dates = f[12:18]            # dates        = YYMMDD
                    year = int(f[10:14])        # year         = YYYY
                    month = int(f[14:16])       # month        = MM
                    day = int(f[16:18])         # day          = DD
                    dn = dt.datetime(year,month,day)
                    doy = dn.timetuple().tm_yday
                print "\n!!! For", name
                # Fix inum for Letters
                if inum[1].isalpha():
                    inum = inum[1]
                if not(inum[0].isalpha()):
                    emails = activeinstruments()[site][instr][inum]['email']


                ##### TEMPLOG CASE: #####
                if instrument in ['cloud', 'templ']:
                    ### Part 1: Sorting Data
                    print "!!! Begin Sorting..."
                    # Create fake checksum for tracker
                    checkname = code[instrument]+'00'+f[5:]
                    os.system('cp '+f+' '+checkname)
                    makeinfo(checkname)
                    os.rename(checkname,'tracking/'+checkname)
                    # Move data into directory
                    dir_data = dir_local + 'templogs/' + code[code[instrument]] + '/' + site + '/'
                    os.rename(f, dir_data + f)
                    os.system('chmod 744 ' + dir_data + f)
                    #os.system('chown airglow.airglow ' + dir_data + r)
                    print "!!! Success Sorting"
                    
                elif instr in ['bwc', 'x3t']:
                    ### Send Error if checkfile
                    print "!!! Begin Sorting..."
                    sortinghat([],f)


                ##### FPI CASE: #####
                elif instr in 'fpi':
                    ### Part 1: Sorting Data
                    dir_data = dir_local + 'fpi/' + code[instr] + inum + '/' + site + '/' + str(year) + '/'
                    result = sortinghat(dir_data,f)
                    if result:
                        # CHMOD all added files
                        for r in result:
                            os.system('chmod u+rwx,go+rX,go-w ' + dir_data + r)
                            os.system('chown airglow.fpi ' + dir_data + r)
                        # Remove files from rx
                        os.system('rm -f ' + name + '*')
                        print "!!! Success Sorting"
                        
                    ### Part 2: Processing Data
                        print "!!! Begin Processing..."
                        # Get correct doy from files
                        for r in result:
                            if r[-4:] in '.img':
                                # Find solar local time and subtract 12 hours. That's the definition of date that we use for FPIs.
                                # This is the only way to ensure that no matter what location and time zone, all files
                                # from a night refer to the same date.
                                ldn = FPI.ReadIMG(dir_data+r).info['LocalTime']
                                site_info = fpiinfo.get_site_info(site, ldn)
                                utdn = ldn.replace(tzinfo=pytz.timezone(site_info['Timezone'])).astimezone(pytz.utc).replace(tzinfo=None)
                                site_lon = np.mod(site_info['Location'][1]+180,360)-180
                                sltdn = utdn + dt.timedelta(hours = 24*site_lon/360.)
                                dn0 = sltdn - dt.timedelta(hours=12) # No matter what time of night, and what location, this will be during the same date
                                doy = dn0.timetuple().tm_yday
                                year = dn0.year
                                break
                        # Run processing script for site
                        try:
                            warning = FPIprocess.process_instr(code[instr] + inum,year,doy)
                            if warning:
                                subject = "!!! Manually inspect (\'" + code[instr]+inum+'\','+str(year)+','+str(doy)+') @ ' + site
                                print subject
                                Emailer.emailerror(emails, subject, warning)
                        except:
                            subject = "!!! Processing error (\'" + code[instr]+inum+'\','+str(year)+','+str(doy)+') @ ' + site
                            print subject
                            Emailer.emailerror(emails, subject, traceback.format_exc())
                        # Run CV processing for project
                        # ?????
                        print "!!! End Processing"


                ##### GPS CASE: #####
                elif instr in ['tec','scn']:
                    ## Part 1: Sorting Data
                    dir_data = dir_local + 'gps/' + code[instr] + inum + '/' + site + '/' + str(year) + '/'
                    # if SCN - Send data to raw folder
                    if instr == 'scn':
                        dir_data = dir_data + '/raw_data/'
                        try:
                            os.makedirs(dir_data)
                            os.system('chmod 755 ' + dir_data)
                        except OSError:
                            print '!!! Raw Folder Exists... moving on'
                    result = sortinghat(dir_data,f)
                    if result:
                        # CHMOD all added files
                        for r in result:
                            os.system('chmod u+rwx,go+rX,go-w ' + dir_data + r)
                            os.system('chown airglow.gps ' + dir_data + r)
                        # Remove files from rx
                        os.system('rm -f ' + name + '*')

                        print "!!! Success Sorting"
                        
                    ### Part 2: Processing Data
                        print "!!! Begin Processing..."
                        # Run processing script for site
                        try:
                            GPSprocess.process_instr(code[instr] + inum,year,doy)
                        except:
                            subject = "!!! Processing error (\'" + code[instr]+inum+'\','+str(year)+','+str(doy)+') @ ' + site
                            print subject
                            Emailer.emailerror(emails, subject, traceback.format_exc())
                        print "!!! End Processing"


                ##### CASES CASE: #####
                elif instr in ['cas']:
                    ### Partr 1: "Sort Data"
                    # Make sure folder exists
                    dir_data = dir_local + 'gps/'+code[instr] + inum + '/' + site + '/' + str(year) + '/'
                    try:
                        os.makedirs(dir_data)
                        os.system('chmod 755 ' + dir_data)
                    except OSError:
                        print '!!! Raw Folder Exists... moving on'
                    # Move info file to tracking
                    os.system('mv ' + f + ' ./tracking')
                    # Remove files from rx (it was a duplicate)
                    os.system('rm -f ' + name + '*')
                    print "!!! Success Sorting"

                    ### Part 2: Processing Data
                    print "!!! Begin Processing..."
                    # Run processing script for site
                    try:
                        GPSprocess.process_instr(code[instr] + inum,year,doy)
                    except:
                        subject = "!!! Processing error (\'" + code[instr]+inum+'\','+str(year)+','+str(doy)+') @ ' + site
                        print subject
                        Emailer.emailerror(emails, subject, traceback.format_exc())                       
                        
                
                ##### IMAGER CASE: #####
                elif instr in ['asi','nfi','pic','sky','swe']:
                    ### Part 1: Sorting Data
                    dir_data = dir_local + 'imaging/' + code[instr] + inum + '/' + site + '/' + str(year) + '/'
                    result = sortinghat(dir_data,f)
                    if result:
                        if asiinfo.get_site_info(site)['share']:
                            # Check that share folder for copy exists
                            dir_copy = dir_share + site + '/' + str(year) + '/' + str(doy) + '/'
                            try:
                                os.makedirs(dir_copy)
                                os.system('chmod 755' + dir_copy)
                                print "!!! Share Folder Created"
                            except OSError:
                                print '!!! Share Folder Exists... moving on'
                        # CHMOD all added files
                        for r in result:
                            os.system('chmod u+rwx,go+rX,go-w ' + dir_data + r)
                            os.system('chown airglow.imaging ' + dir_data + r)
                            #os.system('mv ' + dir_data + r + ' ' + dir_data + str(doy) + '/.')
                            # Copy files if needed
                            if asiinfo.get_site_info(site)['share']:
                                os.system('cp -r ' + dir_data + r + ' ' + dir_copy)
                        # Remove files from rx
                        os.system('rm -f ' + name + '*')
                        print "!!! Success Sorting"
                        
                    ### Part 2: Processing Data
                        print "!!! Begin Processing..."
                        ## Get correct doy from files
                        #for r in result:
                        #    if r[-4:] in '.tif':
                        #        ldn = Image.open(dir_data+r).info['UniversalTime'] # Local is standard, but ASI is JJM's timechoice
                        #        if ldn.hour < 12:
                        #            ldn -= dt.timedelta(days = 1)
                        #        doy = ldn.timetuple().tm_yday
                        #        year = ldn.year
                        #        break
                        # Run processing script for site
                        # TODO: Mimic FPIprocess warnings
                        msg = ASIprocess.process_instr(code[instr] + inum,year,doy)
                        if msg:
                            subject = "!!! Processing Issue (\'" + code[instr]+inum+'\','+str(year)+','+str(doy)+') @ ' + site
                            print subject
                            Emailer.emailerror(emails, subject, msg)
                        print "!!! End Processing"
                        
                        
                ##### BAD INSTR CATCH #####
                else:
                    emails = activeinstruments()['ADMIN']['email'] 
                    subject = "!!! Badly named files: " + name
                    print subject
                    Emailer.emailerror(emails, subject, 'Name is not real instrument...')

    except:
        emails = activeinstruments()['ADMIN']['email'] 
        subject = "!!! Something is wrong..."
        print subject
        Emailer.emailerror(emails, subject+f, traceback.format_exc())
        
    finally:
        print "\n!!! Unpack Complete!"
        os.unlink(pidfile)
    

        
if __name__=="__main__":

    # Parse the command line
    usage = "usage: Sorter -s SITE"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--site", dest="site", help="specific program to run (eg site,other,all)",
                   metavar="SITE", type="str",default='all')
    #parser.add_option("-s", "--split", dest="split", help="Number of split Sorter Processes",metavar="SPLIT", type="int",default=1)
    (options, args) = parser.parse_args()
    pgm = options.site.lower()

    # CASE= if 'all'
    if pgm == 'all':
        print 'All sites to be sorted - no additional calls may be made!'
        sites = fpiinfo.get_all_sites_info().keys() + \
                asiinfo.get_all_sites_info().keys() + \
                gpsinfo.get_all_sites_info().keys()
        sites = np.unique(sites)
    
    # CASE= if 'other' (not fpi sites)
    elif pgm == 'other':
        fpis  = fpiinfo.get_all_sites_info().keys()
        other = asiinfo.get_all_sites_info().keys() + \
                gpsinfo.get_all_sites_info().keys()
        sites = [x for x in other if x not in fpis]
        
    else:
        sites = [pgm]
    
    # Sort by site!
    for san in sites:
        sorter(san,pgm)
    
