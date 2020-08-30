#!/usr/bin/python
'''
Script to Zip, Split and Move recent data to Sending Folder
  -s to add SITE - XXX location name
  -i to add INSTRUMENT - XXX instrument name
  -n to add NUMBER - # instrument number/letter
  -p to add PRIORDAYS - # Send additional days back
  -y to add YEAR - (USE WITH -d) #### Send a specific Year-Doy back 
  -d to add DOY - (USE WITH -y) ### Send a specific Year-Doy back

History: 2 Oct 2012 - initial script written
        17 Jul 2012 - new server update
        12 Feb 2014 - Updated to v3.0 - txtcheck
        31 Mar 2016 - Added email addresses so emails can be sent directly

Written by Daniel J. Fisher (dfisher2@illinois.edu)
'''

# Import required modules
import os
import sys
from glob import glob
import datetime as dt
import urllib
import tarfile
import shutil
from optparse import OptionParser
from collections import defaultdict
import Emailer


def activeinstruments():
    '''
    Summary:
        code = activeinstruments()
        Is the dictionary for use in DATAMANAGEMENT ONLY of all active instruments (site-instr-num) with:
            send_dir:  Sending Directory
            local_dir: Data Directory
            email:     Email list who wants warnings
            split:     Location of split command [for FPI and ASI]
            url:       Address of temp logs [for X3T]
    
    Outputs:
        code =  dictionary of all active site/instrument/num info

    History:
        2/25/14 -- Written by DJF (dfisher2@illinois.edu)
        3/31/16 -- Added emails DJF
    '''
    
    # Set up site dictionary: Comment out or delete non-active instrument/sites
    # This dictionary has 2 functions:
    # 1. for data-taking computer list the path to the data directories and special functions
    # 2. for remote2 server for emailing person of interest.

    # Email list
    UIemail = ['lnav@illinois.edu','grawe2@illinois.edu','jmakela@illinois.edu']
    BZemail = ['rburiti.ufcg@gmail.com']
    AKemail = []
    MOemail = ['zouhair@uca.ac.ma']
    
    # The dictionary!
    code = defaultdict(lambda: defaultdict(dict))

    code['uao']['fpi']['05'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    #code['uao']['sky']['01'] = {'send_dir':'D:/Sending/', 'local_dir':'D:/Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    #code['uao']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Users/MiniME/Documents/ClarityII/','email':UIemail}
    code['uao']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Users/MiniME/Documents/Interactiveastronomy/SkyAlert/','email':UIemail}
    code['uao']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.0.2/log.txt','email':UIemail}

    # EKU is being moved as of 2-17-2016
    #code['eku']['fpi']['07'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    #code['eku']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Documents and Settings/meriwej/My Documents/ClarityII/','email':UIemail}
    #code['eku']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://157.89.43.12/log.txt','email':UIemail}

    # ANN offline ATM
    #code['ann']['fpi']['08'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    #code['ann']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Users/fpi/My Documents/ClarityII/','email':UIemail}
    #code['ann']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.0.102/log.txt','email':UIemail}
    
    # PAR deconstructed and being sent to S Africa 2 Aug 2017 
    #code['par']['fpi']['06'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    #code['par']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Documents and Settings/meriwej/My Documents/ClarityII/','email':UIemail}
    #code['par']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://10.20.1.2/log.txt','email':UIemail}
    #code['par']['fpi']['09'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    
    # SAO site is a copy of PAR
    code['sao']['fpi']['06'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    code['sao']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Documents and Settings/meriwej/My Documents/ClarityII/','email':UIemail}
    code['sao']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.0.110/log.txt','email':UIemail}

    # VTI offline
    #code['vti']['fpi']['09'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    #code['vti']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Documents and Settings/meriwej/My Documents/ClarityII/','email':UIemail}
#    code['vti']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://10.20.1.2/log.txt','email':UIemail}
    
#    code['caj']['fpi']['02'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail+BZemail}
#    code['caj']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Documents and Settings/MiniME/My Documents/ClarityII/','email':UIemail+BZemail}
#    code['caj']['pic']['05'] = {'send_dir':'/data/Sending/', 'local_dir':'/data/', 'split':'/usr/bin/split','email':UIemail+BZemail}
#    code['caj']['scn']['0U'] = {'send_dir':'/home/scintmon/Sending/', 'local_dir':'/home/scintmon/cascade-1.62/', 'split':'/usr/bin/split','email':UIemail+BZemail}
#    code['caj']['tec']['02'] = {'send_dir':'/home/gps/Sending/', 'local_dir':'/home/gps/data/', 'split':'/usr/bin/split','email':UIemail+BZemail}
#    code['caj']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.1.2/log.txt','email':UIemail+BZemail}

    # Car offline due to water damage Feb 5 2018 
    #code['car']['fpi']['01'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail+BZemail}
    #code['car']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Documents and Settings/MiniME/My Documents/ClarityII/','email':UIemail+BZemail}
    # code['car']['scn']['0S'] = {'send_dir':'F:/Sending/', 'local_dir':'F:/Data/ScintmonS/', 'split':'/usr/bin/split','email':UIemail+BZemail}
    # code['car']['scn']['0T'] = {'send_dir':'F:/Sending/', 'local_dir':'F:/Data/ScintmonT/', 'split':'/usr/bin/split','email':UIemail+BZemail}
    #code['car']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.1.2/log.txt','email':UIemail+BZemail}

    # Morocco is back online 17 July 2018
    code['mor']['fpi']['03'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    code['mor']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Users/admin/Documents/ClarityII/','email':UIemail}
#    code['mor']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.1.204/log.txt','email':UIemail}
#    code['mor']['pic']['04'] = {'send_dir':'/data/Sending/', 'local_dir':'/data/', 'split':'/usr/bin/split','email':UIemail}

    #code['bog']['scn']['0O'] = {'send_dir':'/home/gps/Sending/', 'local_dir':'/home/gps/cascade/', 'split':'usr/bin/split','email':UIemail}

    # CTO is broken, CCD overheats and Filterwheel does not work
    #code['cto']['pic']['02'] = {'send_dir':'/home/airglow/Sending/', 'local_dir':'/data/', 'split':'/usr/bin/split','email':UIemail}
    #code['cto']['scn']['0L'] = {'send_dir':'/home/gps/Sending/', 'local_dir':'/home/gps/cascade/', 'split':'/usr/bin/split','email':UIemail}
    #code['cto']['scn']['0M'] = {'send_dir':'/home/gps/Sending/', 'local_dir':'/home/gps2/cascade/', 'split':'/usr/bin/split','email':UIemail}

    # CASI stopped sending data 2/10/2019
#     code['hka']['asi']['01'] = {'send_dir':'C:/Sending/', 'local_dir':'/cygdrive/c/Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    code['hka']['nfi']['01'] = {'send_dir':'/home/airglow/Sending/', 'local_dir':'/data/', 'split':'/usr/bin/split','email':UIemail}
    #code['hka']['scn']['0K'] = {'send_dir':'/home/gps/Sending', 'local_dir':'/home/gps/cascade', 'split':'/usr/bin/split','email':UIemail}
#    code['hka']['cas']['01'] = {'send_dir':'/mnt/data/Sending/', 'local_dir':'mnt/data/','split':'/usr/bin/split','email':UIemail}

    #REMOVED Feb 27, 2017
    #code['tht']['pic']['06'] = {'send_dir':'/home/airglow/Sending/', 'local_dir':'/data/', 'split':'/usr/bin/split','email':UIemail}

    #NSO now in Texas April 10, 2017
    #code['nso']['pic']['03'] = {'send_dir':'/home/airglow/Sending/', 'local_dir':'/data/', 'split':'/usr/bin/split','email':UIemail}

    # pic03 currently being shipped to UTA
#    code['uta']['pic']['03'] = {'send_dir':'/home/airglow/Sending/', 'local_dir':'/data/', 'split':'/usr/bin/split','email':UIemail}
    

    #code['kaf']['fpi']['04'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    #code['kaf']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Users/???/My Documents/ClarityII/','email':UIemail}
   
   
    # BDR has been offline for a while -- originally turned off for rainy season (7 Sep 2017 BJH)
    # It seems BDR came back online for a few days. Turning back on to process these. (16 Jan 2019 BJH)
    # Disabling BDR once more (20 jan 2019)
    #code['bdr']['fpi']['94'] = {'send_dir':'C:/Sending/', 'local_dir':'/cygdrive/c/FPIData/', 'split':'/usr/bin/split','email':UIemail}
    #code['bdr']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Documents and Settings/User/My Documents/ClarityII/','email':UIemail}
    #code['bdr']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.1.204/log.txt','email':UIemail}

#    code['arg']['fpi']['93'] = {'send_dir':'C:/Sending/', 'local_dir':'/cygdrive/c/FPIData/', 'split':'/usr/bin/split','email':UIemail}
#    code['arg']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Users/meriwej/Documents/ClarityII/','email':UIemail}
    #code['arg']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.1.204/log.txt','email':UIemail}
    #code['xxx']['xxx']['00'] = {'local_dir':'.','email':UIemail}

    # Kwaj FPI
    code['kwj']['fpi']['09'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'/usr/bin/split','email':UIemail}
    code['kwj']['bwc']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Users/User/My Documents/ClarityII/','email':UIemail}
    #code['kwj']['x3t']['00'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Scripts/Python/modules/', 'url':'http://192.168.0.2/log.txt','email':UIemail}
   
    # Leo FPI
    code['leo']['fpi']['80'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/FPI_Data/', 'split':'C:/cygwin/bin/split','email':UIemail}
    #code['leo']['bwc']['80'] = {'send_dir':'C:/Sending/', 'local_dir':'C:/Users/meriwej/Documents/ClarityII/','email':UIemail}

    # Last index is for remote2 admin related issues (not site based)
    code['ADMIN']={'email':['jmakela@illinois.edu']}

    return(code)


def doer(site,instr,num,prior=1,pyear=0,pdoy=0):
    '''
    Summary:
        doer(site,instr,num,prior,pyear,pdoy)
        Zips, splits, and moves all data to Sending Folder
    
    Inputs:
        site =  XXX - location name
        instr = XXX - instrument name
        num =   ##  - instrument number or letter
        prior = #   - Send additional days back (optional)
        pyear = ### - Year to send back (optional)
        pdoy =  ### - DOY to send back (optional)

    History:
        7/17/13 -- Written by DJF (dfisher2@illinois.edu)
    '''
    
    # Minimum file size for writing (to find "empty" folders)
    mfs = 1000
    
    code = activeinstruments()

    # Go through all possible days
    for dback in range(prior,0,-1):

        # Get yesterdays day's date
        d = dt.datetime.today()-dt.timedelta(dback)
        doy = d.strftime('%j')
        year = d.strftime('%Y')
        yr = d.strftime('%y')
        month = d.strftime('%m')
        mon = d.strftime('%b')
        day = d.strftime('%d')
        
        # Get today's date
        d = dt.datetime.today()-dt.timedelta(dback-1)
        nday = d.strftime('%d')
        now = dt.datetime.utcnow()
        
        # Get todays day's date
        d = dt.datetime.today()
        tyear = d.strftime('%Y')
        tmonth = d.strftime('%m')
        tday = d.strftime('%d')

        # Do single day (assumes no -p)
        if pyear > 0:
            d = dt.date.fromordinal(dt.datetime(pyear,1,1).toordinal()-1+pdoy)
            doy = d.strftime('%j')
            year = d.strftime('%Y')
            yr = d.strftime('%y')
            month = d.strftime('%m')
            mon = d.strftime('%b')
            day = d.strftime('%d')
        
            # Get prior day's date
            d = dt.date.fromordinal(dt.datetime(pyear,1,1).toordinal()-1+pdoy+1)
            nday = d.strftime('%d')
            
        print instr+': '+day+'-'+mon+'-'+year
        # Go to local directory!
        os.chdir(code[site][instr][num]['local_dir'])
        filename = "%03s%02s_%s_%04s%02s%02s.tar.gz" %(instr,num,site,year,month,day)
        checkname = "%03s%02s_%s_%04s%02s%02s.txt" %(instr,num,site,year,month,day)
  
        
    ########## TO ZIP UP FPI ##########
        if 'fpi' == instr:
            # Grab all files made within 24 hour period
            name = [f for f in glob(os.path.join('*','*.img')) if (dt.datetime(int(year),int(month),int(day),12) < dt.datetime.fromtimestamp(os.path.getmtime(f)) and dt.datetime.fromtimestamp(os.path.getmtime(f)) < (dt.datetime(int(year),int(month),int(day),12)+dt.timedelta(1)))]
            zipper(name,filename)
            splitter(site,instr,num,code,filename,checkname,mfs)
            os.remove(code[site][instr][num]['local_dir']+filename)


    ########## TO ZIP UP NFI/ASI ##########  
        if 'asi' == instr: # or 'nfi' == instr:
            name = "%02s-%02s" %(day,nday)
            os.chdir(code[site][instr][num]['local_dir']+year)
            # Gun-Tar Folder
            tar = tarfile.open(filename,"w:gz")
            try:
                # Change folder name to doy before tarballing & then change back
                os.rename(code[site][instr][num]['local_dir']+year+'/'+mon+name,code[site][instr][num]['local_dir']+year+'/'+doy)
                tar.add(doy);
                os.rename(code[site][instr][num]['local_dir']+year+'/'+doy,code[site][instr][num]['local_dir']+year+'/'+mon+name)
            except:
                print "!!! No Data"
            tar.close()
            splitter(site,instr,num,code,filename,checkname,mfs)
            os.remove(code[site][instr][num]['local_dir']+year+'/'+filename)

                
    ########## TO ZIP UP x300 ##########
        if 'x3t' == instr:
            filename = "TempL_%s_%04s%02s%02s.txt" %(site,year,month,day)
            urllib.urlretrieve(code[site][instr][num]['url'],"X300_temp_log.txt")
            #url = urllib2.urlopen(code[site][instr][num]['url'])
            #log = open(code[site][instr][num]['local_dir']+"X300_temp_log.txt", 'wb')
            #log.write(url.read())
            #log.close()
            sys.path = sys.path + [code[site][instr][num]['local_dir']]
            import X300Sensor
            try:
                X300Sensor.WriteLogFromRaw(dt.datetime(int(year),1,1)+ dt.timedelta(days = int(doy)-1),code[site][instr][num]['local_dir'],code[site][instr][num]['send_dir'],site)
            except:
                try:
                    X300Sensor.WriteLogFromRaw(dt.datetime(int(year),1,1)+ dt.timedelta(days = int(doy)-1),code[site][instr][num]['local_dir'],code[site][instr][num]['send_dir'],site)
                except:
                    print 'X300 Error!'
            # Check filesize for offline X300
            if os.stat(code[site][instr][num]['send_dir']+filename).st_size < 100:
                # Send checkfile that system is down!
                print 'X300 ERROR!!!'
                check = open(checkname, 'w')
                check.write(filename+'\n0\n0\n'+str(now)+'\n999\n')
                # Legend for checkfile + boosts size for Sending
                check.write('1: 1st file name\n2: Number of parts\n3: Size of tar.gz\n4: Time of Creation\n5: GB Disk Free')
                
            # Get as much current day as possible to cover all night
            if prior == 1:
                filename = "TempL_%s_%04s%02s%02s.txt" %(site,tyear,tmonth,tday)
                X300Sensor.WriteLogFromRaw(dt.datetime(int(year),1,1)+ dt.timedelta(days = int(doy)),code[site][instr][num]['local_dir'],code[site][instr][num]['send_dir'],site)
            #os.system('mv '+filename+' '+code[site][instr][num]['send_dir'])
            #os.system('mv log.txt '+code[site][instr][num]['send_dir']+filename)


    ########## TO ZIP UP CLOUD ##########
        if 'bwc' == instr:
            name = "%04s-%02s-%02s.txt" %(year,month,day)
            filename = "Cloud_%s_%04s%02s%02s.txt" %(site,year,month,day)
            # try SkyAlert
            if "SkyAlert" in code[site][instr][num]['local_dir']:
                name = "/Weatherdata Files/%02s-%02s-%04s.txt" %(month,day,year)
            if os.path.isfile(code[site][instr][num]['local_dir']+name):
                #os.rename(code[site][instr][num]['local_dir']+name,code[site][instr][num]['send_dir']+filename)
                shutil.copy2(code[site][instr][num]['local_dir']+name,code[site][instr][num]['send_dir']+filename)
                # If Downloading yesterday get todays partial file too.
                if prior == 1:
                    name = "%04s-%02s-%02s.txt" %(tyear,tmonth,tday)
                    if "SkyAlert" in code[site][instr][num]['local_dir']:
                        name = "weatherdatafiles.txt"
                    filename = "Cloud_%s_%04s%02s%02s.txt" %(site,tyear,tmonth,tday)
                    shutil.copy2(code[site][instr][num]['local_dir']+name,code[site][instr][num]['send_dir']+filename)
            else:
                ## Send checkfile that system is down!
                print 'No Data Error!'
                os.chdir(code[site][instr][num]['send_dir'])
                check = open(checkname, 'w')
                check.write(filename+'\n0\n0\n'+str(now)+'\n999\n')
                # Legend for checkfile + boosts size for Sending
                check.write('1: 1st file name\n2: Number of parts\n3: Size of tar.gz\n4: Time of Creation\n5: GB Disk Free')
                
                
    ########## TO ZIP UP PIC ##########
        if 'pic' == instr or 'nfi' == instr:
            os.chdir(code[site][instr][num]['local_dir']+year)
            # Grab all files made within 24 hour period
            name = [f for f in glob(os.path.join('*','*.tif')) if (dt.datetime(int(year),int(month),int(day),12) < dt.datetime.fromtimestamp(os.path.getmtime(f)) < (dt.datetime(int(year),int(month),int(day),12)+dt.timedelta(1)))]
            zipper(name,filename)
            splitter(site,instr,num,code,filename,checkname,mfs)
            os.remove(code[site][instr][num]['local_dir']+year+'/'+filename)


    ########## TO ZIP UP SKY ##########
        if 'sky' == instr:
            name = "%04s%02s%02s/" %(year,month,day)
            zipper(name,filename)
            splitter(site,instr,num,code,filename,checkname,mfs)
            os.remove(code[site][instr][num]['local_dir']+filename)


    ########## TO ZIP UP SCN ##########  
        if 'scn' == instr:
            name = "%02s%02s%02s*" %(yr,month,day)
            zipper(name,filename)
            splitter(site,instr,num,code,filename,checkname,mfs)
            os.remove(code[site][instr][num]['local_dir']+filename)


    ########## TO ZIP UP TEC ##########
        if 'tec' == instr:
            name = glob("%02s%02s%02s*" %(yr,month,day))
            zipper(name,filename)
            splitter(site,instr,num,code,filename,checkname,mfs)
            os.remove(code[site][instr][num]['local_dir']+filename)
            
        if 'xxx' == instr:
            print '\nPlease check your input:\nsite -s\ninstrumetn -i\nnumber -n\n'


def zipper(name,filename):
    '''
    Summary:
        zipper(name,filename)
        Zips up all data
    
    Inputs:
        name    = X     - name* of files to be zipped up
        filename    = XXX##_XXX_########.tar.gz - tarball filename
        
    History:
        7/17/13 -- Written by DJF (dfisher2@illinois.edu)
    '''
    
    # Guntar folder
    tar = tarfile.open(filename,"w:gz")
    if(type(name).__name__=='str'):
        files = glob(name)
    else:
        files = name
    for oscar in files:
        tar.add(oscar)
    tar.close()


def splitter(site,instr,num,code,filename,checkname,mfs):
    '''
    Summary:
        splitter(site,instr,num,code,filename,checkname,mfs)
        Splits and moves all data to Sending Folder + checkfile
    
    Inputs:
        site    = XXX   - location name
        instr   = XXX   - instrument name
        num     = ##    - instrument number or letter
        code    = DICT  - Dictionary of locations
        filename    = XXX##_XXX_########.tar.gz - tarball filename
        checkname   = XXX##_XXX_########.txt    - checksum filename
        mfs     = #     - minimum file size

    History:
        7/17/13 -- Written by DJF (dfisher2@illinois.edu)
    '''

    # Find size of data
    statinfo = os.stat(filename)
    now = dt.datetime.utcnow()
    # Split and move files to Sending
    if statinfo.st_size > mfs :
        os.system(code[site][instr][num]['split']+' -a 6 -b 5242880 -d '+filename+' '+code[site][instr][num]['send_dir'] + filename)
        # Check Disk Space Free in GB
        s = os.statvfs(code[site][instr][num]['local_dir'])
        df = round(float(s.f_bavail * s.f_frsize) / 1024**3,2)
        # Make checkfile 
        os.chdir(code[site][instr][num]['send_dir'])
        check = open(checkname, 'w')
        pieces = glob(filename + '*')
        check.write(pieces[0]+'\n' + str(len(pieces))+'\n' + str(statinfo.st_size)+'\n'+str(now)+'\n'+str(df)+'\n')
    else:
        print 'COLLECTION ERROR!!!'
        # Make error checkfile 
        os.chdir(code[site][instr][num]['send_dir'])
        check = open(checkname, 'w')
        check.write(filename+'\n0\n0\n'+str(now)+'\n999\n')
    # Legend for checkfile + boosts size for Sending
    check.write('1: 1st file name\n2: Number of parts\n3: Size of tar.gz\n4: Time of Creation\n5: GB Disk Free')
    check.close()
    
    
    
if __name__=="__main__":

    # Main module allows zipper to be run via command line

    # Parse the command line
    usage = "usage: Zipper -s SITE -i INSTRUMENT -n NUMBER -p DAYSPRIOR -y YEAR -d DOY"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--site", dest="site", help="Site to be run",
                   metavar="SITE", type="str",default='XXX')
    parser.add_option("-i", "--instrument", dest="instr", help="Instrument to be run",
                   metavar="INSTRUMENT", type="str",default='XXX')
    parser.add_option("-n", "--number", dest="num", help="Number/Letter of instrument",
                   metavar="NUMBER", type="str",default='00')
    parser.add_option("-p", "--prior", dest="prior", help="Days Prior to be run",
                   metavar="PRIORDAYS", type="int",default=1)
    parser.add_option("-y", "--year", dest="pyear", help="Year to be run",
                   metavar="YEAR", type="int",default=0)
    parser.add_option("-d", "--doy", dest="pdoy", help="Doy to be run",
                   metavar="DOY", type="int",default=0)

    (options, args) = parser.parse_args()
    site = options.site.lower()
    instr = options.instr.lower()
    num = options.num.zfill(2)
    prior = options.prior
    pyear = options.pyear
    pdoy = options.pdoy

    if prior > 1 and pyear > 0:
        print 'Prior Days (-p) and Specific Days (-y & -d) cannot be used together'
        sys.exit()
        
    if len(instr)>3:
        print 'Usable Instrument Keys: fpi,asi,nfi,pic,sky,swe,cas,tec,scn,bwc,x3t'
        
    doer(site,instr,num,prior,pyear,pdoy)

    print 'Zip Complete...'
    
