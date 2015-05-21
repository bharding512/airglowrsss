#!/usr/bin/python
'''
Run program to track what files are being processed.
Email if lacking checkfile for a day or two... 
Email if disc space on a computer is low...

History: 26 Feb 2014 - initial script written

Written by Daniel J. Fisher (dfisher2@illinois.edu)
'''

# Import required modules
import os
from glob import glob
import datetime as dt
import Emailer
from Zipper import activeinstruments
   
    
if __name__=="__main__":

    # Get Active Instruments:
    code = activeinstruments()
    
    # Get current time
    now = dt.datetime.utcnow()
    
    # Set warning limits
    hoursback = 27  # Must be slightly longer than 24 hrs - allows transfer time to be sluggish. Probably less than 30 hrs so warnings don't come after we've left work.
    diskfree = 5    # GB remaining. Some drives are 1TB, others are 20GB. Some systems use 100MB a day, others 1GB. So 5 seems acceptable.  Currently you must manually run code to get space.
    
    #locations
    tracking = '/rdata/airglow/rx/tracking/'
    rx = '/rdata/airglow/rx/'
    scripts = '/home/airglow/DataManagementRO/'
    os.chdir(tracking)
    os.system('chmod 770 *')
    os.chdir(rx)
    os.system('chmod 770 *.txt')
    
    # Loop through each instrument to find most recent file.
    for site in code.keys():
        for instr in code[site].keys():
            for num in code[site][instr].keys():
                name = instr+num+'_'+site
                # Sort files (newest on bottom)
                os.chdir(tracking)
                files = glob(name+'*')
                files.sort()
                # Remove older files
                for f in files[:-1]:
                    os.system('rm -f ' + f)
                if 'x3t' == instr:
                    # Ignore missing x3t warnings so leave loop
                    break
                # Append checksums not yet sorted
                os.chdir(rx)
                files = files + glob(name+'*.txt')
                files.sort()
                # Send email if no files exist
                if files == []:
                    subject = "!!! No data received:" + name
                    print subject
                    Emailer.emailerror(subject, 'There are no checkfiles for this site/instrument.\nIs it active? Did you set up the scripts?')
                else:
                    # Check file stats
                    try:
                        info = open(tracking+files[-1], 'r')
                    except:
                        info = open(rx+files[-1], 'r')
                    info.readline()
                    info.readline()
                    info.readline()
                    time = dt.datetime.strptime(info.readline()[:19],'%Y-%m-%d %H:%M:%S')
                    df = float(info.readline())
                    info.close()
                    # Send email if file is too old
                    age = (now - time).total_seconds()/3600.0
                    if age > hoursback:
                        subject = "!!! No data received:" + site
                        print subject
                        Emailer.emailerror(subject, '%s: this file is %i hours old.\nIs the internet connected?  Is the computer down?\nBad Zip? Try -p %i' % (name,age,age/24))
                    # Check disk space
                    if df < diskfree:
                        subject = "!!! Disk space low:" + name
                        print subject
                        try:
                            os.system('sh %slist_%s.sh' % (scripts,name))
                            Emailer.emailerror(subject, 'There are %3.2f GB remaining on the drive. You may want to free some space soon...\nRun sync on remote computer.' % df)
                        except:
                            Emailer.emailerror(subject, 'There are %3.2f GB remaining on the drive. You may want to free some space soon...\nCreate list script first.' % df)
                    
                
