#!/usr/bin/python
'''
Script to Zip, Split and Move PERU data to Sending Folder
  -s to add SITE - XXX location name
  -i to add INSTRUMENT - XXX instrument name
  -n to add NUMBER - # instrument number/letter
  -p to add PRIORDAYS - # Send additional days back
  -y to add YEAR - (USE WITH -d) #### Send a specific Year-Doy back 
  -d to add DOY - (USE WITH -y) ### Send a specific Year-Doy back

History: 2 Oct 2012 - initial script written
        17 Jul 2012 - new server update

Written by Daniel J. Fisher (dfisher2@illinois.edu)
'''

# Import required modules
import os
import sys
from glob import glob
import commands
import datetime as dt
import urllib
#import urllib2
import tarfile
import shutil
from optparse import OptionParser
import Emailer

def zipper(site,instr,num,prior=1,pyear=0,pdoy=0):
    '''
    Summary
    -------
        zipper(site,instr,num,prior,pyear,pdoy)
        Zips, splits, and moves all data to Sending Folder
    
    Inputs
    ------
        site =  XXX - location name
        instr = XXX - instrument name
        num =   ##  - instrument number or letter
        prior = #   - Send additional days back (optional)
        pyear = ### - Year to send back (optional)
        pdoy =  ### - DOY to send back (optional)

    History
    -------
        7/17/13 -- Written by DJF (dfisher2@illionis.edu)
    '''
    # Minimum file size for writing (to find "empty" folders"
    mfs = 1000
    website = 'http://jro-app.igp.gob.pe/cielodatabase/data/'
    rxdir = '/rdata/airglow/rx/'
    
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
        
        # Get prior day's date
        d = dt.datetime.today()-dt.timedelta(dback+1)
        pday = d.strftime('%d')
        
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
            d = dt.date.fromordinal(dt.datetime(pyear,1,1).toordinal()-1+pdoy-1)
            pday = d.strftime('%d')
            
        print instr+': '+day+'-'+mon+'-'+year
        # Go to local directory!
        os.chdir(rxdir)
        name = "%03s01_%s_%04s%02s%02s.tar.gz" %(instr.upper(),site.upper(),year,month,day)
        path = "%03s/FPI/%04s/%02s/%02s/" %(site.upper(),year,month,day)
        filename = "%03s%02s_%s_%04s%02s%02s.tar.gz" %(instr.lower(),num,site.lower(),year,month,day)
        
        print website+path+name
        urllib.urlretrieve(website+path+name,filename)
        #os.rename(name, filename)
        '''
        # Check filesize for offline X300
        if os.stat(rxdir+filename).st_size < 100:
            ## Emails Warning that system is down!
            print 'FPI ERROR!!!'
            msg = "FPI down at %s!\n" % site
            subject = "Site Error Log %02s-%02s-%02s" %(day,month,yr)
            Emailer.emailerror(subject,msg)
        '''

if __name__=="__main__":

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
        
    zipper(site,instr,num,prior,pyear,pdoy)

    print 'Zip Complete...'
    

