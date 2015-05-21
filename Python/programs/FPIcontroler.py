# -*- coding: utf-8 -*-
# Script to dynamically control the NATION network
#
# History: 10 May 2013 - Look directions read in from text file
#	   15 Oct 2012 - Converted to reward based decision tree
#	   14 Aug 2012 - initial script to control just the UAO FPI
#	   30 Aug 2012 - added control of EKU
#	    6 Sep 2012 - added control of PARI
#         5-9 Aug 2013 - modified to use optparser, update dictionaries from files, and other improvments (jjm)
#
#Laser Delay 10 minutes in in-code backup. Change to match obsfile when/if permanent value decided.
#LINES 798 TEMPORARILY SET FOR 2 FPIs at PARI. REMOVE THE 3 FOLLOWING LINES WHEN DONE.
#
# Written by Tom Gehrels (gehrels2@illinois.edu) based on a template written by Jonathan J. Makela (jmakela@illinois.edu)

# Import required modules
import ephem
import datetime
import time
import sys
import os
import numpy as np
import subprocess
import math
import string
import pdb
import BoltwoodSensor as bs
import pytz
from pytz import timezone
import logging
import socket
import sys, traceback
import fpiinfo
from optparse import OptionParser

def GrabInfo(file_name):
    from datetime import datetime

    fid = open(file_name,'r')

    OBSUID = int(fid.readline().split()[1])

    time = fid.readline().split()[1]
    year = int(time[0:4])
    month = int(time[4:6])
    day = int(time[6:8])
    hour = int(time[8:10])
    minute = int(time[10:12])
    second = int(time[12:14])
    dn = datetime(year,month,day,hour,minute,second)

    try:
        # See if this file has an intensity and integration time in it
	intensity = float(fid.readline().split()[1])
	int_time = float(fid.readline().split()[1])
    except:
	intensity = -999.9
	int_time = -999.9

    fid.close()

##    print file_name, intensity, int_time

    return OBSUID, dn, intensity, int_time

def MoonAngle(site, Az, Ze, obsTime):
        # Function that returns True if the observation is obscured by the moon
        # and False if the moon is not within 37 degrees of the requested look
        # direction
        #
        # INPUT: site - a sites dictionary, as defined below in the main program
        #        Az, Ze - the azimuth and zenith of the reqested observation in degrees
        #        obsTime - a datetime for the requested observation
        #
        # HISTORY: Written on 11 Sep 2012 by Jonathan J. Makela

        # Threshold for the moon/observation angle (37 degrees from original Master scripts)
        coneThresh = 37 * math.pi/180

        # Set the time of the sites ephem structure
        site['ephem'].date = obsTime

        # Find the moon at this time
        moon = ephem.Moon(site['ephem'])

        # Grab the moon az and ze (converted from elevation angle)
        moonAz = moon.az
        moonZe = ephem.degrees(math.pi/2 - moon.alt)

        # Run the math (based on code in original Master scripts for FPI observations
        a = math.cos(Az*math.pi/180)*math.sin(Ze*math.pi/180)
        b = math.sin(Az*math.pi/180)*math.sin(Ze*math.pi/180)
        aMoon = math.cos(moonAz)*math.sin(moonZe)
        bMoon = math.sin(moonAz)*math.sin(moonZe)
        moonAngle = math.acos(a*aMoon + b*bMoon + math.cos(Ze*math.pi/180) * math.cos(moonZe))

        if (moonAngle>coneThresh) or (Ze>90 or Ze<-90):
                return (False, moonAngle)
        else:
                return (True, moonAngle)

def Cloudy(site):
    # Function that returns the sky status of the requested site based on Boltwood Sensor data
    #
    # INPUT: site - a sites dictionary, as defined below in the main program
    #
    # OUTPUT: returns True for cloudy conditions, false otherwise.  Also returns the sky-ambient temperature
    #
    # HISTORY: Written on 26 Nov 2012 by Jonathan J. Makela
    #
    # TODO: CHECK THAT FILE ISNT OUTDATED

    # Read the file
    dns, sky_temp, amb_temp = bs.ReadRawTempLog(site['LocalCloudFile'], site['TimeZone'])

    if (not sky_temp) and (not sky_temp==0):
	if site['State']=='cloudy':
		sky_temp = [site['CloudThresh']+10]
	else:
		sky_temp = [site['CloudThresh']-10]

    if (sky_temp[0]>site['CloudThresh']) or np.isnan(sky_temp[0]):
        return (True, sky_temp)
    else:
        return (False, sky_temp)

def permutations(availt,sitest):
	avail=dict.copy(availt)
	sites=dict.copy(sitest)
	perms={}
	n=0
	#pdb.set_trace()
	allsite=np.unique(avail)
	if not len(allsite):
		return perms
	allmeasures=sites[allsite[0]]['Directions']
	del avail[allsite[0]]
	for cmeasure in np.unique(allmeasures):
		omeasures={}
		cavail=dict.copy(avail)
		breakout=0
		for site in np.unique(sites):
			if site == allsite[0]:
				continue
			if string.find(cmeasure.lower(),site.lower()) > -1:
				if avail.get(site,0):
					omeasures[site]=cmeasure
					del cavail[site]
				else:
					breakout=1
		if breakout:
			continue
		rperms=permutations(cavail,sites)
		if len(rperms):
			for k in xrange(0,len(rperms)):
				perms[n]={}
				perms[n][allsite[0]]=cmeasure
				if len(omeasures):
					perms[n].update(omeasures)
				perms[n].update(rperms[k])
				n=n+1
		else:
			perms[n]={}
			perms[n][allsite[0]]=cmeasure
			if len(omeasures):
				perms[n].update(omeasures)
			n+=1

	return perms

def rewardfunction(perm,t,taction,sites):
	tscale=.2
	dscale=datetime.timedelta(seconds = 20)
	reward=0
	for site in np.unique(perm):
		CVcard=False
		if sites[site]['State']=='CV':
			iscale=2
			cvscale=2
			cscale=0
			lscale=.5
			zscale=.5
			usedelay=1
		elif sites[site]['State']=='cardinal':
			iscale=0
			cvscale=0
			cscale=.6
			lscale=.5
			zscale=.5
			usedelay=1
		else:
			iscale=0
			cvscale=0
			cscale=.6
			lscale=.5
			zscale=.5
			usedelay=0

		if sites[site]['Directions'][perm[site]]['last_exp'] is None:
			lastt=min(sites[site]['Sunset'],taction-datetime.timedelta(minutes = 60))
		else:
			lastt=sites[site]['Directions'][perm[site]]['last_exp']

		delay = (usedelay!=0)*datetime.timedelta(seconds=sites[site]['Directions'][perm[site]]['delay'])

		if string.find(perm[site],'Laser') > -1:
			scale=lscale
		elif string.find(perm[site],'Zenith') > -1:
			scale=zscale
		elif string.find(perm[site],'CV') > -1:
			scale=cvscale
		elif string.find(perm[site],'IN') > -1:
			scale=iscale
		else:
			if sites[site]['State']=='CV':
				CVcard=True
				Cardinal = dict.keys(sites[site]['Directions'])
				if ('Laser' in Cardinal):
					Cardinal.remove('Laser')
				if('Zenith' in Cardinal):
					Cardinal.remove('Zenith')
				for card in Cardinal:
					if string.find(card.lower(),site.lower()) > -1:
						Cardinal.remove(card)
				
				lasttime={}
				for card in Cardinal:
					lasttime[card]=sites[site]['Directions'][card]['last_exp']
	
				cardnone=[k for k, v in lasttime.iteritems() if v==None]
				cardsort=sorted([(value,key) for (key,value) in lasttime.items() if value!=None])
				creward={}
				cval=lscale
				for card in cardnone:
					creward[card]=cval
				for card in cardsort:
                                        cval=cval-.01
					creward[card[1]]=cval

			else:
				scale=cscale
		if CVcard:
			reward = reward + creward[perm[site]]
		else:	
			reward = reward + scale*np.exp(tscale*float((taction-lastt-delay).seconds+86400*(taction-lastt-delay).days)/float((taction-t+dscale).seconds))

	return reward

# Set up the command line parser
usage = "usage: FPIcontroler -n NETWORK -d -s"
parser = OptionParser(usage=usage)
parser.add_option("-n", "--network", dest="network", help="network to use. (renoir, nation)",
                  metavar="NETWORK", type="string", default="nation")
parser.add_option("-d", "--dynamic", dest="dynamicintegration", action="store_true", help="use dynamic integration timing", default=False)
(options, args) = parser.parse_args()
network = options.network.lower()
dynamicintegration = options.dynamicintegration

# Set up the logger
today = datetime.datetime.utcnow()
logname = '/home/master/LOG_' + network + '_' + today.strftime('%Y%m%d') + '.log'
logging.basicConfig(filename=logname, level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.info('Logfile Opened')

#File to read observation information from
obsfile='/home/master/SourceCode/Python/programs/obsfile_%s.txt' % network
sitefile='/home/master/SourceCode/Python/programs/sites_%s.txt' % network

# Flag to calculate dynamic integration timing
if dynamicintegration:
    logging.info('Configured for dynamic integration times')

# Data transfers can take some time and there may be some slack between
# clocks on all of the machines.  Times within TimeSlack seconds
# of one another are treated as if they are the same time
TimeSlack = datetime.timedelta(seconds = 5)

# Get the ephems for the sun
sun = ephem.Sun()

# Read the sites dictionary from the fpiinfo function and update it with the sites to be used as
# defined in the sitesfil
sites = fpiinfo.get_network_info(network)
try:
	f=open(sitefile,'r')
	temp=eval(f.read())
	f.close()

	# Remove any sites not included in sitefile (http://stackoverflow.com/questions/11277432/how-to-remove-a-key-from-dictionary)
        missing_sites = []
        for site in sites:
            if site in temp:
                sites[site].update(temp[site])
            else:
                missing_sites.append(site)

        for site in missing_sites:
            sites.pop(site,None)

	logging.info('Sites dictionary updated')
            
except:
	logging.info('%s read failed', sitefile)
	traceback.print_exc(file=sys.stdout)

# Read the observation dictionary from the obsfile and update the Direction dictionary for each site
# Attempt to read the observation file from the requested location.  If it succeeds,
# use these look directions.  Otherwise, fall back to the backup defined above.
try:
	f=open(obsfile,'r')
	obs=eval(f.read())
	f.close()

        # First, deactivate all directions
	for site in sites:
            for d in sites[site]['Directions']:
                sites[site]['Directions'][d]['active'] = False

        # Second, update the observation dictionary from what was read from obsfile
        for site in obs:
            if site in sites:
                for d in obs[site]:
                    sites[site]['Directions'][d].update(obs[site][d])
            else:
                logging.info('obsfile contained %s which was not in available sites', site)

        logging.info('Directions dictionary updated')

        # Print out active observations
##        for site in sites:
##            for d in sites[site]['Directions']:
##                if sites[site]['Directions'][d]['active']:
##                    logging.info('%s %s is active with exptime=%.1f, delay=%.1f', site, d, sites[site]['Directions'][d]['exptime'], sites[site]['Directions'][d]['delay'])
except:
	logging.info('%s read failed', obsfile)
	traceback.print_exc(file=sys.stdout)



# Set up Observer structure and other info for each site
for site in np.unique(sites):
        #print site
        sites[site]['ephem'].lat = str(sites[site]['Location'][0])
        sites[site]['ephem'].lon = str(sites[site]['Location'][1])
        sites[site]['ephem'].elevation = sites[site]['Location'][2]
        sites[site]['ephem'].pressure = 0
        sites[site]['ephem'].horizon = '-8:0' # for astro twilight
        sites[site]['Sunset'] = sites[site]['ephem'].next_setting(sun).datetime()
        sites[site]['Sunrise'] = sites[site]['ephem'].next_rising(sun).datetime()

        # Check if we are after sunset
        if sites[site]['Sunset'] > sites[site]['Sunrise']:
                sites[site]['Sunset'] = datetime.datetime.utcnow() + datetime.timedelta(seconds = 5)

        logging.info('%s sunset at %s, sunrise at %s', site, sites[site]['Sunset'], sites[site]['Sunrise'])

        # Send sun file to site
        SunFile = open(sites[site]['LocalSunFile'],'w')
        SunFile.write('SUNUID % s\n' % sites[site]['Sunset'].strftime('%Y%m%d'))
        SunFile.write('beginNightTime % s\n' % (sites[site]['Sunset']).strftime('%Y.%m.%d.%H:%M:%S'))
        SunFile.write('endNightTime % s\n' % (sites[site]['Sunrise']).strftime('%Y.%m.%d.%H:%M:%S'))
        SunFile.close()
	
	# SCP file over using a thread
	port = '-P %d' % sites[site]['scpPort']
	source = sites[site]['LocalSunFile']
	destination = '%s@localhost:%s' % (sites[site]['scpUser'], sites[site]['scpDir'])
	subprocess.Popen(['scp', port, source, destination])

        # Check if the site is online
	s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	try:
		s.connect(('localhost',sites[site]['scpPort']))
		s.shutdown(2)
		sites[site]['online']=True
		logging.info('%s is currently online', site)
	except:
		sites[site]['online']=False
		logging.info('%s is currently not online', site)

        # Reset observation counter
        sites[site]['NumObservations'] = 0

        # Set first observation time
	sites[site]['CurrentStartTime']=sites[site]['Sunset']
	sites[site]['NextStartTime']=sites[site]['Sunset']

	# Read in the LastImage file
	sites[site]['remoteLastOBSUID'], sites[site]['remoteLastStartTime'], sites[site]['remoteLastIntensity'], sites[site]['remoteLastIntTime'] = GrabInfo(sites[site]['LocalLastImageFile'])
	logging.info('Read %s and found %s %s %s %s', sites[site]['LocalLastImageFile'], sites[site]['remoteLastOBSUID'], sites[site]['remoteLastStartTime'], sites[site]['remoteLastIntensity'], sites[site]['remoteLastIntTime'])

# Calculate the start time based on the sunsets at each site.
# Also calculate the end time based on the sunrises at each site.
StartTime = datetime.datetime.utcnow() + datetime.timedelta(days = 5)
EndTime = datetime.datetime.utcnow() - datetime.timedelta(days = 5)
MaxExposure = datetime.timedelta(seconds=0)
for site in np.unique(sites):
        if sites[site]['Sunset'] < StartTime:
                # The current site has an earlier sunset than the current overall starttime
                StartTime = sites[site]['Sunset']

        if sites[site]['Sunrise'] > EndTime:
                # The current site has a later sunrise than the current overall endtime
                EndTime = sites[site]['Sunrise']

	temp = datetime.timedelta(seconds = (sites[site]['Directions']['Zenith']['exptime']+sites[site]['BufferTime']))
	if temp > MaxExposure:
		MaxExposure = temp

# NextEventTime contains the time that the next event will start
NextStartTime = StartTime
NextEventTime = NextStartTime - datetime.timedelta(seconds=60)

# Calculate time to wait until observations start
WaitTime = NextEventTime - datetime.datetime.utcnow()

logging.info('Waiting %d seconds', WaitTime.days*86400 + WaitTime.seconds)
logging.info('Start time %s', StartTime)
logging.info('End time %s', EndTime)

# Wait for the observations to start
if(WaitTime.days >= 0):
	time.sleep(WaitTime.days*86400 + WaitTime.seconds)

# Enter the main control loop, which runs until sunrise at the last site
while NextEventTime <= EndTime:
        # Wait for the next event time to occur in a while loop.
        # Eventually, we will check in with the sites during this
        # loop to determine cloud coverage and connectivity
        #print '(', datetime.datetime.utcnow(), ') Waiting for %s (current time: %s)' % (NextEventTime + TimeSlack, datetime.datetime.utcnow())

        # Flag to force reading of the obsfile during the loop              
	toread = True

	# While waiting for the next event to pass, see if changes to the observation have been requested.
        while(datetime.datetime.utcnow() <= NextEventTime + TimeSlack):
##		lastexp={}
##		for site in sites:
##			lastexp[site]={}
##			for observation in sites[site]['Directions']:
##				lastexp[site][observation]=sites[site]['Directions'][observation]['last_exp']

		if toread:
                    # Read the sites dictionary from the fpiinfo function and update it with the sites to be used as
                    # defined in the sitesfil
##                    try:
##                            f=open(sitefile,'r')
##                            temp=eval(f.read())
##                            f.close()
##                            
##                            #Only update the keys in the if statement below
##                            for site in temp:
##                                topop=[]
##                                for key in temp[site]:
##                                    if (key is not "CloudThresh") and (key is not "desiredIntensity") and (key is not "minIntTime") and (key is not "maxIntTime"):
##                                        topop.append(key)
##                                for key in topop:
##                                    temp[site].pop(key,None)
##
##                            # Remove any sites not included in sitefile (http://stackoverflow.com/questions/11277432/how-to-remove-a-key-from-dictionary)
##                            missing_sites = []
##                            for site in sites:
##                                if site in temp:
##                                    sites[site].update(temp[site])
##                                else:
##                                    missing_sites.append(site)
##
##                            for site in missing_sites:
##                                sites.pop(site,None)
##
##                            logging.info('Sites dictionary updated')
##                                
##                    except:
##                            logging.info('%s read failed', sitefile)
##                            traceback.print_exc(file=sys.stdout)

                                        
                    # Read the observation dictionary from the obsfile and update the Direction dictionary for each site
                    # Attempt to read the observation file from the requested location.  If it succeeds,
                    # use these look directions.  Otherwise, fall back to the backup defined above.
                    try:
                            f=open(obsfile,'r')
                            obs=eval(f.read())
                            f.close()

                            # First, deactivate all directions
                            for site in sites:
                                for d in sites[site]['Directions']:
                                    sites[site]['Directions'][d]['active'] = False

                            # Second, update the observation dictionary from what was read from obsfile
                            for site in obs:
                                if site in sites:
                                    for d in obs[site]:
                                        sites[site]['Directions'][d].update(obs[site][d])
                                else:
                                    logging.info('obsfile contained %s which was not in available sites', site)

                            logging.info('Directions dictionary updated')

                            # Print out active observations
##                            for site in sites:
##                                for d in sites[site]['Directions']:
##                                    if sites[site]['Directions'][d]['active']:
##                                        logging.info('%s %s is active with exptime=%.1f, delay=%.1f', site, d, sites[site]['Directions'][d]['exptime'], sites[site]['Directions'][d]['delay'])
                    except:
                            logging.info('%s read failed', obsfile)
                            traceback.print_exc(file=sys.stdout)
                            
		    toread = False

                # Try downloading the LastImage files from the sites 
                for site in np.unique(sites):

			#Check if site's LastImage file is more than 10 minutes old
			#if sites[site]['Sunset']<=datetime.datetime.utcnow() and sites[site]['Sunrise']>datetime.datetime.utcnow():
			    #LastImageDelay=(time.mktime(time.localtime())-os.path.getmtime(sites[site]['LocalLastImageFile']))/60
			    #if LastImageDelay>10:
				#logging.info('%s LastImage file hasn\'t updated in %s m', site, LastImageDelay)


			# Check if we are after the requested next start time
			if sites[site]['CurrentStartTime'] <= datetime.datetime.utcnow():

                            # Check if we have confirmed that this site's observation has begun.  Check if we already are downloading the file as well
                            if sites[site]['remoteObsBegan'] is False and sites[site]['remoteLastImageFileThread'].poll() == 0:
                            
                                # Read the LastImageFile
                                sites[site]['remoteLastOBSUID'], sites[site]['remoteLastStartTime'], sites[site]['remoteLastIntensity'], sites[site]['remoteLastIntTime'] = GrabInfo(sites[site]['LocalLastImageFile'])
                                
                                # Check if this is the correct observation
                                if sites[site]['remoteLastOBSUID'] == sites[site]['LastOBSUID']:
                                    sites[site]['CurrentStartTime'] = sites[site]['remoteLastStartTime']
                                    
                                    # Update the next start time based on the actual start time for this observation
##				    print 'here3'
##				    print sites[site]['Directions'][sites[site]['ObservationHistory'][-1][0]]['exptime']
##				    print sites[site]['ObservationHistory'][-1][2]
				    sites[site]['NextStartTime'] = sites[site]['CurrentStartTime'] + datetime.timedelta(seconds = sites[site]['ObservationHistory'][-1][2]) + datetime.timedelta(seconds = sites[site]['BufferTime'])
##                                    sites[site]['NextStartTime'] = sites[site]['CurrentStartTime'] + datetime.timedelta(seconds = sites[site]['Directions'][sites[site]['ObservationHistory'][-1][0]]['exptime']) + datetime.timedelta(seconds = sites[site]['BufferTime'])
                                    sites[site]['NextStartTime'] = max(sites[site]['NextStartTime'],datetime.datetime.utcnow()+TimeSlack+datetime.timedelta(seconds=sites[site]['BufferTime']))

                                    if sites[site]['NextStartTime'] < NextStartTime:
                                            NextStartTime = max(sites[site]['NextStartTime'],datetime.datetime.utcnow()+TimeSlack+datetime.timedelta(seconds=sites[site]['BufferTime']))

                                    if (sites[site]['NextStartTime']-datetime.timedelta(seconds=sites[site]['BufferTime'])) < NextEventTime:
                                            NextEventTime = sites[site]['NextStartTime'] - datetime.timedelta(seconds=sites[site]['BufferTime'])

                                    # We have begun the requested observation
                                    sites[site]['remoteObsBegan'] = True

                                    logging.info('%s New start: %s, Next: %s', site, sites[site]['CurrentStartTime'], sites[site]['NextStartTime'])
                                #else:
                                    #logging.info('%s Obs mismatch -remote: %s, master: %s', site, sites[site]['remoteLastOBSUID'], sites[site]['LastOBSUID'])
					

		time.sleep(1)

 	
	#Figure out all possible start times
	t_avail=[]
        for site in sites:
		#Check connectivity
		s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		try:
			s.connect(('localhost',sites[site]['scpPort']))
			s.shutdown(2)
			sites[site]['online']=True
		except:
			sites[site]['online']=False
			logging.info('%s is currently not online', site)

		#Check if cloudy
		cloud, st = Cloudy(sites[site])
		logging.info('%s cloud check.  Sky-ambient: %.1f.  Thresh: %.1f', site, st, sites[site]['CloudThresh'])
		if cloud:
			sites[site]['State']='cloudy'
		else:
			sites[site]['State']='cardinal'

		if (sites[site]['NextStartTime']<sites[site]['Sunrise']) and sites[site]['online']:
			t_avail.append(sites[site]['NextStartTime'])

	#Get rid of times that are too far ahead
	t_avail = sorted(t_avail)
	for n in xrange(1,np.size(t_avail)):
		if (t_avail[n]-t_avail[0])>MaxExposure:
			while np.size(t_avail)>n:
				junk=t_avail.pop()
			break

	#Find all possible permutations at each start time
	allperm={}
	allt={}
	n=0
	for times in np.unique(t_avail):
		avail={}
		for site in sites:
			if (sites[site]['NextStartTime']<=times) and (sites[site]['NextStartTime']<sites[site]['Sunrise']) and sites[site]['online']:
				avail[site]=1
	
		perms=permutations(avail,sites)

		#Get rid of permutations that have moon, CV ones that are cloudy, and deactivated ones
		for p_i in perms:
			keep=True
			perm=perms[p_i]
			for site in perm:
				moon,ma=MoonAngle(sites[site],sites[site]['Directions'][perm[site]]['az'],sites[site]['Directions'][perm[site]]['ze'],times)
				if moon:
					keep=False
					#logging.info('%s - %s mooned out',site,perm[site])
					#logging.info('Removing mooned out direction %s %s', site, perm[site])
					break

				# This should contain the directions to observe during cloudy conditions 
				dowhilecloudy = string.find(perm[site],site.lower())==-1
				
				if (sites[site]['State']=='cloudy') and (not dowhilecloudy):
                                        #logging.info('Removing cloudy direction %s %s', site, perm[site])
					keep=False
					break
				if not sites[site]['Directions'][perm[site]]['active']:
                                        #logging.info('Removing deactivated direction %s %s', site, perm[site])
                                        keep=False
                                        break
			if keep:
				allperm[n]=perm
				allt[n]=times
				n+=1

	#Determine if site is cloudy, CV, or cardinal
	for site in sites:
		if sites[site]['State']=='cloudy':
			continue
		for p_i in allperm:
			perm=allperm[p_i]
			if site not in perm:
				continue
			if string.find(perm[site].lower(),site.lower()) > -1:
				sites[site]['State']='CV'
				break

	#Find best permutation
	best=0
	bperm=[]
	tbest=[]
	for p_i in allperm:
		reward=rewardfunction(allperm[p_i],NextStartTime,allt[p_i],sites)
		if reward>best:
				best=reward
				tbest=allt[p_i]
				bperm=allperm[p_i]

	if len(bperm):
		logging.info('%s - %s - %s',bperm, tbest, datetime.datetime.utcnow())
	else:
		logging.info('No sites online')

	startt={}
	for site in np.unique(bperm):
		logging.info('%s - %s', site, sites[site]['State'])
		if sites[site]['NextStartTime']<tbest:
			double=0
			if string.find(bperm[site].lower(),site.lower()) > -1:
				double=1
			if double:
				startt[site]=tbest
			else:
				startt[site]=max(sites[site]['NextStartTime'],tbest-datetime.timedelta(seconds=300))
		else:
			startt[site]=tbest
						



        # Record the current observation in the history at each site for which we are 
	for site in np.unique(bperm):
                # Read in the LastImage file
        	sites[site]['remoteLastOBSUID'], sites[site]['remoteLastStartTime'], sites[site]['remoteLastIntensity'], sites[site]['remoteLastIntTime'] = GrabInfo(sites[site]['LocalLastImageFile'])
                logging.info('Read %s and found %s %s %s %s', sites[site]['LocalLastImageFile'], sites[site]['remoteLastOBSUID'], sites[site]['remoteLastStartTime'], sites[site]['remoteLastIntensity'], sites[site]['remoteLastIntTime'])

                if (dynamicintegration) and (bperm[site] != 'Laser') and (sites[site]['State']!='cloudy'):
                    # Calculate dynamic integration time, if available (TODO: BASED ON LOOK DIRECTION?)
                    logging.info('%s last intensity %.1f (desired: %.1f)', site, sites[site]['remoteLastIntensity'],sites[site]['desiredIntensity'])
                    if (sites[site]['remoteLastIntensity'] > 0.0):
                        # Scale the integration time based on the laser integration time used and intensity achieved 
                        int_time = sites[site]['desiredIntensity']/sites[site]['remoteLastIntensity']*sites[site]['remoteLastIntTime']

                        # Round this up to the nearest 10 seconds
                        int_time = int_time + (10-int_time % 10)

                        # Set the exposure time
                        sites[site]['NextExpTime'] = np.min([np.max([int_time,sites[site]['minIntTime']]),sites[site]['maxIntTime']])
                        logging.info('%s dynamic integration time %.1f s', site, sites[site]['NextExpTime'])
                    else:
                        sites[site]['NextExpTime'] = sites[site]['Directions'][bperm[site]]['exptime']
                        logging.info('%s static integration time %.1f s (no remote intensity information)', site, sites[site]['NextExpTime'])
                else:
                    sites[site]['NextExpTime'] = sites[site]['Directions'][bperm[site]]['exptime']
                    logging.info('%s static integration time %.1f s (%s, %s)', site, sites[site]['NextExpTime'],bperm[site],sites[site]['State'])

	print 'about to set times'
        # Write the observation control files
        for site in np.unique(bperm):   
                # Save information to the site dictionary
		sites[site]['Directions'][bperm[site]]['n_exp'] += 1
		sites[site]['Directions'][bperm[site]]['last_exp'] = startt[site]
		sites[site]['CurrentStartTime'] = startt[site]
		sites[site]['ObservationHistory'].append((bperm[site],sites[site]['CurrentStartTime'],sites[site]['NextExpTime']))
		sites[site]['NumObservations'] += 1
		sites[site]['NextStartTime'] = sites[site]['CurrentStartTime'] + datetime.timedelta(seconds = sites[site]['NextExpTime'] + sites[site]['BufferTime'])
#		sites[site]['NextStartTime'] = sites[site]['CurrentStartTime'] + datetime.timedelta(seconds = sites[site]['Directions'][bperm[site]]['exptime'] + sites[site]['BufferTime'])

		# Send Observation file to site
                ControlFile = open(sites[site]['LocalControlFile'],'w')
                ControlFile.write('OBSUID % 04d\n' % sites[site]['NumObservations'])
                ControlFile.write('StartTime % s\n' % (sites[site]['CurrentStartTime']).strftime('%Y.%m.%d.%H:%M:%S'))
                ControlFile.write('ze % .2f\n' % sites[site]['Directions'][bperm[site]]['ze'])
                ControlFile.write('az % .2f\n' % sites[site]['Directions'][bperm[site]]['az'])
		ControlFile.write('exptime % .2f\n' % sites[site]['NextExpTime'])
#                ControlFile.write('exptime % .2f\n' % sites[site]['Directions'][bperm[site]]['exptime'])
                if bperm[site] == 'Laser':
                        ControlFile.write('laser 1\n')
                else:
                        ControlFile.write('laser 0\n')
                if bperm[site] == 'Dark':
                        ControlFile.write('dark 1\n')
                else:
                        ControlFile.write('dark 0\n')
                ControlFile.close()

                sites[site]['LastOBSUID'] = sites[site]['NumObservations']
		sites[site]['remoteObsBegan'] = False

##		# Check cloud conditions
##		cloud, st = Cloudy(sites[site])
##		logging.info('%s cloud check.  Sky-ambient: %.1f.  Thresh: %.1f', site, st, sites[site]['CloudThresh'])
##		if cloud:
##			sites[site]['State']='cloudy'
##		else:
##                        sites[site]['State']='cardinal'
        
		#logging.info('Sending request to %s', site)

	        # SCP file over using a thread
                port = '-P %d' % sites[site]['scpPort']
                source = sites[site]['LocalControlFile']
    	        destination = '%s@localhost:%s' % (sites[site]['scpUser'], sites[site]['scpDir'])
    	        subprocess.Popen(['scp', port, source, destination])

                logging.info('Sent request to %s to observe %s (%.1f, %.1f) beginning at %s with integration time %.1f', site, bperm[site], sites[site]['Directions'][bperm[site]]['ze'], sites[site]['Directions'][bperm[site]]['az'], sites[site]['CurrentStartTime'], sites[site]['NextExpTime'])
		moon,ma=MoonAngle(sites[site],sites[site]['Directions'][bperm[site]]['az'],sites[site]['Directions'][bperm[site]]['ze'],startt[site])

	logging.info('\n')

	NextStartTime = datetime.datetime.utcnow() + datetime.timedelta(days = 5)
	NextEventTime = datetime.datetime.utcnow() + datetime.timedelta(days = 5)
	for site in sites:
		if sites[site]['NextStartTime'] < datetime.datetime.utcnow() + datetime.timedelta(seconds=sites[site]['BufferTime']) + TimeSlack:
			sites[site]['NextStartTime'] = datetime.datetime.utcnow() + datetime.timedelta(seconds=sites[site]['BufferTime']) + TimeSlack

		if sites[site]['NextStartTime'] < NextStartTime and sites[site]['NextStartTime'] > sites[site]['Sunset'] and sites[site]['NextStartTime'] < sites[site]['Sunrise'] and sites[site]['online']:
			NextStartTime = sites[site]['NextStartTime']

		if (sites[site]['NextStartTime']-datetime.timedelta(seconds=sites[site]['BufferTime'])) < NextStartTime and sites[site]['NextStartTime'] > sites[site]['Sunset'] and sites[site]['NextStartTime'] < sites[site]['Sunrise'] and sites[site]['online']:
			NextEventTime = sites[site]['NextStartTime'] - datetime.timedelta(seconds=sites[site]['BufferTime'])

	if not len(bperm):
		NextEventTime = datatime.datetime.utcnow() + TimeSlack + datetime.timedelta(seconds=30)
		NextStartTime = NextEventTime + datetime.timedelta(seconds=45)
		

	if NextEventTime < datetime.datetime.utcnow() + TimeSlack:
		NextEventTime = datetime.datetime.utcnow() + TimeSlack

logging.info('Sequence Finished')

obsname = '/home/master/obs_' + today.strftime('%Y%m%d') + '.log'
obsid=open(obsname,'w')
obsid.write(str(sites))
obsid.close()


# Set up tomorrow's observations for each site
for site in np.unique(sites):
        sites[site]['Sunset'] = sites[site]['ephem'].next_setting(sun).datetime()
        sites[site]['Sunrise'] = sites[site]['ephem'].next_rising(sun).datetime()

        # Check if we are after sunset
        if sites[site]['Sunset'] > sites[site]['Sunrise']:
                sites[site]['Sunset'] = datetime.datetime.utcnow() + datetime.timedelta(seconds = 5)

        logging.info('%s sunset at %s, sunrise at %s', site, sites[site]['Sunset'], sites[site]['Sunrise'])

        # Send sun file to site
        SunFile = open(sites[site]['LocalSunFile'],'w')
        SunFile.write('SUNUID % s\n' % sites[site]['Sunset'].strftime('%Y%m%d'))
        SunFile.write('beginNightTime % s\n' % (sites[site]['Sunset']).strftime('%Y.%m.%d.%H:%M:%S'))
        SunFile.write('endNightTime % s\n' % (sites[site]['Sunrise']).strftime('%Y.%m.%d.%H:%M:%S'))
        SunFile.close()

        # SCP file over using a thread
        port = '-P %d' % sites[site]['scpPort']
        source = sites[site]['LocalSunFile']
        destination = '%s@localhost:%s' % (sites[site]['scpUser'], sites[site]['scpDir'])
        subprocess.Popen(['scp', port, source, destination])

logging.info('Finished')
