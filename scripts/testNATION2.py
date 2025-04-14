# -*- coding: utf-8 -*-
# Script to dynamically control the NATION network
#
# History: 14 Aug 2012 - initial script to control just the UAO FPI
#
# Written by Jonathan J. Makela (jmakela@illinois.edu)

# Import required modules
import ephem
import datetime
import time
import sys
import os
import numpy as np

# Sequence to follow
Sequence = ['Zenith','North','East','South','West','Laser']

# Data transfers can take some time and there may be some slack between
# clocks on all of the machines.  Times within TimeSlack seconds
# of one another are treated as if they are the same time
TimeSlack = datetime.timedelta(seconds = 5)

# Get the ephems for the sun
sun = ephem.Sun()

# Set up sites dictionary
sites = {'UAO' : {'Location': (40.133, -88.2, 200),
       'Observations': {'Laser': {'ze': 180, 'az': 87, 'exptime': 30, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'Zenith': {'ze': 0, 'az': 0, 'exptime': 180, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'North': {'ze': 45, 'az': 0, 'exptime': 180, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'South': {'ze': -45, 'az': 0, 'exptime': 180, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'East': {'ze': 45, 'az': 90, 'exptime': 180, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'West': {'ze': -45, 'az': 90, 'exptime': 180, 
                                  'n_exp': 0, 'last_exp': None}},
       'BufferTime': 45,
       'Sunset': None,
       'Sunrise': None,
       'ObservationHistory': [],
       'CurrentStartTime': None,
       'NextStartTime': None,
       'LocalSunFile': '/home/master/NATION/Control/UAOSunriseSet.txt',
       'LocalControlFile': '/home/master/NATION/Control/UAOControl.txt',
       'scpUser': 'MiniME',
       'scpPort': 19999,
       'scpDir': '/cygdrive/c/NATIONControl',
       'ephem': ephem.Observer(),
       'SequenceCounter': 0,
       'online': True}}

sites['EKU'] = {'Location': (37.75, -84.29, 300),
       'Observations': {'Laser': {'ze': 180, 'az': 104, 'exptime': 30, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'Zenith': {'ze': 0, 'az': 0, 'exptime': 300, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'North': {'ze': 45, 'az': 0, 'exptime': 300, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'South': {'ze': -45, 'az': 0, 'exptime': 300, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'East': {'ze': 45, 'az': 90, 'exptime': 300, 
                                  'n_exp': 0, 'last_exp': None}, 
                        'West': {'ze': -45, 'az': 90, 'exptime': 300, 
                                  'n_exp': 0, 'last_exp': None}},
       'BufferTime': 45,
       'Sunset': None,
       'Sunrise': None,
       'ObservationHistory': [],
       'CurrentStartTime': None,
       'NextStartTime': None,
       'LocalSunFile': '/home/master/NATION/Control/EKUSunriseSet.txt',
       'LocalControlFile': '/home/master/NATION/Control/EKUControl.txt',
       'scpUser': 'meriwej',
       'scpPort': 19998,
       'scpDir': '/cygdrive/c/NATIONControl',
       'ephem': ephem.Observer(),
       'SequenceCounter': 0,
       'online': True}

# Set up Observer structure and other info for each site
for site in np.unique(sites):
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

        print site, ' sunset at ', sites[site]['Sunset'], ', sunrise at ', sites[site]['Sunrise']

        # Send sun file to site
        SunFile = open(sites[site]['LocalSunFile'],'w')
        SunFile.write('SUNUID % s\n' % sites[site]['Sunset'].strftime('%Y%m%d'))
        SunFile.write('beginNightTime % s\n' % sites[site]['Sunset'].strftime('%Y.%m.%d.%H:%M:%S'))
        SunFile.write('endNightTime % s\n' % sites[site]['Sunrise'].strftime('%Y.%m.%d.%H:%M:%S'))
        SunFile.close()

        # SCP file over to system
        result = os.system('scp -P %d %s %s@localhost:%s' % (sites[site]['scpPort'], sites[site]['LocalSunFile'], sites[site]['scpUser'], sites[site]['scpDir']))
        if result != 0:
                # Copy failed for some reason
                sites[site]['online'] = False
                print site, 'not responding'
        else:
                sites[site]['online'] = True

        # Reset observation counter
        sites[site]['NumObservations'] = 0

        # Set first observation time
        sites[site]['CurrentStartTime'] = sites[site]['Sunset']
        sites[site]['NextStartTime'] = sites[site]['CurrentStartTime'] + datetime.timedelta(seconds = sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['exptime']) + datetime.timedelta(seconds = sites[site]['BufferTime'])
        print site, 'set next observation of ', Sequence[sites[site]['SequenceCounter']], ' to begin at ', sites[site]['CurrentStartTime']

        # Send Observation file to site
        ControlFile = open(sites[site]['LocalControlFile'],'w')
        ControlFile.write('OBSUID % 04d\n' % sites[site]['NumObservations'])
        ControlFile.write('StartTime % s\n' % sites[site]['CurrentStartTime'].strftime('%Y.%m.%d.%H:%M:%S'))
        ControlFile.write('ze % .2f\n' % sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['ze'])
        ControlFile.write('az % .2f\n' % sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['az'])
        ControlFile.write('exptime % .2f\n' % sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['exptime'])
        if Sequence[sites[site]['SequenceCounter']] == 'Laser':
            ControlFile.write('laser 1\n')
        else:
            ControlFile.write('laser 0\n')
        if Sequence[sites[site]['SequenceCounter']] == 'Dark':
                ControlFile.write('dark 1\n')
        else:
                ControlFile.write('dark 0\n')
        ControlFile.close()

        # SCP file over to system
        result = os.system('scp -P %d %s %s@localhost:%s' % (sites[site]['scpPort'], sites[site]['LocalControlFile'], sites[site]['scpUser'], sites[site]['scpDir']))
        if result != 0:
                # Copy failed for some reason
                sites[site]['online'] = False
                print site, 'not responding'
        else:
                sites[site]['online'] = True

# Calculate the start time based on the sunsets at each site.
# Also calculate the end time based on the sunrises at each site.
StartTime = datetime.datetime.utcnow() + datetime.timedelta(days = 5)
EndTime = datetime.datetime.utcnow() - datetime.timedelta(days = 5)
for site in np.unique(sites):
        if sites[site]['Sunset'] < StartTime:
                # The current site has an earlier sunset than the current overall starttime
                StartTime = sites[site]['Sunset']

        if sites[site]['Sunrise'] > EndTime:
                # The current site has a later sunrise than the current overall endtime
                EndTime = sites[site]['Sunrise']

# Calculate time to wait until observations start
WaitTime = StartTime - datetime.datetime.utcnow()

print 'Waiting ', WaitTime.days*86400 + WaitTime.seconds, ' seconds'
print 'Start time ', StartTime
print 'End time ', EndTime

# NextEventTime contains the time that the next event will start
NextEventTime = StartTime

# Wait for the observations to start
if(WaitTime.days >= 0):
	time.sleep(WaitTime.days*86400 + WaitTime.seconds)

# TODO: CHECK MOON LOCATION!!

# Enter the main control loop, which runs until sunrise at the last site
while NextEventTime <= EndTime:
        # Wait for the next event time to occur in a while loop.
        # Eventually, we will check in with the sites during this
        # loop to determine cloud coverage and connectivity
        print 'Waiting for %s (current time: %s)' % (NextEventTime + TimeSlack, datetime.datetime.utcnow())
        while(datetime.datetime.utcnow() <= NextEventTime + TimeSlack):
                time.sleep(1)

        # Check the next start times at each site for the next event that will occur
        NextEventTime = datetime.datetime.utcnow() + datetime.timedelta(days = 5)
        for site in np.unique(sites):
                if sites[site]['NextStartTime'] < NextEventTime and sites[site]['NextStartTime'] > sites[site]['Sunset'] and sites[site]['NextStartTime'] < sites[site]['Sunrise']:
                        NextEventTime = sites[site]['NextStartTime']

        # Find the sites who's next event is within TimeSlack of the NextEventTime
        NextSites = []
        for site in np.unique(sites):
                if abs(sites[site]['NextStartTime'] - NextEventTime) < TimeSlack:
                        NextSites.append(site)

        print 'Next sites:', NextSites, 'at', NextEventTime

        # Record the current observation in the history at each site for which we are 
        for site in NextSites:
                sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['n_exp'] += 1
                sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['last_exp'] = sites[site]['CurrentStartTime']
                sites[site]['ObservationHistory'].append((Sequence[sites[site]['SequenceCounter']], sites[site]['CurrentStartTime']))

                # Increment this sites sequence counter
                sites[site]['SequenceCounter'] += 1
                if sites[site]['SequenceCounter'] >= len(Sequence):
                        sites[site]['SequenceCounter'] = 0

                # Increment the total number of observations at this site
                sites[site]['NumObservations'] += 1

                # Set up the next observation
                sites[site]['CurrentStartTime'] = sites[site]['NextStartTime']
                sites[site]['NextStartTime'] = sites[site]['CurrentStartTime'] + datetime.timedelta(seconds = sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['exptime']) + datetime.timedelta(seconds = sites[site]['BufferTime'])

                # Send Observation file to site
                ControlFile = open(sites[site]['LocalControlFile'],'w')
                ControlFile.write('OBSUID % 04d\n' % sites[site]['NumObservations'])
                ControlFile.write('StartTime % s\n' % sites[site]['CurrentStartTime'].strftime('%Y.%m.%d.%H:%M:%S'))
                ControlFile.write('ze % .2f\n' % sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['ze'])
                ControlFile.write('az % .2f\n' % sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['az'])
                ControlFile.write('exptime % .2f\n' % sites[site]['Observations'][Sequence[sites[site]['SequenceCounter']]]['exptime'])
                if Sequence[sites[site]['SequenceCounter']] == 'Laser':
                        ControlFile.write('laser 1\n')
                else:
                        ControlFile.write('laser 0\n')
                if Sequence[sites[site]['SequenceCounter']] == 'Dark':
                        ControlFile.write('dark 1\n')
                else:
                        ControlFile.write('dark 0\n')
                ControlFile.close()
        
		print '(', datetime.datetime.utcnow(), ') Sending request to', site

                # SCP file over to system
                result = os.system('scp -P %d %s %s@localhost:%s' % (sites[site]['scpPort'], sites[site]['LocalControlFile'], sites[site]['scpUser'], sites[site]['scpDir']))
                if result != 0:
                       # Copy failed for some reason
                       sites[site]['online'] = False
                       print site, 'not responding'
                else:
                       sites[site]['online'] = True

                print '(', datetime.datetime.utcnow(), ') Sent request to', site, 'to observe', Sequence[sites[site]['SequenceCounter']], 'beginning at', sites[site]['CurrentStartTime']

##        # Figure out the next event time
##        NextEventTime = datetime.datetime.utcnow() + datetime.timedelta(days = 5)
##        for site in np.unique(sites):
##                if sites[site]['NextStartTime'] < NextEventTime:
##                        NextEventTime = sites[site]['NextStartTime']

print 'Sequence Finished'
