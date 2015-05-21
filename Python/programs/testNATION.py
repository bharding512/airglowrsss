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

# Get the ephems for the sun
sun = ephem.Sun()

# Set up UAO location dictionary
UAO = {'Location': (40.133, -88.2, 200),
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
       'DropboxSunFile': '/home/master/NATION/Control/UAOSunriseSet.txt',
       'DropboxFile': '/home/master/NATION/Control/UAOControl.txt',
       'scpUser': 'MiniME',
       'scpPort': 19999,
       'scpDir': '/cygdrive/d/NATIONControl'}

# Set up Observer structure for UAO
uao = ephem.Observer()
uao.lat = str(UAO['Location'][0])
uao.lon = str(UAO['Location'][1])
uao.elevation = 200
uao.pressure = 0
uao.horizon = '-12:0' # for astro twilight

# Calculate and save next sunrise and sunset times
UAO['Sunset'] = uao.next_setting(sun).datetime()
UAO['Sunrise'] = uao.next_rising(sun).datetime()

# Check if we are after sunset
if UAO['Sunset'] > UAO['Sunrise']:
	UAO['Sunset'] = datetime.datetime.utcnow() + datetime.timedelta(seconds = 5)

print 'Sunset: ', UAO['Sunset'], 'Sunrise: ', UAO['Sunrise']

# Send sun file to UAO
UAO_sun = open(UAO['DropboxSunFile'],'w')
UAO_sun.write('SUNUID % s\n' % UAO['Sunset'].strftime('%Y%m%d'))
UAO_sun.write('beginNightTime % s\n' % UAO['Sunset'].strftime('%Y.%m.%d.%H:%M:%S'))
UAO_sun.write('endNightTime % s\n' % UAO['Sunrise'].strftime('%Y.%m.%d.%H:%M:%S'))
UAO_sun.close()

# SCP file over to system
os.system('scp -P %d %s %s@localhost:%s' % (UAO['scpPort'], UAO['DropboxSunFile'], UAO['scpUser'], UAO['scpDir']))

# Sequence to follow
Sequence = ['Zenith','North','East','South','West','Laser']

# Reset the observation number value
UAO['NumObservations'] = 0

# Reset sequence counter
SequenceCounter = 0

# Set first observation time
UAO['CurrentStartTime'] = UAO['Sunset']
UAO['NextStartTime'] = UAO['CurrentStartTime'] + datetime.timedelta(seconds = UAO['Observations'][Sequence[SequenceCounter]]['exptime']) + datetime.timedelta(seconds = UAO['BufferTime'])
print 'Set next observation of ', Sequence[SequenceCounter], ' to begin at ', UAO['CurrentStartTime']

# Send Obervation file to UAO
UAO_control = open(UAO['DropboxFile'],'w')
UAO_control.write('OBSUID % 04d\n' % UAO['NumObservations'])
UAO_control.write('StartTime % s\n' % UAO['CurrentStartTime'].strftime('%Y.%m.%d.%H:%M:%S'))
UAO_control.write('ze % .2f\n' % UAO['Observations'][Sequence[SequenceCounter]]['ze'])
UAO_control.write('az % .2f\n' % UAO['Observations'][Sequence[SequenceCounter]]['az'])
UAO_control.write('exptime % .2f\n' % UAO['Observations'][Sequence[SequenceCounter]]['exptime'])
if Sequence[SequenceCounter] == 'Laser':
    UAO_control.write('laser 1\n')
else:
    UAO_control.write('laser 0\n')
UAO_control.write('dark 0')
UAO_control.close()

# SCP file over to system
os.system('scp -P %d %s %s@localhost:%s' % (UAO['scpPort'], UAO['DropboxFile'], UAO['scpUser'], UAO['scpDir']))

# Calculate time to wait until observations start
WaitTime = UAO['Sunset'] - datetime.datetime.utcnow()

if (WaitTime.days*86400 + WaitTime.seconds) > 600:
        HouseKeepingSeconds = datetime.timedelta(seconds = 600)
        StartHouseKeeping = WaitTime - HouseKeepingSeconds
else:
        StartHouseKeeping = WaitTime

print 'Waiting ', StartHouseKeeping.days*86400 + StartHouseKeeping.seconds, ' seconds'
print 'Start time ', UAO['Sunset']

if(StartHouseKeeping.days >= 0):
	time.sleep(StartHouseKeeping.days*86400 + StartHouseKeeping.seconds)

while datetime.datetime.utcnow() <= UAO['Sunrise']:

    # Wait until the current observation start time (in 15 second steps)
    print 'Waiting for %s (current time: %s)' % (UAO['CurrentStartTime'], datetime.datetime.utcnow())
    while (datetime.datetime.utcnow() <= (UAO['CurrentStartTime'])) :
        time.sleep(1)

    # We are now after the CurrentStartTime.  Pause for 10 seconds (must be small enough such that the
    # shortest exposure time + sleep time in the wait loop won't get the remote computer out of step)
    time.sleep(10)

    # Record the current observation in the history
    UAO['Observations'][Sequence[SequenceCounter]]['n_exp'] += 1
    UAO['Observations'][Sequence[SequenceCounter]]['last_exp'] = UAO['CurrentStartTime']
    UAO['ObservationHistory'].append((Sequence[SequenceCounter], UAO['CurrentStartTime']))

    # Increment the sequence counter
    SequenceCounter += 1
    if SequenceCounter >= len(Sequence):
        SequenceCounter = 0

    # Increment the total number of observations at this site
    UAO['NumObservations'] += 1

    # Set up the next observation
    UAO['CurrentStartTime'] = UAO['NextStartTime']
    UAO['NextStartTime'] = UAO['CurrentStartTime'] + datetime.timedelta(seconds = UAO['Observations'][Sequence[SequenceCounter]]['exptime']) + datetime.timedelta(seconds = UAO['BufferTime'])

    # Send Obervation file to UAO
    UAO_control = open(UAO['DropboxFile'],'w')
    UAO_control.write('OBSUID % 04d\n' % UAO['NumObservations'])
    UAO_control.write('StartTime % s\n' % UAO['CurrentStartTime'].strftime('%Y.%m.%d.%H:%M:%S'))
    UAO_control.write('ze % .2f\n' % UAO['Observations'][Sequence[SequenceCounter]]['ze'])
    UAO_control.write('az % .2f\n' % UAO['Observations'][Sequence[SequenceCounter]]['az'])
    UAO_control.write('exptime % .2f\n' % UAO['Observations'][Sequence[SequenceCounter]]['exptime'])
    if Sequence[SequenceCounter] == 'Laser':
    	UAO_control.write('laser 1\n')
    else:
        UAO_control.write('laser 0\n')
    UAO_control.write('dark 0')
    UAO_control.close()

    # SCP file over to system
    os.system('scp -P %d %s %s@localhost:%s' % (UAO['scpPort'], UAO['DropboxFile'], UAO['scpUser'], UAO['scpDir']))

    print '(', datetime.datetime.utcnow(), ') Sent request to observe', Sequence[SequenceCounter], 'beginning at', UAO['CurrentStartTime']

print 'Sequence Finished'
