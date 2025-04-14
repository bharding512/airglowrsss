#!/usr/bin/python
#
# Script to read in and process temperatures given by the BoltWood
# and X300 sensors.
#
# Creates plot of ambient temperature, room temperature, and dome temperature
#
# Run as process_UAO_temperatures -y YEAR -d DOY
#
# If no options given, run for yesterday's data
#
# Written by Timothy Duly (duly2@illinois.edu) on 23 Aug 2012

import matplotlib
matplotlib.use('AGG')

import datetime
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import sys
sys.path.append("../modules")
import BoltwoodSensor
import X300Sensor

import MySQLdb as mdb
import pytz
from pytz import timezone

# Change these parameters depending
# on the machine you're running in:
# =================================
path_templogs = '/mnt/FPIData/NATION/uao/TempLogs/' # airglow
# path_templogs = "/Users/duly/data/FPIData/temps/" # tim's
# switch_save = 0
# =================================

# Database information:
site_id = 6
inst_id = 79
summary_stub = '/data/Urbana/auxdata/'
web_stub = 'SummaryImages/Urbana/'

# Timezones
utc = pytz.utc
central = pytz.timezone('US/Central')

''' 
get command line arguments for year and doy.
Default is yesterdays date.
'''
usage = "usage: process_UAO_temperatures -y YEAR -d DOY"
parser = OptionParser(usage=usage)

# Get yesterday's date information as a default
dn = datetime.date.today()-datetime.timedelta(days = 1)
yyyy = dn.year
ddd = dn.timetuple().tm_yday

# Parse the command line
parser.add_option("-y", "--year", dest="year", help="year to be processed",
		metavar="YEAR", type="int",default=yyyy)
parser.add_option("-d", "--doy", dest="doy", help="day of year to be processed",
		metavar="DOY", type="int",default=ddd)

(options, args) = parser.parse_args()
year = options.year
doy = options.doy

# Grab the datetime of our input
dn = datetime.datetime(year,1,1) + datetime.timedelta(days=doy-1)

# call modules from BoltWood and X300 sensor
dnsB, sky_temp, amb_temp = BoltwoodSensor.ReadTempLog(path_templogs + "Cloud_UAO_" + dn.strftime("%Y%m%d") + ".txt")
dnsX, temp1, temp2 = X300Sensor.ReadTempLog(path_templogs + "TempL_UAO_" + dn.strftime("%Y%m%d") + ".txt")

# Plot the temperature data
fig = plt.figure(1); plt.clf()
ax = fig.add_subplot(111)
ax.plot(dnsB,amb_temp,'.',dnsX,temp1,'.',dnsX,temp2,'.')
plt.grid('on')
plt.axis([dn, dn+datetime.timedelta(days=1), 0., 100.])
plt.title("UAO " + dn.strftime("%Y-%m-%d"))
plt.xlabel("LT"); plt.ylabel("Temperature [C]")
plt.legend(("amb [BoltWood]", "dome [X300]","room [X300]"),loc=2)
ax.xaxis.set_major_formatter( mdates.DateFormatter("%H:%M") )
fig.autofmt_xdate()
plt.draw(); plt.show()

# Save image file
summary_png = 'UAOaux_' + dn.strftime('%Y%m%d') + '.png'
fig.savefig(summary_stub + summary_png)

# Get the start and stop UT
d0 = central.localize(min(dnsB))
d1 = central.localize(max(dnsB))
startut = d0.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')
stoput = d1.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')

# Open and populate the database (see http://zetcode.com/databases/mysqlpythontutorial/ for help)
con = mdb.connect(host='localhost',user='WebUser',passwd='webpass',db='webdatabase')
cur = con.cursor()

# First find out if the entry is in there (i.e., we are just updating the png and avi file)
sql_cmd = 'SELECT id FROM DataSet WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (site_id, inst_id, startut)
cur.execute(sql_cmd)
rows = cur.fetchall()
if len(rows) == 0:
    sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime,SummaryImage,SummaryMovie) VALUES(%d, %d, \"%s\", \"%s\", \"%sauxdata/%s\", \"\")' % (site_id, inst_id, startut, stoput, web_stub, summary_png)
    cur.execute(sql_cmd)
