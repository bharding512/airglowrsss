#!/usr/bin/python
#
# Script to read in and process images from the UAO SkyMonitor camera
# system (a white light imager saving data in the .img format
#
# Reads in Boltwood Cloud Sensor data and prints the ambient-sky
# temperature data on that images as well.
#
# Inserts data into the database
#
# Run as process_SkyMon -y YEAR -d DOY
#
# If no options given, run for yesterday's data
#
# Written by Jonathan J. Makela (jmakela@illinois.edu) on 22 Aug 2012

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import pylab
import numpy
import Image
import glob
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import BoltwoodSensor
import datetime
from optparse import OptionParser
import MySQLdb as mdb
import pytz
from pytz import timezone

usage = "usage: process_SkyMon -y YEAR -d DOY"

parser = OptionParser(usage=usage)

# Get yesterday's date information as a default
dn = datetime.date.today()-datetime.timedelta(days = 1)
yyyy = dn.year
ddd = dn.timetuple().tm_yday

# Database IDs
site_id = 6
inst_id = 78

# Timezones
utc = pytz.utc
central = pytz.timezone('US/Central')

# Directories
img_stub = '/mnt/ImagingData/NATION/uao/2012/'
bw_stub = '/mnt/FPIData/NATION/uao/TempLogs/'
movie_stub = '/data/Urbana/movies/'
summary_stub = '/data/Urbana/skymon/'
web_stub = 'SummaryImages/Urbana/'

# Parse the command line
parser.add_option("-y", "--year", dest="year", help="year to be processed",
		metavar="YEAR", type="int",default=yyyy)
parser.add_option("-d", "--doy", dest="doy", help="day of year to be processed",
		metavar="DOY", type="int",default=ddd)

(options, args) = parser.parse_args()

# Create the image directory name to use
dn = datetime.datetime(options.year,1,1)+datetime.timedelta(days=options.doy-1)
img_dir = dn.strftime('%Y%m%d')

dn1 = dn - datetime.timedelta(days=1)
bw_file1 = 'Cloud_UAO_%s.txt' % dn1.strftime('%Y%m%d')
bw_file2 = 'Cloud_UAO_%s.txt' % dn.strftime('%Y%m%d')

# Read in Boltwood sensor data
bw_date1, bw_sky1, bw_amb1 = BoltwoodSensor.ReadTempLog(bw_stub+bw_file1)
bw_date2, bw_sky2, bw_amb2 = BoltwoodSensor.ReadTempLog(bw_stub+bw_file2)
bw_date = numpy.hstack((bw_date1,bw_date2))
bw_sky = numpy.hstack((bw_sky1, bw_sky2))
bw_amb = numpy.hstack((bw_amb1, bw_amb2))

# Change directories to the Image directory
os.chdir(img_stub+img_dir)

# Find files to be used
sky = glob.glob('S*.img')

sky.sort()

# List for times
myLT = []

# Process each image
for fname in sky:
    # Open the image
    im = Image.open(fname)
    im = im.convert('I')
    data = numpy.array(im.getdata())
    data.resize(im.size[1],im.size[0])

    # Add time
    myLT.append(im.info['LocalTime'])
    
    # Create the image figure
    f = plt.figure()
    i = plt.imshow(data, cmap=plt.cm.gray)
    
    # Find the Boltwood data corresponding to this image time
    closest = sorted(bw_date, key=lambda d:abs(im.info['LocalTime']-d))[0]

    # Format the plot
    plt.clim(1000,10000)
    i.axes.get_xaxis().set_visible(False)
    i.axes.get_yaxis().set_visible(False)
    plt.text(10,500,'%s LT: %.1f C' % (closest, bw_sky[(bw_date == closest).nonzero()][0]),color='white',ha='left',va='bottom',bbox=dict(boxstyle="Round", fc="gray", ec="0.5", alpha="0.5"))
    plt.title('Sky Monitor: %s LT' % im.info['LocalTime'])

    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(cax=cax)

    # Save the png file
    outname = '_' + fname[0:-4] + '.png'
    f.savefig(outname)
    
    plt.close('all')

# Create the movie in a two-pass mpeg4 (see http://matplotlib.sourceforge.net/rc/v1.1.1rc2/examples/old_animation/movie_demo.html and
# http://mariovalle.name/mencoder/mencoder.html)
command = ('mencoder',
           'mf://_S*.png',
           '-mf',
           'type=png:w=800:h=600:fps=5',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4:vpass=1:vbitrate=2343750:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq',
           '-oac',
           'copy',
           '-nosound',
           '-o',
           '/dev/null')

os.spawnvp(os.P_WAIT, 'mencoder', command)

outname = 'SkyMon_' + fname[1:-7] + '.avi'

command = ('mencoder',
           'mf://_S*.png',
           '-mf',
           'type=png:w=800:h=600:fps=5',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4:vpass=2:vbitrate=2343750:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq',
           '-oac',
           'copy',
           '-nosound',
           '-o',
           outname)

os.spawnvp(os.P_WAIT, 'mencoder', command)

# Move files to where they can be seen by the database
skypng = glob.glob('_S*.png')
summary_png = dn.strftime('%Y%m%d.png')
os.system("mv " + skypng[0] + " " + summary_stub + summary_png)
os.system("mv " + outname + " " + movie_stub)
d0 = central.localize(min(myLT))
d1 = central.localize(max(myLT))
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
    sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime,SummaryImage,SummaryMovie) VALUES(%d, %d, \"%s\", \"%s\", \"%sskymon/%s\", \"%smovies/%s\")' % (site_id, inst_id, startut, stoput, web_stub, summary_png, web_stub, outname)
    cur.execute(sql_cmd)

# Cleanup
skypng = glob.glob('_S*.png')
for fname in skypng:
    os.remove(fname)
    
os.remove('divx2pass.log')
