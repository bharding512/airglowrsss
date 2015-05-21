#
# Script to reprocess the UAO NATION wind and temperature files
# and create a plot that incorporates the Boltwood cloud sensor
# data
#
# Run at reprocess_UAO_plots.py -y YEAR -d DOY
#
# If not options given, run for yesterday's data
#
# Written by Jonathan Makela (jmakela@illinois.edu) on 18 Oct 2012

import matplotlib
matplotlib.use('AGG')

import scipy.io
import scipy.interpolate as sp
import matplotlib.pyplot as plt
import BoltwoodSensor
from matplotlib.font_manager import FontProperties
import glob
import numpy
import math
import os

import datetime
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sys
sys.path.append("../modules")

# Directory containing the MATLAB save files to work with
path_matfiles = '/mnt/FPIData/Results/'
path_templogs = '/mnt/FPIData/NATION/uao/TempLogs/'
path_png = '/data/NATION/natuaopngs/'

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
dn1 = datetime.datetime(year,1,1) + datetime.timedelta(days=doy)

# Find the requested MATLAB save file
mat_files = glob.glob(path_matfiles + 'NATUAO*_%04d%03d.mat' % (year, doy))

for mat_file in mat_files:
	mat = scipy.io.loadmat(mat_file)
	FPIResults = mat['FPIResults']
	time = mat['t'].flatten()

	# call modules from BoltWood for the two days required and combine the data
	bw_date1, bw_sky1, bw_amb1 = BoltwoodSensor.ReadTempLog(path_templogs + "Cloud_UAO_" + dn.strftime("%Y%m%d") + ".txt")
	bw_date2, bw_sky2, bw_amb2 = BoltwoodSensor.ReadTempLog(path_templogs + "Cloud_UAO_" + dn1.strftime("%Y%m%d") + ".txt")
	bw_date = numpy.hstack((bw_date1,bw_date2))
	bw_sky = numpy.hstack((bw_sky1, bw_sky2))
	bw_amb = numpy.hstack((bw_amb1, bw_amb2))

	# Files to save resultant plots to
	temp_fname = mat_file[len(path_matfiles):-4] + '_temperatures.png'
	wind_fname = mat_file[len(path_matfiles):-4] + '_winds.png'

	# Color scheme to use
	colors = ['b','r','g','y','m']
	
	# Empty array to store directions used
	direction = []

	# Range of sky temperatures considered for fully clear and fully cloudy
	full_clear = -30
	full_cloud = -15

	# Plotting setup
	fontP = FontProperties()
	fontP.set_size('small')

	Doppler_Fig = plt.figure()
	Doppler_Graph = Doppler_Fig.add_subplot(111)
	Doppler_Graph.hold(True)

	Temperature_Fig = plt.figure()
	Temperature_Graph = Temperature_Fig.add_subplot(111)
	Temperature_Graph.hold(True)

	# Loop through the directions
	for i in range(0,5):
	    # Read values from the matlab FPIResults structure
	    direction.append(str(FPIResults[:,i][0][0])[3:-2])
	    t = FPIResults[:,i][0][8].flatten()
	    sign = FPIResults[:,i][0][3].flatten()
	    ze = FPIResults[:,i][0][1].flatten()
	    Doppler_Laser = FPIResults[:,i][0][22].flatten()
	    e_Doppler_Laser = FPIResults[:,i][0][26].flatten()
	    Temperature = FPIResults[:,i][0][24].flatten()
	    e_Temperature = FPIResults[:,i][0][28].flatten()
    
	    # Determine if this is zenith or horiztonal measurement
	    if i == 0:
	        # This is the zenith measurement.  Setup the interpolation so that the
	        # horizontal look directions can compensate for the vertical wind
	        w = Doppler_Laser
	        time_w = t
	        w_int = sp.interp1d(time_w,w,kind='linear',bounds_error = False, fill_value = 0)
	        horizontal_wind = Doppler_Laser*sign
	    else:
	        # Use the zenith measurement to correct for vertical wind
	        my_w = w_int(t)
	        horizontal_wind = (Doppler_Laser-my_w*math.cos(abs(ze*math.pi/180)))*sign/math.sin(abs(ze*math.pi/180))
	    
	    # Initial plot.  For winds, plot a light line.  For both winds and temperatures,
	    # plot a single point at a bogus value (-999,-999) to set the legent correctly
	    Doppler_Graph.plot(t,horizontal_wind,color=colors[i],alpha=0.3,marker=None)
	    Doppler_Graph.plot(-999,-999,color=colors[i],marker='o',label=direction[i])
    
	    Temperature_Graph.plot(-999,-999,color=colors[i],marker='o',label=direction[i])
    
	    # Loop through each point (needed to be done this way to allow setting the 
	    # alpha value for each individual point
	    for (x,y,ey,z,ez) in zip(t,horizontal_wind,e_Doppler_Laser,Temperature,e_Temperature):
	        # Convert from matlab datenum to datetime
	        trash, d = divmod(x,1)
	        tt = datetime.datetime.fromordinal(int(x))-datetime.timedelta(days = 366)
	        mt = tt + datetime.timedelta(seconds = d*86400)
	        
	        # Find the time in the cloud data that corresponds to this point
	        closest = sorted(bw_date,key=lambda d:abs(mt-d))[0]
	        sky_temp = bw_sky[(bw_date == closest).nonzero()][0]
        
	        # Calculate the alpha value to be used for this point.  Scaled linearly between
	        # [.1,1] corresponding to a sky temp range of [full_cloud, full_clear]
	        aval = 1-max(0,min(((sky_temp-full_clear)/(full_cloud-full_clear)),.9))
	        
	        # Plot the points and error bars on the graphs
	        Doppler_Graph.scatter(x,y,alpha=aval,color=colors[i])
	        Doppler_Graph.plot([x,x],[y-ey,y+ey],alpha=aval,color=colors[i])
        
	        Temperature_Graph.scatter(x,z,alpha=aval,color=colors[i])
	        Temperature_Graph.plot([x,x],[z-ez,z+ez],alpha=aval,color=colors[i])

	# Add a reference line at 0 m/s and format labels        
	Doppler_Graph.plot(time,numpy.zeros(len(time)),'k--')
	Doppler_Graph.set_ylim([-200,200])
	Doppler_Graph.set_xlim([time[0],time[-1]])
	Doppler_Graph.legend(ncol=5,prop=fontP)
	Doppler_Graph.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
	Doppler_Graph.xaxis_date()
	Doppler_Graph.axes.set_xlabel('LT')
	Doppler_Graph.axes.set_ylabel('m/s')

	# Format labels
	Temperature_Graph.set_ylim([400,1400])
	Temperature_Graph.set_xlim([time[0],time[-1]])
	Temperature_Graph.legend(ncol=5,prop=fontP)
	Temperature_Graph.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
	Temperature_Graph.xaxis_date()
	Temperature_Graph.axes.set_xlabel('LT')
	Temperature_Graph.axes.set_ylabel('K')

	# Get the start and end time for the title
	i, d = divmod(time[0],1)
	tt = datetime.datetime.fromordinal(int(i))-datetime.timedelta(days = 366)
	t0 = tt + datetime.timedelta(seconds = d*86400)

	i, d = divmod(time[-1],1)
	tt = datetime.datetime.fromordinal(int(i))-datetime.timedelta(days = 366)
	t1 = tt + datetime.timedelta(seconds = d*86400)

	if t0.day != t1.day:
	    Doppler_Graph.axes.set_title("UAONATION Winds: %02d-%02d %s" % (t0.day, t1.day, t1.strftime('%b %Y')))
	    Temperature_Graph.axes.set_title("UAONATION Temperatures: %02d-%02d %s" % (t0.day, t1.day, t1.strftime('%b %Y')))
	else:
	    Doppler_Graph.axes.set_title("UAONATION Winds: %s" % (t1.strftime('%d %b %Y')))
	    Temperature_Graph.axes.set_title("UAONATION Tempratures: %s" % (t1.strftime('%d %b %Y')))

	# Move the original png files and then save the new png files
	os.system('mv %s%s %sorig' % (path_png, temp_fname, path_png))
	os.system('mv %s%s %sorig' % (path_png, wind_fname, path_png))

	Temperature_Fig.savefig(path_png + temp_fname)
	Doppler_Fig.savefig(path_png + wind_fname)
