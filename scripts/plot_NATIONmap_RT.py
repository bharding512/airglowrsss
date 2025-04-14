import matplotlib
matplotlib.use('AGG')

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import pytz
from pytz import timezone
import os
import BoltwoodSensor as bs
import glob as glob
import datetime

def convert_azze_latlon(lat, lon, az, ze, ht=250):
    # Converts a sites az/ze look directions to lat/lon pierce points at the requested altitude
    #
    # INPUTS: site - site dictionary to use in the conversion
    # OPTIONAL INPUTS: ht - altitude to perform the conversion (default = 250 km)
    #
    # HISTORY: Written on 4 Jan 2013 by Jonathan J. Makela based on MATLAB code
    
    # Define the radius of the Earth in km
    Re = 6371.2
    
    if ze < 0:
        ze = abs(ze)
        az = az + 180.
    
    # Convert angles and locations to radian
    el_r = (90.0 - ze) * np.pi/180.
    az_r = np.mod(az * np.pi/180.,2*np.pi)
    lat_r = lat * np.pi/180.
    lon_r = lon * np.pi/180.
    
    # Calculate the differential angle, alpha
    temp = np.cos(el_r)/(1+(ht/Re))
    alpha = np.arccos(temp) - el_r
    
    # Calculate the pierce point latitude
    temp = np.sin(lat_r) * np.cos(alpha) + np.cos(lat_r)*np.cos(az_r)*np.sin(alpha)
    lat_r = np.arcsin(temp)
    
    # Calculate the pierce point longitude
    temp = np.sin(alpha) * np.sin(az_r) / np.cos(lat_r)
    lon_r = np.arcsin(temp) + lon_r
    
    # Convert back to degrees
    return lat_r * 180.0/np.pi, lon_r * 180.0/np.pi

sites = {'UAO': {'lat': 40.0985,
                 'lon': -88.2291,
                 'color': 'b',
                 'tz': 'US/Central'}}

sites['PAR'] = {'lat': 35.2040,
                 'lon': -82.8792,
                 'color': 'm',
                 'tz': 'US/Eastern'}

sites['ANN'] = {'lat': 42.4028,
               'lon': -83.9247,
                 'color': 'g',
               'tz': 'US/Eastern'}

sites['EKU'] = {'lat': 37.7267,
                'lon': -84.3001,
                 'color': 'r',
                'tz': 'US/Eastern'}

cloud_temp = -25.0
local = pytz.timezone('US/Central')
utc = pytz.timezone('UTC')

plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')

# setup Lambert Conformal basemap.
m = Basemap(llcrnrlon=-95.,llcrnrlat=30.,urcrnrlon=-75.,urcrnrlat=45.,
            projection='mill',resolution='i')
#            projection='merc', area_thresh=1000,
#            resolution='i',lat_1=45.,lat_2=55,lat_0=40,lon_0=-85.)

# draw coastlines and fill the continents (alpha value is set to partial transparancy because for some
# reason the fill is over top the quivers used to denote the wind vectors
m.drawcoastlines()
m.drawstates()
m.fillcontinents(color='0.8')
m.drawparallels(np.arange(30,45,5),labels=[1,1,0,0])
m.drawmeridians(np.arange(-95,-75,5),labels=[0,0,0,1])

# draw coastlines and fill the continents (alpha value is set to partial transparancy because for some
# reason the fill is over top the quivers used to denote the wind vectors
m.drawcoastlines()
m.drawstates()

m.plot(0,0,'g^',label='Site Online')
m.plot(0,0,'r^',label='Site Cloudy')
m.plot(0,0,'y^',label='No Cloud Data')
m.plot(0,0,'rx',label='Site Offline')
    
for s in sites:
    # Parse the cloud data
    dns, sky_temp, amb_temp = bs.ReadRawTempLog('/home/master/NATION/Control/'+s+'Cloud.txt', sites[s]['tz'])

    # convert to map projection coords.
    # Note that lon,lat can be scalars, lists or numpy arrays.
    xpt,ypt = m(sites[s]['lon'],sites[s]['lat'])
    
    # put some text next to the dot, offset a little bit
    # (the offset is in map projection coordinates)
    plt.text(xpt+50000,ypt+50000,'%s' % (s))
    
    # Parse the LastImage file
    f = glob.glob('/home/master/NATION/Control/'+ s + 'Control.txt')
    
    try:
        fid = open(f[0],'r')
    except IOError as e:
        print f, 'does not exist'
        
    data = []
    for line in fid:
        single_line = line.split()
        data.append(single_line)
        
    fid.close()
    
    # Convert the ze/az to lat lon
    lat, lon = convert_azze_latlon(sites[s]['lat'], sites[s]['lon'], float(data[3][1]), float(data[2][1]), ht=250)
    xpt_obs,ypt_obs = m(lon,lat)

    print s, lat, lon

    # Check the time
    tmp = data[1][1].replace(':','.').split('.')
    a = datetime.datetime(int(tmp[0]),int(tmp[1]),int(tmp[2]),int(tmp[3]),int(tmp[4]),int(tmp[5]))

    # convert to map projection coords.
    # Note that lon,lat can be scalars, lists or numpy arrays.
    xpt,ypt = m(sites[s]['lon'],sites[s]['lat'])

    if abs(dns - local.localize(datetime.datetime.now())) > datetime.timedelta(minutes = 5):
        valid_cloud = False
    else:
        valid_cloud = True

    if abs(utc.localize(a) - local.localize(datetime.datetime.now())) > datetime.timedelta(minutes = 15):
        valid_obs = False
    else:
        valid_obs = True

    # Plot the site locations
#    if abs(utc.localize(a) - local.localize(datetime.datetime.now())) > datetime.timedelta(minutes = 15):
    if (not valid_obs):
        m.plot(xpt,ypt,'rx')
    elif (not valid_cloud):
        m.plot(xpt,ypt,'y^')  # yellow means no cloud data for site
    elif sky_temp < -25.0:
        m.plot(xpt,ypt,'g^')  # plot a blue dot for valid site with cloud data
    else:
        m.plot(xpt,ypt,'r^')  # plot a red dot for clouded over site

    # put some text next to the dot, offset a little bit
    # (the offset is in map projection coordinates)
    plt.text(xpt+50000,ypt+50000,'%s' % (s))

    # Plot the observation locations
#    if abs(utc.localize(a) - local.localize(datetime.datetime.now())) > datetime.timedelta(minutes = 15):
    if (not valid_obs):
        m.plot(xpt_obs,ypt_obs,color=sites[s]['color'],marker='o',alpha=0.25)
    else:
        m.plot(xpt_obs,ypt_obs,color=sites[s]['color'],marker='o',alpha=0.75)
    m.plot(0,0,color=sites[s]['color'],marker='o',label=s)

    plt.legend(numpoints=1,loc=2,prop={'size':8})
    
    plt.title(local.localize(datetime.datetime.now()).astimezone(utc).strftime('NATION %Y-%m-%d %H:%M:%S UT'))

    plt.savefig('/usr/share/drupal7/images/NATIONmap_RT.png')
