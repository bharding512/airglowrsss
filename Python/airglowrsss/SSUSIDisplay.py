from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
from scipy import mean
from pyglow import pyglow
from pyglow.pyglow import Point
from datetime import datetime, timedelta
import math
import numpy as np
from numpy import array as arr
from pyglow import coord as coord
from datetime import datetime, timedelta
from pyglow.pyglow import Line   
import multiprocessing
from multiprocessing import Pool
import time
import Image
import TifImagePlugin
from datetime import datetime, timedelta
import asiinfo
import ASI
from pyglow.pyglow import Point
import glob as glob
import matplotlib.pyplot as plt

def bin_SSUSI(ssusi_fname, width, height):
    # Read in the data
    f = netcdf.netcdf_file(ssusi_fname)
    
    # Pull variables from requested file
    t = f.variables['TIME']
    jd = f.variables['JULDAY']
    intensity = f.variables['DISK_RADIANCEDATA_INTENSITY']
    i1356 = intensity.data[:,:,:,2]
    error = f.variables['DISK_COUNTERROR_TOTAL']
    e1356 = error.data[:,:,:,2]
    angles = f.variables['DISK_SCAN_ANGLES'].data
    pp_alt = f.variables['PIERCEPOINT_NIGHT_ALTITUDE'].data
    pp_lat = f.variables['PIERCEPOINT_NIGHT_LATITUDE'].data#N, 132 eastwest lon, 16 updown lat
    pp_lon = f.variables['PIERCEPOINT_NIGHT_LONGITUDE'].data#N, 132, 16
    dmsp_coords_time = f.variables['DMSP_COORDS_TIME'].data
    starting_time = datetime.strptime(f.STARTING_TIME, '%Y%j%H%M%S')
    time = f.variables['TIME'].data
    angles = f.variables['DISK_SCAN_ANGLES'].data[0,:] # Why are there 3 in here?
    
    # Close the file
    f.close()

    # Used to correct for slant-path.  Assumes optically thin emission.
    cosan = np.cos(np.deg2rad(angles))**0.91 # Why 0.91?


    # cosan is just the angles for one strip.  Replicate this so it can easily be used against the i1356 value
    sz = np.shape(i1356)
    xx = np.tile(cosan,(sz[2],1))
    yy = xx.transpose()
    cosan = np.tile(yy,(sz[0],1,1))

    # The grid we will bin onto
    latedges = np.arange(-90,90,height)
    lonedges = np.arange(0,360,width)

    # Flatten arrays
    a = pp_lon
    lon = np.array(a.flatten())

    a = pp_lat
    lat = np.array(a.flatten())

    # Correct intensity by the slant-path factor
    data = np.array(i1356.flatten())*np.array(cosan.flatten())
    
    # Get times
    # Need to 
    time_data = np.array(time.flatten())
    time_data = np.repeat(time_data, np.shape(lat)[0]/np.shape(time_data)[0])

    # Count the number of pixels in each bin
    grid_N, xedges, yedges = np.histogram2d(lat,lon,bins=(latedges,lonedges))

    # Sum the intensities in each bin
    grid_bias_sum, xedges, yedges = np.histogram2d(lat,lon,weights=data,bins=(latedges,lonedges))

    # Calcualte the average intensity in the bin
    binned_i1356 = grid_bias_sum/grid_N

    # Calculate the average time in the bin
    time_bias_sum, exedges, yedges = np.histogram2d(lat,lon,weights=time_data, bins=(latedges,lonedges))

    binned_time = time_bias_sum/grid_N
    
    # The histograms above use the edges of the bins for each histogram bin.  We want to know the 
    # center of each bin for projection to apex coordinates
    lat_centers = 0.5 * np.diff(latedges) + latedges[:-1]
    lon_centers = 0.5 * np.diff(lonedges) + lonedges[:-1]
    
    binned_dn = [[0 for i in range(len(lon_centers))] for j in range(len(lat_centers))]#np.empty(shape(binned_times), dtype='datetime64[us]')
    for i in range(len(lat_centers)):
        for j in range(len(lon_centers)):
            if not np.isnan(binned_time[i][j]):
                binned_dn[i][j] = starting_time.replace(hour=0, minute=0, second=0, microsecond=0)+timedelta(0,binned_time[i][j])
            else:
                binned_dn[i][j] = np.nan
    
    return lat_centers, lon_centers, binned_i1356, binned_dn



def Line_star(a_b):
# Helper function to parallelize the Line function in pyglow
    if np.isfinite(a_b[1]) and np.isfinite(a_b[2]):
        pts = Line(a_b[0],a_b[1],a_b[2],a_b[4],target_ht=a_b[4]-1, step=a_b[5]);
        
        # Grab the individual points
        la = []
        lo = []
        al = []
        for p in pts:
            la.append(p.lat)
            lo.append(p.lon)
            al.append(p.alt)

        # Find the index of the maximum altitude, which is the location of the apex
        apex_i = np.argmax(al)

        # Wrap longitudes into the range of [0,360]
        lo = np.mod(np.array(lo),360.)

        # Save out the apex values
        apex_lat=la[apex_i]
        apex_lon=lo[apex_i]
        apex_ht=al[apex_i]
        apex_i1356=a_b[3]

        # Figure out if this is in the northern (+1) or southern hemisphere (-1)
        if np.sqrt((a_b[2]-lo[0])**2+(a_b[1]-la[0])**2) < np.sqrt((a_b[2]-lo[-1])**2+(a_b[1]-la[-1])**2):
            hemisphere = 1
        else:
            hemisphere = -1
    
        return (apex_lat, apex_lon, apex_ht, apex_i1356, hemisphere)
    else:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)

def project_SSUSI(ssusi_fname, (lat0, lat1, lon0, lon1), bin_resolution = 0.5, emission_ht = 350., trace_step = 10., num_processors=8):
# Loads in a SSUSI L1B file and projects the low-latitude data up to the magnetic equator (apex altitude).
#
# INPUTS:
#   ssusi_fname - the full path to the SSUSI file to read
#   (lat0, lat1, lon0, lon1) - bounding box of the area to project to the magnetic equator
#   bin_resolution - the resolution in degress of the latitude/longitude to bin data onto (default = 0.5 deg)
#   emission_ht - assumed emission height of the 1356 layer (default = 350 km)
#   trace_step - the step size to use while tracing the magnetic field lines (default = 10 km)
#   num_processors - the number of processors to use in the tracing function (default = 8)
#
# OUTPUTS:
#

    # Bin the SSUSI 1356 data onto a uniform grid
    lat_centers, lon_centers, binned_i1356, binned_dn = bin_SSUSI(ssusi_fname, bin_resolution, bin_resolution)
    
    # Find the indicies for the nearest points in the binned lat/lon grid
    ilat0 = np.abs(lat_centers-lat0).argmin()
    ilat1 = np.abs(lat_centers-lat1).argmin()
    ilon0 = np.abs(lon_centers-lon0).argmin()
    ilon1 = np.abs(lon_centers-lon1).argmin()

    # Perform the projection
    projection_dn = None
    trace_step = 25.

    job_args = []
    for i in range(ilat0,ilat1):
        for j in range(ilon0,ilon1):
            if not np.isnan(binned_i1356[i][j]):
                if projection_dn is None:
                    # This is the first projection.  Set the time to use for all projections
                    projection_dn = binned_dn[i][j].replace(hour=0, minute=0, second=0, microsecond=0)

                job_args.append((projection_dn,lat_centers[i],lon_centers[j],binned_i1356[i][j],emission_ht,trace_step))
            else:
                job_args.append((projection_dn,np.nan,np.nan,np.nan,emission_ht,trace_step))
                
    # Trace along field lines
    N = multiprocessing.cpu_count()

    # Create the pool.  Be nice.  Don't use all the cores!
    pool = Pool(processes=min(N/4, num_processors))

    results = pool.map(Line_star,job_args)

    pool.close()
    pool.join()
                                                 
    # Unpack results
    la, lo, ht, i1356, hemi = zip(*results)
    apex_i1356 = np.array(i1356).reshape((ilat1-ilat0, ilon1-ilon0))
    apex_lat = np.array(la).reshape((ilat1-ilat0, ilon1-ilon0))
    apex_lon = np.array(lo).reshape((ilat1-ilat0, ilon1-ilon0))
    apex_ht = np.array(ht).reshape((ilat1-ilat0, ilon1-ilon0))
    hemisphere = np.array(hemi).reshape((ilat1-ilat0, ilon1-ilon0))
                    
    return lat_centers, lon_centers, binned_i1356, binned_dn, apex_lat, apex_lon, apex_ht, apex_i1356, hemisphere



def get_cnfi_dn(files, dn):
    # Array to store times
    dt = []
    for fn in files:
        d = Image.open(fn)
        dt.append(d.info['UniversalTime'])
        
    # Find the difference between the times in all images and the requested time
    dt = np.array(dt)
    del_t = dt - dn
    helper = np.vectorize(lambda x: np.abs(x.total_seconds()))
    dt_sec = helper(del_t)
    
    # The minimum time difference
    a = dt_sec.min()
    i = dt_sec.argmin()   
    
    return a,i

# Specify the filename to read 
# Location is in the following format
# /rdata/airglow/imaging/<instr_name>/<site_name>/<year>/<day number of year>/<file: 6300_xxxx.tif>
# (instr_name is either casi01 or cnfi01, not sure which 
#  one he wants you to use. probably cnfi)
def get_cnfi_image(dn, ht=250):
    doy = dn.timetuple().tm_yday
    year = dn.timetuple().tm_year
    
    #print doy
    #print year
    
    # Grab all files (need to extend to search doy-1 to doy+1)
    files = glob.glob('/rdata/airglow/imaging/cnfi01/hka/%i/%03i/6300_*.tif' % (year, doy))
    files.sort()
    
    a,i = get_cnfi_dn(files, dn)
    
    # If time difference is too big, try previous/next day
    if a > 30*60.:
        if a > 0:
            # Try previous day
            files = glob.glob('/rdata/airglow/imaging/cnfi01/hka/%i/%03i/6300_*.tif' % (year, doy-1))
            files.sort()
        else:
            # Try next day
            files = glob.glob('/rdata/airglow/imaging/cnfi01/hka/%i/%03i/6300_*.tif' % (year, doy+1))
            files.sort()
        
        a,i = get_cnfi_dn(files, dn)

        # Exit out if nothing is within 30 minutes
        if a > 30*60.:
            print a
            return np.nan, np.nan, np.nan, np.nan
        
    # Load the closest image
    fn = files[i]
    d = Image.open(fn)
    
    im = np.asarray(d)
    dt = d.info['UniversalTime']
    print 'SSUSI:', dn, 'CNFI:', dt
    
    # Load appropriate calibration image
    instr_info = asiinfo.get_instr_info('cnfi01', dt)
    cal_file = np.load(instr_info['cal_file'])
    az = cal_file['az']
    el = cal_file['el']
    rx_lat = cal_file['rx_lat']
    rx_lon = cal_file['rx_lon']
    lat, lon = ASI.ConvertAzEl2LatLon(az, el, ht, rx_lat, rx_lon)
    
    apex_ht, apex_lon = load_cnfi_proj()
    
    return lat, lon, apex_ht, apex_lon, im, dt

def load_cnfi_proj():
    apex_ht = np.loadtxt('/usr/local/share/airglowrsss/Python/notebooks/ugrad/mhoch2/apex_ht.txt')
    apex_lon = np.loadtxt('/usr/local/share/airglowrsss/Python/notebooks/ugrad/mhoch2/apex_lon.txt')
    
    return apex_ht, apex_lon

def create_plots(target_lat, target_lon,lat_centers, lon_centers, binned_i1356, binned_dn, apex_lat, apex_lon, apex_ht, apex_i1356, hemisphere, emission_ht = 350., SSUSIbubbles = None):
    # Cast longitudes into [-180 180] range to match how Basemap works
    lon_centers[np.where(lon_centers > 170)] -= 360.
    apex_lon[np.where(apex_lon > 170)] -= 360.
    if target_lon > 170:
        target_lon -= 360.

   # Find the indicies for the nearest points in the binned lat/lon grid
    ilat = np.abs(lat_centers-target_lat).argmin()
    ilon = np.abs(lon_centers-target_lon).argmin()

    ssusi_north_dn = binned_dn[ilat][ilon]

    # Southern time
    pts = Line(ssusi_north_dn,target_lat,target_lon,emission_ht,target_ht=emission_ht-1,step=25.)
    target_lat = pts[-1].lat
    target_lon = pts[-1].lon
#    if(target_lon < 0): target_lon+=360.
    if(target_lon > 170.): target_lon-=360.

    grid_lon, grid_lat = np.meshgrid(lon_centers,lat_centers)

    X = np.ma.masked_where(np.isnan(binned_i1356),grid_lon)
    Y = np.ma.masked_where(np.isnan(binned_i1356),grid_lat)

    dist = np.sqrt((X-target_lon)**2+(Y-target_lat)**2)
    idx = np.argmin(dist)
    ssusi_south_dn = np.array(binned_dn).flatten()[idx]

    map_lat = 30

    map_min_lon = np.nanmin(apex_lon.flatten())
    map_max_lon = np.nanmax(apex_lon.flatten())

    ##
    plt.figure()
    m = Basemap(llcrnrlon=map_min_lon-0.5,
                llcrnrlat=-map_lat-0.5,
                urcrnrlon=map_max_lon+0.5,
                urcrnrlat=map_lat+0.5,
                projection='merc',
                area_thresh=1000,
                resolution='i')

    lon_grid, lat_grid = np.meshgrid(lon_centers,lat_centers)
    xpt, ypt = m(lon_grid,lat_grid)

    I = np.ma.masked_where(np.isnan(binned_i1356),binned_i1356)
    
    m.pcolormesh(xpt,ypt,I,vmin=0,vmax=200)

    m.drawparallels(np.arange(-80.,81.,20.),labels=[True,False,False,False])
    m.drawmeridians(np.arange(-180.,181.,20.),labels=[False,False,False,True])
    plt.title('SSUSI 135.6-nm binned data')
    c = plt.colorbar();
    c.set_label('[R]')

    # draw coastlines and fill the continents (alpha value is set to partial transparancy because for some
    # reason the fill is over top of the quivers used to denote the wind vectors)
    m.drawcoastlines(color=[.6,.6,.6])
    m.drawstates(color=[.6,.6,.6])
    m.drawcountries(color=[.6,.6,.6])

    # plot the magnetic equator
    fid = open('/home/jmakela/MagEq.txt')
    m_lon = []
    m_lat = []
    for line in fid:
        temp = line.split()
        m_lon.append(float(temp[2]))
        m_lat.append(float(temp[3]))
    fid.close()

    m_lon = np.array(m_lon)
    m_lat = np.array(m_lat)

    xpt,ypt = m(m_lon,m_lat)

    pts = np.array(sorted(zip(xpt,ypt)))

    m.plot(pts[:,0],pts[:,1],'k',linewidth=2)
    
    apex_lon = np.mod(apex_lon,360.)
        
    ##
    plt.figure()
    I = np.ma.masked_where(hemisphere != 1,apex_i1356)
    X = np.ma.masked_where(hemisphere != 1,apex_lon)
    Y = np.ma.masked_where(hemisphere != 1,apex_ht)
    plt.pcolormesh(X,Y,I,vmin=0,vmax=200)
    plt.xlabel('Geog. Longitude [deg]')
    plt.ylabel('Apex Altitude [km]')
    plt.title('SSUSI 135.6 (north)')
    c = plt.colorbar();
    c.set_label('[R]')
    plt.xlim([np.nanmin(apex_lon.flatten()),np.nanmax(apex_lon.flatten())])
    plt.ylim([200,1400]);

    # If we are asked to plot SSUSI bubble locations, do it
    if SSUSIbubbles is not None:
        apex_bub = []

        print np.shape(SSUSIbubbles)
        for b in SSUSIbubbles:
            lon = b[0]
            lat = b[1]
            alt = b[2]
            ndep = b[3]

            lo = []
            la = []
            al = []
            if alt > 360.:

                # Trace the field
                pts = Line(ssusi_north_dn,lat,np.mod(lon,360),alt,target_ht=alt-1,step=25.)

                for p in pts:
                    lo.append(np.mod(p.lon,360.)) # Force to 0...360
                    la.append(p.lat)
                    al.append(p.alt)

                i = np.argmax(al)
                apex_bub.append([lo[i], la[i], al[i], ndep])

        apex_bub = np.array(apex_bub)
        print np.shape(apex_bub)

        plt.scatter(apex_bub[:,0],apex_bub[:,2],c=apex_bub[:,3],marker='.',linewidth=0,s=25)

    ##
    plt.figure()
    I = np.ma.masked_where(hemisphere != -1,apex_i1356)
    X = np.ma.masked_where(hemisphere != -1,apex_lon)
    Y = np.ma.masked_where(hemisphere != -1,apex_ht)
    plt.pcolormesh(X,Y,I,vmin=0,vmax=200)
    plt.xlabel('Geog. Longitude [deg]')
    plt.ylabel('Apex Altitude [km]')
    plt.title('SSUSI 135.6 (south)')
    c = plt.colorbar();
    c.set_label('[R]')
    plt.xlim([np.nanmin(apex_lon.flatten()),np.nanmax(apex_lon.flatten())])
    plt.ylim([200,1400]);
    
#    # Based on http://matplotlib.1069221.n5.nabble.com/colormap-from-blue-to-transparent-to-red-td26440.html
#    theCM = plt.cm.get_cmap('gray')
#    theCM._init()
#    #alphas = linspace(0,0.3,259)
#    alphas = ones([259])
#    theCM._lut[:,-1] = alphas
    
    # Get the closest CNFI image to the north pass

    plt.figure()
    cnfi_lat,cnfi_lon,cnfi_apex_ht,cnfi_apex_lon,im,cnfi_dn = get_cnfi_image(ssusi_north_dn)
    #cnfi_apex_lon[np.where(cnfi_apex_lon > 170)] -= 360.
    #apex_lon[np.where(apex_lon > 170)] -= 360.
   
    im_min = im.mean()-im.std()
    im_max = im.mean()+2*im.std()

    I = np.ma.masked_where(hemisphere != 1,apex_i1356)
    X = np.ma.masked_where(hemisphere != 1,apex_lon)
    Y = np.ma.masked_where(hemisphere != 1,apex_ht)
    plt.pcolormesh(X,Y,I,vmin=0,vmax=200);
    c = plt.colorbar();
    c.set_label('[R]')

    I = np.ma.masked_where(np.isnan(im),im)
    X = np.ma.masked_where(np.isnan(cnfi_apex_lon),cnfi_apex_lon)
    Y = np.ma.masked_where(np.isnan(cnfi_apex_ht),cnfi_apex_ht)
    plt.pcolormesh(X,Y,I,cmap='gray',vmin=im_min,vmax=im_max);

    plt.xlabel('Geog. Longitude [deg]')
    plt.ylabel('Apex Altitude [km]')
    plt.title('SSUSI (north): %s UT\n CNFI: %s UT' % (ssusi_north_dn.strftime('%d %b, %Y %H:%M:%S'),cnfi_dn.strftime('%d %b, %Y %H:%M:%S')))
    plt.xlim([np.nanmin(apex_lon.flatten()),np.nanmax(apex_lon.flatten())])
    plt.ylim([200,1400]);
    
    if SSUSIbubbles is not None:
        plt.scatter(apex_bub[:,0],apex_bub[:,2],c=apex_bub[:,3],marker='.',linewidth=0,s=25)

    # Get the closest CNFI image to the south pass
    plt.figure()
    cnfi_lat,cnfi_lon,cnfi_apex_ht,cnfi_apex_lon,im,cnfi_dn = get_cnfi_image(ssusi_south_dn)
    #cnfi_apex_lon[np.where(cnfi_apex_lon > 170)] -= 360.

    im_min = im.mean()-im.std()
    im_max = im.mean()+2*im.std()

    I = np.ma.masked_where(hemisphere != -1,apex_i1356)
    X = np.ma.masked_where(hemisphere != -1,apex_lon)
    Y = np.ma.masked_where(hemisphere != -1,apex_ht)
    plt.pcolormesh(X,Y,I,vmin=0,vmax=200);
    c = plt.colorbar();
    c.set_label('[R]')

    I = np.ma.masked_where(np.isnan(im),im)
    X = np.ma.masked_where(np.isnan(cnfi_apex_lon),cnfi_apex_lon)
    Y = np.ma.masked_where(np.isnan(cnfi_apex_ht),cnfi_apex_ht)
    plt.pcolormesh(X,Y,I,cmap='gray',vmin=im_min,vmax=im_max);

    plt.xlabel('Geog. Longitude [deg]')
    plt.ylabel('Apex Altitude [km]')
    plt.title('SSUSI (south): %s UT\n CNFI: %s UT' % (ssusi_south_dn.strftime('%d %b, %Y %H:%M:%S'),cnfi_dn.strftime('%d %b, %Y %H:%M:%S')))
    plt.xlim([np.nanmin(apex_lon.flatten()),np.nanmax(apex_lon.flatten())])
    plt.ylim([200,1400]);

    # If we are asked to plot SSUSI bubble locations, do it
    if SSUSIbubbles is not None:
        plt.scatter(apex_bub[:,0],apex_bub[:,2],c=apex_bub[:,3],marker='.',linewidth=0,s=25)
