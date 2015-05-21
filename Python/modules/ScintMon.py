import numpy as np
import datetime as datetime
#import math as math

def parse_sum(fname):
    # Function to parse a .sum file created by lsum4
    # INPUT:
    #     fname - full path to the file to be parsed
    # OUTPUT:
    #     t - times of measurements
    #     s4 - s4 index of measurements
    #     rx_pos - the receiver ECEF position
    #     sx - the x positions of the satellites
    #     sy - the y positions of the satellites
    #     sz - the z position of the satellites
    #     site - the name of the site saved in the header
    # HISTORY
    #     Written by Jonathan J. Makela on 4 Feb 2014 based on GPS_ParseSum.m
    
    try:
        fid = open(fname,'r')
    except IOError as e:
        print file, 'does not exist'
        return [], [], [], [], [], [], []
    
    fid = open(fname,'r')
    data = []
    for line in fid:
        single_line = line.split()
        data.append(single_line)
    fid.close()
    
    N = len(data) # Number of lines in file
    
    site = data[0]
    out = ''
    for s in site:
        out = out + s + ' '
    site = out[0:-1]
    
    comment = data[1]
    rx_pos = np.array([float(data[2][0])*1e3, float(data[2][1])*1e3, float(data[2][2])*1e3])
    
    time_info = datetime.datetime(int(data[3][0]), int(data[3][1]), int(data[3][2]), int(data[3][3]), int(data[3][4]))

    maxPRN = 40
    l = 5
    t = []
    sx = []
    sy = []
    sz = []
    s4 = []
    while l < N:
        # Grab the time.  It extends past 23:59, so check this case and modify the date and mod the hour value as appropraite
        if int(data[l][0][0:2]) > 23:
            t.append(datetime.datetime(time_info.year, time_info.month, time_info.day, np.mod(int(data[l][0][0:2]),24), int(data[l][0][2:4]))+datetime.timedelta(days=1))
        else:
            t.append(datetime.datetime(time_info.year, time_info.month, time_info.day, int(data[l][0][0:2]), int(data[l][0][2:4])))
        num_sats = int(data[l][1])
        
        # Create dummy arrays of values
        temp_sx = np.empty((maxPRN,))
        temp_sy = np.empty((maxPRN,))
        temp_sz = np.empty((maxPRN,))
        temp_s4 = np.empty((maxPRN,))
        temp_sx[:] = np.nan
        temp_sy[:] = np.nan
        temp_sz[:] = np.nan
        temp_s4[:] = np.nan
    
        # Read in the data record
        for k in range(l+1,l+num_sats+1):
            PRN = int(data[k][0])-1
            temp_sx[PRN] = float(data[k][1])*1e3
            temp_sy[PRN] = float(data[k][2])*1e3
            temp_sz[PRN] = float(data[k][3])*1e3
            temp_s4[PRN] = float(data[k][7])
            
            l = l+1
        
        sx.append(temp_sx)
        sy.append(temp_sy)
        sz.append(temp_sz)
        s4.append(temp_s4)
        
        l = l+1
        
    sx = np.array(sx)
    sy = np.array(sy)
    sz = np.array(sz)
    s4 = np.array(s4)
    
    return t, s4, rx_pos, sx, sy, sz, site

def ecef2wgs84(xyz):
# Function to convert from ECEF to WGS84.
# INPUT:
#     xyz - the ECEF coordinate to convert [m]
# OUTPUT:
#     (lat, lon, alt) - the WGS84 LLA coordinate
# HISTORY:
#     Written by Jonathan J. Makela on 5 Feb 2014 based on MATLAB code

    a = 6378137;                  # semi-major axis of the earth [m]
    b = 6356752.3145;             # semi-minor axis of the earth [m]
    e = 0.08181919035596          # eccentricity of the earth = sqrt(1-(b^2)/(a^2))
    lat_accuracy_thresh = 1.57e-6 # 10 meter latitude accuracy
    
    # Parse the xyz coordinates sent in
    x = xyz[0]
    y = xyz[1];
    z = xyz[2];
    
    run = True;    # Flag used to stop iteration
    
    # Compute longitude
    lon = np.arctan2(y,x) * (180./np.pi);
    
    # Guess initial latitude ( assume you are on the surface (h=0) ) 
    p = np.sqrt(x**2 + y**2);
    lat0 = np.arctan(z/p * (1-e**2)**-1);  #Result in radians
    
    while run:
        # Use initial latitude to estimate N
        N = a**2 / np.sqrt(a**2 * (np.cos(lat0))**2 + b**2 * (np.sin(lat0))**2);
    
        # Estimate altitude
        h = p/np.cos(lat0) - N;
    
        # Estimate new latitude using new height
        lat1 = np.arctan(z/p * (1 - ((e**2*N)/(N+h)))**-1);
    
        if abs(lat1 - lat0) < lat_accuracy_thresh:
            run = False;
        
        # Replace our guess latitude with most recent estimate
        lat0 = lat1;
    
    lat = lat1 * (180./np.pi)
    
    return (lat, lon, h)

def elaz(rx_pos, sats):
# Function to calcuate the elevation and azimuth from one point to another
# INPUTS:
#     rx_pos - ECEF position of the point from which the el/az is calcuated [m]
#     sats - the ECEF position of the point to which the el/az is calculated [m].
#            This can be an [3,N] array of points
# OUTPUTS:
#     (el, az) - the elevation and azimuths
# HISTORY:
#     Written by Jonathan J. Makela on 5 Feb 2014 based on MATLAB code
    
    # Get the size of the incoming array
    (M,N) = np.shape(sats)
    
    # Convert the receiver location to WGS84
    (lat,lon,alt) = ecef2wgs84(rx_pos)
    
    # Create variables with the latitude and longitude in radians
    lat = lat * (np.pi/180.);
    lon = lon * (np.pi/180.);
    
    # Create the 3 x 3 transform matrix from ECEF to VEN
    VEN = np.array([ [np.cos(lat)*np.cos(lon),     np.cos(lat)*np.sin(lon),     np.sin(lat)],
            [        -np.sin(lon),              np.cos(lon),            0],
           [-np.sin(lat)*np.cos(lon),    -np.sin(lat)*np.sin(lon),     np.cos(lat)]]);
    
    # Replicate the Rx array to be the same size as the satellite array
    rx_array = np.tile(rx_pos,(N,1)).T;
    
    # Calculate the pseudorange for each satellite
    p = sats - rx_array;
    
    # Calculate the length of this vector
    n = np.sqrt(p[0,:]**2 + p[1,:]**2 + p[2,:]**2);
    
    # Create the normalized unit vector
    p = p / np.tile(n,(3,1));
    
    # Perform the transform of the normalized psueodrange from ECEF to VEN
    p_VEN = np.dot(VEN,p)
    
    # Calculate elevation and azimuth in degrees
    el = (np.pi/2 - np.arccos(p_VEN[0,:])) * 180./np.pi;
    az = np.arctan2(p_VEN[1,:],p_VEN[2,:]) * 180/np.pi;
    
    return (el, az)


