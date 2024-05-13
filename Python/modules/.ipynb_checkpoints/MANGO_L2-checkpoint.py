import matplotlib.pyplot as plt
import xarray as xr
import h5py
import pandas as pd
import pytz
from scipy.signal import medfilt2d
import subprocess
import fpiinfo
import BoltwoodSensor
import asiinfo
import ephem
import gabor
import multiprocessing as mp
import numpy.fft as fft
from skimage.feature import peak_local_max
import numpy as np
import datetime

def runcmd(cmd, verbose = False, *args, **kwargs):
# Runs a commond line command (useful for downloading data)
# From https://www.scrapingbee.com/blog/python-wget/

    process = subprocess.Popen(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        text = True,
        shell = True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)
    pass

def load_hdf5_to_xarray(file_path, kernel_size=7):
# Function to load data from HDF5 file into xarray DataArrays
# 
# Input:
#   file_path: complete path to an MANGO HDF5 file
#   kernel_size: (optional) size of 2D median filter to apply to each image
#
# Output:
#   ds: xarray data structure containing information from the HDF5 file
#       The data structure has the following format:
#       Coordinates:
#           time (datetime64) - time stamps of each image in the data set
#           east (float) - east km distance from the center of the image. This is a uniform step.
#           north (float) - north km distance from the center of the image. This is a uniform step.
#       Data Variable:
#           integration_time (float) - integration time for each image in sec (dim: time)
#           ImageData (int) - intensity values in the image (dims: time, north, east)
#                             the image data has had a 7x7 median filter applied to remove stars
#           Azimuth (float) - azimuth for each pixel, in deg (dims: north, east)
#           Elevation (float) - elevation for each pixel, in deg (dims: north, east)
#           Latitude (float) - latitude for each pixel, in deg N (dims: north, east)
#           Longitude (float) - longtitude for each pixel, in deg E (dims: north, east)
#           CCDTemperature (float) - CCD temperature, in C (dim: time)

    # Read in the HDF5 file into a temporary array
    data_arrays = {}
    with h5py.File(file_path, 'r') as hdf:
        for key in hdf.keys():
            # Check if the item is a dataset
            if isinstance(hdf[key], h5py.Dataset):
                # Load the dataset into an xarray DataArray
                data_arrays[key] = xr.DataArray(hdf[key][...], name=key)

    # Convert Unix times to datetime objects
    unix_times_start = data_arrays['UnixTime'][0].values
    start_times = pd.to_datetime(unix_times_start, unit='s', utc=True)
    
    unix_times_stop = data_arrays['UnixTime'][1].values
    stop_times = pd.to_datetime(unix_times_stop, unit='s', utc=True)
    
    # Calculate integration time and ensure it's a float
    integration_time_values = (stop_times - start_times).total_seconds()
    integration_time_values = integration_time_values.astype(float)
    
    # Ensure 'time' is numpy.datetime64
    start_times_np = np.array(start_times, dtype='datetime64[ns]')

    # Rename the dimensions of the ImageData
    image_data = data_arrays['ImageData'].rename({'dim_0': 'time', 'dim_1': 'north', 'dim_2': 'east'})
    image_data.attrs['description'] = 'Pixel values of image at each time step'

    # Process other variables
    azimuth = data_arrays['Azimuth'].rename({'dim_0': 'north', 'dim_1': 'east'})
    azimuth.attrs['units'] = 'Degrees'
    azimuth.attrs['description'] = 'Azimuth of each pixel'
    elevation = data_arrays['Elevation'].rename({'dim_0': 'north', 'dim_1': 'east'})
    elevation.attrs['units'] = 'Degrees'
    elevation.attrs['description'] = 'Elevation of each pixel'
    latitude = data_arrays['Latitude'].rename({'dim_0': 'north', 'dim_1': 'east'})
    latitude.attrs['units'] = 'Degrees N'
    latitude.attrs['description'] = 'Geodetic latitude of each pixel'
    longitude = data_arrays['Longitude'].rename({'dim_0': 'north', 'dim_1': 'east'})
    longitude.attrs['units'] = 'Degrees E'
    longitude.attrs['description'] = 'Geodetic longitude of each pixel'
    CCDTemperature = data_arrays['CCDTemperature'].rename({'dim_0': 'time'})
    CCDTemperature.attrs['units'] = 'Celcius'
    CCDTemperature.attrs['description'] = 'Temperature of CCD'

    north_array = xr.DataArray(data_arrays['PixelCoordinates'][0][0,:], dims=('north'))
    north_array.attrs['units'] = 'km'
    north_array.attrs['description'] = 'Distance in north direction from site at airglow layer'
    east_array = xr.DataArray(data_arrays['PixelCoordinates'][1][:,0], dims=('east'))
    east_array.attrs['units'] = 'km'
    east_array.attrs['description'] = 'Distance in east direction from site at airglow layer'

    # Apply median filter to each image
    im_ds = []
    for id in image_data:
        temp = xr.DataArray(medfilt2d(id,kernel_size=kernel_size), dims=('north','east'))
        im_ds.append(temp)

    combined = xr.concat(im_ds, dim=('time'))
    combined.attrs['description'] = 'Pixel values of image at each time step after 7x7 median filter'

    integration_time = xr.DataArray(integration_time_values, coords={'time': start_times_np}, dims=('time'))
    integration_time.attrs['units'] = 's'
    integration_time.attrs['description'] = 'Integration time of the image'

    # Create xarray Dataset
    ds = xr.Dataset({'integration_time': integration_time,
                    'ImageData': combined,
                    'Azimuth': azimuth,
                    'Elevation': elevation,
                    'Latitude': latitude,
                    'Longitude': longitude,
                    'CCDTemperature': CCDTemperature,
                    'east': east_array,
                    'north': north_array})
    
    return ds

def load_cloud_and_moon(analysis_parameters, times):
# Load in the cloud data, if it exists, and calculate the moon angle for requested times.
#
# INPUTS:
#   analysis_parameters - dictionary containing at least the "site" entry with the site abbreviation to process and "date" entry with the datenum to process
#   times - list of times, generally corresponding to the times of images for which we want to check if the moon is up
#
# OUTPUTS:
#   results_df: Pandas dataframe structure containing times, cloud temperature, and moon angle
#       The data structure has the following format:
#           UT (datetime64) - universal time of each estimate
#           Cloud Temperature (float) - temperature from the cloud sensor in C
#           Moon Altitude (float) - angle above the horizon of the moon (deg)
    
    if analysis_parameters['site'] in ['cvo','low','blo']:
        my_site = fpiinfo.get_site_info(analysis_parameters['site'],analysis_parameters['date'])

        # Because of timezones and how the cloud senor data files are broken up, 
        # read on day earlier and later just to be sure we capture the data needed
        dns0, sky_temp0, amb_temp0 = BoltwoodSensor.ReadTempLog('/rdata/airglow/templogs/cloudsensor/%s/Cloud_%s_%s.txt' % (analysis_parameters['site'], analysis_parameters['site'], (analysis_parameters['date']-datetime.timedelta(days=1)).strftime('%Y%m%d')),my_site['Timezone'])
        dns1, sky_temp1, amb_temp1 = BoltwoodSensor.ReadTempLog('/rdata/airglow/templogs/cloudsensor/%s/Cloud_%s_%s.txt' % (analysis_parameters['site'], analysis_parameters['site'], analysis_parameters['date'].strftime('%Y%m%d')),my_site['Timezone'])
        dns2, sky_temp2, amb_temp2 = BoltwoodSensor.ReadTempLog('/rdata/airglow/templogs/cloudsensor/%s/Cloud_%s_%s.txt' % (analysis_parameters['site'], analysis_parameters['site'], (analysis_parameters['date']+datetime.timedelta(days=1)).strftime('%Y%m%d')),my_site['Timezone'])
        dns = np.concatenate((dns0,dns1,dns2))
        sky_temp = np.concatenate((sky_temp0, sky_temp1, sky_temp2))
        amb_temp = np.concatenate((amb_temp0, amb_temp1, amb_temp2))
    else:
        # No cloud data
        dns = times.tz_localize(pytz.utc)
        sky_temp = np.full(times.shape, np.nan)
        amb_temp = np.full(times.shape, np.nan)

    # find the moon
    s_info = asiinfo.get_site_info(analysis_parameters['site'])
    moon = ephem.Moon
    obs = ephem.Observer()
    obs.lat = str(s_info['Location'][0])
    obs.lon = str(s_info['Location'][1])

    moon_alt = []
    for t in times:
        obs.date = t
        moon = ephem.Moon(obs)
        moon_alt.append(moon.alt.real * 180./np.pi)

    # Create dataframes to merge data and interpolate easily
    image_df = pd.DataFrame(times, columns=['UT'])
    image_df.index = pd.to_datetime(times)
    image_df.index = image_df.index.tz_localize(pytz.utc)

    # Create the timzone aware cloud data frame
    cloud_df = pd.DataFrame(np.array([dns,sky_temp]).T, columns=['LT','Cloud Temperature'])
    cloud_df.index = pd.to_datetime(dns)
    cloud_df.index = cloud_df.index.tz_convert(pytz.utc)

    # Merge the cloud data onto the image data. This is done by finding the previous measurement
    # in cloud_df. May want to look into some sort of interpolation scheme?
    results_df = pd.merge_asof(image_df, cloud_df, left_index=True, right_index=True,
                                    tolerance=pd.Timedelta("60s"),)
    results_df = results_df.drop(['LT'],axis=1)
    results_df['Moon Altitude'] = moon_alt

    return results_df

def run_gabor(subset, real_lam, theta, n_stds, Bf, Btheta, NUM_PROCESSORS = 32):
# Run the Gabor filter on a dataset
# INPUTS:
#   subset - xarray data structure containing information on the dataset to be processed with at least the following format:
#       Coordinates:
#           time (datetime64) - time stamps of each image in the data set
#           east (float) - east km distance from the center of the image. This is a uniform step.
#           north (float) - north km distance from the center of the image. This is a uniform step.
#       Data Variable:
#           FilteredImageData (float) - intensity values in the image (dims: time, north, east)
#                                       the image data has had a 7x7 median filter applied to remove stars and has been temporally filtered 
#   real_lam - array of wavelengths to run the Gabor filtering through (meters)
#   theta - array of orientations to run the Gabor filtering through (radians)
#   n_stds - width of the Gabor filter kernel (see Ch 3.2 of Grawe thesis)
#   Bf - frequency bandwidth for the Gabor filter kernel (see Ch 3.2 of Grawe thesis)
#   Btheta - orientation bandwidth for the Gabor filter kernel (see Ch 3.2 of Grawe thesis)
#   NUM_PROCESSORS - (optional) number of cores to use in the processing
# OUTPUTS:
#   energies_data - size (time, len(real_lam), len(theta)) array containing the energy in a specified Gabor filter kernel at each time step
#   LAM - size (len(real_lam), len(theta)) array containing the meshgrid of wavelengths
#   THETA - size (len(real_lam), len(theta)) array containing the meshgrid of orientations

    # Set up the Gabor filter
    lam        = real_lam/subset['dr'].values # This dr must be what gets the spatial information!
    u          = 1./lam
    g          = gabor.Create_FilterBank_Gabor(n_stds, u, theta, Bf, Btheta)

    gabor.g_kernels = g

    # Need to run filter on data with the background removed
    background_subtract = subset.FilteredImageData - subset.FilteredImageData.mean(dim=['east','north'])
    gabor.data_frames = background_subtract.transpose('time','north','east')

    print("largest kernel size: " + str(g[-1][0].real.shape[1]) + " x " + str(g[-1][0].real.shape[0]) + " pixels")

    # Setup the grid of wavelengths and orientations
    THETA, LAM = np.meshgrid(np.degrees(theta), real_lam/1000)

    energies_data = []

    # Setup multiprocessing and run the analysis
    pool = mp.Pool(processes = NUM_PROCESSORS)
    energies_data = pool.map(gabor.Energies_Direct_PLL, np.arange(0,len(background_subtract.time)))

    return energies_data, LAM, THETA

def run_rcp(subset, mask):
# Run the rolling cross periodogram process (see Ch 3.3 of Grawe thesis)
# INPUTS:
#   subset - xarray data structure containing information on the dataset to be processed with at least the following format:
#       Coordinates:
#           time (datetime64) - time stamps of each image in the data set
#           east (float) - east km distance from the center of the image. This is a uniform step.
#           north (float) - north km distance from the center of the image. This is a uniform step.
#       Data Variable:
#           GaborFilteredImageData (float) - intensity values in the image (dims: time, north, east)
#                                       which has been run through specific Gabor filters to isolate a wave of interest
#           dr (float) - resolution in meters of the projected image data
#   mask - array (of length corresponding to time in subset) containing True/False for whether a specific image should be analyzed
# OUTPUTS:
#   periods_output - array (of length corresponding to time in subset) of estimated periods (minutes)
#   wavelengths_output - array (of length corresponding to time in subset) of estimated wavelengths (km)
#   orientations_output - array (of length corresponding to time in subset) of estimate orientations (deg)
#   wavespeeds_output - array (of length corresponding to time in subset) of estimate wavespeeds (m/s)

    dx = subset['dr'].values
    dy = subset['dr'].values


    periods_output = []
    wavelengths_output = []
    orientations_output = []
    wavespeeds_output = []

    # First entry is nan
    periods_output.append(np.nan)
    wavelengths_output.append(np.nan)
    orientations_output.append(np.nan)
    wavespeeds_output.append(np.nan)

    # Now run the cross periodogram on consecutive images
    for i in np.arange(0,len(subset['time'])-1):
        if mask[i]:
            # The two frames we are using
            frame1 = subset['GaborFilteredImageData'].isel(time=i)
            frame2 = subset['GaborFilteredImageData'].isel(time=i+1)
        
            # Time step
            dt    = pd.to_timedelta((subset['time'].isel(time=i+1) - subset['time'].isel(time=i)).values)
        
            # Windowing
            sx, sy = frame1.shape[1], frame1.shape[0]
            windowx, windowy = np.hamming(sx), np.hamming(sy) ## hamming window
            win2D  = np.outer(windowy, windowx) # this may be backwards
            
            frame1 *= win2D ## smooth with hamming window
            frame2 *= win2D ## size = 599*602
        
            P2 = fft.fftshift(fft.fft2(frame2, (1024, 1024))) ## Y2[kx,ky]
            P1 = fft.fftshift(fft.fft2(frame1, (1024, 1024))) ## Y1[kx,ky]
            nx, ny = P1.shape[1], P1.shape[0] ## n-point DFT, nx=ny=1024
            ns = float(nx*ny)
            
            C  = P1*P2.conj()/ns ## Cross-spectural density SY1Y2[k] (3.20), size = 1024*1024
            
            kx = np.linspace(-np.pi/(dx), np.pi/(dx), nx) ## normalize in angular space
            ky = np.linspace(-np.pi/(dy), np.pi/(dy), ny)
            KX, KY = np.meshgrid(kx*1000, ky*1000)
            absC = np.abs(C)
        
            # kx and ky index of the peak in the cross-periodogram
            kym, kxm = peak_local_max(absC[512:,:], min_distance = 1, num_peaks = 1)[0] # was [512:,:512]
        
            # kx and ky values of the peak in the cross-periodogram
            kym=kym+512
            kxmax, kymax = KX[kym, kxm], KY[kym, kxm]
        
            # Get the mirror image location
            kym2 = ny - kym
            kxm2 = nx - kxm
            kxmax2, kymax2 = KX[kym2-1, kxm2-1], KY[kym2-1, kxm2-1]
        
            phase = np.angle(C)[kym, kxm]
            k    = np.sqrt(kxmax**2 + kymax**2)
            ori  = np.degrees(np.arctan2(kymax,kxmax))
            if phase < 0:
                phase = np.angle(C)[kym2, kxm2]
                k    = np.sqrt(kxmax2**2 + kymax2**2)
                ori  = np.degrees(np.arctan2(kymax2,kxmax2))
                # Need to use the other peak
                kxmax = kxmax2
                kymax = kymax2
        
            lamb = 2*np.pi/k*1000
        
            delr = lamb*phase/(2*np.pi) ## distance @ peak
            ws   = delr/dt.total_seconds() ## wavespeed @ peak
        
            # note: potential issues with 360 ambiguity resolution
            # and the use of pi phase unwrapping
        
            T = lamb/ws/60. ## period @ peak

            periods_output.append(T)
            wavelengths_output.append(lamb/1000.)
            orientations_output.append(np.degrees(np.arctan2(kymax,kxmax)))
            wavespeeds_output.append(ws)
        else:
            periods_output.append(np.nan)
            wavelengths_output.append(np.nan)
            orientations_output.append(np.nan)
            wavespeeds_output.append(np.nan)

    return periods_output, wavelengths_output, orientations_output, wavespeeds_output