#!/usr/bin/env python3
# MANGODisplay_OSN.py
# Script to create MANGO display movies with FPI data overplotted and store results on OSN

from datetime import datetime, timedelta
import os
import shutil
from pandas import date_range
from glob import glob
from matplotlib import ticker, gridspec, dates
from cartopy import crs, feature
import gc
import numpy as np
import prepare_agimages
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import calendar
import pytz
import plotmangodasi as pmd
import fpiinfo
from os.path import exists
import shutil
import matplotlib.pyplot as plt
from optparse import OptionParser
import BoltwoodSensor
import subprocess
import asiinfo
import ephem
import MANGO_L2
import prepare_agimages
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D
import logging
from pathlib import Path
from cloud_storage import CloudStorage, Configuration
import cv2


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


def read_FPI(fname, reference='zenith'):
    """Read FPI data from file."""
    # Read FPI data
    npzfile = np.load(fname, allow_pickle='False', encoding='latin1')
    FPI_Results = npzfile['FPI_Results']
    FPI_Results = FPI_Results.reshape(-1)[0]
    del npzfile.f  # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
    npzfile.close()

    # Process the data for each direction
    FPI_ut = {}
    FPI_wind = {}
    FPI_error = {}
    FPI_cloud = {}
    FPI_wq = {}

    FPI_ut['North'], FPI_wind['North'], FPI_error['North'], FPI_cloud['North'], FPI_wq['North'] = Process_FPI(FPI_Results, 'North', reference=reference)
    FPI_ut['East'], FPI_wind['East'], FPI_error['East'], FPI_cloud['East'], FPI_wq['East'] = Process_FPI(FPI_Results, 'East', reference=reference)
    FPI_ut['South'], FPI_wind['South'], FPI_error['South'], FPI_cloud['South'], FPI_wq['South'] = Process_FPI(FPI_Results, 'South', reference=reference)
    FPI_ut['West'], FPI_wind['West'], FPI_error['West'], FPI_cloud['West'], FPI_wq['West'] = Process_FPI(FPI_Results, 'West', reference=reference)
    FPI_ut['Zenith'], FPI_wind['Zenith'], FPI_error['Zenith'], FPI_cloud['Zenith'], FPI_wq['Zenith'] = Process_FPI(FPI_Results, 'Zenith', reference=reference)

    return FPI_ut, FPI_wind, FPI_error, FPI_cloud, FPI_wq


def Process_FPI(FPI_Results, desired_dir, reference='zenith'):
    """Process FPI data for a specific direction."""
    import FPI
    from scipy import interpolate
    # desired_dir is the direction to process
    
    (ref_Dop, e_ref_Dop) = FPI.DopplerReference(FPI_Results, reference=reference)
    
    # Calculate the vertical wind and interpolate it
    ind = FPI.all_indices('Zenith', FPI_Results['direction'])
    w = (FPI_Results['LOSwind'][ind] - ref_Dop[ind])  # LOS is away from instrument
    sigma_w = FPI_Results['sigma_LOSwind'][ind]
    dt = []
    for x in FPI_Results['sky_times'][ind]:
        diff = (x - FPI_Results['sky_times'][0])
        dt.append(diff.seconds + diff.days * 86400.)
    dt = np.array(dt)

    # Remove outliers
    ind = abs(w) < 200.

    if sum(ind) <= 1:
        # No good data, just use all ind
        ind = np.array([True for i in range(len(w))])  # There has to be a clearer way to do this...

    if len(ind) == 0:
        raise Exception('No Zenith look directions')

    # Interpolate
    w2 = interpolate.interp1d(dt[ind], w[ind], bounds_error=False, fill_value=0.0)
    sigma_w2 = interpolate.interp1d(dt[ind], sigma_w[ind], bounds_error=False, fill_value=0.0)
    dt = []

    for x in FPI_Results['sky_times']:
        diff = (x - FPI_Results['sky_times'][0])
        dt.append(diff.seconds + diff.days * 86400.)
    w = w2(dt)
    sigma_w = sigma_w2(dt)
    
    ind = FPI.all_indices(desired_dir, FPI_Results['direction'])

    if desired_dir == 'Zenith':
        Doppler_Wind = (FPI_Results['LOSwind'][ind] - ref_Dop[ind])
        Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2)
    else:
        Doppler_Wind = (FPI_Results['LOSwind'][ind] - ref_Dop[ind] - w[ind] * np.cos(FPI_Results['ze'][ind] * np.pi / 180.)) / np.sin(FPI_Results['ze'][ind] * np.pi / 180.)
        Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2 + sigma_w[ind]**2)
    if desired_dir == 'South' or desired_dir == 'West':
        Doppler_Wind = -Doppler_Wind
        
    return FPI_Results['sky_times'][ind], Doppler_Wind, Doppler_Error, FPI_Results['Clouds']['mean'][ind], FPI_Results['wind_quality_flag'][ind]


def toTimestamp(d):
    """Convert datetime to timestamp."""
    return calendar.timegm(d.timetuple())


# Define parameters based on sky_line_tag
def get_parameters(sky_line_tag, analysis_parameters):
    """Get plotting parameters based on sky_line_tag."""
    params = {}
    if sky_line_tag == 'XG':
        params['figsize'] = (8, 4)
        params['width_ratios'] = [2., 2]
        params['extent'] = [-122, -95, 25, 55]
        params['ymin'] = -200
        params['ymax'] = 200
        params['title'] = 'Green Line Winds'
        params['scale'] = 25.
        params['xoff'], params['yoff'] = 0.06, 0.05
        params['cloc'] = [0.7, -0.0, 0.4, 0.02]
        if (np.isnan(analysis_parameters['Thi'])):
            params['cmap'] = 'viridis'
        else:
            params['cmap'] = 'Greens_r'
        
    elif sky_line_tag == 'XR':
        params['figsize'] = (9, 4)
        params['width_ratios'] = [2.2, 2]
        params['extent'] = [-125, -75, 25, 55]
        params['ymin'] = -300
        params['ymax'] = 300
        params['title'] = 'Red Line Winds'
        params['scale'] = 100.
        params['xoff'], params['yoff'] = 0.06, 0.05
        params['cloc'] = [0.55, -0.02, 0.5, 0.02]
        if (np.isnan(analysis_parameters['Thi'])):
            params['cmap'] = 'inferno'
        else:
            params['cmap'] = 'Reds'
    return params


def download_fpi_data(cloud_storage, config, site, date, sky_line_tag, fpi_dir):
    """Download FPI data from OSN."""
    logger = logging.getLogger(__name__)
    created_files = []
    
    # Get the instrument name for this site
    instr_name = fpiinfo.get_instr_at(site, date)[0]
    
    # Create the filename and file path
    fpi_datestr = date.strftime('%Y%m%d')
    fpi_filename = f"{instr_name}_{site}_{fpi_datestr}_{sky_line_tag.lower()}.npz"
    
    # Handle the case for older files that might not have the tag
    cloud_key = f"{config.aws_results_prefix}{fpi_filename}"
    local_fpi_path = os.path.join(fpi_dir, fpi_filename)
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(local_fpi_path), exist_ok=True)
    
    # Download the file
    if cloud_storage.download_file(cloud_key, local_fpi_path):
        logger.info(f"Downloaded FPI file: {fpi_filename}")
        created_files.append(local_fpi_path)
    else:
        # Try without the tag for older files
        if sky_line_tag == 'XR':
            fpi_filename = f"{instr_name}_{site}_{fpi_datestr}.npz"
            cloud_key = f"{config.aws_results_prefix}{fpi_filename}"
            local_fpi_path = os.path.join(fpi_dir, fpi_filename)
            
            if cloud_storage.download_file(cloud_key, local_fpi_path):
                logger.info(f"Downloaded FPI file (without tag): {fpi_filename}")
                created_files.append(local_fpi_path)
            else:
                logger.warning(f"Failed to download FPI file for {site} on {fpi_datestr}")
        else:
            logger.warning(f"Failed to download FPI file for {site} on {fpi_datestr}")
    
    return created_files, local_fpi_path


def delete_files_and_empty_dirs(file_list):
    """Delete all files in the provided list and remove empty directories."""
    logger = logging.getLogger(__name__)
    
    # First delete all files
    for file_path in file_list:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                logger.info(f"Successfully deleted: {file_path}")
            else:
                logger.warning(f"File not found: {file_path}")
        except Exception as e:
            logger.error(f"Error deleting {file_path}: {e}")
    
    # Track directories to check for deletion
    dirs_to_check = set(os.path.dirname(file_path) for file_path in file_list)
    
    # Recursively delete empty directories
    for directory in sorted(dirs_to_check, key=len, reverse=True):
        try:
            if os.path.exists(directory) and not os.listdir(directory):
                os.rmdir(directory)
                logger.info(f"Deleted empty directory: {directory}")
                
                # Check parent directory
                parent_dir = os.path.dirname(directory)
                if parent_dir and parent_dir != directory:
                    dirs_to_check.add(parent_dir)
        except Exception as e:
            logger.error(f"Error removing directory {directory}: {e}")

def create_video_with_opencv(png_pattern, mp4file, framerate=15):
    # Get all matching PNG files and sort them
    png_files = sorted(glob(png_pattern))
    
    if not png_files:
        print(f"No PNG files found matching pattern: {png_pattern}")
        return False
    
    # Read the first image to get dimensions
    img = cv2.imread(png_files[0])
    if img is None:
        print(f"Failed to read image: {png_files[0]}")
        return False
        
    height, width, layers = img.shape
    
    # Create video writer
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # MPEG-4 codec
    video = cv2.VideoWriter(mp4file, fourcc, framerate, (width, height))
    
    # Add each image to the video
    for i, png_file in enumerate(png_files):
        img = cv2.imread(png_file)
        if img is not None:
            video.write(img)
        else:
            print(f"Warning: Failed to read image {i+1}/{len(png_files)}: {png_file}")
    
    # Release the video writer
    video.release()
    print(f"Video saved to {mp4file}")
    return True

# Example usage to replace your original code
def generate_video(png_pattern, mp4file):
    logger.info(f"Generating video from pattern: {png_pattern}")
    result = create_video_with_opencv(png_pattern, mp4file, framerate=15)
    if result:
        logger.info(f"Successfully generated video: {mp4file}")
    else:
        logger.error(f"Failed to generate video: {mp4file}")

def MakeSummaryMovies(system_parameters, analysis_parameters, cloud_storage, config, delete_working_files=True):
    """
    Runs the full process of creating a fusion movie with MANGO imaging data and
    FPI data overplotted. Results in a mp4 file. All data downloaded and created
    are removed upon completion and the movie is uploaded to OSN.
    """
    logger = logging.getLogger(__name__)
    
    if analysis_parameters['emission'] == 'green':
        sky_line_tag = 'XG'
    elif analysis_parameters['emission'] == 'red':
        sky_line_tag = 'XR'
    else:
        sky_line_tag = 'X'

    downloaded_files = []
    created_files = []

    # Download data
    if analysis_parameters['download_data']:
        for site in analysis_parameters['sites_asi']:
            try:
                if sky_line_tag == 'XG':
                    cmd = '/usr/bin/wget -r -nH --cut-dirs=8 --no-parent -P %s/%s/%s https://data.mangonetwork.org/data/transport/mango/archive/%s/greenline/level1/%04d/%s/' % \
                    (system_parameters['ASI_directory'], site, analysis_parameters['date'].year, site, analysis_parameters['date'].year, analysis_parameters['date'].strftime('%j'))
                elif sky_line_tag == 'XR':
                    cmd = '/usr/bin/wget -r -nH --cut-dirs=8 --no-parent -P %s/%s/%s https://data.mangonetwork.org/data/transport/mango/archive/%s/redline/level1/%04d/%s/' % \
                    (system_parameters['ASI_directory'], site, analysis_parameters['date'].year, site, analysis_parameters['date'].year, analysis_parameters['date'].strftime('%j'))
                logger.info(f"Downloading ASI data for {site}: {cmd}")
                MANGO_L2.runcmd(cmd)
            except Exception as e:
                logger.error(f"Failure to download {site} on {analysis_parameters['date']}: {str(e)}")

            # This is where the data should have been downloaded to
            target_directory = '%s/%s/%s/%s' % (system_parameters['ASI_directory'], site, analysis_parameters['date'].year, analysis_parameters['date'].strftime('%j'))
            # Add downloaded files to the list
            for root, dirs, files in os.walk(target_directory):
                for file in files:
                    downloaded_files.append(os.path.join(root, file))

    # Download FPI data
    fpi_dir = os.path.join(system_parameters['FPI_directory'])
    os.makedirs(fpi_dir, exist_ok=True)
    
    # Create a day earlier date to match FPI file naming convention
    fpi_dt = analysis_parameters['date'] + timedelta(days=-1)
    
    # Download FPI data for each site
    for fpi_site in analysis_parameters['sites_fpi']:
        files, _ = download_fpi_data(cloud_storage, config, fpi_site, fpi_dt, sky_line_tag.lower(), fpi_dir)
        downloaded_files.extend(files)

    # Convert to xarray
    ds = {}

    # Load each site
    for site in analysis_parameters['sites_asi']:
        files = glob('%s/%s/%s/*%s*.hdf5' % (system_parameters['ASI_directory'], site, analysis_parameters['date'].strftime('%Y/%j'), analysis_parameters['emission']))
        
        if len(files) > 0:
            file_path = files[0]

            # Load the data
            ds[site] = MANGO_L2.load_hdf5_to_xarray(file_path)

    # Process all images
    for site in ds.keys():
        logger.info(f"Processing images for site: {site}")
        # Apply a temporal filter to the data on this night
        IM3D = ds[site].ImageData.values
        
        # Swap dimensions in image data to get it in the right format for the processing
        IM3Dt = np.transpose(IM3D, (1, 2, 0))
        times = pd.to_datetime(ds[site].time.values)
        
        # Create the filter and implement it
        if (np.isnan(analysis_parameters['Tlo'])) and (np.isnan(analysis_parameters['Thi'])):
            # No filtering
            IM3Dfilt = IM3Dt
        else:
            b = prepare_agimages.initialize_airglow_filter(analysis_parameters['ntaps'], analysis_parameters['Tlo'], analysis_parameters['Thi'], times)
            IM3Dfilt = prepare_agimages.filter_airglow(IM3Dt, b, analysis_parameters['ntaps'])
        
        # Add this to the dataset, rearranging the dimensions back to (time, north, east)
        filt_array = xr.DataArray(np.transpose(IM3Dfilt, (2, 0, 1)), dims=('time', 'north', 'east'))
        ds[site]['FilteredImageData'] = filt_array

        # Get cloud and moon information
        analysis_parameters['site'] = site
        results_df = MANGO_L2.load_cloud_and_moon(analysis_parameters, times)

        # Create a new variable (e.g., temperature) with the same length as times
        cloud_temperature = xr.DataArray(results_df['Cloud Temperature'], dims=('time'), coords={'time': times})
        moon_angle = xr.DataArray(results_df['Moon Altitude'], dims=('time'), coords={'time': times})
        
        # Add the temperature variable to the dataset
        ds[site]['cloud_temperature'] = cloud_temperature
        ds[site]['moon_angle'] = moon_angle
        
        # You can add attributes or metadata to the temperature variable if needed
        ds[site]['cloud_temperature'].attrs['units'] = 'Celsius'
        ds[site]['moon_angle'].attrs['units'] = 'Degrees'

    # Process FPI Data
    FPI_ut = {}
    FPI_wind = {}
    FPI_error = {}
    FPI_cloud = {}
    FPI_wq = {}
    FPI_tt = {}
    FPI_walpha = {}

    for fpi in analysis_parameters['sites_fpi']:
        fname = os.path.join(fpi_dir, f"{fpiinfo.get_instr_at(fpi, fpi_dt)[0]}_{fpi}_{fpi_dt.strftime('%Y%m%d')}_{sky_line_tag.lower()}.npz")
        if (not exists(fname)) and (sky_line_tag == 'XR'):
            # Some of the older npz files for redline don't have a tag
            fname = os.path.join(fpi_dir, f"{fpiinfo.get_instr_at(fpi, fpi_dt)[0]}_{fpi}_{fpi_dt.strftime('%Y%m%d')}.npz")

        if exists(fname):
            FPI_ut[fpi], FPI_wind[fpi], FPI_error[fpi], FPI_cloud[fpi], FPI_wq[fpi] = read_FPI(fname)
            FPI_tt[fpi] = {}
            FPI_walpha[fpi] = {}
            for d in FPI_ut[fpi].keys():
                FPI_tt[fpi][d] = np.array([toTimestamp(d.astimezone(pytz.utc)) for d in FPI_ut[fpi][d]])

                FPI_walpha[fpi][d] = FPI_wq[fpi][d].copy()
                FPI_walpha[fpi][d][FPI_walpha[fpi][d] == 2] = 0.3
                FPI_walpha[fpi][d][FPI_walpha[fpi][d] == 1] = 0.7
                FPI_walpha[fpi][d][FPI_walpha[fpi][d] == 0] = 1.0

    # Find the unique times to work with
    all_times = []
    for k in ds.keys():
        if hasattr(ds[k], 'time') and ds[k].time is not None:
            all_times += list(ds[k].time.values)

    if not all_times:
        logger.error("No time data found in datasets")
        return

    all_times.sort()

    all_fpi_times = []
    for k in FPI_ut.keys():
        for d in FPI_ut[k].keys():
            all_fpi_times += list(FPI_ut[k][d])

    all_fpi_times.sort()

    # Remove values from alltimes that are too close in time (e.g., if one imager is a couple of seconds out of sync)
    all_unique_times = np.sort(np.unique(all_times))

    unique_times = [all_unique_times[0]]
    for i in range(1, len(all_unique_times)):
        if all_unique_times[i] - all_unique_times[i-1] > np.timedelta64(10, 's'):
            unique_times.append(all_unique_times[i])

    # Load parameters
    params = get_parameters(sky_line_tag, analysis_parameters)

    my_colors = {}
    my_colors['cvo'] = '#1f77b4'
    my_colors['blo'] = '#ff7f0e'
    my_colors['low'] = '#9467bd'
    my_colors['uao'] = '#8c564b'

    # The largest gap size [in minutes] to allow in plotting
    max_gap = 120.
    max_fpi_gap = 30.

    # Ensure output directory exists
    os.makedirs(system_parameters['output_directory'], exist_ok=True)

    # THIS WOULD BE THE LOOP START
    for target_time in unique_times:
        fig = plt.figure(figsize=params['figsize'])
        spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig,
                                left=0.04, right=0.94, bottom=0.1, top=0.90,
                                wspace=0.05, hspace=0.25,
                                width_ratios=params['width_ratios'])
        legend_elements = []

        axes00 = fig.add_subplot(spec[:, 0], projection=crs.Orthographic(np.nanmean(params['extent'][:2]), np.nanmean(params['extent'][2:])))
        axes00.add_feature(feature.COASTLINE)
        axes00.add_feature(feature.STATES, alpha=0.2, zorder=1)

        for k in ds.keys():
            if hasattr(ds[k], 'time') and ds[k].time is not None:
                nearest_time = ds[k].sel(time=target_time, method='nearest')
                time_difference = abs(pd.to_timedelta(nearest_time.time.values - target_time).total_seconds())

                if time_difference < 10:
                    data = ds[k].sel(time=target_time, method='nearest')

                    latitude = data.Latitude
                    longitude = data.Longitude
                    mask = data.Elevation > analysis_parameters['el_cutoff']

                    pc = axes00.pcolormesh(longitude, latitude, data.FilteredImageData.where(mask), transform=crs.PlateCarree(), cmap=params['cmap'])

        # Set the title and limits of the map
        axes00.set_title('%s UT' % pd.to_datetime(str(target_time)))
        axes00.set_extent(params['extent'], crs=crs.PlateCarree())
        axes00.set_facecolor('lightgray')

        pos = axes00.get_position()
        cax = fig.add_axes([pos.x0 + params['cloc'][0] * pos.width, pos.y0 + params['cloc'][1], params['cloc'][2] * pos.width - 0.05, + params['cloc'][3]])
        cbar = fig.colorbar(pc, cax=cax, orientation='horizontal')
        cbar.set_label('Intensity [units]', fontsize=8)
        cbar.set_ticks([])
        cbar.set_ticklabels([])
        cbar.ax.tick_params(labelsize=8)

        axes01 = fig.add_subplot(spec[0, 1])
        axes02 = fig.add_subplot(spec[1, 1])

        for fpi in FPI_ut.keys():
            for i in range(len(FPI_ut[fpi]['North']) - 1):
                if (FPI_ut[fpi]['North'][i + 1] - FPI_ut[fpi]['North'][i]).total_seconds() < (max_gap * 60.):
                    axes01.plot(FPI_ut[fpi]['North'][i:i + 2], FPI_wind[fpi]['North'][i:i + 2], color=my_colors[fpi], label='_nolegend_', alpha=FPI_walpha[fpi]['North'][i + 1])
            axes01.scatter(FPI_ut[fpi]['North'], FPI_wind[fpi]['North'], ec=None, color=my_colors[fpi], label='_nolegend_', alpha=FPI_walpha[fpi]['North'], s=10)

            for i in range(len(FPI_ut[fpi]['West']) - 1):
                if (FPI_ut[fpi]['West'][i + 1] - FPI_ut[fpi]['West'][i]).total_seconds() < (max_gap * 60.):
                    axes02.plot(FPI_ut[fpi]['West'][i:i + 2], FPI_wind[fpi]['West'][i:i + 2], color=my_colors[fpi], label='_nolegend_', alpha=FPI_walpha[fpi]['West'][i + 1])
            axes02.scatter(FPI_ut[fpi]['West'], FPI_wind[fpi]['West'], ec=None, color=my_colors[fpi], label='_nolegend_', alpha=FPI_walpha[fpi]['West'], s=10)

            legend_elements.append(Line2D([0], [0], marker='o', linestyle='None', color=my_colors[fpi], label=fpi, markerfacecolor=my_colors[fpi], markersize=np.sqrt(10)))

        axes01.set_ylim([params['ymin'], params['ymax']])
        axes01.axhline(y=0, color='k', linewidth=1)
        axes01.axvline(x=target_time, color='r', linewidth=1)
        axes01.xaxis.set_major_locator(dates.HourLocator(interval=2))
        axes01.xaxis.set_minor_locator(dates.HourLocator())
        axes01.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        axes01.set_title(params['title'])
        axes01.tick_params(axis='both', which='major', size=6, direction='in', right=True, top=True, labelright=True, labelleft=False)
        axes01.tick_params(axis='both', which='minor', size=3, direction='in', right=True, top=True, labelright=True, labelleft=False)
        axes01.set_ylabel('Northward Winds [m/s]')
        axes01.legend(handles=legend_elements, loc='upper right', ncol=2, fontsize=8)
        if len(all_fpi_times) > 0:
            axes01.set_xlim([min(all_fpi_times), max(all_fpi_times)])

        axes02.set_ylim([params['ymin'], params['ymax']])
        axes02.axhline(y=0, color='k', linewidth=1)
        axes02.axvline(x=target_time, color='r', linewidth=1)
        axes02.xaxis.set_major_locator(dates.HourLocator(interval=2))
        axes02.xaxis.set_minor_locator(dates.HourLocator())
        axes02.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        axes02.tick_params(axis='both', which='major', size=6, direction='in', right=True, top=True, labelright=True, labelleft=False)
        axes02.tick_params(axis='both', which='minor', size=3, direction='in', right=True, top=True, labelright=True, labelleft=False)
        axes02.yaxis.tick_right()
        axes02.yaxis.set_label_position("left")
        axes02.set_ylabel('Eastward Wind [m/s]')
        axes02.set_xlabel('UT [hrs]')
        if len(all_fpi_times) > 0:
            axes02.set_xlim([min(all_fpi_times), max(all_fpi_times)])

        low_quality = Line2D([], [], color=my_colors['cvo'], alpha=0.3, label='Low')
        med_quality = Line2D([], [], color=my_colors['cvo'], alpha=0.7, label='Medium')
        high_quality = Line2D([], [], color=my_colors['cvo'], alpha=1.0, label='High Quality')
        axes02.legend(handles=[low_quality, med_quality, high_quality], loc='upper right', ncols=3, fontsize=6)

        width = 0.015
        headwidth = 4
        headlength = 4
        headaxislength = headlength - 1
        minshaft = 2
        sc = axes00.bbox.width / fig.dpi / (10. * params['scale'])

        obj = None
        for fpi in FPI_ut.keys():
            valid_t = False

            if (np.min(abs(toTimestamp(pd.to_datetime(str(target_time))) - FPI_tt[fpi]['North'])) < max_fpi_gap * 60) & (
                    np.min(abs(toTimestamp(pd.to_datetime(str(target_time))) - FPI_tt[fpi]['West'])) < max_fpi_gap * 60):
                valid_t = True

            if valid_t:
                u = np.interp(toTimestamp(pd.to_datetime(str(target_time))), FPI_tt[fpi]['West'], FPI_wind[fpi]['West']) * sc
                v = np.interp(toTimestamp(pd.to_datetime(str(target_time))), FPI_tt[fpi]['North'], FPI_wind[fpi]['North']) * sc
                a = 1
            else:
                u = np.nan
                v = np.nan
                a = 0

            if a < 0.5:
                a = 0

            glat, glon, x = fpiinfo.get_site_info(fpi)['Location']

            obj = axes00.quiver(np.array([glon]), np.array([glat]), np.array([u]), np.array([v]),
                                angles='uv', scale_units='inches', scale=1, width=width,
                                pivot='tail', headwidth=headwidth, headlength=headlength, alpha=a,
                                minshaft=minshaft, headaxislength=headaxislength, color=my_colors[fpi], transform=crs.PlateCarree(), zorder=1000)
            
        x0, x1, y0, y1 = axes00.get_extent()
        if obj is not None:
            quiverobj1 = axes00.quiverkey(obj, params['xoff'], params['yoff'], sc * params['scale'], r'$%i\,\frac{m}{s}$' % params['scale'],
                                        labelsep=0.05, color="black", alpha=1, coordinates='axes', labelpos='N',
                                        transform=axes00.transAxes, fontproperties={'size': 8})
        
        today_date = pd.to_datetime('today').strftime('%m-%d-%Y')
        if np.isnan(analysis_parameters['Tlo']) and np.isnan(analysis_parameters['Thi']):
            text1 = 'No filtering'
        elif np.isnan(analysis_parameters['Tlo']):
            text1 = 'Highpass filter cutoff at %d min' % analysis_parameters['Thi']
        elif np.isnan(analysis_parameters['Thi']):
            text1 = 'Lowpass filter cutoff at %d min' % analysis_parameters['Tlo']
        else:
            text1 = 'Temporal filter: %d-%d min' % (analysis_parameters['Tlo'], analysis_parameters['Thi'])
        text2 = 'Plotted on ' + today_date
        text_x = 0
        axes00.annotate(text1, xy=(text_x, -0.05), xycoords='axes fraction', ha='left', va='bottom', fontsize=8)
        axes00.annotate(text2, xy=(text_x, -0.05 - 0.05), xycoords='axes fraction', ha='left', va='bottom', fontsize=8)

        spec.tight_layout(fig)

        # Create output directory if it doesn't exist
        os.makedirs(system_parameters['output_directory'], exist_ok=True)
        
        # Save the PNG
        savename = os.path.join(system_parameters['output_directory'], 
                               f"MANGO_{pd.to_datetime(str(target_time)).strftime('%Y%m%d_%H%M%S')}_{sky_line_tag}.png")
        plt.savefig(savename, dpi=300)
        created_files.append(savename)

        plt.close(fig)

    # Generating video from PNG frames
    linetag = "red" if sky_line_tag == 'XR' else "green"
    mp4file = os.path.join(system_parameters['output_directory'], 
                           f"MANGO_{analysis_parameters['date'].strftime('%Y%m%d')}_{linetag}.mp4")
    
    # Create ffmpeg command to create MP4 from PNGs
    png_pattern = os.path.join(system_parameters['output_directory'], f"*MANGO_{pd.to_datetime(str(target_time)).strftime('%Y%m%d')}*_{sky_line_tag}.png")
 #   cmd = f'/usr/bin/ffmpeg -framerate 15 -pattern_type glob -i "{png_pattern}" -c:v mpeg4 -q:v 1 -y {mp4file}'
 #   
 #   logger.info(f"Generating video with command: {cmd}")
 #   os.system(cmd)
    generate_video(png_pattern, mp4file)
    
    # Upload the movie to OSN
    if os.path.exists(mp4file):
        logger.info(f"Uploading movie {mp4file} to OSN")
        cloud_key = f"{config.aws_mango_movies_prefix}MANGO_{analysis_parameters['date'].strftime('%Y%m%d')}_{linetag}.mp4"
        cloud_storage.upload_file(mp4file, cloud_key)
        logger.info(f"Movie uploaded to {cloud_key}")
    else:
        logger.error(f"Failed to create movie {mp4file}")

    # Delete all tracked files
    if delete_working_files:
        logger.info("Cleaning up downloaded and created files")
        delete_files_and_empty_dirs(downloaded_files + created_files)

if __name__=="__main__":
    # Main module allows this to be run via command line
    logger = setup_logging()

    # Parse the command line
    usage = "usage: MANGODisplay_OSN.py -y YEAR -d DOY -t EMISSION_TAG [-e ENV_FILE]"
    parser = OptionParser(usage=usage)
    parser.add_option("-y", "--year", dest="year", help="Year to be run", metavar="YEAR", type="int", default=0)
    parser.add_option("-d", "--doy", dest="doy", help="Day of year to be run", metavar="DOY", type="int", default=0)
    parser.add_option("-t", "--tag", dest="sky_line_tag", help="XG for greenline, X or XR for redline", metavar="TAG", type="str", default="XG")
    parser.add_option("-l", "--download", dest="download", help="Download data, False to not download", metavar="DOWNLOAD", type="str", default="True")
    parser.add_option("-e", "--env", dest="env_file", help="Path to .env file", metavar="ENV_FILE", type="str", default=".env")

    (options, args) = parser.parse_args()
    year = options.year
    doy = options.doy
    sky_line_tag = options.sky_line_tag
    download = options.download
    env_file = options.env_file

    # Defaults if no date is given
    if (doy == 0) or (year == 0):
        dt = datetime.now()
        # DOY offset 2 days in the past to allow for data transfer latency
        doy = int(dt.strftime("%j")) - 2
        year = int(dt.strftime("%Y"))
        download_data = True

    if download == "False":
        download_data = False
    else:
        download_data = True

    # Load configuration from .env file
    config = Configuration(env_file)
    cloud_storage = CloudStorage(config)

    # System parameters that need to be set
    system_parameters = {
        'ASI_directory': '/home/jmakela/tmp/mango/asi/',
        'FPI_directory': '/home/jmakela/tmp/mango/fpi/',
        'output_directory': '/home/jmakela/tmp/mango/movies/'
    }

    # Create directories if they don't exist
    for directory in system_parameters.values():
        os.makedirs(directory, exist_ok=True)
    
    analysis_parameters = {
        'date': datetime(year, 1, 1) + timedelta(days=doy-1),
        'ntaps': 13,
        'download_data': download_data,
        'el_cutoff': 20.
    }

    # Run the code
    try:
        if sky_line_tag == 'XG':
            analysis_parameters['sites_asi'] = ['cvo', 'low', 'blo', 'cfs', 'mro', 'bdr', 'new']
            analysis_parameters['sites_fpi'] = ['cvo', 'low', 'blo']
            analysis_parameters['emission'] = 'green'
            analysis_parameters['Tlo'] = 2
            analysis_parameters['Thi'] = 20
        elif sky_line_tag == 'XR':
            analysis_parameters['sites_asi'] = ['cfs', 'cvo', 'eio', 'mdk', 'mto', 'par']
            analysis_parameters['sites_fpi'] = ['cvo', 'low', 'blo', 'uao']
            analysis_parameters['emission'] = 'red'
            analysis_parameters['Tlo'] = 8
            analysis_parameters['Thi'] = 30

        logger.info(f"Starting processing for {year} DOY {doy} with tag {sky_line_tag}")
        MakeSummaryMovies(system_parameters, analysis_parameters, cloud_storage, config)
        logger.info(f"Processing complete for {year} DOY {doy} with tag {sky_line_tag}")
    except Exception as error:
        logger.error(f"Error {doy} {year} {sky_line_tag}: {str(error)}")
        logger.error(f"Error type: {type(error).__name__}")    
