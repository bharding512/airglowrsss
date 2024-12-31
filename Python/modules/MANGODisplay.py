from datetime import datetime,timedelta
import os,shutil
from pandas import date_range
from glob import glob
from matplotlib import ticker,gridspec,dates
import cv2
from cartopy import crs,feature
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

outfolder = "/home/airglow/scratch_data/DASI_Data/"
fpi_repo = "/rdata/airglow/fpi/results/"
repo_ASI = "/home/airglow/scratch_data/MANGO_Data"

def read_FPI(fname, reference='zenith'):

    # Read FPI data
    npzfile = np.load(fname,allow_pickle='False', encoding='latin1')
    FPI_Results = npzfile['FPI_Results']
    FPI_Results = FPI_Results.reshape(-1)[0]
    del npzfile.f # http://stackoverflow.com/questions/9244397/memory-overflow-when-using-numpy-load-in-a-loop
    npzfile.close()

    # Process the data for each direction
    FPI_ut = {}
    FPI_wind = {}
    FPI_error = {}
    FPI_cloud = {}
    FPI_wq = {}

    FPI_ut['North'],FPI_wind['North'],FPI_error['North'], FPI_cloud['North'], FPI_wq['North'] = Process_FPI(FPI_Results,'North',reference=reference)
    FPI_ut['East'],FPI_wind['East'],FPI_error['East'], FPI_cloud['East'], FPI_wq['East'] = Process_FPI(FPI_Results,'East',reference=reference)
    FPI_ut['South'],FPI_wind['South'],FPI_error['South'], FPI_cloud['South'], FPI_wq['South'] = Process_FPI(FPI_Results,'South',reference=reference)
    FPI_ut['West'],FPI_wind['West'],FPI_error['West'], FPI_cloud['West'], FPI_wq['West'] = Process_FPI(FPI_Results,'West',reference=reference)
    FPI_ut['Zenith'],FPI_wind['Zenith'],FPI_error['Zenith'], FPI_cloud['Zenith'], FPI_wq['Zenith'] = Process_FPI(FPI_Results,'Zenith',reference=reference)

    return FPI_ut, FPI_wind, FPI_error, FPI_cloud, FPI_wq

def Process_FPI(FPI_Results, desired_dir, reference='zenith'):
    import FPI
    from scipy import interpolate
    # desired_dir is the direction to process
    
    (ref_Dop, e_ref_Dop) = FPI.DopplerReference(FPI_Results,reference=reference)
    
    # Calculate the vertical wind and interpolate it
    ind = FPI.all_indices('Zenith',FPI_Results['direction'])
    w = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]) # LOS is away from instrument
    sigma_w = FPI_Results['sigma_LOSwind'][ind]
    dt = []
    for x in FPI_Results['sky_times'][ind]:
        diff = (x - FPI_Results['sky_times'][0])
        dt.append(diff.seconds+diff.days*86400.)
    dt = np.array(dt)

    # Remove outliers
    ind = abs(w) < 200.

    if sum(ind) <= 1:
    # No good data, just use all ind
        ind = np.array([True for i in range(len(w))]) # There has to be a clearer way to do this...

    if len(ind) == 0:
       raise Exception('%s: No Zenith look directions' % f)

    # Interpolate
    w2 = interpolate.interp1d(dt[ind],w[ind],bounds_error=False,fill_value=0.0)
    sigma_w2 = interpolate.interp1d(dt[ind],sigma_w[ind],bounds_error=False,fill_value=0.0)
    dt = []

    for x in FPI_Results['sky_times']:
        diff = (x - FPI_Results['sky_times'][0])
        dt.append(diff.seconds+diff.days*86400.)
    w = w2(dt)
    sigma_w = sigma_w2(dt)
    
    ind = FPI.all_indices(desired_dir,FPI_Results['direction'])

    if desired_dir == 'Zenith':
        Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind])
        Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2)
    else:
        Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
        Doppler_Error = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)
    if desired_dir == 'South' or desired_dir == 'West':
        Doppler_Wind = -Doppler_Wind
        
    return FPI_Results['sky_times'][ind], Doppler_Wind, Doppler_Error, FPI_Results['Clouds']['mean'][ind], FPI_Results['wind_quality_flag'][ind]

def toTimestamp(d):
  return calendar.timegm(d.timetuple())

#def runcmd(cmd, verbose = False, *args, **kwargs):
## From https://www.scrapingbee.com/blog/python-wget/
#    process = subprocess.Popen(
#        cmd,
#        stdout = subprocess.PIPE,
#        stderr = subprocess.PIPE,
#        text = True,
#        shell = True
#    )
#    std_out, std_err = process.communicate()
#    if verbose:
#        print(std_out.strip(), std_err)
#    pass

# Define parameters based on sky_line_tag
def get_parameters(sky_line_tag, analysis_parameters):
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

def MakeSummaryMovies(system_parameters, analysis_parameters, delete_working_files=True):
    # Runs the full process of creating a fusion movie with MANGO imaging data and
    # FPI data overplotted. Results in a mp4 file. All data downloaded and created
    # are removed upon completion
    #
    # INPUTS:
    #   system_parameters: dictionary containing the location of data
    #       'ASI_directory': location to store ASI images downloaded,
    #       'FPI_directory': location where FPI data exists,
    #       'output_directory': location to store the output,
    #   analysis_parameters: dictionary containing information about processing
    #       'date': datetime of the day to process
    #       'sites_asi': array of MANGO site codes, e.g., ['cvo','low','blo','cfs','mro','bdr'],
    #       'sites_fpi': array of FPI site codes, e.g., ['cvo','low','blo'],
    #       'emission': emission to use, either 'green' or 'red,
    #       'ntaps': number of taps in the temporal filter (e.g., 13),
    #       'Tlo': minimum period in the temporal filter in minutes (e.g., 2 for GL or 8 for RL
    #       'Thi': maximum period in the temporal filter in minutes (e.g., 20 for GL or 30 for RL
    #       'download_data': should MANGO be downloaded? True/False
    #       'el_cutoff': elevation angle to mask below in degrees (e.g., 20).

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
                print(cmd)
                MANGO_L2.runcmd(cmd)
            except:
                print('!!! Failure to download %s on %s' % (site, analysis_parameters['date']))

            # This is where the data should have been downlaeded to
            target_directory = '%s/%s/%s/%s' % (system_parameters['ASI_directory'], site, analysis_parameters['date'].year, analysis_parameters['date'].strftime('%j'))
            # Add downloaded files to the list
            for root, dirs, files in os.walk(target_directory):
                for file in files:
                    downloaded_files.append(os.path.join(root, file))

    # Convert to xarray
    ds = {}

    # Load each site
    for site in analysis_parameters['sites_asi']:
        files = glob('%s/%s/%s/*%s*.hdf5'%(system_parameters['ASI_directory'], site, analysis_parameters['date'].strftime('%Y/%j'),analysis_parameters['emission']))
        
        if len(files) > 0:
            file_path = files[0]

            # Load the data
            ds[site] = MANGO_L2.load_hdf5_to_xarray(file_path)
            
            # The data_arrays dictionary now contains xarray DataArrays for each dataset in the HDF5 file

    # Process all images
    for site in ds.keys():
        print(site)
        # Apply a temporal filter to the data on this night
        IM3D = ds[site].ImageData.values
        
        # Swap dimensions in image data to get it in the right format for the processing
        IM3Dt = np.transpose(IM3D, (1, 2, 0))
        times = pd.to_datetime(ds[site].time.values)
        
        # Create the filter and implement it
        if (np.isnan(analysis_parameters['Tlo'])) and (np.isnan(analysis_parameters['Thi'] )):
            # No filtering
            IM3Dfilt = IM3Dt
        else:
            b = prepare_agimages.initialize_airglow_filter(analysis_parameters['ntaps'],analysis_parameters['Tlo'],analysis_parameters['Thi'],times)
            IM3Dfilt = prepare_agimages.filter_airglow(IM3Dt,b,analysis_parameters['ntaps'])
        
        # Add this to the dataset, rearranging the dimensions back to (time, north, east)
        filt_array = xr.DataArray(np.transpose(IM3Dfilt, (2, 0, 1)), dims=('time','north','east'))
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

    fpi_dt = analysis_parameters['date'] + timedelta(days=-1)

    for fpi in analysis_parameters['sites_fpi']:
        fname = system_parameters['FPI_directory'] + fpiinfo.get_instr_at(fpi,fpi_dt)[0] + '_' + fpi + '_' + fpi_dt.strftime('%Y%m%d') + '_' + sky_line_tag.lower() + '.npz'
        if (exists(fname) == False) and (sky_line_tag == 'XR'):
            # Some of the older npz files for redline don't have a tag
            fname = system_parameters['FPI_directory'] + fpiinfo.get_instr_at(fpi,fpi_dt)[0] + '_' + fpi + '_' + fpi_dt.strftime('%Y%m%d') + '.npz'

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
        if ds[k].time is None:
            continue
        all_times += list(ds[k].time.values)

    all_times.sort()

    all_fpi_times = []
    for k in FPI_ut.keys():
        for d in FPI_ut[k].keys():
            all_fpi_times += list(FPI_ut[k][d])

    all_fpi_times.sort()

    # Remove values from alltimes that are too close in time (e.g., if one imager is a couple of seconds out of sync)
    all_unique_times = np.sort(np.unique(all_times))

    unique_times = [all_unique_times[0]]
    for i in range(1,len(all_unique_times)):
        if all_unique_times[i] - all_unique_times[i-1] > np.timedelta64(10,'s'):
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

    # THIS WOULD BE THE LOOP START
    for target_time in unique_times:
        fig = plt.figure(figsize=params['figsize'])
        spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig,
                                left=0.04, right=0.94, bottom=0.1, top=0.90,
                                wspace=0.05, hspace=0.25,
                                width_ratios=params['width_ratios'])
        #fig.clf()
        legend_elements = []

        axes00 = fig.add_subplot(spec[:, 0], projection=crs.Orthographic(np.nanmean(params['extent'][:2]), np.nanmean(params['extent'][2:])))
        axes00.add_feature(feature.COASTLINE)
        axes00.add_feature(feature.STATES, alpha=0.2, zorder=1)

        for k in ds.keys():
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
        cbar.set_label('Intensity [units]',fontsize=8)
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
            axes01.set_xlim([np.unique(all_fpi_times)[0], np.unique(all_fpi_times)[-1]])

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
            axes02.set_xlim([np.unique(all_fpi_times)[0], np.unique(all_fpi_times)[-1]])

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

        folderpngs = system_parameters['output_directory']
        savename = folderpngs + 'MANGO_%s_%s.png' % (pd.to_datetime(str(target_time)).strftime('%Y%m%d_%H%M%S'), sky_line_tag)
        plt.savefig(savename, dpi=300)
        created_files.append(savename)

        plt.close(fig)

    # Generating video from PNG frames
    linetag="red" if sky_line_tag=='XR' else "green"
    mp4file='%s/MANGO_%s_%s.mp4'%(system_parameters['output_directory'], analysis_parameters['date'].strftime('%Y%m%d'),linetag)
    cmd = '/usr/bin/ffmpeg -framerate 15 -pattern_type glob -i \"%s*MANGO_%s*_%s.png\" -c:v mpeg4 -q:v 1 -y %s'%(folderpngs, pd.to_datetime(str(target_time)).strftime('%Y%m%d'), sky_line_tag, mp4file)
    os.system(cmd)

    # Delete all tracked files
    if delete_working_files:
        for file in downloaded_files + created_files:
            try:
                os.remove(file)
            except OSError:
                pass # File does not exist

        # Remove empty directories
        for directory in set(os.path.dirname(file) for file in downloaded_files + created_files):
            try:
                os.removedirs(directory)
            except OSError:
                pass  # Directory not empty or already deleted

if __name__=="__main__":
    # Main module allows this to be run via command line

    # Parse the command line
    usage = "usage: MANGODisplay -y YEAR -d DOY -t EMISSION_TAG"
    parser = OptionParser(usage=usage)
    parser.add_option("-y", "--year", dest="year", help="Year to be run", metavar="YEAR", type="int", default=0)
    parser.add_option("-d", "--doy", dest="doy", help="Day of year to be run", metavar="DOY", type="int", default=0)
    parser.add_option("-t", "--tag", dest="sky_line_tag", help="XG for greenline, X or XR for redline", metavar="TAG", type="str", default="XG")
    parser.add_option("-l", "--download", dest="download", help="Download data, False to not download", metavar="DOWNLOAD", type="str", default="True")

    (options, args) = parser.parse_args()
    year = options.year
    doy = options.doy
    sky_line_tag = options.sky_line_tag
    download = options.download

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

    # System parameters that need to be set
    system_parameters = {
        'ASI_directory': '/home/airglow/scratch_data/MANGO_Data/',
        'FPI_directory': '/rdata/airglow/fpi/results/',
        'output_directory': '/home/airglow/scratch_data/DASI_Data/'
    }

    analysis_parameters = {
        'date': datetime(year,1,1) + timedelta(days=doy-1),
        'ntaps': 13,
        'download_data': download_data,
        'el_cutoff': 20.
    }

    # Run the code
    try:
        if sky_line_tag == 'XG':
            analysis_parameters['sites_asi'] = ['cvo','low','blo','cfs','mro','bdr','new']
            analysis_parameters['sites_fpi'] = ['cvo','low','blo']
            analysis_parameters['emission'] = 'green'
            analysis_parameters['Tlo'] = 2
            analysis_parameters['Thi'] = 20
        elif sky_line_tag == 'XR':
            analysis_parameters['sites_asi'] = ['cfs','cvo','eio','mdk','mto','par']
            analysis_parameters['sites_fpi'] = ['cvo','low','blo','uao']
            analysis_parameters['emission'] = 'red'
            analysis_parameters['Tlo'] = 8
            analysis_parameters['Thi'] = 30

        MakeSummaryMovies(system_parameters,analysis_parameters)
    except Exception as error:
        print("Error %d %d %s" % (doy, year, sky_line_tag))
        print(error)
        print(type(error).__name__)
#    MakeSummaryMovies(year, doy, sky_line_tag, download_data = False, sites_asi = ['mro'])
#    MakeSummaryMovies(year, doy, sky_line_tag, download_data = False)
