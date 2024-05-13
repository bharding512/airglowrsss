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

def runcmd(cmd, verbose = False, *args, **kwargs):
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

def MakeSummaryMovies(year, doy, sky_line_tag,sites_asi = ['cvo','low','blo','cfs','mro','bdr'], sites_fpi = ['cvo','low','blo'], ntaps = 13, Tlo   = 2, Thi   = 20, download_data = True, outfolder = "/home/airglow/scratch_data/DASI_Data/", fpi_repo = "/rdata/airglow/fpi/results/", repo_ASI = "/home/airglow/scratch_data/MANGO_Data"):

    try:
        # Date to run
        dt=datetime(year,1,1) + timedelta(days=doy-1)

        # If restarting a day, we don't need to download data and probably start at a non-zero index
        istart=0

        # Make sure the output png folder exists
        folderpngs=outfolder+"/%s/"%(dt.strftime("%Y%m%d")+"_%s"%sky_line_tag)
        if not os.path.exists(folderpngs):
            os.makedirs(folderpngs)

        # Download ASI data
        if download_data:
            for site in sites_asi:
                if sky_line_tag == 'XG':
                    cmd = '/usr/bin/wget -r -nH --cut-dirs=8 --no-parent -P %s/%s/%s https://data.mangonetwork.org/data/transport/mango/archive/%s/greenline/raw/%s/%s/' % (repo_ASI, site, dt.strftime('%Y'), site, dt.strftime('%Y'), dt.strftime('%j'))
                elif sky_line_tag == 'XR':
                    cmd = '/usr/bin/wget -r -nH --cut-dirs=8 --no-parent -P %s/%s/%s https://data.mangonetwork.org/data/transport/mango/archive/%s/redline/raw/%s/%s/' % (repo_ASI, site, dt.strftime('%Y'), site, dt.strftime('%Y'), dt.strftime('%j'))
                print(cmd)
#                os.system(cmd)
                runcmd(cmd)
    except:
        print('!!!! FAILURE DOWNLOADING DATA %s' % dt)

    # Process Imaging Data
    if sky_line_tag == 'XG':
        height = 95
        emission = '5577'
        el_cutoff = 20.
    elif sky_line_tag == 'XR':
        height = 250
        emission = '6300'
        el_cutoff = 20.

    # Define dictionaries
    IM3Dfilt = {}
    lat = {}
    lon = {}
    times = {}
    dns = {}
    amb_temp = {}
    sky_temp = {}

    for instr_name in sites_asi:
        # Code to set up and run the temporal filter
        IM3D, lat[instr_name], lon[instr_name], times[instr_name], emissionHeight, observation_point, success = prepare_agimages.prepare_airglow_dasi(instr_name, year, doy, emission, el_cutoff, coord_nan=False)
  #      print(success)
  #      print(len(times[instr_name]))
  #      print(np.size(IM3D))
        if success == False or len(times[instr_name]) == 0:
            times.pop(instr_name)
            continue
        b = prepare_agimages.initialize_airglow_filter(ntaps,Tlo,Thi,times[instr_name])
        IM3Dfilt[instr_name] = prepare_agimages.filter_airglow(IM3D, b, ntaps)

        # Load cloud sensor data, if it exists. Only works for ['cvo','low','blo']
        if instr_name in ['cvo','low','blo']:
            my_site = fpiinfo.get_site_info(instr_name,dt)

            # Because of timezones and how the cloud senor data files are broken up, 
            # read on day earlier and later just to be sure we capture the data needed
            dns0, sky_temp0, amb_temp0 = BoltwoodSensor.ReadTempLog('/rdata/airglow/templogs/cloudsensor/%s/Cloud_%s_%s.txt' % (instr_name, instr_name, (dt-timedelta(days=1)).strftime('%Y%m%d')),my_site['Timezone'])
            dns1, sky_temp1, amb_temp1 = BoltwoodSensor.ReadTempLog('/rdata/airglow/templogs/cloudsensor/%s/Cloud_%s_%s.txt' % (instr_name, instr_name, dt.strftime('%Y%m%d')),my_site['Timezone'])
            dns2, sky_temp2, amb_temp2 = BoltwoodSensor.ReadTempLog('/rdata/airglow/templogs/cloudsensor/%s/Cloud_%s_%s.txt' % (instr_name, instr_name, (dt+timedelta(days=1)).strftime('%Y%m%d')),my_site['Timezone'])
            dns[instr_name] = np.concatenate((dns0,dns1,dns2))
            sky_temp[instr_name] = np.concatenate((sky_temp0, sky_temp1, sky_temp2))
            amb_temp[instr_name] = np.concatenate((amb_temp0, amb_temp1, amb_temp2))

    # Process FPI Data
    FPI_ut = {}
    FPI_wind = {}
    FPI_error = {}
    FPI_cloud = {}
    FPI_wq = {}
    FPI_tt = {}
    FPI_walpha = {}

    fpi_dt = dt + timedelta(days=-1)

    for fpi in sites_fpi:
        fname = fpi_repo + fpiinfo.get_instr_at(fpi,fpi_dt)[0] + '_' + fpi + '_' + fpi_dt.strftime('%Y%m%d') + '_' + sky_line_tag.lower() + '.npz'
        if (exists(fname) == False) and (sky_line_tag == 'XR'):
            # Some of the older npz files for redline don't have a tag
            fname = fpi_repo + fpiinfo.get_instr_at(fpi,fpi_dt)[0] + '_' + fpi + '_' + fpi_dt.strftime('%Y%m%d') + '.npz'

        if exists(fname):
            FPI_ut[fpi], FPI_wind[fpi], FPI_error[fpi], FPI_cloud[fpi], FPI_wq[fpi] = read_FPI(fname)
            FPI_tt[fpi] = {}
            FPI_walpha[fpi] = {}
            for d in FPI_ut[fpi].keys():
                FPI_tt[fpi][d] = np.array([toTimestamp(d.astimezone(pytz.utc)) for d in FPI_ut[fpi][d]])

                FPI_walpha[fpi][d] = FPI_wq[fpi][d].copy()
                FPI_walpha[fpi][d][FPI_walpha[fpi][d] == 2] = 0.2
                FPI_walpha[fpi][d][FPI_walpha[fpi][d] == 1] = 0.5
                FPI_walpha[fpi][d][FPI_walpha[fpi][d] == 0] = 1.0

    all_times = []
    for k in times.keys():
        if times[k] is None:
            continue
        all_times += list(times[k])

    all_times.sort()

    all_fpi_times = []
    for k in FPI_ut.keys():
        for d in FPI_ut[k].keys():
            all_fpi_times += list(FPI_ut[k][d])

    all_fpi_times.sort()

    if sky_line_tag == 'XR':
        # Normalize images
        for s in IM3Dfilt.keys():
#            n98 = np.nanpercentile(IM3Dfilt[s],98)
#            n2 = np.nanpercentile(IM3Dfilt[s],2)
#            temp = (IM3Dfilt[s]-n2)/n98
#            IM3Dfilt[s] = temp    
            for i, t in enumerate(times[s]):
#                print(np.shape(IM3Dfilt[s]))
#                print(i,t)
                im = IM3Dfilt[s][:,:,i]
                n98 = np.nanpercentile(im,98)
                n2 = np.nanpercentile(im,2)
                IM3Dfilt[s][:,:,i] = (im-n2)/n98

    # Values for the colorbar
    p98 = []
    p2 = []

    for s in IM3Dfilt.keys():
        # Find the index of images with moon down and use those to find autocontrast settings
        s_info = asiinfo.get_site_info(s)
        moon = ephem.Moon
        obs = ephem.Observer()
        obs.lat = str(s_info['Location'][0])
        obs.lon = str(s_info['Location'][1])
        
        moon_down = []
        for t in times[s]:
            obs.date = t
            moon = ephem.Moon(obs)
            moon_down.append(moon.alt.real < (el_cutoff)*np.pi/180.)

        # if all times are moon up (False) use all data
        if np.sum(moon_down) == 0:
            moon_down = np.full(np.size(moon_down), True)

        print(s)
        p98.append(np.nanpercentile(IM3Dfilt[s][:,:,moon_down],98))
        p2.append(np.nanpercentile(IM3Dfilt[s][:,:,moon_down],2))

    cbar_min = np.max([np.median(p2),0]) # Was Min # Force a 0 floor
    cbar_max = np.median(p98) # Was Max

    # Plotting
    fig = plt.figure(figsize=(8,4))
    fig.patch.set_facecolor('white')
    if sky_line_tag == 'XG':
        spec=gridspec.GridSpec(ncols=2,nrows=2,figure=fig,
                                left=0.04,right=0.94,bottom=0.1,top=0.90,
                                wspace=0.05,hspace=0.25,
                                width_ratios=[2.,2,])
    elif sky_line_tag == 'XR':
        spec=gridspec.GridSpec(ncols=2,nrows=2,figure=fig,
                                left=0.04,right=0.94,bottom=0.1,top=0.90,
                                wspace=0.05,hspace=0.25,
                                width_ratios=[1.5,2,])

    # The lon, lat range to plot in the map
    if sky_line_tag == 'XG':
        if 'bdr' in sites_asi:
            extent=[-122,-95,25,48]
        else:
            extent=[-124,-102,30,48]
    elif sky_line_tag == 'XR':
        extent=[-125, -80, 25, 55]

    my_colors = {}
    my_colors['cvo'] = '#1f77b4'
    my_colors['blo'] = '#ff7f0e'
    my_colors['low'] = '#9467bd'
    my_colors['uao'] = '#8c564b'

    for t in np.unique(all_times):
        fig.clf()

        # Plot images on the map
        axes00 = fig.add_subplot(spec[:,0], projection=crs.Orthographic(np.nanmean(extent[:2]),np.nanmean(extent[2:])))    
        axes00.add_feature(feature.COASTLINE)
        axes00.add_feature(feature.STATES,alpha=0.2,zorder=1)

        for instr_name in times.keys():
            i = np.argwhere(np.array(times[instr_name]) == t)
            if len(i) == 1:
                # We have an image, see if we also have cloud data
                my_alpha = 1
                if instr_name in dns.keys():
                    if len(dns[instr_name]) > 0:
                        closest_cloud_dt = min(dns[instr_name],key = lambda x: abs(x - t.astimezone(pytz.utc)))
                        cloud_i = np.where(dns[instr_name] == closest_cloud_dt)
                        xx1 = -20.
                        yy1 = 1.
                        xx2 = -10.
                        yy2 = 0.2

                        my_a = (yy1-yy2)/(xx1-xx2)
                        my_b = yy1 - xx1*(yy1-yy2)/(xx1-xx2)

                        my_alpha = np.clip(my_a*sky_temp[instr_name][cloud_i]+my_b,yy2,yy1)[0]
                        print(t,instr_name,my_alpha)
#                pc = axes00.pcolormesh(lon[instr_name],lat[instr_name],IM3Dfilt[instr_name][:,:,i[0][0]],transform=crs.PlateCarree(),vmin=cbar_min,vmax=cbar_max,cmap='gray',alpha=my_alpha)
                pc = axes00.pcolormesh(lon[instr_name],lat[instr_name],IM3Dfilt[instr_name][:,:,i[0][0]],transform=crs.PlateCarree(),vmin=cbar_min,vmax=cbar_max,cmap='gray')

        # Set the title and limits of the map
        axes00.set_title(t)
        axes00.set_extent(extent,crs=crs.PlateCarree())

        # Plot meridional winds
        axes01 = fig.add_subplot(spec[0,1])
        axes02 = fig.add_subplot(spec[1,1])

        for fpi in FPI_ut.keys():
#            axes01.plot(FPI_ut[fpi]['North'],FPI_wind[fpi]['North'], alpha=0.3,color=my_colors[fpi], label='_nolegend_')
#            axes01.scatter(FPI_ut[fpi]['North'],FPI_wind[fpi]['North'], alpha=FPI_walpha[fpi]['North'],ec=None, color=my_colors[fpi],label=fpi)
#            axes02.plot(FPI_ut[fpi]['West'],FPI_wind[fpi]['West'], alpha=0.3,color=my_colors[fpi], label='_nolegend_')
#            axes02.scatter(FPI_ut[fpi]['West'],FPI_wind[fpi]['West'], alpha=FPI_walpha[fpi]['West'],ec=None, color=my_colors[fpi],label=fpi)
            axes01.plot(FPI_ut[fpi]['North'],FPI_wind[fpi]['North'],color=my_colors[fpi], label='_nolegend_')
            axes01.scatter(FPI_ut[fpi]['North'],FPI_wind[fpi]['North'],ec=None, color=my_colors[fpi],label=fpi)
            axes02.plot(FPI_ut[fpi]['West'],FPI_wind[fpi]['West'],color=my_colors[fpi], label='_nolegend_')
            axes02.scatter(FPI_ut[fpi]['West'],FPI_wind[fpi]['West'],ec=None, color=my_colors[fpi],label=fpi)

        if sky_line_tag == 'XR':
            ymin = -300
            ymax = 300
        elif sky_line_tag == 'XG':
            ymin = -200
            ymax = 200

        axes01.set_ylim([ymin,ymax])
        axes01.axhline(y=0, color='k',linewidth=1)
        axes01.axvline(x=t, color='r',linewidth=1)
        axes01.xaxis.set_major_locator(dates.HourLocator(interval = 2))
        axes01.xaxis.set_minor_locator(dates.HourLocator())
        axes01.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        if sky_line_tag == 'XG':
            axes01.set_title('Green Line Winds')
        elif sky_line_tag == 'XR':
            axes01.set_title('Red Line Winds')
        axes01.tick_params(axis='both',which='major',size=6,direction='in',right=True,top=True,labelright=True,labelleft=False)
        axes01.tick_params(axis='both',which='minor',size=3,direction='in',right=True,top=True,labelright=True,labelleft=False)
        axes01.set_ylabel('Northward Winds [m/s]')
        axes01.legend(loc='upper right',ncol=2, fontsize=8)
        if len(all_fpi_times) > 0:
            axes01.set_xlim([np.unique(all_fpi_times)[0],np.unique(all_fpi_times)[-1]])

        axes02.set_ylim([ymin,ymax])    
        axes02.axhline(y=0, color='k',linewidth=1)
        axes02.axvline(x=t, color='r',linewidth=1)
        axes02.xaxis.set_major_locator(dates.HourLocator(interval = 2))
        axes02.xaxis.set_minor_locator(dates.HourLocator())
        axes02.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        axes02.tick_params(axis='both',which='major',size=6,direction='in',right=True,top=True,labelright=True,labelleft=False)
        axes02.tick_params(axis='both',which='minor',size=3,direction='in',right=True,top=True,labelright=True,labelleft=False)
        axes02.yaxis.tick_right()
        axes02.yaxis.set_label_position("left")
        axes02.set_ylabel('Eastward Wind [m/s]')
        axes02.set_xlabel('UT [hrs]')
        if len(all_fpi_times) > 0:   
            axes02.set_xlim([np.unique(all_fpi_times)[0],np.unique(all_fpi_times)[-1]])

        # Plot quiver
        # Scaling factor
        if sky_line_tag == 'XG':
            scale=25.
        elif sky_line_tag == 'XR':
            scale=100.

        xoff,yoff=8e4, -1.3e5
        width=0.015
        headwidth = 4
        headlength = 4
        headaxislength= headlength-1
        minshaft = 2
        sc = axes00.bbox.width/fig.dpi/(10.*scale)

        obj = None
        for fpi in FPI_ut.keys():
            valid_t = False

            for d in FPI_tt[fpi].keys():
                if np.min(abs(toTimestamp(t)-FPI_tt[fpi][d])) < 60*60:
                    valid_t = True

            if valid_t:
                # Get interpolated u, v
                u = np.interp(toTimestamp(t),FPI_tt[fpi]['West'],FPI_wind[fpi]['West'])*sc
                v = np.interp(toTimestamp(t),FPI_tt[fpi]['North'],FPI_wind[fpi]['North'])*sc
#                a = np.interp(toTimestamp(t),FPI_tt[fpi]['North'],FPI_walpha[fpi]['North'])
                a = 1

                glat,glon, x = fpiinfo.get_site_info(fpi)['Location']

                obj = axes00.quiver(np.array([glon]),np.array([glat]),np.array([u]),np.array([v]),
                                        angles='uv', scale_units='inches', scale=1,width=width,
                                        pivot='tail', headwidth=headwidth, headlength=headlength, alpha=a,
                                        minshaft=minshaft, headaxislength=headaxislength, color=my_colors[fpi], transform=crs.PlateCarree(),zorder=1000)

        x0,x1,y0,y1=axes00.get_extent()
        if obj is not None:
            quiverobj1=axes00.quiverkey(obj, x0+xoff, y0+yoff, sc*scale, r'$%i\,\frac{m}{s}$'%scale,labelsep=0,color="black", coordinates='data',labelpos='N', transform=crs.PlateCarree(),fontproperties={'size': 8})

        # Plot the colorbar
        cax00 = inset_axes(
            axes00,
            width="50%",
            height="3%",
            loc="lower left",
        )
        fig.colorbar(pc, cax=cax00, orientation="horizontal", ticks=[cbar_min,(cbar_min+cbar_max)/2., cbar_max])
        cax00.xaxis.tick_top()
        cax00.xaxis.tick_top()
        cax00.set_xticklabels(cax00.get_xticklabels(),fontsize=8)

        spec.tight_layout(fig)

        savename = folderpngs + 'MANGO_%s.png' % t.strftime('%Y%m%d_%H%M%S')
        print(savename)
        plt.savefig(savename,dpi=300)
        print(t)

    # Generating video from PNG frames
    linetag="red" if sky_line_tag=='XR' else "green"
    mp4file='%s/DASI_%s_%s.mp4'%(outfolder, dt.strftime('%Y%m%d'),linetag)
    cmd = '/usr/bin/ffmpeg -framerate 15 -pattern_type glob -i \"%s*.png\" -c:v mpeg4 -q:v 1 -y %s'%(folderpngs, mp4file)
    os.system(cmd)

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

    # Run the code
    try:
        MakeSummaryMovies(year, doy, sky_line_tag, download_data = download_data)
    except:
        print("Error %d %d %s" % (doy, year, sky_line_tag))
#    MakeSummaryMovies(year, doy, sky_line_tag, download_data = False, sites_asi = ['mro'])
#    MakeSummaryMovies(year, doy, sky_line_tag, download_data = False)
