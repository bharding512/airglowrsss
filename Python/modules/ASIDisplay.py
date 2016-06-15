import glob
import Image
import ASI
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import numpy.ma as ma
from scipy import ndimage
from scipy import interpolate
from matplotlib import dates
import bottleneck as bn
import TifImagePlugin

def DisplayRaw(f, cmin=None, cmax=None, dark=None, flips=None, info=True, sitename=None, filt=None, kernel_size = 1):
    # Function to display information about a raw ASI TIF image and
    # plot the raw data
    #
    # INPUTS:
    #   f - full path name to image to display
    #
    # OPTION INPUTS:
    #   cmin, cmax - color limits (min, max) to scale the resultant image by
    #   dark - a dark image to subtract from the requested image
    #   flips - list of flips to get the orientation of the images right.
    #           lr : fliplr,  ud : flipud
    #   info - flag on whether to print info or not
    #   sitename, filt - sitename and filter to put in image title.  If not given, reads from header
    #   kernel_size - size of the median filter kernel used (default = 5)
    #
    # OUTPUT:
    #   f - a reference to the figure created
    #
    # HISTORY:
    #   Written by Jonathan J. Makela on 2 July 2013

    # Load the image
    d = Image.open(f)

    # Should we plot out header info?
    if info:
        print 'Filename: %s' % f
        print 'Exposure time: %.1f s' % d.info['ExposureTime']
        print 'CCD Temperature: %.1f C' % d.info['CCDTemperature']

    # Define the figure and plot
    f = plt.figure(num=None,figsize=(6,6),dpi=80)
    p = f.add_subplot(111)

    # Create the image array from the data and dark image (if requested)
    if dark is None:
#        im = np.reshape(d.getdata(), d.size)
	im = np.array(d.getdata(), np.uint16).reshape(d.size)

        # No Dark image used.  Run a small median filter to remove any
        # noise spikes from the CCD
        im = ndimage.filters.median_filter(im,size=kernel_size)
    else:
#        im = np.reshape(d.getdata(), d.size)-dark
	im = np.array(d.getdata(), np.uint16).reshape(d.size) - dark

    # Do we need to flip the image?
    if flips is not None:
        for flip in flips:
            if flip == 'lr':
                im = np.fliplr(im)
            if flip == 'ud':
                im = np.flipud(im)

    # Create the image display
    a = p.matshow(im, cmap=cm.gray)

    # Set color limits
    if cmin is not None and cmax is not None:
	a.set_clim([cmin,cmax])
    else:
	a.set_clim(np.percentile(im,(5,95)))


#    if dark is None:
#        if cmin is not None:
#            a.set_clim([cmin,cmax])
#        else:
##            a.set_clim(d.info['pmin'],d.info['pmax'])
#	    a.set_clim(np.percentile(im,(5,95)))
#    else:
#        if cmin is not None:
#            a.set_clim([cmin-np.median(dark),cmax-np.median(dark)])
#        else:
#            a.set_clim([d.info['pmin']-np.median(dark),d.info['pmax']-np.median(dark)])
#	    a.set_clim(np.percentile(im,5)-np.median(dark),np.percentile(im,95)-np.median(dark))

    # Create the title and colorbar
    if sitename is None:
        sitename = d.info['Observatory'][0]
    if filt is None:
        filt = d.info['Filter']
        
    p.set_title('%s (%s): %s' % (sitename, filt, d.info['LocalTime'].strftime('%d %b %Y %H:%M:%S LT')))
    f.colorbar(a,orientation='vertical',shrink=0.8)

    # Turn off the axis
    p.axes.get_xaxis().set_visible(False)
    p.axes.get_yaxis().set_visible(False)

    plt.close()

    return f

def DisplayMap(f, m, lat, lon, cmin=None, cmax=None, dark=None, sitename=None, filt=None, kernel_size = 5, displayUT=False, displayColorbar=True):
    # Function to display an image on a provide map projection.
    #
    # INPUTS:
    #   f - full path name to the image to display
    #   m - a fully defined Basemap onto which the image will be displayed.
    #       m = Basemap(llcrnrlon=275.,llcrnrlat=-30.,urcrnrlon=300.,urcrnrlat=-13.,
    #            projection='merc', area_thresh=1000,
    #            resolution='i',lat_1=35.,lat_2=55,lat_0=40,lon_0=-85.)
    #   lat, lon - the lat, lon for each pixel in the image (must be the same size as the
    #              imaging data defined in f.  These values should have NaN in pixels outside
    #              of the field of view of the imager.
    #
    # OPTIONAL INPUTS:
    #   cmin, cmax - color limits (min, max) to scale the resultant image by
    #   dark - a dark image to subtract from the requested image
    #   sitename, filt - sitename and filter to put in image title.  If not given, reads from header
    #   kernel_size - size of median filter kernel used (default = 5)
    #
    # HISTORY:
    #   Written by Jonahtan J. Makela on 9 July 2013

    # Load the image
    d = Image.open(f)

    # Create the image array from the data and dark image (if requested)
    if dark is None:
#        im = np.reshape(d.getdata(), d.size)
	im = np.array(d.getdata(), np.uint16).reshape(d.size)

        # No Dark image used.  Run a small median filter to remove any
        # noise spikes from the CCD
    #    im = ndimage.filters.median_filter(im,size=kernel_size)
    else:
#        im = np.reshape(d.getdata(), d.size)-dark
	im = np.array(d.getdata(), np.uint16).reshape(d.size) - dark

    im = ndimage.filters.median_filter(im,size=kernel_size)

    # Find NaN in the projection arrays and mask the image data
    zdata = ma.masked_where(np.isnan(lat),im)

    # Convert to the map's coordinates
    xpt,ypt = m(lon,lat)

    # Create the plot
    a = m.pcolormesh(xpt,ypt,zdata,shading='flat',cmap=cm.gray)

    # Set the colorlimits 
    if cmin is not None and cmax is not None:
	a.set_clim([cmin,cmax])
    else:
	a.set_clim(np.percentile(im,(5,95)))

#    if dark is None:
#        if cmin is not None:
#            a.set_clim([cmin,cmax])
#        else:
##            a.set_clim(d.info['pmin'],d.info['pmax'])
#	    a.set_clim(np.percentile(im,(5,95)))
#    else:
#        if cmin is not None:
#            a.set_clim([cmin,cmax])
#        else:
#            a.set_clim([d.info['pmin']-np.median(dark),d.info['pmax']-np.median(dark)])
#	    a.set_clim(np.percentile(im,5)-np.median(dark),np.percentile(im,95)-np.median(dark))

    # Create the title and colorbar
    if sitename is None:
        sitename = d.info['Observatory'][0]
    if filt is None:
        filt = d.info['Filter']
        
    if displayUT:
       plt.title('%s (%s): %s' % (sitename, filt, d.info['UniversalTime'].strftime('%d %b %Y %H:%M:%S UT')))
    else:
       plt.title('%s (%s): %s' % (sitename, filt, d.info['LocalTime'].strftime('%d %b %Y %H:%M:%S LT')))

    if displayColorbar:
       plt.colorbar(a,orientation='vertical',shrink=0.9)

def Keogram(files, lat, lon, target_lat, target_lon, darks=None, sitename=None, filt=None, cmin=None, cmax=None, kernel_size = 5):
    # Function to generate keograms through the requested images at the requested lat/lon.
    #
    # INPUTS:
    #   files - list of files to generate keogram from (Assumed a glob list)
    #   lat, lon - the lat, lon for each pixel in the image (must be the same size as the imaging
    #             data defined in files.  These values should have NaN in pixels outside of the
    #             field of view of the imager
    #   target_lat, target_lon - the lat/lon at which the cuts will be made through the images to
    #                            generate the keograms.
    #   cmin, cmax - the min/max values for the color map
    #
    # OPTIONAL INPUTS:
    #   darks - list of files to be considered for dark image subtraction.
    #   sitename, filt - sitename and filter to put in image title.  If not given, reads from header
    #   kernel_size - size of the kernel for the median filter applied to the data (default = 5)
    #
    # OUTPUT:
    #   f - a reference to the figure generated
    #
    # HISTORY:
    #   Written by Jonathan J. Makela on 10 July 2013

    # Create a master dark image if requested
    if darks is not None:
        # Find out the size of the image
        d = Image.open(darks[0])
        all_dark = np.zeros_like(np.reshape(d.getdata(), d.size))

        # Sum up all of the dark images
        for dark in darks:
            d = Image.open(dark)
#            all_dark = all_dark + np.reshape(d.getdata(), d.size)
	    all_dark = all_dark + np.array(d.getdata(), np.uint16).reshape(d.size)

        # Divide by the number of dark images to create a mean image
        darks = all_dark/len(darks)

    # Variables to hold the generated keograms
    keo_lat = []
    keo_lon = []
    ut = []
    all_min = []
    all_max = []

    for f in files:
        # Open the file
        d = Image.open(f)

        # Grab the sitename and filter if it wasn't provided
        if sitename is None:
            sitename = d.info['Observatory'][0]
        if filt is None:
            filt = d.info['Filter']

        if cmin is None:
            all_min.append(d.info['pmin'])
            if d.info['pmax'] < 0:
                # DFJ Moon saturation fix...
                all_max.append(32768)
            else:
                all_max.append(d.info['pmax'])

        # Get the Universal Time of image
        ut.append(d.info['UniversalTime'])

        # Create the image array from the data and dark image (if requested)
        if darks is None:
#            im = np.reshape(d.getdata(), d.size)
	    im = np.array(d.getdata(), np.uint16).reshape(d.size)
        else:
#            im = np.reshape(d.getdata(), d.size)-darks*1.0
	    im = np.array(d.getdata(), np.uint16).reshape(d.size) - darks*1.0
            
        # Median filter the data
        d = ndimage.filters.median_filter(im,size=kernel_size)

        # Cut through the image at the requested latitude
        zlat = lat[np.where(abs(lat-target_lat) < 0.5)]
        zlon = lon[np.where(abs(lat-target_lat) < 0.5)]
        zdata = d[np.where(abs(lat-target_lat) < 0.5)]
        lo = np.arange(zlon.min(),zlon.max(),.1)
        yi = np.array([target_lat])
        #GD = interpolate.griddata(zip(zlon.flatten(),zlat.flatten()),zdata.flatten(),np.meshgrid(lo,yi),method='linear')
        GD = interpolate.griddata(zip(zlon.flatten(),zlat.flatten()),zdata.flatten(),tuple(np.meshgrid(lo,yi)),method='linear')
        keo_lon.append(GD.flatten())

        # Cut through the image at the requested longitude
        zlat = lat[np.where(abs(lon-target_lon) < 0.5)]
        zlon = lon[np.where(abs(lon-target_lon) < 0.5)]
        zdata = d[np.where(abs(lon-target_lon) < 0.5)]
        xi = np.array([target_lon])
        la = np.arange(zlat.min(),zlat.max(),0.1)
        #GD = interpolate.griddata(zip(zlon.flatten(),zlat.flatten()),zdata.flatten(),np.meshgrid(xi,la),method='linear')
        GD = interpolate.griddata(zip(zlon.flatten(),zlat.flatten()),zdata.flatten(),tuple(np.meshgrid(xi,la)),method='linear')
        keo_lat.append(GD.flatten())

    # Calculate a guess for cmin/cmax
    if cmin is None:
        # Cast to an array for easy statistics
        all_min = np.array(all_min)
        all_max = np.array(all_max)
        if darks is None:
            cmin = np.median(all_min)
            cmax = np.median(all_max)
        else:
            cmin = np.median(all_min)-np.median(darks)
            cmax = np.median(all_max)-np.median(darks)
        m = bn.nanmedian(keo_lat)
        s = bn.nanstd(keo_lat)
        cmin = np.max([cmin,m-2.*s,0])
        cmax = np.min([cmax,m+2.*s])

    # Cast the outputs as arrays
    keo_lon = np.array(keo_lon)
    keo_lat = np.array(keo_lat)
    ut = np.array(ut)

    # Generate the image and plot the first keogram
    f = plt.figure(num=None)
    p = f.add_subplot(211)
    a = p.imshow(keo_lon.transpose(), aspect='auto',extent=(dates.date2num(ut[0]),dates.date2num(ut[-1]),lo[-1],lo[0]),cmap=cm.gray)
    ax = a.get_axes()
    ax.xaxis_date()

    # Format the axis
    hfmt = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_locator(dates.HourLocator(interval=1))
    ax.xaxis.set_major_formatter(hfmt)
    plt.ylabel('Longitude')

    # Format the title
    if ut[0].date() == ut[-1].date():
        # Data does not span a day
        tstr = '%s (%s): %s' % (sitename, filt, ut[0].strftime('%d %b %Y'))
    elif ut[0].date().month == ut[-1].date().month:
        # Date spans a day, but not a month:
        tstr = '%s (%s): %s-%s' % (sitename,filt, ut[0].strftime('%d'), ut[-1].strftime('%d %b %Y'))
    elif ut[0].date().year == ut[-1].date().year:
        # Date spans a month, but not a year:
        tstr = '%s (%s): %s-%s' % (sitename,filt, ut[0].strftime('%d %b'), ut[-1].strftime('%d %b %Y'))
    else:
        # Date spans a year
        tstr = '%s (%s): %s-%s' % (sitename,filt, ut[0].strftime('%d %b %Y'), ut[-1].strftime('%d %b %Y'))
    
    plt.title(tstr)
    a.set_clim([cmin,cmax])

    # Generate the second keogram
    p = f.add_subplot(212)
    a = p.imshow(keo_lat.transpose(), aspect='auto',extent=(dates.date2num(ut[0]),dates.date2num(ut[-1]),la[-1],la[0]),cmap=cm.gray)
    ax = a.get_axes()
    ax.xaxis_date()

    # Format the axis
    hfmt = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_locator(dates.HourLocator(interval=1))
    ax.xaxis.set_major_formatter(hfmt)
    plt.xlabel('UT')
    plt.ylabel('Latitude')
    a.set_clim([cmin,cmax])
    plt.gca().invert_yaxis()
    
    return f
