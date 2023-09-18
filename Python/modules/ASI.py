import matplotlib
matplotlib.use('AGG')

import ASIDisplay
import numpy as np
import os
try:
    import Image
except ImportError:
    from PIL import Image
import matplotlib.pyplot as plt
import string
import random

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    # Generate a random string of uppercase, lowercase, and numbers.
    # Taken from http://stackoverflow.com/questions/2257441/python-random-string-generation-with-upper-case-letters-and-digits
    return ''.join(random.choice(chars) for x in range(size))

def RawMovie(files, movie_name=None, cmin=None, cmax=None, darks=None, flips=None,sitename=None,filt=None):
    # Function to generate a movie of the requested ASI images.
    #
    # INPUTS:
    #   files - list of files to generate movie from (assumed a glob list)
    #
    # OPTION INPUTS:
    #   movie_name - file name of movie to be generated.  If None provided
    #                we will generate a name based on the information in
    #                the header file of the first image passed
    #   cmin, cmax - color limits (min, max) to scale all images by. If
    #                None provided, we will generate these based on the
    #                median values of the autoscale values in each
    #                image in the files list.
    #   darks - list of files to be considered for dark image subtraction
    #   flips - list of flips to get the orientation of the images right.
    #           lr : fliplr,  ud : flipud
    #   sitename, filt - sitename and filter to put in image title.  If not given, reads from header
    #
    # OUTPUT:
    #   generates a movie saved to movie_name
    #
    # HISTORY:
    #   Written by Jonathan J. Makela on 2 July 2013

    # Create a master dark image if requested
    if darks:
        # Find out the size of the image
        d = Image.open(darks[0])
        all_dark = np.zeros_like(np.reshape(d.getdata(), d.size))

        # Sum up all of the dark images
        for dark in darks:
            d = Image.open(dark)
#            all_dark = all_dark + np.reshape(d.getdata(), d.size)
            all_dark = all_dark + np.array(d.getdata(),np.uint16).reshape(d.size)

        # Divide by the number of dark images to create a mean image
        darks = all_dark/len(darks)
    else:
        darks = None

    # If no cmin, cmax given, calculate best values
    if cmin is None:
        all_min = []
        all_max = []

        for f in files:
            d = Image.open(f)
            im = np.array(d.getdata(),np.uint16).reshape(d.size)

#            all_min.append(d.info['pmin'])
#            all_max.append(d.info['pmax'])
            all_min.append(np.percentile(im,5))
            all_max.append(np.percentile(im,95))

        # Cast to an array for easy statistics
        all_min = np.array(all_min)
        all_max = np.array(all_max)

        cmin = np.median(all_min)
        cmax = np.median(all_max)

    # Now generate the individual png files
    png_names = []

    # Generate a random key (needed if multiple threads are generating png files
    id_str = id_generator(6)

    for f in files:
        # Grab the path and filename given
        pa,fi = os.path.split(f)

        # Create a temporay file name
        png_name = pa + 'tmp_' + fi[0:-4] + '_' + id_str + '.png'

        # Plot the image and save it
        p = ASIDisplay.DisplayRaw(f,cmin=cmin,cmax=cmax,dark=darks, flips=flips, info=False,sitename=sitename,filt=filt)
        p.savefig(png_name, bbox_inches='tight', pad_inches=0.1)

        # Created this image, append the name to the list to work with
        png_names.append(png_name)

    # Now that the frames are generated, create a movie
    # Create the movie in a two-pass mpeg4 (see http://matplotlib.sourceforge.net/rc/v1.1.1rc2/examples/old_animation/movie_demo.html and
    # http://mariovalle.name/mencoder/mencoder.html)
    command = ('mencoder',
           'mf://tmp_*.png',
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

#    command = ('mencoder',
#               'mf://' + pa + '/tmp_*_' + id_str + '.png',
#               '-mf',
#               'fps=5',
#               '-o',
#               '/dev/null',
#               '-ovc',
#               'x264',
#               '-x264encopts',
#               'bitrate=500:pass=1')

    os.spawnvp(os.P_WAIT, 'mencoder', command)

    command = ('mencoder',
           'mf://tmp_*.png',
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
           movie_name)

#    command = ('mencoder',
#               'mf://' + pa + '/tmp_*_' + id_str + '.png',
#               '-mf',
#               'fps=5',
#               '-o',
#               movie_name,
#               '-ovc',
#               'x264',
#               '-x264encopts',
#               'bitrate=500:pass=2')

    os.spawnvp(os.P_WAIT, 'mencoder', command)

    # Remove the temporary png files
    for f in png_names:
        os.remove(f)

    os.remove('divx2pass.log')
    os.remove('divx2pass.log.mbtree')

def MapMovie(files, m, lat, lon, mask, movie_name=None, cmin=None, cmax=None, darks=None, sitename=None, filt=None, kernel_size=5, displayCountries=True,displayUT=False):
    # Function to generate a movie of the requested ASI images.
    #
    # INPUTS:
    #   files - list of files to generate movie from (assumed a glob list)
    #   m - a fully defined Basemap onto which the image will be displayed.
    #       m = Basemap(llcrnrlon=275.,llcrnrlat=-30.,urcrnrlon=300.,urcrnrlat=-13.,
    #            projection='merc', area_thresh=1000,
    #            resolution='i',lat_1=35.,lat_2=55,lat_0=40,lon_0=-85.)
    #   lat, lon - the lat, lon for each pixel in the image (must be the same size as the
    #              imaging data defined in f.  These values should have NaN in pixels outside
    #              of the field of view of the imager.
    #
    # OPTION INPUTS:
    #   movie_name - file name of movie to be generated.  If None provided
    #                we will generate a name based on the information in
    #                the header file of the first image passed
    #   cmin, cmax - color limits (min, max) to scale all images by. If
    #                None provided, we will generate these based on the
    #                median values of the autoscale values in each
    #                image in the files list.
    #   darks - list of files to be considered for dark image subtraction.
    #   sitename, filt - sitename and filter to put in image title.  If not given, reads from header
    #   kernel_size - size of median filter to apply to data (assume square)
    #   displayCountries - draw country borders
    #   displayUT - display UT rather than LT
    #
    # OUTPUT:
    #   generates a movie saved to movie_name
    #
    # HISTORY:
    #   Written by Jonathan J. Makela on 9 July 2013

    # Create a master dark image if requested
    if darks:
        # Find out the size of the image
        d = Image.open(darks[0])
        all_dark = np.zeros_like(np.reshape(d.getdata(), d.size))

        # Sum up all of the dark images
        for dark in darks:
            d = Image.open(dark)
#            all_dark = all_dark + np.reshape(d.getdata(), d.size)
            all_dark = all_dark + np.array(d.getdata(),np.uint16).reshape(d.size)

        # Divide by the number of dark images to create a mean image
        darks = all_dark/len(darks)
    else:
        darks = None

    # If no cmin, cmax given, calculate best values
    if cmin is None:
        all_min = []
        all_max = []

        for f in files:
            d = Image.open(f)
            im = np.array(d.getdata(),np.uint16).reshape(d.size)

            if darks is None:
                print('No Dark Subtraction')
            else:
                im = im - darks

#            if d.info['pmax'] > 0:
#                all_min.append(d.info['pmin'])
#                all_max.append(d.info['pmax'])
            all_min.append(np.percentile(im,5))
            all_max.append(np.percentile(im,95))

        # Cast to an array for easy statistics
        all_min = np.array(all_min)
        all_max = np.array(all_max)

        cmin = np.median(all_min)
        cmax = np.median(all_max)

    # Now generate the individual png files
    png_names = []

    # Generate a random key (needed if multiple threads are generating png files
    id_str = id_generator(6)

    pa = None

    for f in files:
        # Clear the figure
        plt.clf();

        # Grab the path and filename
        pa,fi = os.path.split(f)
        print(pa)

        # Create a temporary file name
        png_name = pa + 'tmp_' + fi[0:-4] + '_' + id_str + '.png'

        # Plot the image and save it
        ASIDisplay.DisplayMap(f,m,lat,lon,mask,cmin=cmin,cmax=cmax,dark=darks,sitename=sitename,filt=filt,kernel_size=kernel_size, displayUT=displayUT)

        # Display coasts, etc
        m.drawcoastlines()
        #m.drawstates()
        if displayCountries:
            m.drawcountries()

        # Draw parallels and meridians
        parallels = np.arange(-85.,85.,5.);
        m.drawparallels(parallels,labels=[True,False,False,True])
        meridians = np.arange(0.,360.,5.)
        m.drawmeridians(meridians,labels=[True,False,False,True])

        # Save the file
        plt.savefig(png_name, bbox_inches='tight', pad_inches=0.1)
        

        # Created this image, append the name to the list to work with
        png_names.append(png_name)

    # Now that the frames are generated, create a movie
    # Create the movie in a two-pass mpeg4 (see http://matplotlib.sourceforge.net/rc/v1.1.1rc2/examples/old_animation/movie_demo.html and
    # http://mariovalle.name/mencoder/mencoder.html)
#    command = ('mencoder',
#           'mf://' + pa + '/tmp_*.png',
#           '-mf',
#           'type=png:w=800:h=600:fps=5',
#           '-ovc',
#           'lavc',
#           '-lavcopts',
#           'vcodec=mpeg4:vpass=1:vbitrate=2343750:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq',
#           '-oac',
#           'copy',
#           '-nosound',
#           '-o',
#           '/dev/null')

    if pa is not None:

        # Generate a unique logfile name (so we don't get clashes when running multiple encoding processes
        passfile = pa + '/divx2pass' + id_str + '.log'

        command = ('mencoder',
                   'mf://' + pa + 'tmp_*_' + id_str + '.png',
                   '-mf',
                   'fps=5',
                   '-o',
                   '/dev/null',
                   '-ovc',
                   'x264',
                   '-x264encopts',
                   'bitrate=500:pass=1',
                   '-passlogfile',
                   passfile)

        os.spawnvp(os.P_WAIT, 'mencoder', command)

#    command = ('mencoder',
#           'mf://' + pa + '/tmp_*.png',
#           '-mf',
#           'type=png:w=800:h=600:fps=5',
#           '-ovc',
#           'lavc',
#           '-lavcopts',
#           'vcodec=mpeg4:vpass=2:vbitrate=2343750:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq',
#           '-oac',
#           'copy',
#           '-nosound',
#           '-o',
#           movie_name)
        command = ('mencoder',
                   'mf://' + pa + 'tmp_*_' + id_str + '.png',
                   '-mf',
                   'fps=5',
                   '-o',
                   movie_name,
                   '-ovc',
                   'x264',
                   '-x264encopts',
                   'bitrate=500:pass=2',
                   '-passlogfile',
                   passfile)

        os.spawnvp(os.P_WAIT, 'mencoder', command)

        os.remove(passfile)
        os.remove(passfile + '.mbtree')

    # Remove the temporary png files
    for f in png_names:
        os.remove(f)

def ConvertAzEl2LatLon(az,el,ht,rx_lat,rx_lon, horizon=-1.):
    # Function to convert from an azimuth/elevation grid to a
    # latitude/longitude grid given an unwarping height and a site location.
    # For elevations below the horizon, returns NaN.
    #
    # INPUTS:
    #   az - M x N array of azimuths to be converted [degrees]
    #   el - M x N array of elevation to be converted [degrees]
    #   ht - height to be used in the conversion [km]
    #   rx_lat, rx_lon - [latitude, longitude] of the site [degrees]
    #
    # OPTION INPUTS:
    #   horizon - elevation angles below this will be masked out.  Defaults to -1 degree
    #
    # OUTPUTS:
    #   lat - M x N array of latitudes [degrees]
    #   lon - M x N array of longitudes [degrees]
    #
    # HISTORY:
    #   Converted from MATLAB by Jonathan J. Makela on 9 July 2013
    #   Updated by BJH to allow scalar values of az,el on 5 Feb 2013

    # Earth radius in km
    Re = 6371.2

    # Convert all input angles from degrees to radians
    el_r = el*np.pi/180.
    az_r = az*np.pi/180.
    lat_r = rx_lat*np.pi/180.
    lon_r = rx_lon*np.pi/180.

    # Calculate the differential angle, alpha
    temp = np.cos(el_r)/(1+(ht/Re))
    alpha = np.arccos(temp) - el_r

    # Calculate the pierce point latitude
    temp = np.sin(lat_r) * np.cos(alpha) + np.cos(lat_r)*np.cos(az_r)*np.sin(alpha)
    lat_r = np.arcsin(temp)

    # Calculate the pierce point longitude
    temp = np.sin(alpha) * np.sin(az_r)/np.cos(lat_r)
    lon_r = np.arcsin(temp) + lon_r

    # Convert back to radian measurements
    lat = lat_r*180./np.pi
    lon = lon_r*180./np.pi

    # NaN out values below the horizon
    # Need to check if value is a scalar or a list/array/whatever.
    # Why is there not a cleaner way of checking this?
#    import collections
#    if isinstance(lat, collections.abc.Iterable): # it's an array/list
#        lat[el < horizon] = np.nan
#        lon[el < horizon] = np.nan
#    else: # it's a scalar
#        if el < horizon:
#            lat = np.nan
#            lon = np.nan


    return lat, lon
