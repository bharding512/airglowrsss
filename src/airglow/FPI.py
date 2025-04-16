#!/usr/bin/python
# Filename: FPI.py

import matplotlib
#matplotlib.use('AGG')

import subprocess
import re
try:
    import Image
except:
    from PIL import Image
import matplotlib as mpl
from math import pi, floor, sqrt
import numpy as np
from lmfit import Minimizer, Parameters, Parameter, minimize
import scipy
from scipy import interpolate
import glob as glob
import time
import matplotlib.pyplot as pyplot
#from matplotlib.mlab import prctile
import matplotlib.dates as dates
from scipy import ndimage
import mahotas
import pytz
from pytz import timezone
import datetime
from scipy import ndimage
import mahotas
from matplotlib.font_manager import FontProperties
import ephem
#import logging
# Since FPI.py is intended to be general purpose, this shouldn't be necessary:
import airglow.fpiinfo as fpiinfo
# but it's the easiest way to deal with the chain-vs-direct drive issue in Peru.
# Is there a better way?
import airglow.ImgImagePlugin
import h5py
from airglow.exceptions import *
from airglow.warning_log import WarningLog
from os.path import basename

def sort_look_directions(valid_az, valid_ze, az, ze, tol):
    '''
    Sort the actual look directions (az,ze) into categories defined by valid_az and
    valid_ze. Directions are said to match if they are within a great-circle distance of
    tol degrees.

    INPUTS:
        valid_az - with valid_ze, an array that specifies the valid look directions [deg]
        valid_ze - with valid_az, an array that specifies the valid look directions [deg]
        az - with ze, an array of actual look directions over the night [deg]
        ze - with az, an array of actual look directions over the night [deg]
        tol - the deviation in degrees between actual and expected directions, above
              which the actual direction will not be recognized as valid

    OUTPUTS:
        direc - a vector of length len(az). The entries are lists containing the indicies of
                the directions which matched the observation, where indices refer to valid_az, valid_ze.
                For example, if direc == [ [4], [7], [2,3], [] ], then (az[0], ze[0]) matched
                (valid_az[4], valid_ze[4]), (az[1], ze[1]) matched (valid_az[7], valid_ze[7]),
                (az[2], ze[2]) matched both 2 and 3, and (az[3], ze[3]) didn't match anything.


    HISTORY:
        2013 Jul 15 - created by Brian Harding
    '''
    ze = np.array(ze)
    az = np.array(az)
    az_rad = az * pi/180.0
    valid_ze = np.array(valid_ze)
    valid_az = np.array(valid_az)
    valid_az_rad = valid_az * pi/180.0

    N = len(az)
    M = len(valid_az)

    # sin/cos backwards because azimuth is deg E of N
    # (this doesn't really matter here though).
    c_x = valid_ze * np.sin(valid_az_rad)
    c_y = valid_ze * np.cos(valid_az_rad)
    x   = ze       * np.sin(az_rad)
    y   = ze       * np.cos(az_rad)

    direc = []
    dists = []
    for j in range(N):
        dist = np.sqrt((c_x - x[j])**2 + (c_y - y[j])**2)
        idx_list = (dist == dist.min()).nonzero()[0]
        direc.append(idx_list)
        dists.append(dist.min())

    return direc, dists

def circle_fit(x,y,method='geometric'):
# Fit a circle to the points (x,y) using either the algebraic or geometric method
# INPUTS:
#     x: vector of x-coordinates
#     y: vector of y-coordinates
#     method: Either 'algebraic' or 'geometric'
#          'algebraic': Fast, robust, and stable, but may be biased,
#                 especially with uneven sampling around the circle.
#          'geometric': Truly minimizes sum-of-squared-distance-error.
#                  The algebraic solution will be used as initial guess
#                  for this nonlinear method.
# OUTPUTS:
#     (x0,y0,R) = center location and radius of fitted circle
    if method == 'algebraic':
        A = np.vstack((2.*x, 2.*y, np.ones(np.shape(x)))).T
        b = x**2 + y**2
        p = np.linalg.lstsq(A,b, rcond=-1)[0]
        x0 = p[0]
        y0 = p[1]
        R = np.sqrt(p[2] + x0**2 + y0**2)
        return x0,y0,R

    elif method == 'geometric':
        # use the linear solution as the initial guess
        x0,y0,R = circle_fit(x, y, 'algebraic')

        # define residual
        def residual(params,x,y):
            x0 = params['x0'].value
            y0 = params['y0'].value
            R  = params['R'].value
            resid = R - np.sqrt( (x-x0)**2 + (y-y0)**2 )
            return resid
        # do inversion using lmfit
        params = Parameters()
        params.add('x0',value=x0)
        params.add('y0',value=y0)
        params.add('R',value=R)
        out = minimize(residual, params, args=(x,y))
        x0fit = out.params['x0'].value
        y0fit = out.params['y0'].value
        Rfit = out.params['R'].value
        return x0fit,y0fit,Rfit

    else:
        raise Exception('Invalid Input: %s' % method)

def FindCenter(img,  circle_fit_method = 'geometric', thresh=None, max_r=None):
    # Function to estimate the center of a ring pattern.
    # INPUTS:
    #      img - The laser image to work with (i.e. fringe pattern).  This is a
    #       2D array.
    #      (optional):
    #      circle_fit_method - string, 'algebraic' or 'geometric'.
    #              Algebraic is faster but can be biased.
    #      thresh - scalar. The threshold used for isolating fringes [pixel counts]
    #      max_r - scalar. Indirectly controls how many fringes are used.
    #              Maximum number of pixels from center that are used [pixels]
    # OUTPUTS:
    #     (cx,cy) - The inferred center.
    # ALGORITHM:
    #     Threshold the image and label contiguous regions with greater than 100 pixels.
    #     Fit a circle to the contiguous regions.  These should correspond to different
    #     fringes.  Return the median center value, because of outliers (e.g.,
    #     partial fringes).  Avoid using fringes near the edge of the image.
    # HISTORY:
    #     01 Apr 2013 - written by Brian Harding, using FPI.FindCenterLabel
    #                   as a guide (which was written by Jonathan J. Makela)
    #     10 Jul 2013 - (BJH) Take array as input instead of Image.
    #                   Don't use label=0 since this is the background.
    #                   Optimize to decrease runtime.
    #

    # Median-filter preprocessing is not necessary. Labeling accounts for noise spikes.

    # Parameters
    pthresh = 70 # percentile for thresholding
    pmax_r = 0.8 # how much of the image to use, if max_r is not given
    min_label_size = 100 # only use labeled regions with more pixels than this

    sx,sy = np.shape(img)
    gx = int(floor(sx/2))
    gy = int(floor(sy/2))

    # If we are not given a threshold, estimate one
    if not thresh:
        thresh = np.percentile(np.reshape(img,np.size(img)),pthresh)

    # If we are not given a max_r, estimate one
    if not max_r:
        max_r = (sx - gx)*pmax_r

    # Create a mask and a thresholded image
    mask = img.astype(bool)
    imgth = img.astype(bool)

    # Calculate radial distances
    mask[gx,gy] = False
    f = np.sqrt(mahotas.distance(mask))
    mask = f < max_r

    # Threshold, mask and label the image
    imgth = (img > thresh) & mask
    labels1, numlabels1 = ndimage.label(imgth)
    labels1 = labels1.astype('int32') # added 2016-03-29 per Fasil Tesema's suggestion

    # Get rid of labels that only have a few points
    labels2 = np.zeros(np.shape(labels1)) # new, processed labels image
    sizes = np.bincount(labels1.ravel())
    mask_sizes = sizes > min_label_size
    fringeidx = mask_sizes.nonzero()[0]
    for j,i in enumerate(fringeidx):
        idx = labels1==i
        labels2[idx] = j
    numlabels2 = len(fringeidx)

    # Fit a circle to the points in each labeled region
    centers = np.zeros((numlabels2,2))
    Rs = np.zeros((numlabels2,1))
    for idx in range(numlabels2):
        label = labels2 == idx
        ty,tx = label.nonzero()
        x0,y0,R = circle_fit(tx, ty, circle_fit_method)
        centers[idx,:] = [x0,y0]
        Rs[idx] = R

    # Remove outliers by using median.
    # Ignore label=0 since this is the background.
    (cx,cy) = np.median(centers[1:,:],0)
    return cx,cy

def ReadIMG(fname):
    # Function to read in an FPI IMG file.  Converts the data from I;16L to I so that the
    # data can easily be plotted and manipulated.
    #
    # INPUT:
    #	fname - the filename to be read
    #
    # OUTPUT:
    #	im - the image read in.  This structure also includes the header information
    #
    # HISTORY:
    #	16 May 2012 - written by Jonathan J. Makela (jmakela@illinois.edu) based on MATLAB code

    if fname.split('.')[-1] == 'img':
        # Load the image in
        im = Image.open(fname)

        # Convert the image from I;16L to I
        im = im.convert('I')
    elif fname.split('.')[-1] == 'hdf5':
        # New HDF5 image formate
        temp = h5py.File(fname,'r')

        # Convert to im structure
        im = temp['image']
        im.info = {'ExposureTime': temp['image'].attrs['ExposureTime'],
                   'zeAngle': temp['image'].attrs['zeAngle'],
                   'azAngle': temp['image'].attrs['azAngle'],
                   'XBinning': 2, ### JJM NEED TO CHANGE!
                   'YBinning': 2,
                   'CCDTemperature': temp['image'].attrs['CCDTemperature'],
                   'Pressure': temp['image'].attrs['Pressure (Pa)'],
                   'OutsideTemperature': temp['image'].attrs['OutsideTemperature (C)'],
                   'LocalTime': datetime.datetime.strptime(temp['image'].attrs['LocalTime'],'%Y-%m-%d %H:%M:%S.%f'),}

    else:
        im = None

    return im

def FindEqualAreas(img,cx,cy,N):
        # Function to find the equal area annular regions required for annular summation
        #
        # INPUT:
        #       img - the image to work with.  2D array
        #       cx - the center pixel of the image in the x direction (0 index!)
        #       cy - the center pixel of the image in the y direction (0 index!)
        #       N - the number of annular regions to create
        #
        # OUTPUT:
        #       annuli - a dictionary containing the start ('inner') and stop ('outer') pixel
        #                        of each annular region as well as the indecies ('ind') into a sorted
        #                        array of the pixel radii from the input center pixel
        #
        # HISTORY:
        #       17 May 2012 - written by Jonathan J. Makela (jmakela@illinois.edu) based on MATLAB code
        #       02 Apr 2013 - For portability, take the 2D array as a parameter instead of an Image instance

        # Get the size of the image we are working with
        (sx,sy) = np.shape(img)

        # Raise an exception if the specified center is outside the image
        if (cx < 0) or (cx > sx) or (cy < 0) or (cy > sy):
            raise Exception('Center Location (%.1f, %.1f) is outside the allowable range' % (cx,cy))

        # Calculate the area of the largest circle inscribed within the image centered on the requested
        # center pixel. (cx,cy) is zero indexed.
        area = min((sx-1-cx)**2*pi/N,(cx)**2*pi/N,(sy-1-cy)**2*pi/N,(cy)**2*pi/N)

        # Calculate the radius of each pixel from the requested center pixel
        [y,x] = np.mgrid[0:sy,0:sx]
        r = np.sqrt((x-cx)**2+(y-cy)**2)

        # Sort the radius array and save the indecies of the sorted array
        I = np.argsort(r,axis=None)
        rs = r.flatten()
        Y = rs[I]

        # Define output vectors
        inner = np.zeros((N))
        outer = np.zeros((N))
        ri = np.zeros((N))

        # Loop through the requested number of annular regions.  Find the indecies within
        # the sorted array bounding each region
        rinner = 0
        for i in range(0,N):
                # The outer radius of the current annular region being considered
                #ri[i] = np.sqrt(area/pi+rinner**2)
                router = np.sqrt(area/pi + rinner**2)
                ri[i] = (router+rinner)/2 # radius that matters

                # The indecies of pixels in the current annular region
                temp = np.logical_and(rinner<=Y,Y<router).nonzero()
                ind = temp[0]

                # Save the first
                inner[i] = ind[0]
                outer[i] = ind[-1]

                # The current outer radius becomes the next regions inner radius
                rinner = router

        # create the dictionary
        annuli = {'inner':inner, 'outer': outer, 'ind': I, 'r': ri}

        # return the annuli dictionary
        return annuli

def AnnularSum(img,annuli,bg = None):
    # Function to perform the annular summation of a 2D FPI interferogram.
    #
    # INPUT:
    #	img - the image to work with.  This is a 2D array
    #	annuli - a dictionary containing information on the annular regions.  This is returned by
    #			FPI.FindEqualAreas
    #	bg - a constant background value to remove from the image data.  If not supplied, the mean
    #		value of the image is used
    #
    # OUTPUT:
    #	spectra - the summed and normalized 1D spectra
    #
    # HISTORY:
    #	17 May 2012 - written by Jonathan J. Makela (jmakela@illinois.edu) based on MATLAB code
    #   29 Jun 2012 - added None default parameter for bg

    # Get the size of the image we are working with
    (sx,sy) = np.shape(img)
    data = img.ravel()


    # Check if a background value was provided.
    if bg is None:
        # No, it wasn't.  Use the mean value of the image
        bg = data.mean()

    data = data - bg

    # Loop through each annular region
    N = int(annuli['inner'].size)
    spectra = np.zeros(N)
    sigma = np.zeros(N)
    for i in range(0,N):
        # The indicies to work with
        ind = annuli['ind'][int(annuli['inner'][i]):int(annuli['outer'][i]+1)]

        # The indicies of pixels that are within 3 std of the mean of the annular region
        ind2 = (abs(data[ind] - data[ind].mean()) < (3*data[ind].std())).nonzero()

        # Check if we have any values to remove and remove them.
        while(np.not_equal(np.size(ind), np.size(ind2))):
            # swap the indecies
            ind = ind[ind2]

            # The indecies of pixels that are within 3 std of the mean of the annular region
            ind2 = (abs(data[ind] - data[ind].mean()) < (3*data[ind].std())).nonzero()

        # swap the indecies
        ind = ind[ind2]

        # Perform the summation
        spectra[i] = data[ind].mean()

        # Calculate the uncertanty in the estimate which is the std divided by the
        # sqrt of the number of pixels used
        sigma[i] = data[ind].std()/np.sqrt(len(ind))

    # return the summed spectra
    return spectra, sigma


def weighted_avg_and_std(values, weights):
    """
    Returns the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.

    from http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    """
    ind = (np.isfinite(values).nonzero() and np.isfinite(weights).nonzero())[0]

    if len(ind) > 0 and sum(weights[ind]) != 0.:
        average = np.average(values[ind], weights=weights[ind])
        variance = np.dot(weights[ind], (values[ind]-average)**2)/weights[ind].sum()  # Fast and numerically precise
    else:
        average = np.nan
        variance = np.nan

    return (average, sqrt(variance))

def all_indices(value, qlist):
# Function to return all of the indicies in a list that match a requested value
#
# INPUTS:
#       value - value to be matched
#       qlist - list to search through
#
# OUTPUTS:
#       indices - indecies in qlist that match value
#
# HISTORY:
#       10 Aug 2012 - written by Jonathan J. Makela (jmakela@illinois.edu) based on
#       http://stackoverflow.com/questions/176918/in-python-how-do-i-find-the-index-of-an-item-given-a-list-containing-it

    indices = []
    idx = -1
    while True:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices

def Laser_FringeModel(params, r):
# Function to simulate a 1D fringe pattern
#
# INPUTS:
#       r - indecies running from 1:n_elements
#       params - a dictionary contaiing the fit parameters.
#
# OUTPUTS:
#       fringe - the 1-D simulated fringe
#
# HISTORY:
#       25 Oct 2012 - written by Jonathan J. Makela based on FPI_FringePattern.m
#                     which was written by Brian Harding.
#       06 Aug 2013 - added blur kernel estimation parameters b0,b1,b2 and removed
#                     Ra1 and Ra2.

    # unpack the parameter dictionary
    n = params['n'].value
    t = params['t'].value
    lam = params['lam'].value
    R = params['R'].value
    alpha = params['alpha'].value
    I = params['I'].value
    B = params['B'].value
    a1 = params['a1'].value
    a2 = params['a2'].value
    b0 = params['b0'].value
    b1 = params['b1'].value
    b2 = params['b2'].value


    x = np.pi * r/r.max()
    FF = (1 + a1*x + a2*x**2) # intensity falloff
    refl = R

    # Finesse
    F = 4*refl/(1-refl)**2

    # m
    th = np.arctan(alpha * r)
    m = 2*n*t/lam*np.cos(th)

    # full expression
    fringe = FF*I/(1+F*np.sin(np.pi*m)**2)+B

    # apply gaussian blur
    b = b0 + b1*np.cos(x) + b2*np.sin(x)
    if b0:
        fringeblur = np.zeros(len(fringe))
        for i in range(len(fringeblur)):
            z = (r - r[i])/b[i]
            w = np.exp(-z**2)
            w = w / w.sum() # normalize
            fringeblur[i] = np.dot(w, fringe)
    else:
        fringeblur = fringe

    return fringeblur

def Laser_Residual(pars, x, sigma=None, data=None):
# Simple return function for the use of the lmfit Minimizer
#
# INPUT:
#      pars - parameters to be fit
#      x - x values
#
# OUTPUT:
#      returns the residuals between the model and the data
#
# HISTORY:
#      25 Oct 2012 - written by Jonathan J. Makela based on https://github.com/newville/lmfit-py/blob/master/tests/test_covar.py
#				  and http://newville.github.com/lmfit-py/parameters.html#using-parameters-instead-of-variables

    # Run the model
    model = Laser_FringeModel(pars, x)

    # If no data was passed, just return the model (e.g., use to simulate a fringe)
    if data is None:
        return model
    if sigma is None:
        return (model-data)

    return (model - data) / sigma

def Sky_FringeModel(params, r, lamvec, A):
# Function to simulate a 1D fringe pattern
#
# INPUTS:
#       r - indicies running from 1:n_elements
#       params - a dictionary contaiing the fit parameters (see below for format)
#
# OUTPUTS:
#       fringe - the 1-D simulated fringe
#
# HISTORY:
#       25 Oct 2012 - written by Jonathan J. Makela based on FPI_FringePattern.m
#                     which was written by Brian Harding.
#       20 May 2013 - modified by Brian Harding to take v,T as parameters and to
#                     scale skyI and skyB to have units of (approximately) detector counts
#       03 Jun 2013 - reverted slightly by Brian Harding to take lamc as parameter instead
#                     of v, because terms like (1+v/c) might be within machine precision
#                     of (1+ (v+delta_v)/c), yielding problems with estimating
#                     the Jacobian during the inversion.
#       09 Jul 2013 - Brian Harding added skym parameter in order to simulate background
#                     continuum slope in the sky spectrum and removed skya1 and skya2
#       04 Dec 2013 - Brian Harding correctly normalized the sky spectrum Gaussian
#                     so that skyI does not contain a contribution from temperature.

    # Constants
    c = 299792458.
    k = 1.3806503e-23
    m = 16/6.0221367e26

    # unpack the parameter dictionary
    T = params['T'].value # Doppler temperature
    skyI = params['skyI'].value # line intensity
    skyB = params['skyB'].value # continuum background
    ccdB = params['ccdB'].value # CCD bias
    lam0 = params['lam0'].value # nominal line center
    skym = params['skym'].value # slope of continuum background, relative to line intensity (i.e. -1 and 1 are extreme values)
    lam_c = params['lamc'].value # Doppler-shifted center wavelength

    # TODO: something about this
    if T < 0.0:
        T = 0.1
        print('Zero T!')

    # Adjust scaling of I and B
    intA = np.mean(np.sum(A,1))
    skyI = skyI / intA
    skyB = skyB / intA

    # Parameterized sky spectrum:
    broad = lam0/c * np.sqrt(k*T/(2*m)) # width of spectrum
    broad0 = lam0/c * np.sqrt(k*1000/(2*m)) # width at 1000 K, for scaling
    skyspectrum = skyB + skyI * (broad0/broad) * np.exp(-((lam_c-lamvec)/(2*broad))**2)
    # Add background slope
    bg = skym * skyI * (lamvec - lamvec.min())/(lamvec.max() - lamvec.min())
    skyspectrum = skyspectrum + bg

    # Calculate fringe pattern
    fringe = A.dot(skyspectrum)

    # Apply ccd bias
    fringe = fringe + ccdB

    return fringe



def Sky_Residual(pars, x, lamvec, A, sigma=None, data=None):
# Simple return function for the use of the lmfit Minimizer
#
# INPUT:
#      pars - parameters to be fit
#      x - x values
#
# OUTPUT:
#      returns the residuals between the model and the data
#
# HISTORY:
#      25 Oct 2012 - written by Jonathan J. Makela based on https://github.com/newville/lmfit-py/blob/master/tests/test_covar.py
#				  and http://newville.github.com/lmfit-py/parameters.html#using-parameters-instead-of-variables

    # Run the model
    model = Sky_FringeModel(pars, x, lamvec, A)

    # If no data was passed, just return the model (e.g., use to simulate a fringe)
    if data is None:
        return model
    if sigma is None:
        return (model-data)

    return (model - data) / sigma

def get_conv_matrix_1D(instr_params, r, L, lam0, nFSRs = 4.):
    ''' (A, lamvec) = get_conv_matrix_1D(instr_params, M, L, lam0)
    Precomputes the matrix that implements the convolution between the instrument
    function and the sky spectrum
    INPUTS:
          instr_params - a Parameters object containing the instrument parameters
                         (see Laser_FringeModel for format)
          r - the radii of the annular bins
          L - The number of wavelength bins to use (sets the spectral resolution
              of the forward model)
          lam0 - The emission wavelength in meters (e.g., 630.0e-9). The modeled spectrum
                 will be a couple FSRs around lam0
          nFSRs - The size of the domain for the spectrum in the forward model, measured in
                  units of FSR. (default 4.)

    OUTPUTS:
          A - the shape-(M,L) matrix that implements the convolution as:
                     sky_fringes = A.dot(sky_spectrum)
          lamvec - the wavelength coordinates over which "sky_spectrum" needs to be defined


    HISTORY:
          17 May 2013 - Written by Brian Harding '''

    B = instr_params['B'].value

    # Set up the matrix that will implement the convolution ahead of time.
    FSR = lam0**2/(2*instr_params['t'].value) # free spectral range
    # Q: does the source of "t" matter? A: probably not, as long as we use a wide-enough range

    # Find a reasonable range to add wavelength contributions over that includes the whole emission
    lamvec = np.linspace(lam0-nFSRs/2*FSR, lam0+nFSRs/2*FSR, L+1)
    lamvec = lamvec[0:-1] # the last wavelength in this list is indistinguishable from the first
    dellam = lamvec[1]-lamvec[0]
    A = np.zeros((len(r),L))
    for i in range(len(lamvec)):
        # Create a "laser image" at this wavelength and add this to the cumulative sky image
        ptemp = Parameters()
        for p in list(instr_params.values()):
            ptemp.add(p.name,value=p.value,vary=p.vary)
        ptemp['lam'].value = lamvec[i]
        single_wavelength_image = Laser_FringeModel(ptemp,r) - B
        # put image in one column of A
        A[:,i] = single_wavelength_image*dellam

    return A, lamvec


def ParameterFit(instrument, site, laser_fns, sky_fns, sky_line_tag='X',direc_tol = 10.0, N=500, N0=0, N1=500, 
                 logfile=None, diagnostic_fig=None, horizon_cutoff = -6, reference='laser',
                 warning_log=None):
    '''
     Function that solves for temperatures and Doppler velocities for all data
     in the present directory.  Uses a parameter fit algorithm to first estimate
     the instrument parameters of the FPI
     and then uses these instrument parameters to fit the sky images. Only uses
     data that were taken when the sun was down (defined by horizon_cutoff).

     INPUT:
           instrument - an instrument dictionary containing necessary information about the instrument:
                        ('focal_length', 'pix_size', 'lam_laser', 'lam0', 'nominal_t') all in meters
           site - a site dictionary containing information about the site:
                 ('Name', 'Abbreviation', 'Timezone', 'Location', 'Directions')
           laser_fns - a list of strings containing the paths to the laser images to be analyzed
           sky_fns   - a list of strings containing the paths to the sky images to be analyzed

     OPTIONAL INPUT:
           direc_tol - [degrees] The tolerance used to match expected look directions to actual directions.
                        (TODO: Is there a better place for this parameter to live?)
           N, N0, N1 - number of annular regions to use as well as the innermost and outermost
                       rings to consider in the parameter fit
           logfile - the txt file to write logs to.  If None, it will be created.
           diagnostic_fig - the figure that diagnostic plots will be drawn to. If None, no plots
                            will be produced.
           horizon_cutoff - int - the angle above which the sun is considered "up". (e.g., -6)
           reference - str, 'laser' or 'zenith'. Whether or not to use the laser images. For example,
                       if the laser images are not usable, then specify 'Zenith', and default instrument
                       parameters (specified in the instruments dictionary) will be used in the
                       analysis.
            warning_log - WarningLog to log the text associated with any "notify_the_human" warnings.


     OUTPUT:
           FPI_Results - a dictionary with many parameters pertaining to the inversion (see below)
           notify_the_humans - a boolean, True if the inversion has a problem that is
                               best solved by a human, such as if the laser shutter doesn't
                               appear to be open. This is a softer warning than raising an
                               Exception, and can be ignored in most cases.

     HISTORY:
           14 Nov 2012 - written by Jonathan J. Makela based on MATLAB code written by
                         Brian Harding.
           15 Nov 2012 - modifications to use instrument_params and site variables
           12 Dec 2012 - allowed center pixel to change with time
           15 Apr 2013 - two-pass smoothing of gap parameter (jjm)
           10 Jul 2013 - (bjh) clean up, change calling convention. Made inversion more robust
                         by rolling in changes from Monte Carlo simulation. Removed second-pass
                         of laser fit. Added diagnostic plots. Changed spline fit to laser params.
           06 Aug 2013 - (bjh) Replace "reflectivity falloff" with "blur kernel estimation"

    '''


    # Tunable parameters
    CENTER_VARIATION_MAX = 1.5 # pixels, variation throughout night. Above this, outliers >1 stdev will be ignored.
    # Set threshold below which, the solved-for lamc
    # will be used to center the grid search for the next iteration.
    # This has a small effect; it is used to discourage FSR-jumps in wind.
    SIGMA_V_THRESH = 30.
    # How much of the fringe to use around the peak.
    # Between 0.0 and 1.0
    FRINGEFACTOR = 1.0
    MOON_THRESH = 37. # deg. Don't use samples that are closer than this angle to the moon.
    # What temporal interpolation to use for instrument parameters
    LASERPARAMINTERP = 'linear' # 'linear' or 'spline'
    # Whether to estimate blur parameters. This was found necessary
    # because the Morocco FPI appears to have something wrong with it, and
    # there are fewer artifacts in the estimated temperature if this is set
    # to False.
    ESTIMATE_BLUR = True
    if instrument['name'] in ['minime03','minime11',]:
        ESTIMATE_BLUR = False
    # The threshold for CCD temperature, above which a human should be
    # notified, because there is probably something wrong with the CCD.
    MAX_CCD_TEMP = -55.0 # C
    LAS_I_THRESH = 13. # counts. If the laser intensity is less than this, ignore it.


    ############## TEMPORARY ################
    if 'Abbreviation' in list(site.keys()):
        if   site['Abbreviation'] == 'a3o':
            FRINGEFACTOR = 0.6
        elif site['Abbreviation'] == 'mrh':
            FRINGEFACTOR = 0.82
        elif site['Abbreviation'] == 'nzk':
            FRINGEFACTOR = 0.9
    #########################################

    # Parse inputs
    uselaser = False
    if reference == 'laser':
        uselaser = True
    elif reference == 'zenith':
        uselaser = False
    else:
        raise Exception('Invalid option ("%s") for "reference" parameter. Try "zenith" or "laser"' % reference)

    # Copy the filenames and sort
    lasers_all = laser_fns[:]
    sky_all =    sky_fns[:]
    lasers_all.sort()
    sky_all.sort()


    if len(sky_all)==0:
        raise Exception('No sky images specified.')

    notify_the_humans = False

    # Check if a logfile was sent in
    created_logfile = False
    if logfile is None:
        # No file was given, generate a default logfile in the current directory
        logfile = open('analyze.log','a')
        created_logfile = True

    if warning_log is None:
        # No warning log provided, open a new one
        warning_log = WarningLog()

    # Only use images during sundown conditions.
    local = pytz.timezone(site['Timezone']) # pytz timezone of the site
    sun = ephem.Sun()
    obs = ephem.Observer()
    obs.lat =       str(site['Location'][0])
    obs.lon =       str(site['Location'][1])
    obs.elevation = site['Location'][2]
    obs.pressure =  0.
    obs.horizon =  '%i:0' % horizon_cutoff # for astro twilight
    # Read the first file, subtract 12 hours, then ask for next sunset time
    # Read the first file, and ask for next sunrise time.
    d = ReadIMG(sky_all[0])
    dt0 = local.localize(d.info['LocalTime']).astimezone(pytz.utc)
    refdt = dt0 - datetime.timedelta(hours=12)
    obs.date = refdt
    # Calculate sunset time on this day.
    # Make datetime instance "aware" by specifying utc.
    sunset = obs.next_setting(sun).datetime().replace(tzinfo = pytz.utc)
    obs.date = dt0
    # Same for sunrise.
    sunrise = obs.next_rising(sun).datetime().replace(tzinfo = pytz.utc)

    lasers = []
    for fn in lasers_all:
        d = ReadIMG(fn)
        dt = local.localize(d.info['LocalTime'])
        if dt > sunset and dt < sunrise:
            lasers.append(fn)
        else:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Sun up during %s. Ignoring this image...\n' % fn )


    # Initialize data holder
    laser_times_center = []
    laser_intT = []
    laser_temperature = []
    dt_laser = []
    center = None
    fnames_to_remove = []
    laser_fringes = []
    laser_annuli = []
    lt0 = local.localize(datetime.datetime(1989,4,18)) # arbitrary default
    if uselaser: # analyze the laser images for the center
        for fname in lasers:
            # Read in the laser image and estimate the center location
            d = ReadIMG(fname)
            img = np.asarray(d)
            (cx,cy) = FindCenter(img)

            appx_I = np.percentile(img,95) - np.percentile(img,5)
            if appx_I < LAS_I_THRESH:
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                                '%s laser fringes too dim (%.1f < %.1f). Ignoring this image...\n' % (fname, appx_I, LAS_I_THRESH))
                fnames_to_remove.append(fname)

            # Make sure valid data returned
            elif cx is not None and np.isfinite(cx):
                if center is None:
                    center = [cx, cy]
                else:
                    center = np.vstack((center, [cx, cy]))

                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + '%s center: %.2f, %.2f\n' % (fname, cx, cy))

                # Append the time of the laser image
                laser_times_center.append(local.localize(d.info['LocalTime']))
                laser_intT.append(d.info['ExposureTime'])
                laser_temperature.append(d.info['CCDTemperature'])

                dt_laser.append((laser_times_center[-1]-laser_times_center[0]).seconds)

            else:
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                                '%s center is None or nan. Ignoring this image...\n' % (fname))
                fnames_to_remove.append(fname)

        for fname in fnames_to_remove:
            lasers.remove(fname)

        dt_laser = np.array(dt_laser)

        # Try to do some initial quality control by judging the centerfinding
        # If there are no usable center locations, raise an exception. We can't analyze
        # the data without any usable laser images.
        if center is None or len(center)==0:
            mid = np.ceil(np.shape(d)[0]/2.0)
            cx = lambda t: mid
            cy = lambda t: mid
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                    'Centerfinding failed for all laser images. <BADLASER>\n')
#            raise Exception('Centerfinding failed for all laser images. This day was tagged with <BADLASER> flag.')
            raise BadLaserError('Centerfinding failed for all laser images. This day was tagged with <BADLASER> flag.')

        # If there is only one usable point, throw a warning but use this point.
        elif np.shape(center)==(2,):
            cx = lambda t: center[0]
            cy = lambda t: center[1]
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                    'Centerfinding failed for all but 1 laser image. <BADLASER>\n')
#            raise Exception('Centerfinding failed for all but 1 laser image. This day was tagged with <BADLASER> flag.')
            raise BadLaserError('Centerfinding failed for all but 1 laser image. This day was tagged with <BADLASER> flag.')

        else: # There are > 1 center points. Use them.
            # If we have too much variation in the center locations, remove outliers
            if (np.std(center,0) > CENTER_VARIATION_MAX).any(): # if variation is more than a few pixels, there's probably a bad point.
                cenx = center[:,0]
                ceny = center[:,1]
                good = (abs(cenx - np.mean(cenx)) < 2*np.std(cenx)) & (abs(ceny - np.mean(ceny)) < 2*np.std(ceny))
                center = center[good,:]
                dt_laser = dt_laser[good]
                # Should we ignore the laser images in the subsequent analysis? For now,
                # ignore them since these are probably invalid data, and will cause problems
                # with laser parameter interpolation.
                lasers = [fn for (fn, valid) in zip(lasers, good) if valid]
                laser_times_center = [t for (t, valid) in zip(laser_times_center, good) if valid]
                nignored = len(good) - sum(good)
                if nignored > 0:
                    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                            'WARNING: %03i laser images ignored because the center is an outlier.\n' % nignored)
                if nignored > 1: # Only warn humans if there's more than one problem image
                    notify_the_humans = True
                    warning_log.add(message='%03i laser images ignored because the center is an outlier.' % nignored,
                                    title='Bad laser center',
                                    label='warning:laser_center')

            # Find large jumps and remove the files that contributed to them
            cenx = center[:,0]
            ceny = center[:,1]
            center_diff = np.sqrt((cenx[1:]-cenx[:-1])**2 + (ceny[1:]-ceny[:-1])**2)
            jumpy = center_diff > CENTER_VARIATION_MAX
            if any(jumpy):
                # Remove all points before and after the jump
                good = ~ (np.concatenate(([False],jumpy)) | np.concatenate((jumpy,[False])))
                center = center[good,:]
                dt_laser = dt_laser[good]
                lasers = [fn for (fn, valid) in zip(lasers, good) if valid]
                laser_times_center = [t for (t, valid) in zip(laser_times_center, good) if valid]
                nignored = len(good) - sum(good)
                if nignored > 0:
                    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                            'WARNING: %03i laser images ignored because the center jumped by too much.\n' % nignored)
                if nignored > 2: # Only warn humans if there's more than one problem jump
                    notify_the_humans = True
                    warning_log.add(message='%03i laser images ignored because the center jumped by too much.' % nignored,
                                    title='Bad laser center',
                                    label='warning:laser_center')

            # If there are enough points after this, fit a poly to the center positions
            if len(dt_laser) > 0:
                lt0 = laser_times_center[0]
                npoly = np.floor(len(dt_laser)/10) # use an adaptive degree polynomial to fit
                # Limit the maximum degree, because of numerical sensitivity for large orders.
                if npoly > 10:
                    npoly = 10
                pf_cx = np.polyfit(dt_laser,center[:,0],npoly)
                cx = np.poly1d(pf_cx)
                pf_cy = np.polyfit(dt_laser,center[:,1],npoly)
                cy = np.poly1d(pf_cy)

                #uncomment for linear interpolation
                #from scipy import interpolate as sinterpolate
                #cx=sinterpolate.interp1d(dt_laser,center[:,0],kind='linear',fill_value='extrapolate')
                #cy=sinterpolate.interp1d(dt_laser,center[:,1],kind='linear',fill_value='extrapolate')

                if len(dt_laser) < 6:
                    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                                'WARNING: Very few (%i) laser images for center trending. Consider using Zenith reference. <BADLASER>\n' % len(dt_laser))
                    notify_the_humans = True
                    warning_log.add(message='Very few (%i) laser images for center trending. Consider using Zenith reference. <BADLASER>' % len(dt_laser),
                                    title='Too few laser centers',
                                    label='warning:laser_center')

            else: # I don't think we can ever get to this point, but just in case:
                mid = np.ceil(np.shape(d)[0]/2.0)
                cx = lambda t: mid
                cy = lambda t: mid
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                        'WARNING: No usable center points for polyfit. Using midpoint of CCD, ' + \
                        'but suggested solution is to use Zenith reference. <BADLASER>\n')
                notify_the_humans = True
                warning_log.add(message='No usable center points for polyfit. Using midpoint of CCD, but suggested solution is to use Zenith reference. <BADLASER>',
                                title='Not enough laser centers to fit',
                                label='warning:laser_center')

    else: # uselaser is False. Use the center locations from the instrument dictionary
        # Create a dummy function of time that always returns the same value.
        cx = lambda t: instrument['default_params']['center'][0]
        cy = lambda t: instrument['default_params']['center'][1]
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                                'Using zenith reference. Using hard-coded value of (%.1f, %.1f) as center pixel. \n' % (cx(0), cy(0)))


    laser_spfits = {} # laser_spfits['t'] contains the interpolation function for t
    if uselaser: # Analyze the laser images to obtain the spline fits to the parameters
        # Initialize values
        output = []
        last_chi = -999
        laser_times = []
        laser_fns = []
        laser_redchi = []
        last_t = instrument['nominal_t']
        lam_laser = instrument['lam_laser']

        # Loop through all of the lasers
        for fname in lasers:


            # Read in the laser image and estimate the center location
            d = ReadIMG(fname)
            img = np.asarray(d)

            # Calculate the annuli to use for this time
            dt = (local.localize(d.info['LocalTime'])-lt0).seconds
            annuli = FindEqualAreas(img,cx(dt),cy(dt),N)

            # Perform annular summation
            laser_spectra, laser_sigma = AnnularSum(img,annuli,0)

            if not np.isnan(laser_spectra).any(): # Do the inversion. Otherwise, ignore it.

                ####### Find good initial guesses for the parameters ######

                # Magnification parameter by using geometry (assuming square binning)
                alpha = instrument['pix_size']/instrument['focal_length'] * d.info['XBinning']
                # Intensity and background by looking at fringes
                I = laser_spectra.max() - laser_spectra.min()
                B = laser_spectra.min()

                laser_params = Parameters()
                laser_params.add('n',     value = 1.0,       vary = False)
                laser_params.add('t',     value = None,      vary = False) # We will search for this
                laser_params.add('lam',   value = lam_laser, vary = False)
                laser_params.add('R',     value = 0.5,        vary = False)
                laser_params.add('alpha', value = alpha,     vary = False)
                laser_params.add('I',     value = I,         vary = False)
                laser_params.add('B',     value = B,         vary = False)
                laser_params.add('a1',    value = -0.1,      vary = False)
                laser_params.add('a2',    value = 0.0,       vary = False)
                laser_params.add('b0',    value = 0.5,       vary = False)
                laser_params.add('b1',    value = 0.0,       vary = False)
                laser_params.add('b2',    value = 0.0,       vary = False)

                # To find a good initial guess for "t", we'll need to do a grid search.  Search
                # over 1 FSR around the last solved-for value (or the nominal value, if this is first trial).
                # TODO: make this grid search a separate general function.
                def goodness(t):
                    # return the correlation between the fringe with this t and the data
                    laser_params['t'].value = t
                    fringe = Laser_FringeModel(laser_params, annuli['r'])
                    return np.dot(fringe, laser_spectra)/(np.linalg.norm(fringe)*np.linalg.norm(laser_spectra))
                t_candidates = np.linspace(last_t - lam_laser/4, last_t + lam_laser/4, 50)
                corrvec = np.array([goodness(t) for t in t_candidates])
                best_t = t_candidates[corrvec.argmax()]
                laser_params['t'].value = best_t

                ####### Inversion of laser image ##########

                # Now do least-squares fit, but in stages, varying only certain parameters at a time, according to:
                order = [
                         ['alpha'],\
                         ['t','alpha'],\
                         ['B','I'],\
                         ['R'],\
                         ['t','alpha','B','I','R','a1','a2'], \
                         ]
                if ESTIMATE_BLUR:
                    order.append(['t','alpha','B','I','R','a1','a2','b0','b1','b2'])

                data = laser_spectra[N0:N1]
                sigma = laser_sigma[N0:N1]

                for group in order: # set all the params in this group to "vary=True", and then run inversion
                    for param in list(laser_params.keys()):
                        if param in group:
                            laser_params[param].vary = True
                        else:
                            laser_params[param].vary = False
                    # Set falloff terms to false, if this instrument only has a couple fringes
                    if not instrument['many_fringes']:
                        for param in ['a2', 'b1', 'b2']:
                            laser_params[param].vary=False
                            laser_params[param].value=0.0

                    laser_fit = Minimizer(Laser_Residual,laser_params, \
                                        fcn_args=(annuli['r'][N0:N1],), fcn_kws={'data': data, 'sigma': sigma}, \
                                        scale_covar = True)
                    result = laser_fit.leastsq()
                    laser_params=result.params

                #if (laser_fit.redchi/laser_fit.params['I'].value < 1.): # TODO: try to find a better way to do this
                if True: # TODO: just keep this? Let's just do quality control at a later stage.
                    # This is a valid image, store time
                    laser_times.append(local.localize(d.info['LocalTime']))
                    laser_fns.append(fname)
                    # Save the results
                    output.append(result)

                    last_chi = result.redchi
                    last_t = best_t
                    laser_redchi.append(last_chi)
                    laser_fringes.append(data)
                    laser_annuli.append(annuli['r'][N0:N1])

                    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + '%s reduced chisqr: %.4f\n' % (fname, last_chi))

            else:
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + '%s is invalid and not analyzed.\n' % fname)

        #Remove outliers in laser intensity

        # Check if any laser images were actually analyzed.
        if len(output)==0:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No usable laser images found. <BADLASER>\n')
            if created_logfile:
                logfile.close()
 #           raise Exception('No usable laser images found.')
            raise BadLaserError('No usable laser images found. <BADLASER>')

        # Convert laser_params to array
        n = len(output)
        laser_value = {}
        laser_stderr = {}
        for p in output[-1].params:
            laser_value[p] = np.zeros(n)
            laser_stderr[p] = np.zeros(n)
            for i in range(0,n):
                laser_value[p][i] = output[i].params[p].value
                laser_stderr[p][i] = output[i].params[p].stderr


        ######### Spline fit to the parameters ###########

        # Grab the seconds from the first laser observation (required for spline interpolation)
        dt = []
        for x in laser_times:
            diff = (x - lt0)
            dt.append(diff.seconds+diff.days*86400.)
        dt = np.array(dt)

        # Spline fit to the parameters.
        # Take care for times when there are so few measurements the spline fit fails.
        # (This will work for both the fitted and constant params)
        for param in laser_fit.params:
            # Grab the estimates and the associated (normalized) uncertanty
            p = []
            ep = []
            for o in output:
                p.append(o.params[param].value)
                ep.append(o.params[param].stderr)

            p = np.array(p)
            ep = np.array(ep)
            m = len(ep)
            s = m  # Expected value of chi^2

            if param == 't': # force 't' to fit more smoothly
                # i.e., we don't trust those small error bars, so make them bigger
                ep = ep * 3 # TODO: optimal value?  1 seems too small, 10 seems too big

            # set weights
            if laser_fit.params[param].vary:
                # use the error bars from the inversion, but
                # replace the values that are infinite,
                # (i.e., where the errorbars couldn't be calculated)
                # with very large numbers.
                w = [1/stderr if stderr > 0 else 1/(1e6*paramval) for paramval,stderr in zip(p,ep) ]
            else:
                # use arbitrary error bars (doesn't matter - it will fit exactly)
                w = np.ones(np.shape(ep))

            # Try spline, then linear, then zeroth-order interpolation

            try: # spline interpolation
                if LASERPARAMINTERP == 'spline':
                    sfit = interpolate.UnivariateSpline(np.array(dt),p, w=w, s=s)
                elif LASERPARAMINTERP == 'linear':
                    sfit = interpolate.interp1d(np.array(dt), p)
                else:
                    raise Exception('Unrecognized LASERPARAMINTERP: %s' % LASERPARAMINTERP)
            except: # Try linear.
                notify_the_humans = True
                warning_log.add(message='Spline interpolation failed for laser param "%s". Defaulting to linear interpolation.' % (param),
                                title='Parameter fitting failed',
                                label='warning:parameter_fitting')
                try: # if there's only one laser image, then even linear won't work.
                    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                        'Spline interpolation failed for laser param "%s". Defaulting to linear interpolation.\n' % (param))
                    sfit = interpolate.interp1d(np.array(dt), p)
                except: # all else failed. There's probably only one laser image.
                    logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                        'Spline and linear interpolation both failed for laser param "%s". Defaulting to zeroth-order interpolation.\n' % (param))
                    sfit = lambda t,m=p[0] : m # whatever the input, just return the parameter.


            laser_spfits[param] = sfit

    else: # uselaser is False.
        # Create dummy interpolation function for laser parameters using default values.

        # Use default values from the instrument dictionary.
        laser_spfits['n']       = lambda t: 1.0
        laser_spfits['t']       = lambda t: instrument['nominal_t']
        laser_spfits['lam']     = lambda t: instrument['lam_laser']
        laser_spfits['R']       = lambda t: instrument['default_params']['R']
        laser_spfits['alpha']   = lambda t: instrument['default_params']['alpha']
        laser_spfits['I']       = lambda t: instrument['default_params']['I']
        laser_spfits['B']       = lambda t: instrument['default_params']['B']
        laser_spfits['a1']      = lambda t: instrument['default_params']['a1']
        laser_spfits['a2']      = lambda t: instrument['default_params']['a2']
        laser_spfits['b0']      = lambda t: instrument['default_params']['b0']
        laser_spfits['b1']      = lambda t: instrument['default_params']['b1']
        laser_spfits['b2']      = lambda t: instrument['default_params']['b2']

        # Create placeholders for the output.
        laser_times = []
        laser_temperature = []
        laser_redchi = []
        laser_value = {}
        laser_stderr = {}


    ############# Invert the sky images ###############
    # Define function to check moon angle
    def get_moon_angle(t, lat, lon, az, ze):
        '''
        Calculate the great-circle angle between the direciton of the moon
        and the direction given by (az,ze), for the datetime t and location (lat,lon)
        '''
        obs = ephem.Observer()
        obs.lat = str(lat)
        obs.lon = str(lon)
        obs.date = t.astimezone(pytz.utc)

        moon = ephem.Moon(obs)
        moonAz = moon.az.real
        moonZe = pi/2 - moon.alt.real

        # Code copied from original Master scripts for FPI observations
        a = np.cos(az*pi/180)*np.sin(ze*pi/180)
        b = np.sin(az*pi/180)*np.sin(ze*pi/180)
        aMoon = np.cos(moonAz)*np.sin(moonZe)
        bMoon = np.sin(moonAz)*np.sin(moonZe)
        moonAngle = np.arccos(a*aMoon + b*bMoon + np.cos(ze*pi/180) * np.cos(moonZe))

        return moonAngle*180./pi


    # Define function to characterize the wind error arising from etalon
    # gap fluctuations
    def calc_wind_calib_error(dt, gaps, my_dt):
        '''
        Calculate the error induced by the inaccurate (not imprecise!) knowledge
        of the etalon gap, at the time indicated.
        INPUTS:
            dt: array of floats, number of seconds since first laser exposure,
                for each laser exposure. The first element should be 0.
            gaps: array of floats, the estimated gap at each laser time
            my_dt: number of seconds since first laser exposure: the time you want
                   to calculate the calibration error at.
        OUTPUTS:
            wind_cal_error: m/s. This is nan if it's not straddled by laser
                            exposures.
        HISTORY:
            2015 Feb 02: Written by Brian J. Harding
            2015 Jun 30: Rewritten by Brian J. Harding to use gap temporal curvature
                         instead of gradient as proxy for calibration error.
        '''

        c = 3e8
        #dt = np.array([(laser_time - laser_times[0]).total_seconds() for laser_time in laser_times])

        if my_dt < dt[0] or my_dt > dt[-1]: # not bounded by laser images
            return np.nan
        elif my_dt < dt[1]: # extrapolate from second laser image
            my_dt = dt[1]
        elif my_dt > dt[-2]: # extrapolate from second-to-last laser image
            my_dt = dt[-2]

        # Calculate the cal_error at each laser image and linearly interpolate to my_dt
        cal_error_vec = np.zeros(len(dt))
        for k in range(1,len(dt)-1):
            # Extract gap before, during, and after this image
            gt0 = dt[k-1]
            gt1 = dt[k]
            gt2 = dt[k+1]
            g0 = gaps[k-1]
            g1 = gaps[k]
            g2 = gaps[k+1]
            # Estimate uncertainty due to curvature
            g1interp = (gt1-gt0) * (g2-g0)/(gt2-gt0) + g0
            g1diff = g1 - g1interp
            # Convert to an uncertainty in wind
            cal_error = abs(c/g1*g1diff)
            cal_error_vec[k] = cal_error

        try:
            # linearly interpolate to my_dt
            sfit = interpolate.interp1d(dt[1:-1], cal_error_vec[1:-1])
            wind_cal_error = sfit(my_dt)
        except:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                    'WARNING: Interpolation failed in calc_wind_calib_error. Returning NaN.\n')

            wind_cal_error = np.nan 
        
        return wind_cal_error



    # Grab all sky images. Don't use sky images when sun is up,
    # or before (after) the first (last) laser image (if reference=='laser'),
    # or when the moon is up
    sky = []

    for fn in sky_all:

        d = ReadIMG(fn)
        dt = local.localize(d.info['LocalTime'])
        sundown = dt > sunset and dt < sunrise
        moonangle = get_moon_angle(dt, site['Location'][0], site['Location'][1], d.info['azAngle'], d.info['zeAngle'])
        if not sundown: # Write to log and don't use this image
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'Sun up during %s. Ignoring this image...\n' % fn )
        elif uselaser and not (dt > laser_times[0] and dt < laser_times[-1]):
            # Write to log and don't use this image
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + '%s not bounded by laser images. Ignoring this image...\n' % fn )
        elif moonangle < MOON_THRESH:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + '%s is looking too close to the moon (%.1f deg). Ignoring this image...\n' % (fn,moonangle))
        else:
            sky.append(fn)


    if len(sky) == 0:
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No usable sky images found.\n')
        if created_logfile:
            logfile.close()
        raise Exception('No usable sky images found.')

    # Constants
    c = 299792458.
    k = 1.3806503e-23
    m = 16/6.0221367e26

    # Initialize data holder
    sky_times = []
    az = []
    ze = []
    LOSwind = []
    sigma_LOSwind = []
    sigma_fit_LOSwind = []
    sigma_cal_LOSwind = []
    T = []
    sigma_T = []
    direction = []
    sky_redchi = []
    skyI = []
    skyB = []
    ccdB = []
    sigma_skyI = []
    sigma_skyB = []
    sigma_ccdB = []
    sky_out = []
    sky_intT = []
    sky_temperature = []
    sky_fns = []
    sky_fringes = []
    sky_annuli = []


    # Categorize the look directions in the directory
    valid_az = np.array([site['Directions'][direc]['az'] for direc in site['Directions']])
    valid_ze = np.array([site['Directions'][direc]['ze'] for direc in site['Directions']])
    direc_names = list(site['Directions'].keys())
    all_az = []
    all_ze = []
    for fname in sky:
        d = ReadIMG(fname)
        all_az.append(d.info['azAngle'])
        all_ze.append(d.info['zeAngle'])
    all_az = np.array(all_az)
    all_ze = np.array(all_ze)
    # Apply a correction to the reported angles.
    # Unless there was a problem, this function does nothing.
    dn = ReadIMG(sky_all[0]).info['LocalTime']
    all_az, all_ze = fpiinfo.angle_correction(all_az, all_ze, instrument['name'], dn)

    direc_idx, dists = sort_look_directions(valid_az, valid_ze, all_az, all_ze, direc_tol)

    # Determine the filter used
    lam0=instrument['lam0']
    if 'G' in sky_line_tag:
        lam0 = instrument['lam0_gl']
    last_lamc = lam0 # nominal value of line center
    lamc_start_set = False # Becomes true when a reasonable search interval center is found

    my_params = dict()
    Novertemp = 0

    for i, fname in enumerate(sky):

        ignore_this_one = False

        # Read in the sky image
        d = ReadIMG(fname)
        img = np.asarray(d)

        # Check the CCD temperature and send an email if necessary
        if d.info['CCDTemperature'] > MAX_CCD_TEMP:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                    'WARNING: %s was taken with a CCD temperature of %.1f C. It will be ignored.\n' % (fname, d.info['CCDTemperature']))
            notify_the_humans = True
            warning_log.add(message='%s was taken with a CCD temperature of %.1f C. It will be ignored.' % (fname, d.info['CCDTemperature']),
                            title='Image with high CCD temperature',
                            label='warning:ccd_temperature')
            ignore_this_one = False # We'll still analyze these images, but they should get a quality flag of 1.
            Novertemp += 1

        # Calculate the annuli to use for this time
        dt = (local.localize(d.info['LocalTime'])-lt0).seconds
        annuli = FindEqualAreas(img,cx(dt),cy(dt),N)

        # Perform the annular summation
        sky_spectra, sky_sigma = AnnularSum(img,annuli,0)

        # See if there are any NaNs in the spectra.  If so, skip this one.
        if np.isnan(sky_spectra).sum() == 0 and not ignore_this_one:

            # Append the time of the sky image
            sky_times.append(local.localize(d.info['LocalTime']))
            intT = d.info['ExposureTime']
            sky_intT.append(intT)
            sky_temperature.append(d.info['CCDTemperature'])
            sky_fns.append(fname)

            # Record the look direction.
            # For the record, use "human-friendly" azimuth and zenith:
            temp_az = d.info['azAngle']
            temp_ze = d.info['zeAngle']
            if temp_ze < 0.0:
                temp_ze = abs(temp_ze)
                temp_az = np.mod(temp_az+180.0,360)
            az.append(temp_az)
            ze.append(temp_ze)
            # Change direction index into string
            idx_match = list(direc_idx[i]) # This is a list like [1] or [1,3] or []
            # How many matches were there? Do something different for each case.
            direcstr = ''
            if len(idx_match)==0: # No match was found
                direcstr = 'Unknown'
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                    'WARNING: %s has a look direction that is not recognized: Az=%.1f, Ze=%.1f.\n' % (fname, temp_az, temp_ze))
                notify_the_humans = True
#                warning_log.add(message='%s has a look direction that is not recognized: Az=%.1f, Ze=%.1f.' % (basename(fname), temp_az, temp_ze),
#                                title='Image with unknown direction',
#                                label='warning:unknown_direction')
            elif len(idx_match)==1: # unique recognized look direction
                direcstr = direc_names[idx_match[0]]
            else: # more than one match was found
                # If "Zenith" is one of them, use it, otherwise, use any of them.
                # Write a warning to the log.
                direcstr_list = [direc_names[i] for i in idx_match]
                if 'Zenith' in direcstr_list:
                    direcstr = 'Zenith'
                else:
                    direcstr = direcstr_list[0]
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') +\
                    '%s has a look direction with > 1 match: Az=%.1f, Ze=%.1f. Matches: %s. Choosing: %s\n' % (fname, temp_az, temp_ze, str(direcstr_list), direcstr))

            # If the match is too far from any known directions, warn the user
            if dists[i] > direc_tol:
                notify_the_humans = True
                logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') +\
                    'WARNING: %s has a look direction (az,ze) = (%.1f,%.1f) that is too far (%.1f degrees) from known directions. Labeling as "Unknown".\n' % (fname, temp_az, temp_ze, dists[i]))
#                warning_log.add(message='%s has a look direction (az,ze) = (%.1f,%.1f) that is too far (%.1f degrees) from known directions. Labeling as "Unknown".' % (basename(fname), temp_az, temp_ze, dists[i]),
#                                title='Image with unknown direction',
#                                label='warning:unknown_direction')
                direcstr = 'Unknown'

            direction.append(direcstr)

            # Grab the sky times
            diff = sky_times[-1] - lt0
            my_dt = diff.seconds+diff.days*86400.

            # Set the fitting parameter guesses
            sky_params = Parameters()


            # Use the spline-interpolated values for the instrument function at this time
            for param in laser_spfits.keys():
                sky_params.add(param, value = laser_spfits[param](my_dt), vary = False)

            # But change the sky wavelength to lam0
            sky_params['lam'].value = lam0

            # Set up the forward model
            L = 301
            A_1D, lamvec = get_conv_matrix_1D(sky_params, annuli['r'][N0:N1], L, lam0)

            # "Ignore" first and last few samples due to the way that the blurring in the
            # forward model works along the edges of the image
            sky_sigma[N0:N0+5] = sky_sigma[N0:N0+5]*10
            sky_sigma[N1-5:N1] = sky_sigma[N1-5:N1]*10

            # Come up with good initial guesses for the sky parameters
            skyI_guess = sky_spectra[N0:N1].max() - sky_spectra[N0:N1].min()
            skyB_guess = 0.0
            ccdB_guess = sky_spectra[N0:N1].min()

            sky_params.add('lamc', value = lam0,       vary = True ) # We'll do a grid search to find a good value for lamc
            sky_params.add('T',    value = 1000,       vary = True , min = 20., max = 5000.0)
            sky_params.add('skyI', value = skyI_guess, vary = True , min = 0.0)
            sky_params.add('skyB', value = skyB_guess, vary = True )
            sky_params.add('skym', value = 0.0,        vary = False) # Don't try to estimate skym (or should we?)
            sky_params.add('ccdB', value = ccdB_guess, vary = True )
            sky_params.add('lam0', value = lam0,       vary = False) # This is determined by chemistry so it can't change.

            # Do a grid search to find a good starting value for lamc.
            # The instrument bias might mean that lam0 is a really bad value to start with.
            # Search over one FSR around the last solved-for value (or the nominal value, if it's the first run).
            # It's not really important to sort by direction, since the instrument bias dominates.
            def goodness(lamc):
                # return the correlation between the fringe with this lamc and the data
                sky_params['lamc'].value = lamc
                fringe = Sky_FringeModel(sky_params, annuli['r'][N0:N1], lamvec, A_1D)
                return np.dot(fringe, sky_spectra[N0:N1])/(np.linalg.norm(fringe)*np.linalg.norm(sky_spectra[N0:N1]))
            FSR = lam0**2/(2*sky_params['t'].value)

            lamc_candidates = np.linspace(last_lamc - FSR/2, last_lamc + FSR/2, 50)
            corrvec = np.array([goodness(lamc) for lamc in lamc_candidates])
            best_lamc = lamc_candidates[corrvec.argmax()]
            sky_params['lamc'].value = best_lamc

            # Take the inversion in steps
            order = [
                     ['skyI'],\
                     ['skyI','ccdB'],\
                     ['skyI', 'ccdB', 'skyB'], \
                     ['skyI', 'ccdB', 'skyB', 'lamc', 'T'], \
                     ]

            for group in order: # set all the params in this group to "vary=True", and then run inversion
                for param in list(sky_params.keys()):
                    if param in group:
                        sky_params[param].vary = True
                    else:
                        sky_params[param].vary = False

                sky_fit = Minimizer(Sky_Residual,sky_params,fcn_args=(annuli['r'][N0:N1],lamvec,A_1D), \
                          fcn_kws={'data': sky_spectra[N0:N1], 'sigma': sky_sigma[N0:N1]}, scale_covar=True)
                sky_fit.prepare_fit()
                result = sky_fit.leastsq()
                sky_params=result.params

            # Redo the fit using only points near the spectral peak (determined by fringefactor)
            if FRINGEFACTOR < 1.0:
                alpha = sky_params['alpha'].value
                t =     sky_params['t'].value
                n =     sky_params['n'].value
                lamc =  sky_params['lamc'].value

                r = annuli['r'][N0:N1]
                m = 2*n*t/lamc*np.cos(alpha*r)

                idx = abs(m - m.round()) < FRINGEFACTOR/2

                ####### TEMPORARY ########
                # only use one fringe for a3o
                if site['Abbreviation'] == 'a3o':
                    # Find the nearest peak to the middle
                    m0 = round(m[round(len(m)/2)])
                    idx = abs(m - m0) < FRINGEFACTOR/2
                ##########################

                r2 = r[idx]
                temp = sky_spectra[N0:N1]
                sky_spectra2 = temp[idx]
                temp = sky_sigma[N0:N1]
                sky_sigma2 = temp[idx]

                A_1D2, lamvec2 = get_conv_matrix_1D(sky_params, r2, L, lam0)

                sky_fit = Minimizer(Sky_Residual,sky_params,fcn_args=(r2,lamvec2,A_1D2), \
                          fcn_kws={'data': sky_spectra2, 'sigma': sky_sigma2})
                sky_fit.prepare_fit()
                sky_fit.scale_covar = True
                result = sky_fit.leastsq()
                sky_params=result.params

             # if the inversion failed, replace values with nans and set error bars to inf
            if not result.success or not result.errorbars or np.isnan(result.params['lamc'].stderr):
                for p in result.params:
                    result.params[p].value = np.nan
                    result.params[p].stderr = np.inf

            sky_redchi.append(result.redchi)
            sky_out.append(result)
            sky_fringes.append(sky_spectra[N0:N1])
            sky_annuli.append(annuli['r'][N0:N1])

            # Calculate calibration wind error
            sigma_cal_v = 0.0 # default to 0 for zenith reference
            if uselaser: # if laser reference, calculate calibration error
                caltref = laser_times[0]
                calt = np.array([(last - caltref).total_seconds() for last in laser_times])
                skyt = (sky_times[-1] - caltref).total_seconds()
                sigma_cal_v = calc_wind_calib_error(calt, laser_value['t'], skyt)

                if np.isnan(sigma_cal_v):
                    notify_the_humans = True
                    warning_log.add(message='calc_wind_calib_error returned NaN. Sky images are not bounded by laser images.',
                                    title='Interpolation error in wind calibration',
                                    label='warning:wind_calibration')

            # Transform from lamc to v and collect all parameters
            lamc = sky_params['lamc'].value
            sigma_lamc = sky_params['lamc'].stderr
            # Incorporate error from laser images (etalon gap estimate) if using laser reference.
            # Assume error in gap is constant and adds in quadrature.
            if reference=='laser':
                sigma_lamc_gap = np.median(laser_stderr['t']) * lam0/np.median(laser_value['t'])
                sigma_lamc = sqrt(sigma_lamc**2 + sigma_lamc_gap**2)
            # convert Doppler shift to velocity
            v = c * (lamc/lam0 - 1)
            sigma_fit_v = c * sigma_lamc / lam0
            LOSwind.append(v) # positive LOSwind is away from instrument
            sigma_fit_LOSwind.append(sigma_fit_v) # fit error
            sigma_cal_LOSwind.append(sigma_cal_v)
            sigma_LOSwind.append(sqrt(sigma_fit_v**2 + sigma_cal_v**2))
            T.append(result.params['T'].value)
            sigma_T.append(result.params['T'].stderr)
            # Normalize skyI and skyB by integration time
            skyI.append(sky_params['skyI'].value/intT)
            sigma_skyI.append(sky_params['skyI'].stderr/intT)
            skyB.append(sky_params['skyB'].value/intT)
            sigma_skyB.append(sky_params['skyB'].stderr/intT)
            ccdB.append(sky_params['ccdB'].value)
            sigma_ccdB.append(sky_params['ccdB'].stderr)
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + '%s reduced chisqr: %.4f.\n\t\t\tMessage:"%s"\n' % (fname, result.redchi, result.message))


            # use the lamc value to center the next grid search, if the error bars are small enough
            if sigma_fit_v < SIGMA_V_THRESH and sigma_fit_v > 0 and not lamc_start_set: # sigma_v==0 means it didn't converge
                last_lamc = lamc # save for next time
                lamc_start_set = True


        else:
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + '%s is invalid and not analyzed\n' % fname)

    # Warn again if the CCD was too hot
    if Novertemp > 0:
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'WARNING: %i/%i sky exposures were taken with a CCD temperature that was too hot. They were flagged with the quality flag.\n' % (Novertemp, len(sky)))

    # Warn if dynamic exposure time wasn't working:
    if np.std(sky_intT) == 0.0:
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'WARNING: Exposure time is not changing.\n')
        notify_the_humans = True
        warning_log.add(message='Exposure time is not changing.',
                        title='Exposure time is not changing',
                        label='warning:static_exposure_time')

    # Convert sky_params to array
    n = len(sky_out)
    if len(sky_out)==0:
        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + 'No usable sky images found.\n')
        if created_logfile:
            logfile.close()
        raise Exception('No usable sky images found.')
    sky_value = {}
    sky_stderr = {}
    for p in sky_out[-1].params:
        sky_value[p] = np.zeros(n)
        sky_stderr[p] = np.zeros(n)
        for i in range(0,n):
            sky_value[p][i] = sky_out[i].params[p].value
            sky_stderr[p][i] = sky_out[i].params[p].stderr


    # Close the logfile if we created it
    if created_logfile:
            logfile.close()

    # If the user specified a figure, draw some diagnostic plots to it
    if diagnostic_fig:

        fig = diagnostic_fig # easier to type

        fontP = FontProperties()
        fontP.set_size('small')
        #mpl.rcParams.update({'font.size': 5})
        #fontsize = 8

        ################ Plot center location ######################
        # only plot this one if uselaser and centerfinding succeeded for > 1 point
        if uselaser and center is not None and len(center)>0 and np.shape(center)!=(2,):
            tdiffvec = [laser_time - lt0 for laser_time in laser_times_center]
            dt = np.array([tdiff.seconds + tdiff.days*86400. for tdiff in tdiffvec])
            t = np.linspace(0, dt.max(), 500)
            datet = [lt0 + datetime.timedelta(seconds=ts) for ts in t]

            xdev = center[:,0] - center[0,0]
            ydev = center[:,1] - center[0,1]
            xdevi = cx(t) - center[0,0]
            ydevi = cy(t) - center[0,1]

            ax = fig.add_subplot(421)
            ax.plot(datet, xdevi, 'b', label='x')
            ax.plot(laser_times_center, xdev, 'b.')
            ax.plot(datet, ydevi, 'r', label='y')
            ax.plot(laser_times_center, ydev, 'r.')
            ax.legend(loc='best', prop=fontP)
            ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
            ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
            m = 0.5 # maximum deviation to show on plot
            if all(xdev < m) and all(xdev > -m) and all(ydev < m) and all(ydev > -m): # Set the x and y lims at +/- m
                ax.set_ylim([-m, m])

            ax.set_ylabel('Center deviation\n[pixels]')
            #ax.set_xlabel('Universal Time, [hours]')
            ax.set_title(site['Abbreviation'] + ':' + \
                laser_times_center[0].strftime(' %d %b, %Y %H:%M LT') + ' - ' + laser_times_center[-1].strftime('%H:%M LT') )

            ax.grid(True)



        ####################### Laser Fit Chi^2 #######################

        if uselaser:
            ax = fig.add_subplot(422)
            ax.plot(laser_times, laser_redchi,'k.-')
            ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
            ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
            #ax.set_xlabel('Universal Time')
            ax.set_ylabel('Laser Fit\nReduced Chi^2')
            ax.set_title(site['Abbreviation'] + ':' + \
                    laser_times[0].strftime(' %d %b, %Y %H:%M LT') + ' - ' + laser_times[-1].strftime('%H:%M LT') )
            ax.grid(True)

        ####################### Spline fit for I #######################
        # Show the laser intensity halfway out in the spectrum
        if uselaser:
            param = 'I'
            r = 0.5 # r/rmax, where r is the radius of the radial bin at which to measure the intensity
            I =  np.array([o.params['I'].value for o in output])
            a1 = np.array([o.params['a1'].value for o in output])
            a2 = np.array([o.params['a2'].value for o in output])
            Ir = I * ( 1 + a1*r + a2*r**2 )
            # assume errorbar at r=rmax/2 is approximately equal to that at r=0
            Ie = np.array([o.params[param].stderr for o in output])
            ax = fig.add_subplot(423)
            ax.errorbar(laser_times,Ir,yerr=Ie,fmt='k.-')
            ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
            ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
            ax.set_ylabel('Laser Intensity\n[counts]')
            #ax.set_xlabel('Universal Time')
            #title = ax.set_title(site['Abbreviation'] + ':' + \
            #         laser_times[0].strftime(' %d %b, %Y %H:%M LT') + ' - ' + laser_times[-1].strftime('%H:%M LT') )
            #title.set_y(1.09)
            ax.grid(True)


        ####################### Spline fit for t #######################
        # Show the spline fit for a certain parameter
        if uselaser:
            param = 't'
            dt = []
            for x in laser_times:
                diff = (x - lt0)
                dt.append(diff.seconds+diff.days*86400.)
            dt_vec = np.linspace(0, max(dt), 500)
            p_vec = laser_spfits[param](dt_vec)
            p = [o.params[param].value for o in output]
            pe = [o.params[param].stderr for o in output]
            datet_vec = [laser_times[0] + datetime.timedelta(seconds=ts) for ts in dt_vec]

            ax = fig.add_subplot(425)
            ax.plot(datet_vec,p_vec,'k')
            ax.errorbar(laser_times,p,yerr=pe,fmt='k.')
            ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
            ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
            ax.set_ylabel('Etalon Gap, [m]')
            #ax.set_xlabel('Universal Time')
            #title = ax.set_title(site['Abbreviation'] + ':' + \
            #         laser_times[0].strftime(' %d %b, %Y %H:%M LT') + ' - ' + laser_times[-1].strftime('%H:%M LT') )
            #title.set_y(1.09)
            ax.grid(True)



        ##################### Look Direction Validation #######################

        # Plot the look directions
        all_ze = np.array(all_ze)
        all_az = np.array(all_az)
        valid_ze = np.array(valid_ze)
        valid_az = np.array(valid_az)

        # Flip az for negative ze angles
        idx = all_ze < 0
        all_ze[idx] = -all_ze[idx]
        all_az[idx] = all_az[idx] + 180
        idx = valid_ze < 0
        valid_ze[idx] = -valid_ze[idx]
        valid_az[idx] = valid_az[idx] + 180

        az_rad = all_az * pi/180.0
        valid_az_rad = valid_az * pi/180.0

        ax = fig.add_subplot(428, projection='polar')
        ax.plot(valid_az_rad, valid_ze, 'kx', label = 'valid')
        valid_idx = [d != 'Unknown' for d in direction]
        invalid_idx = [not i for i in valid_idx]
        ax.plot(az_rad[np.where(valid_idx)], all_ze[np.where(valid_idx)], 'k.', label = 'actual')
        ax.plot(az_rad[np.where(invalid_idx)], all_ze[np.where(invalid_idx)], 'r.', label = 'unrecognized')
        # Now make it look like a cardinal plot, not a math plot
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi/2)
        ax.set_rmax(90.)
        ax.set_rgrids([30,60])
        if sky_line_tag=='X':
            ax.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.,numpoints=1,prop=fontP)
        elif 'G' in sky_line_tag:
            ax.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.,numpoints=1,prop=fontP,title='Green')
        #title = ax.set_title(site['Abbreviation'] + ':' +\
        #        laser_times[0].strftime(' %d %b, %Y %H:%M LT') + ' - ' + laser_times[-1].strftime('%H:%M LT') )
        #title.set_y(1.09)
        # Indicate any non-displayed points (zenith > 90)
        nnshown = sum(abs(all_ze) > 90.)
        if nnshown > 0:
            ax.set_xlabel('%03i points not shown (|ze| > 90)' % nnshown)
        # notify humans if there are Unknown look directions
        if sum(invalid_idx) > 0.:
            notify_the_humans = True # The skyscanner might be acting up.
            logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
                'WARNING: %03i unrecognized look directions. Is the skyscanner acting up?\n' % sum(invalid_idx))
            warning_log.add(message='%03i unrecognized look directions. Is the skyscanner acting up?' % sum(invalid_idx),
                            title='Images with unknown direction',
                            label='warning:unknown_direction')

        ##################### Sky Fit Chi^2 #######################

        ax = fig.add_subplot(424)
        ax.plot(sky_times, sky_redchi,'k.-')
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        ax.set_xlim([sky_times[0] - datetime.timedelta(hours=0.5), sky_times[-1] + datetime.timedelta(hours=0.5)])
        #ax.set_xlabel('Universal Time')
        if sky_line_tag=='X':
            ax.set_ylabel('Sky Fit\nReduced Chi^2')
        elif 'G' in sky_line_tag:
            ax.set_ylabel('Sky Fit (Green Images)\nReduced Chi^2')
        #title = ax.set_title(site['Abbreviation'] + ':' + \
        #        laser_times[0].strftime(' %d %b, %Y %H:%M LT') + ' - ' + laser_times[-1].strftime('%H:%M LT') )
        ax.grid(True)

        ####################### Plot of skyI #######################
        ax = fig.add_subplot(426)
        for direc in list(set(direction)):
            # account for different exposure times
            I = np.array([si for (si,d) in zip(skyI, direction) if d == direc])
            t = np.array([si for (si,d) in zip(sky_times, direction) if d == direc])
            ax.semilogy(t, I, '.-', label=direc)
        tp0 = sky_times[0] - datetime.timedelta(hours=0.5)
        tp1 = sky_times[-1] + datetime.timedelta(hours=0.5)
        sky_thresh = instrument['skyI_quality_thresh']
        ax.semilogy([tp0, tp1],[sky_thresh[0], sky_thresh[0]],'k--',lw=0.5,label='qual thresh (q=1)')
        ax.semilogy([tp0, tp1],[sky_thresh[1], sky_thresh[1]],'k--',lw=0.5,label='qual thresh (q=2)')
        ax.set_xlim([tp0, tp1])
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        if sky_line_tag=='X':
            ax.set_ylabel('Line Intensity\n[arbitrary]')
        elif 'G' in sky_line_tag:
            ax.set_ylabel('Green Line Intensity\n[arbitrary]')
        ax.set_xlabel('Universal Time')
        ax.legend(loc='best', prop={'size':6}, numpoints=1, ncol=5, framealpha=0.5)
        ax.grid(True)


        # Do some aesthetic things
        #fig.subplots_adjust(hspace = 0.35, wspace = 0.3) # Is this still necessary if we use tight_layout?

        logfile.write(datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S %p: ') + \
            'Created diagnostic plots.\n')



    # Pack and return (TODO: IMPROVE THIS!)
    FPI_Results = { 'sky_times': np.array(sky_times),\
                    'sky_chisqr': np.array(sky_redchi),\
                    'skyI': np.array(skyI),\
                    'skyB': np.array(skyB),\
                    'ccdB': np.array(ccdB),\
                    'sigma_skyI': np.array(sigma_skyI),\
                    'sigma_skyB': np.array(sigma_skyB),\
                    'sigma_ccdB': np.array(sigma_ccdB),\
                    'LOSwind': np.array(LOSwind),\
                    'sigma_fit_LOSwind': np.array(sigma_fit_LOSwind),\
                    'sigma_cal_LOSwind': np.array(sigma_cal_LOSwind),\
                    'sigma_LOSwind': np.array(sigma_LOSwind),\
                    'az': np.array(az), 'ze': np.array(ze),\
                    'T': np.array(T),\
                    'sigma_T': np.array(sigma_T),\
                    'direction': direction,\
                    'laser_times': np.array(laser_times),\
                    'laser_chisqr': np.array(laser_redchi),\
                    'laser_value': laser_value,\
                    'laser_stderr': laser_stderr,\
                    'sky_value': sky_value,\
                    'sky_stderr': sky_stderr,\
                    'laser_intT': np.array(laser_intT),\
                    'laser_ccd_temperature': np.array(laser_temperature),\
                    'sky_intT': np.array(sky_intT),\
                    'sky_ccd_temperature': np.array(sky_temperature),\
                    'reference': reference,\
                    'lam0': lam0,\
                    'sky_fns': sky_fns,\
                    'laser_fns':laser_fns,\
                    'center_pixel':center,\
                    'laser_fringes':laser_fringes,\
                    'laser_annuli':laser_annuli,\
                    'sky_fringes':sky_fringes,\
                    'sky_annuli':sky_annuli,\
                    }

    # Apply Doppler reference
    dref,drefe = DopplerReference(FPI_Results, reference=reference)
    FPI_Results['LOSwind'] = FPI_Results['LOSwind'] - dref
    FPI_Results['doppler_reference'] = dref
    # For now, do not add extra error bar Doppler reference,
    # because it is small and it's not clear how it should
    # be calculated.

    return (FPI_Results, notify_the_humans)




def dt2h(dt, tz=None):
# Function to return the decimal hour (in LT, by default) of the provided datetime
# Following http://datadebrief.blogspot.com/2010/10/plotting-sunrise-sunset-times-in-python.html
#
# INPUTS
#       dt - datetime array to convert to decimal hour
# OPTION INPUTS:
#       tz - if given, will convert to the requested timezone
# Written by Jonathan J. Makela on 25 Nov 2012

        out = np.zeros(len(dt))

        count = 0

        for d in dt:
                if tz is None:
                        # Use the timezone for each measurement
                        if d.hour > 12:
                                out[count] = d.hour + d.minute/60. + d.second/3600.
                        else:
                                out[count] = d.hour+24 + d.minute/60. + d.second/3600.
                else:
                        # Use the requested tz
                        if d.astimezone(tz).hour > 12:
                                out[count] = d.astimezone(tz).hour + d.astimezone(tz).minute/60. + d.astimezone(tz).second/3600.
                        else:
                                out[count] = d.astimezone(tz).hour+24 + d.astimezone(tz).minute/60. + d.astimezone(tz).second/3600.
                count += 1

        return (out)

def bin_and_mean(dates, values, errors, bins=np.arange(17,32,0.25)):
# Function to bin up values into bins defined by hmin to hmax with a step of
# dh (all times in local time hours).  Returns the weigthed average in each
# bin as well as the weigthed standard deviation
#
# INPUTS:
#       dates - array of datetimes (localized)
#       values - values to be binned
#       errors - errors of values to be binned (1/errors**2 used as weights)
#
# OPTIONAL INPUTS:
#       bins - localtime for the binning (default np.arange(17,32,0.25))
#
# OUTPUT:
#       bin_means, bin_stds - the mean and standard deviation of values in each
#                             bin
#
# HISTORY:
#       Written by Jonathan J. Makela on 25 Nov 2012

        # Create bins to use
##        bins = np.arange(hmin, hmax, dh)

        # Transform the datetime vector into decimal hour
        h = dt2h(dates)

        # Find which bin each balue belongs to
        digitized = np.digitize(h, bins)

        # Perform the binning
        bin_vals = [weighted_avg_and_std(values[digitized == i], 1/(errors[digitized == i]**2)) for i in range(1, len(bins))]
        bin_means = np.array(bin_vals)[:,0]
        bin_stds = np.array(bin_vals)[:,1]

        return (bin_means, bin_stds)

def average_data(files, bins=np.arange(17,32,0.25),
                 Tmin=500, Tmax=1500, Dmin=-200, Dmax=200,
                 cloudy_temperature = -15.0, reference='laser'):
# Function to create average and standard deviations for a desired
# set of days.  Returns temperatures, zonal, and meridional winds
#
# INPUTS:
#       files - files to be used in creating the mean/std
#
# OPTION INPUTS:
#       bins - localtime for the binngin (Default np.arange(17,32,0.25))
#       Tmin, Tmax - min and max values for "valid" temperature data
#       Dmin, Dmax - min and max values for "valid" Doppler data
#
# OUTPUT:
#       center_time - time of the center of each bin
#       (bin_T, bin_eT) - mean and standard deviation of temperautres
#       (bin_U, bin_eU) - mean and standard deviation of zonal winds
#       (bin_V, bin_eV) - mean and standard deviation of meridional winds
#
# HISTORY:
#       Written on 6 Dec 2012 by Jonathan J. Makela

        # Create empty arrays
        st_u = []
        st_v = []
        st_u2 = []
        st_v2 = []
        st_T = []
        all_T = []
        all_eT = []
        u_dop = []
        eu_dop = []
        u2_dop = []
        eu2_dop = []
        v_dop = []
        ev_dop = []
        v2_dop = []
        ev2_dop = []
        center_time = []
        bin_T = []
        bin_eT = []
        bin_U = []
        bin_eU = []
        bin_U2 = []
        bin_eU2 = []
        bin_V = []
        bin_eV = []
        bin_V2 = []
        bin_eV2 = []

        # Load in each file
        for f in files:
            try:
                # Load in the data
                npzfile = np.load(f,allow_pickle=True)
                FPI_Results = npzfile['FPI_Results']
                FPI_Results = FPI_Results.reshape(-1)[0]
                site = npzfile['site']
                site = site.reshape(-1)[0]
                npzfile.close()

                # Calculate the day of year
                doy = FPI_Results['sky_times'][0].timetuple().tm_yday

                # Find the zero offset of the Doppler shift
                (ref_Dop, e_ref_Dop) = DopplerReference(FPI_Results,reference=reference)

                # Calculate the vertical wind and interpolate it
                ind = all_indices('Zenith',FPI_Results['direction'])
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
                    raise Exception('%s: No Zenith look direction' %f)

                # Interpolate
                w2 = interpolate.interp1d(dt[ind],w[ind],bounds_error=False,fill_value=0.0)
                sigma_w2 = interpolate.interp1d(dt[ind],sigma_w[ind],bounds_error=False,fill_value=0.0)
                dt = []

                for x in FPI_Results['sky_times']:
                    diff = (x - FPI_Results['sky_times'][0])
                    dt.append(diff.seconds+diff.days*86400.)
                w = w2(dt)
                sigma_w = sigma_w2(dt)

                for x in np.unique(FPI_Results['direction']):
                    if x == 'Zenith':
                        ind = all_indices(x, FPI_Results['direction'])

                        # Grab the data
                        st = FPI_Results['sky_times'][ind]
                        T = FPI_Results['T'][ind]
                        eT = FPI_Results['sigma_T'][ind]
                        dop = (FPI_Results['LOSwind'][ind]-ref_Dop[ind])
                        # TODO: ERROR BASED ON BOTH LOS AND ZENITH REFERENCE
                        edop = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)

                        # Check if clouds are provided
                        if 'Clouds' in list(FPI_Results.keys()):
                            if FPI_Results['Clouds'] is not None:
                                clouds = FPI_Results['Clouds']['mean'][ind]
                                idx = clouds < cloudy_temperature
                                st = st[idx]
                                T = T[idx]
                                eT = eT[idx]

                        # Find bad data points
                        idx = (eT < 100) & (eT > 0)
                        st = st[idx]
                        T = T[idx]
                        eT = eT[idx]

                        if len(st) > 0:
                            ## Bin the data
                            #(bin_T, bin_eT) = FPI.bin_and_mean(st,T,eT,bins)

                            # Bin the data
                            st_T = np.hstack((st_T,st))
                            all_T = np.hstack((all_T,T))
                            all_eT = np.hstack((all_eT,eT))

                    # Use east for zonal wind
                    # TODO: GENERALIZE THIS TO ALLOW USER TO DETERMINE WHICH DIRECTION
                    if x == 'East':
                        ind = all_indices(x,FPI_Results['direction'])

                        # Grab the data
                        st = FPI_Results['sky_times'][ind]
                        dop = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
                        # TODO: ERROR BASED ON BOTH LOS AND ZENITH REFERENCE
                        edop = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)

                        # Check if clouds are provided
                        if 'Clouds' in list(FPI_Results.keys()):
                            if FPI_Results['Clouds'] is not None:
                                clouds = FPI_Results['Clouds']['mean'][ind]
                                idx = clouds < cloudy_temperature
                                st = st[idx]
                                dop = dop[idx]
                                edop = edop[idx]

                        # Find bad data points
                        idx = (edop < 50) & (edop > 0)
                        st = st[idx]
                        dop = dop[idx]
                        edop = edop[idx]

                        # Bin the data
                        if len(st) > 0:
                            #(bin_dop, bin_edop) = FPI.bin_and_mean(st,dop,edop,bins)

                            # Bin the data
                            st_u = np.hstack((st_u,st))
                            u_dop = np.hstack((u_dop,dop))
                            eu_dop = np.hstack((eu_dop,edop))

                    if x == 'West':
                        ind = all_indices(x,FPI_Results['direction'])

                        # Grab the data
                        st = FPI_Results['sky_times'][ind]
                        dop = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
                        # TODO: ERROR BASED ON BOTH LOS AND ZENITH REFERENCE
                        edop = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)

                        # Check if clouds are provided
                        if 'Clouds' in list(FPI_Results.keys()):
                            if FPI_Results['Clouds'] is not None:
                                clouds = FPI_Results['Clouds']['mean'][ind]
                                idx = clouds < cloudy_temperature
                                st = st[idx]
                                dop = dop[idx]
                                edop = edop[idx]

                        # Find bad data points
                        idx = (edop < 50) & (edop > 0)
                        st = st[idx]
                        dop = -dop[idx] # - because west
                        edop = edop[idx]

                        # Bin the data
                        if len(st) > 0:
                            #(bin_dop, bin_edop) = FPI.bin_and_mean(st,dop,edop,bins)

                            # Bin the data
                            st_u2 = np.hstack((st_u2,st))
                            u2_dop = np.hstack((u2_dop,dop))
                            eu2_dop = np.hstack((eu2_dop,edop))

                    # Use north for meridional wind
                    # TODO: GENERALIZE THIS TO ALLOW USER TO DETERMINE WHICH DIRECTION
                    if x == 'North':
                        ind = all_indices(x,FPI_Results['direction'])

                        # Grab the data
                        st = FPI_Results['sky_times'][ind]
                        dop = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
                        # TODO: ERROR BASED ON BOTH LOS AND ZENITH REFERENCE
                        edop = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)

                        # Check if clouds are provided
                        if 'Clouds' in list(FPI_Results.keys()):
                            if FPI_Results['Clouds'] is not None:
                                clouds = FPI_Results['Clouds']['mean'][ind]
                                idx = clouds < cloudy_temperature
                                st = st[idx]
                                dop = dop[idx]
                                edop = edop[idx]

                        # Find bad data points
                        idx = (edop < 50) & (edop > 0)
                        st = st[idx]
                        dop = dop[idx]
                        edop = edop[idx]

                        # Bin the data
                        if len(st) > 0:
                            #(bin_dop, bin_edop) = FPI.bin_and_mean(st,dop,edop,bins)

                            # Bin the data
                            st_v = np.hstack((st_v,st))
                            v_dop = np.hstack((v_dop,dop))
                            ev_dop = np.hstack((ev_dop,edop))

                    if x == 'South':
                        ind = all_indices(x,FPI_Results['direction'])

                        # Grab the data
                        st = FPI_Results['sky_times'][ind]
                        dop = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
                        # TODO: ERROR BASED ON BOTH LOS AND ZENITH REFERENCE
                        edop = np.sqrt(FPI_Results['sigma_LOSwind'][ind]**2+sigma_w[ind]**2)

                        # Check if clouds are provided
                        if 'Clouds' in list(FPI_Results.keys()):
                            if FPI_Results['Clouds'] is not None:
                                clouds = FPI_Results['Clouds']['mean'][ind]
                                idx = clouds < cloudy_temperature
                                st = st[idx]
                                dop = dop[idx]
                                edop = edop[idx]

                        # Find bad data points
                        idx = (edop < 50) & (edop > 0)
                        st = st[idx]
                        dop = -dop[idx] # Minus sign because of south
                        edop = edop[idx]

                        # Bin the data
                        if len(st) > 0:
                            #(bin_dop, bin_edop) = FPI.bin_and_mean(st,dop,edop,bins)

                            # Bin the data
                            st_v2 = np.hstack((st_v2,st))
                            v2_dop = np.hstack((v2_dop,dop))
                            ev2_dop = np.hstack((ev2_dop,edop))

            except:
                print((f + 'does not exist'))


        center_time = (bins[0:-1]+bins[1:])/2.

        if(len(st_T) > 0):
            # Perform the binning
            mask_T = (np.isfinite(all_T)) & (all_eT < 100) & (all_T < Tmax) & (all_T > Tmin)
            (bin_T, bin_eT) = bin_and_mean(st_T[mask_T],all_T[mask_T],all_eT[mask_T],bins=bins)

        if(len(st_u) > 0):
            mask_u = (np.isfinite(u_dop)) & (eu_dop < 100) & (u_dop < Dmax) & (u_dop > Dmin)
            (bin_U, bin_eU) = bin_and_mean(st_u[mask_u],u_dop[mask_u],eu_dop[mask_u],bins=bins)

        if(len(st_u2) > 0):
            mask_u2 = (np.isfinite(u2_dop)) & (eu2_dop < 100) & (u2_dop < Dmax) & (u2_dop > Dmin)
            (bin_U2, bin_eU2) = bin_and_mean(st_u2[mask_u2],u2_dop[mask_u2],eu2_dop[mask_u2],bins=bins)

        if(len(st_v) > 0):
            mask_v = (np.isfinite(v_dop)) & (ev_dop < 100) & (v_dop < Dmax) & (v_dop > Dmin)
            (bin_V, bin_eV) = bin_and_mean(st_v[mask_v],v_dop[mask_v],ev_dop[mask_v],bins=bins)

        if(len(st_v2) > 0):
            mask_v2 = (np.isfinite(v2_dop)) & (ev2_dop < 100) & (v2_dop < Dmax) & (v2_dop > Dmin)
            (bin_V2, bin_eV2) = bin_and_mean(st_v2[mask_v2],v2_dop[mask_v2],ev2_dop[mask_v2],bins=bins)

        return center_time, (bin_T, bin_eT), (bin_U, bin_eU), (bin_U2, bin_eU2), (bin_V, bin_eV), (bin_V2, bin_eV2)

def DopplerReference(FPI_Results, reference='zenith', statistic='mode'):
#
# Function to calculate the Doppler reference for a night of observations.
# By default, assumes zero vertical wind and uses the zenith measurements to
# form this reference, linearly interpolating between zenith measurement points
# to the individual look directions
#
# INPUTS:
#       FPI_Results - an FPI_Results dictionary structre defined in FPI.py
#       reference - a string indicating whether to use the laser or zenith measurements
#                   'laser': by choosing this, you trust the stability of the laser, and thus
#                            trust the vertical wind data, to within a constant offset.
##                           This offset will be estimated.
#                   'zenith': by choosing this, you do not trust, or do not have, a laser. The
#                             vertical wind will be assumed to be zero at all times.
#       statistic - only used if reference='laser'. Specifies the method used to estimate the offset:
#                   'mean': assume the mean of the vertical wind is zero (minimize L2 norm of the residual)
#                   'median': assume the median of the vertical wind is zero (minimize L1 norm of the residual)
#                   'mode': assume the "mode" of the vertical wind is zero (minimize small-L norm of the residual)
#
# OUTPUTS:
#       (ref_Dop, e_ref_Dop) - the reference and estimated uncertainty in the estimate of the Doppler reference
# HISTORY:
#       Written by Jonathan J. Makela on 10 December 2012
#       Updated by Brian J. Harding on 20 August 2014
# TODO:
#       1) Weight interpolation by uncertanties and take cloud coverage info
#               into account

        sigma_v_thresh = 150 # don't use samples with an error bar greater than this

        if reference == 'zenith':
                # First, find all zenith measurements and linearly interpolate to each time
                ind = all_indices('Zenith', FPI_Results['direction'])

                # Grab times, convert to decimal hour
                zt_h = dt2h(FPI_Results['sky_times'][ind])
                all_h = dt2h(FPI_Results['sky_times'])

                # Grab the zenith Doppler values
                z = FPI_Results['LOSwind'][ind]
                ez = FPI_Results['sigma_LOSwind'][ind]

                # Only use valid inversions
                sigma_z = FPI_Results['sigma_LOSwind'][ind]
                ind_good = sigma_z < sigma_v_thresh
                z = z[ind_good]
                ez = ez[ind_good]
                zt_h = zt_h[ind_good]

                # Perform interpolation
                if len(z) > 0:
                    ref_Dop = np.interp(all_h, zt_h, z)
                    e_ref_Dop = np.interp(all_h,zt_h,ez)
                else:
                    ref_Dop = np.nan*all_h
                    e_ref_Dop = np.nan*all_h
        elif reference == 'laser':
                # Find the zero offset
                ind = all_indices('Zenith',FPI_Results['direction'])
                if len(ind)==0:
                    raise BadDopplerReferenceError('Cannot establish Doppler reference: No zenith samples')

                # Eliminate samples with large error bars
                sigma_v = FPI_Results['sigma_LOSwind'][ind]
                ind_good = sigma_v < sigma_v_thresh
                ind = np.array(ind)[ind_good]

                if len(ind)==0:
                    raise BadDopplerReferenceError('Cannot establish Doppler reference: No trustworthy zenith samples')
                if len(ind) == 1:
                    raise BadDopplerReferenceError('Cannot establish Doppler reference: Only 1 trustworthy zenith sample')

                if np.array(FPI_Results['sigma_LOSwind'])[ind].sum() == 0:
                    ref_Dop = np.zeros(len(FPI_Results['LOSwind']))
                    e_ref_Dop = np.zeros(len(FPI_Results['LOSwind']))
                else:
                        (offset, e_offset) = weighted_avg_and_std(FPI_Results['LOSwind'][ind],1/FPI_Results['sigma_LOSwind'][ind])
                        if statistic == 'mean':
                            offset = np.mean(FPI_Results['LOSwind'][ind])
                            ref_Dop = offset*np.ones(len(FPI_Results['LOSwind']))
                            e_ref_Dop = e_offset*np.ones(len(FPI_Results['LOSwind']))
                        elif statistic == 'median':
                            offset = np.median(FPI_Results['LOSwind'][ind])
                            ref_Dop = offset*np.ones(len(FPI_Results['LOSwind']))
                            e_ref_Dop = e_offset*np.ones(len(FPI_Results['LOSwind'])) # TODO: better way?
                        elif statistic == 'mode': # find value which minimizes Lp-norm of residual, for small p
                            res = 0.1 # m/s, desired resolution
                            p = 0.1 # norm
                            rawzenith = FPI_Results['LOSwind'][ind]
                            rawsigma_zenith = FPI_Results['sigma_LOSwind'][ind]
                            # Find minimum weighted norm dopp ref.
                            mxz = np.percentile(rawzenith, 90)
                            mnz = np.percentile(rawzenith, 10)
                            if mxz < mnz: # There aren't many samples. Just use actual min/max
                                mxz = max(rawzenith)
                                mnz = min(rawzenith)
                            N = int((mxz-mnz)/res)
                            if N > 1e6: # This is probably a low quality day. Try median instead of mode
                                return DopplerReference(FPI_Results, reference=reference, statistic='median')
                            dtest = np.linspace(mnz, mxz, N)
                            cost = np.zeros(N)
                            for j in range(N):
                                d = dtest[j]
                                # doppler reference
                                zcorr = rawzenith - d
                                # evaluate metric
                                cost[j] = sum(abs(zcorr/rawsigma_zenith)**p)
                            jstar = np.argmin(cost)
                            offset = dtest[jstar]
                            ref_Dop = offset*np.ones(len(FPI_Results['LOSwind']))
                            e_ref_Dop = e_offset*np.ones(len(FPI_Results['LOSwind'])) # TODO: better way?
                        else: # invalid argument
                            raise Exception('Unrecognized parameter: "%s". Try "mean", "median", or "mode".' % statistic)
        else: # invalid argument
            raise Exception('Unrecognized parameter: "%s". Try "zenith" or "laser".' % reference)

        return (ref_Dop, e_ref_Dop/sqrt(len(ind))) # sqrt because many samples are used to determine it


