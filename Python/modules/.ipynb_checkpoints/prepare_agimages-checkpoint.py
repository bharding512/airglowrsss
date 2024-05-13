import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import asiinfo

try:
    import Image, ImageFile
except ImportError:
    from PIL import Image, ImageFile

import glob
import ASI

try:
    from skimage import filters
except ImportError:
    from skimage import filter as filters

from scipy import ndimage
from scipy import interpolate
from scipy import signal
import datetime
import multiprocessing as mp
import matplotlib.pyplot as plt

import bu_process

NUM_PROCESSORS = 16

#repo_ASI = "/Users/jonathanmakela/Data/MANGO_Data/ALLSKY/"
#daily_dir="/Users/jonathanmakela/Documents/GitHub/Analysis/"

class ASIimage:
    """
    HISTORY:
        2013 - written by L. Navarro (luis.navarro.dominguez@gmail.com)
    """
    def __init__(self,path='',calibration_params=None):
        self.path=path
        self.calibration_params=calibration_params
        self.calibration_map=False #flag to indicate if the calibration indexes were already calculated
    def execute_calibration_command(self,data,commands):
        """
        Method to apply a list of flips/rotations on data

        INPUT
            data: image data being a numpy object

        OUTPUT
            data: data after the different flips and rotations.

        HISTORY:
            2013 - written by L. Navarro (luis.navarro.dominguez@gmail.com)
        """
        for command in commands:
            command = str(command)
            if 'flipud' in command:
                data=np.flipud(data)
            elif 'fliplr' in command:
                data=np.fliplr(data)
            elif 'rot90ccw' in command:
                data=np.rot90(data, k=1)
            elif 'rot90cw' in command:
                data=np.rot90(data,k=-1)
        return data
    def rotateCoord(self,x,y,theta):
        """
        Method to rotate a set of coordinates x,y over a theta angle. Theta is positive counterclowise and ) pointing to the right.

        INPUT
            x,y: numpy objects representing x,y coordinates
            theta: angle to rotate. In degrees

        OUTPUT
            xx,yy: rotated coordinates

        HISTORY:
            2013 - written by L. Navarro (luis.navarro.dominguez@gmail.com)
        """
        theta=np.radians(theta)
        xx=x*np.cos(theta)-y*np.sin(theta)
        yy=x*np.sin(theta)+y*np.cos(theta)
        return xx,yy
    def geometricalCalibration(self,cobermax=256,H=95,latlon0=None):
        """
        Method to spatially calibrate the image using a known lens function.
        Lens funcions can be a 4th or 6th degree polynomium. It takes a zenith value in radians and returns the distance from
        the center of the image in kms.
        Based on  Garcia, Taylor, Kelley[1997] Two-dimensional spectral analysis of mesospheric airglow image data

        INPUT
            cobermax[default=256]: is the maximum distance in km to be mapped in the final image
            self.data: image on pixel coordinates
            H[default=95]: height of airglow emission. Used in the geometric conversion
            latlon0[default=None]: 2-items array of lat,lon of center point. Used to calcuilate all map in lat, lon coords

        OUTPUT
            self.lens_function: function object that take a zenith angle in radians and returns a distance in kms
            self.resolution: resolution in pixel/km. Usually 1km/pixel
            self.coberture: numpy array of coberture. Example: [-256...255]
            self.setpoint_coberture: saving coberture for future calculations
            self.geocal0: I indexes that map the image from pixels coordinates into geographic coordinates. Calculated only once for night.
            self.geocal1: J indexes that map the image from pixels coordinates into geographic coordinates. Calculated only once for night.
            self.data: image data in geographical coordinates where top border is North and resolution self.resolution

        HISTORY:
            2013 - written by L. Navarro (luis.navarro.dominguez@gmail.com)
        """
        data=np.copy(self.data)
        R=6371.2   #earth radius
        #using calibration_params
        p=self.calibration_params['coefficients']
        if len(p)==4:
            C4,C3,C2,C1=p
            lens_calibration=lambda ze:C4*ze**4+C3*ze**3+C2*ze**2+C1*ze
        else:
            C6,C5,C4,C3,C2,C1=p
            lens_calibration=lambda ze:C6*ze**6+C5*ze**5+C4*ze**4+C3*ze**3+C2*ze**2+C1*ze
        self.lens_function=lens_calibration
        ci,cj,rotation=self.calibration_params['coords']
        commands=self.calibration_params['commands']
        data=self.execute_calibration_command(data,commands)
        newdata=np.zeros((cobermax*2,cobermax*2))
        newdata[:]=np.nan
        self.calibration_map=True
        coberx=np.arange(-cobermax,cobermax)*1.
        cobery=np.arange(cobermax,-cobermax,step=-1)*1.
        self.coberture=coberx
        self.setpoint_coberture=self.coberture
        self.resolution=self.coberture[1]-self.coberture[0]#km/pixel
        self.setpoint_resolution=self.resolution
        coberx,cobery=np.meshgrid(coberx,cobery)
        self.coberx=np.copy(coberx)
        self.cobery=np.copy(cobery)
        r=np.hypot(coberx,cobery)
        psi=r/(R+H)
        alpha=(np.pi-psi)/2.
        a=2.*(R+H)*np.sin(psi/2.)
        c=np.sqrt(H**2+a**2-2.*H*a*np.cos(alpha))
        el=np.arccos(a*np.sin(alpha)/c)
        az=np.arctan2(coberx,cobery)
        ind=az<0
        az[ind]=az[ind]+2*np.pi
        self.geoelev=np.copy(np.degrees(el))
        self.geoazimuth=np.copy(np.degrees(az))
        #===================================================================
        # #getting lat,lons matrix just for display purposes
        # Calculate the differential angle, alpha
        if latlon0 is not None:
            lat_r,lon_r=np.radians(latlon0)

            temp = np.cos(el)/(1+(H/R))
            alpha = np.arccos(temp) - el

            # Calculate the pierce point latitude
            temp = np.sin(lat_r) * np.cos(alpha) + np.cos(lat_r)*np.cos(az)*np.sin(alpha)
            lat_r = np.arcsin(temp)

            #
            temp = np.sin(alpha) * np.sin(az)/np.cos(lat_r)
            lon_r = np.arcsin(temp) + lon_r

            lat = np.degrees(lat_r)
            lon = np.degrees(lon_r)
            self.geolat=np.copy(lat)
            self.geolon=np.copy(lon)
        #=====================================================================
        ze=np.pi/2-el
        d=lens_calibration(ze)
        dx=d*np.sin(az)
        dy=d*np.cos(az)
        dx,dy=self.rotateCoord(dx, dy,rotation)
        I=ci-dy
        J=cj+dx
        #interpolation
        ind1=np.logical_and(I>=0,I<=data.shape[0]-1)
        ind2=np.logical_and(J>=0,J<=data.shape[1]-1)
        ind=np.logical_and(ind1,ind2)
        I_floor=np.floor(I[ind]).astype(int)
        I_ceil=np.ceil(I[ind]).astype(int)
        J_floor=np.floor(J[ind]).astype(int)
        J_ceil=np.ceil(J[ind]).astype(int)
        ul=data[I_floor,J_floor].astype(float)
        dl=data[I_ceil,J_floor].astype(float)
        ur=data[I_floor,J_ceil].astype(float)
        dr=data[I_ceil,J_ceil].astype(float)
        tempu=(ur-ul)*(J[ind]-J_floor)/(J_ceil-J_floor)+ul
        tempd=(dr-dl)*(J[ind]-J_floor)/(J_ceil-J_floor)+dl
        newdata[ind]=(tempu-tempd)*(I[ind]-I_ceil)/(I_floor-I_ceil)+tempd
        ci=int(newdata.shape[0]/2)#for center of image
        newdata[ci,ci]=np.nanmean(newdata[ci-5:ci+6,ci-5:ci+6])#for center of image
        self.geocal0=I
        self.geocal1=J
        self.data=newdata

def readfile(ipath,data=False):
    import h5py
    fg = h5py.File(ipath,'r')
    gdata = fg['image']
    latitude=fg['image'].attrs['latitude']
    longitude=fg['image'].attrs['longitude']
    start_time=datetime.datetime(1970,1,1)+datetime.timedelta(seconds=int(fg['image'].attrs['start_time']))
    raw0=None
    if data:
        raw0=gdata[:,:]
    fg.close()
    return start_time,latitude,longitude,raw0

def readfile(ipath,data=False):
    import h5py
    fg = h5py.File(ipath,'r')
    gdata = fg['image']
    latitude=fg['image'].attrs['latitude']
    longitude=fg['image'].attrs['longitude']
    start_time=datetime.datetime(1970,1,1)+datetime.timedelta(seconds=int(fg['image'].attrs['start_time']))
    raw0=None
    if data:
        raw0=gdata[:,:]
    fg.close()
    return start_time,latitude,longitude,raw0

def latlon2xyz(lat,lon,alt):
    '''
    Convert a lat, lon, alt point to an ECEF x,y,z location
    INPUTS:
        lat - latitude of point [deg]
        lon - longitude of point [deg]
        alt - altitude of point [km]
    OUTPUTS:
        x, y, z - ECEF location [m]
    HISTORY:
        11-Mar-2014: Create by Jonathan J. Makela (jmakela@illinois.edu)
    '''

    # Convert to radians and m
    lat = lat * np.pi/180.
    lon = lon * np.pi/180.
    alt = alt * 1000.

    # Constants describing the WGS-84 ellipsoid
    a = 6378137
    b = 6356752.3145
    e = np.sqrt(1-(b**2)/(a**2))

    # Conversion
    x = a * np.cos(lon) / np.sqrt(1+(1-e**2) * (np.tan(lat))**2) + alt * np.cos(lon) * np.cos(lat)
    y = a * np.sin(lon) / np.sqrt(1+(1-e**2) * (np.tan(lat))**2) + alt * np.sin(lon) * np.cos(lat)
    z = a * (1-e**2) * np.sin(lat) / np.sqrt(1 - e**2 * (np.sin(lat))**2) + alt * np.sin(lat)

    return x,y,z

def GPS_VEN(c1, c2, c3, ECEF = False):
    """
    Revision history:
      11-Jun-2006: Initial template created by Jonathan J. Makela
      (jmakela@uiuc.edu)

      26-Apr-2014: Ported to Python by Matthew Grawe (mgrawe@gmail.com)

      if ECEF = False:  c1,c2,c3 interpreted as lat,lon,alt
      if ECEF = True: c1,c2,c3 interpreted as x,y,z in ECEF coordinates on WGS84 ellipsoid
    """

    phi  = c1
    lamb = c2
    alt  = c3

    if ECEF is True:
        phi, lamb, alt = GPS_WGS84(c1,c2,c3)

    phi = np.radians(phi)
    lamb = np.radians(lamb)

    xRow = np.array([np.cos(phi)*np.cos(lamb), np.cos(phi)*np.sin(lamb), np.sin(phi)])
    yRow = np.array([-np.sin(lamb), np.cos(lamb), 0])
    zRow = np.array([-np.sin(phi)*np.cos(lamb), -np.sin(phi)*np.sin(lamb), np.cos(phi)])
    rot  = np.array([xRow, yRow, zRow])

    return rot # VEN

def process_image(filename):
    d          = Image.open(filename)
    im         = np.reshape(d.getdata(), d.size)
    im_nonans  = im.copy().astype('float')

    nans            = np.where(np.isnan(im))
    im_nonans[nans] = 0.0

    print('Processing %s' % filename)

    im_mf = ndimage.filters.median_filter(im_nonans, size=[10,10])
    im_mf[nans] = np.nan

    return im_mf, d.info

def process_image_bu(filename):

    im = bu_process.load_bu(filename)
    im_nonans  = im.copy().astype('float')

    im_nonans[BU_NANS] = 0.0

    print('Processing %s' % filename)

    im_mf = ndimage.filters.median_filter(im_nonans, size=[10,10])
    #im_mf = im_nonans.copy() # "disables" median filtering
    im_mf[BU_NANS] = np.nan

    return im_mf

def process_image_dasi(filename, calibration_params, cobermax=300, height=95, daily_dir="/home/airglow/scratch_data/"):
    dt,lat0,lon0,currentdata=readfile(filename,data=True)

    ## TODO CHECK IF WE ACTUALLY WANT TO CHANGE THE HISTOGRAM. IS THIS MODIFYING THE ACUTAL DATA?
#    cleared=changeHistogram2(currentdata,nbins=12)
    img = ASIimage(daily_dir,calibration_params=calibration_params)
    img.data = currentdata #cleared
    img.geometricalCalibration(cobermax=cobermax,H=height,latlon0=[lat0,lon0])
#    F=(np.nanmax(img.data)-np.nanmin(img.data))/256
#    img.data=(img.data-np.nanmin(img.data))/F

#    _,mask=cv2.threshold(img.data.astype('uint8'),0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
#    mask=mask.astype(np.float)
#    idx=mask==0
#    mask[idx]=np.nan
#    idx=img.geoelev<15.
#    mask[idx]=np.nan

#    img=ASIimage(daily_dir,calibration_params=calibration_params)
#    img.data=cleared
#    img.geometricalCalibration(cobermax=cobermax,H=height,latlon0=[lat0,lon0])
#    img.data = img.data*mask

    #apply 7x7 median filter
    from scipy import signal
    img.data=signal.medfilt2d(img.data, kernel_size=7)

    return img, dt

def projection_to_regular(args):

    latij, lonij, emissionHeight, \
    localOriginX, localOriginY, localOriginZ, VEN = args

    # construct vectors from observation point to each grid point
    lx, ly, lz = latlon2xyz(latij, lonij, emissionHeight)
    lineX = lx-localOriginX
    lineY = ly-localOriginY
    lineZ = lz-localOriginZ

    RANGE_VECTOR             = np.array([lineX, lineY, lineZ])
    p                        = np.dot(VEN, RANGE_VECTOR)

    return np.array([p[1],p[2],p[0]])

def get_times(fileName):
    d  = Image.open(fileName)
    time = d.info['UniversalTime']
    return time

def prepare_airglow(instr_name, year, doy, start, stop, emission, el_cutoff, site_name_arg = None):

    success = True
    # Obtain datetime object for the imaging data
    # doy-1 so that the datetime starts at the BEGINNING of the day
    dn = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(days=doy-1)

    # Get the info for the instrument on this date
    site         = asiinfo.get_site_of(instr_name, dn)
    instr        = asiinfo.get_instr_info(instr_name, dn)
    site_info    = asiinfo.get_site_info(site)
    filter_names = instr['filter_names']
    unwarp_ht    = instr['unwarp_ht']
    t_lat        = instr['t_lat']
    t_lon        = instr['t_lon']

    loc =  site_info['Location']
    observation_point = (loc[0], loc[1])

    #if site_name_arg is None:
    #   site_name = site_info['Name']
    #else:
    #   site_name = site_name_arg

    # Height to use for unwarping
    ht = unwarp_ht[0] # kilometers
    #ht = 246
    #ht = 90

    print('projection height: %0.2f km' % ht)

    # Load the calibrations for this imager
    npzfile      = np.load(instr['cal_file'])
    el           = npzfile['el']
    az           = npzfile['az']
    rx_lat       = npzfile['rx_lat']
    rx_lon       = npzfile['rx_lon']

    print("rx_lat", rx_lat)
    print("rx_lon", rx_lon)

    npzfile.close()

    print('emission: %s' % emission)
    print('site: %s' % site)

    # Unwarp
    lat, lon = ASI.ConvertAzEl2LatLon(az,el,ht,rx_lat,rx_lon,horizon=el_cutoff)

    # Always keep the true doy as the first element in the doys list. This is assumed.
    doys        = [doy, doy-1, doy+1]
    dns         = []
    directories = []
    all_files   = []

    for doyi in doys:
        directories.append('/rdata/airglow/imaging/%s/%s/%i/%03i/' % (instr_name, site, year, doyi,))
        dns.append(datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(days=doyi-1))

    for i, directory in enumerate(directories):

        if instr_name.upper() == 'CNFI01':
            chk = len(sorted(glob.glob('%s/CNF_*.tif' % (directory))))
            if chk > 0:
                site = 'CNF'

        files = sorted(glob.glob('%s/%s_*.tif' % (directory, emission)))
        files = np.array(files)

        K = len(files)

        if emission == 'Dark' and K == 0:
            files = sorted(glob.glob('%s/%s_*.tif' % (directory, 'DARK')))
            files = np.array(files)
            K = len(files)

        if K == 0:
            #print '%s/%s_%s_%d%0.2d%0.2d_*' % (directory, site.upper(), emission, dn.year, dn.month, dns[i].day)
            files = sorted(glob.glob('%s/%s_%s_%d%0.2d%0.2d_*' % (directory, site.upper(), emission, dn.year, dn.month, dns[i].day)))
            files = np.array(files)
            K = len(files)

        if emission == 'Dark' and K == 0:
            #print '%s/%s_%s_%d%0.2d%0.2d_*' % (directory, site.upper(), emission, dn.year, dn.month, dns[i].day)
            files = sorted(glob.glob('%s/%s_%s_%d%0.2d%0.2d_*' % (directory, site.upper(), 'DARK', dn.year, dn.month, dns[i].day)))
            files = np.array(files)
            K = len(files)

        all_files.extend(files)

    all_files = np.array(sorted(all_files))
    K = len(all_files)
    if K == 0:
        if emission == 'BGND' or 'Dark' or 'DARK':
            return None, None, None, None, None, None, False
        else:
            raise Exception("Unable to locate any correctly named images for this insturment and date.")

    N    = 512

    print("Processing images..")

    # need to build the times vector first in order to determine size for IM3D
    image_range = []
    times_range = []
    pool   = mp.Pool(processes = NUM_PROCESSORS)
    times  = pool.map(get_times, all_files[range(K)])
    pool.close()
    pool.join()

    # sift through all of the times and return only the ones on the correct day
    for i in range(K):
        time = times[i]
        print(time, all_files[i])
        if time.year == dn.year and time.month == dn.month and time.day == dn.day:
            image_range.append(i)
            times_range.append(time)
            #print 'frame %d, %s' % (i, time)

    # create the container for the image sequence (space and time)
    IM3D = np.zeros((N,N,len(image_range))) # 3D Image (space and time)
    info  = []

    pool   = mp.Pool(processes = NUM_PROCESSORS)
    output = pool.map(process_image, all_files[image_range])
    pool.close()
    pool.join()

    #print '%s/%s_*.tif' % (folder, emission)

    for ki, k in enumerate(output):

        im_mf, im_info  = k

        IM3D[:,:,ki] = im_mf

    #   import matplotlib.pyplot as plt
    #   print ki,
    #   plt.imshow(im_mf)
    #   plt.savefig('NSO_%d' % ki)
    #   plt.close()
    #   info.append(im_info)

    print("Processing finished.")
    return IM3D, lat, lon, times_range, ht, observation_point, success

BU_NANS = None

def prepare_airglow_bu(station, year, doy, height, el_cutoff):

    success = True
    N = 512

    global BU_NANS

    # calibrations are hardcoded for now
    if station == 'mho':
        fn_az_cal = '/rdata/airglow/grawe/bu_star_calibrations/mho_az1.mat'
        fn_el_cal = '/rdata/airglow/grawe/bu_star_calibrations/mho_el1.mat'
        latlon_imager = (42.62, -71.5)
    elif station == 'sao':
        fn_az_cal = '/rdata/airglow/grawe/bu_star_calibrations/sao_az2.mat'
        fn_el_cal = '/rdata/airglow/grawe/bu_star_calibrations/sao_el2.mat'
        latlon_imager = (-32.37, 20.81)

    # get lat/lon calibration
    lat, lon = bu_process.bu_latlon_cal(fn_az_cal, fn_el_cal, latlon_imager, projection_height = height, horizon = el_cutoff)

    # build date from doy
    dn = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(days=doy-1)

    # get list of filenames (images, darks, backgrounds)
    # get list of timestamps for the files
    Bfiles, Cfiles, Dfiles, times_B, times_C, times_D = bu_process.get_datetime_paths(dn, station = station)
    BU_NANS = np.where(np.isnan(lat))

    if station == 'mho':
        files = Cfiles
        times = times_C
    elif station == 'sao':
        files = Dfiles
        times = times_D

    # run process_image_bu on each of the files in parallel
    pool = mp.Pool(processes = NUM_PROCESSORS)
    output = pool.map(process_image_bu, files)
    pool.close()
    pool.join()

    # create the container for the image sequence (space and time)
    IM3D = np.zeros((N, N, len(files))) # 3D Image (space and time)

    # populate
    for ki, k in enumerate(output):
       im_mf = k
       IM3D[:,:,ki] = im_mf

    return IM3D, lat, lon, times, height, latlon_imager, success

def prepare_airglow_dasi(instr_name, year, doy, emission, el_cutoff,cobermax=300, daily_dir="/home/airglow/scratch_data/", repo_ASI="/home/airglow/scratch_data/MANGO_Data",coord_nan=True):

    success = True

    # Check if we need to use the other prepare_airglow_dasi function because the 
    # calibration npz file is different
    if instr_name in ['mro','eio','cvo','low','blo','cfs','bdr']:
        return prepare_airglow_dasi2(instr_name, year, doy, emission, el_cutoff, daily_dir=daily_dir, repo_ASI=repo_ASI,coord_nan=coord_nan)

    # Obtain datetime object for the imaging data
    # doy-1 so that the datetime starts at the BEGINNING of the day
    dn = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(days=doy-1)

    if emission == '5577':
        tag = 'greenline'
        ht = 95.
        cobermax=300
    elif emission == '6300':
        tag = 'redline'
        ht = 250.
        cobermax=1250

    if instr_name == 'low':
        observation_point = [34.752,-111.423]
        calibration_path=daily_dir+'calibration_%s_%s_*.npz'%(instr_name, tag,)
        calibration_path=glob.glob(calibration_path)[0]
        cal = np.load(calibration_path,allow_pickle=True,encoding='bytes')
        calibration_params = {key:cal[key] for key in cal.keys()}
        cal.close()
    elif instr_name == 'blo':
        observation_point = [41.934,-111.421]
        calibration_path=daily_dir+'calibration_%s_%s_*.npz'%(instr_name, tag,)
        calibration_path=glob.glob(calibration_path)[0]
        cal = np.load(calibration_path,allow_pickle=True,encoding='bytes')
        calibration_params = {key:cal[key] for key in cal.keys()}
        cal.close()
    elif instr_name == 'cvo':
        observation_point = [43.233,-120.683]
        calibration_path=daily_dir+'calibration_%s_%s_*.npz'%(instr_name, tag,)
        print(calibration_path)
        calibration_path=glob.glob(calibration_path)[0]
        cal = np.load(calibration_path,allow_pickle=True,encoding='bytes')
        calibration_params = {key:cal[key] for key in cal.keys()}
        cal.close()
    elif instr_name == 'cfs':
        observation_point = [38.185,-111.179]
        calibration_path=daily_dir+'calibration_%s_%s_*.npz'%(instr_name, tag,)
        calibration_path=glob.glob(calibration_path)[0]
        cal = np.load(calibration_path,allow_pickle=True,encoding='bytes')
        calibration_params = {key:cal[key] for key in cal.keys()}
        cal.close()

    print(calibration_path)
    print('projection height: %0.2f km' % ht)

    # All of the files
    files = glob.glob('%s/%s/%s/*/*%s*.hdf5'%(repo_ASI,instr_name,dn.strftime('%Y/%j'),tag))
    files.sort()
#    print('%s/%s/%s/*/*greenline*.hdf5'%(repo_ASI,instr_name,dn.strftime('%Y/%j')))

    # Load the calibrations for this imager
    # Need an image to do this; use the first one in files
    if len(files) == 0:
        return None, None, None, None, None, None, False
    dt,lat0,lon0,currentdata=readfile(files[0],data=True)
#    cleared=changeHistogram2(currentdata,nbins=12)
    img = ASIimage(files[0],calibration_params=calibration_params)
    img.data = currentdata#cleared
    img.geometricalCalibration(cobermax=cobermax,H=ht,latlon0=observation_point)
    lat = img.geolat
    lon = img.geolon

    print("Processing images..")

    times = []
    N = 2*cobermax #600 # size of image (Need to understand why this is 600!!)
    IM3D = np.zeros((N, N, len(files)))

    for i, f in enumerate(files):
#        print(i)
        img, dt = process_image_dasi(f, calibration_params, cobermax=cobermax)
        times.append(dt)
        # Apply horizon
        ind = np.where(img.geoelev < el_cutoff)
        img.data[ind] = np.nan
        IM3D[:,:,i] = img.data

    print("Processing finished.")
    return IM3D, lat, lon, times, ht, observation_point, success

def prepare_airglow_dasi2(instr_name,year,doy,emission,el_cutoff,
                          daily_dir="/home/airglow/scratch_data/", repo_ASI="/home/airglow/scratch_data/MANGO_Data",coord_nan=True):
    
    success = True
    # Obtain datetime object for the imaging data
    # doy-1 so that the datetime starts at the BEGINNING of the day
    dn = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(days=doy-1)

    if emission == '5577':
        tag = 'greenline'
    elif emission == '6300':
        tag = 'redline'
        
    # Load the calibration file
    calibration_path=daily_dir+'calibration_%s_%s_*.npz'%(instr_name, tag,)
    calibration_path=glob.glob(calibration_path)[0]
    npzfile = np.load(calibration_path,allow_pickle='False',encoding='latin1')
    lat = npzfile['lat']
    lon = npzfile['lon']
    observation_point = npzfile['observation_point']
    ht = npzfile['ht']
    el = npzfile['el']
    az = npzfile['az']
    npzfile.close()
    
    # All of the files
    files = glob.glob('%s/%s/%s/*/*%s*.hdf5'%(repo_ASI,instr_name,dn.strftime('%Y/%j'),tag))
    files.sort()
#    print('%s/%s/%s/*/*greenline*.hdf5'%(repo_ASI,instr_name,dn.strftime('%Y/%j')))

    print("Processing images..")
    
    times = []
    IM3D = np.zeros((lat.shape[0],lat.shape[1],len(files)))
    
    for i, f in enumerate(files):
        dt,lat0,lon0,currentdata=readfile(f,data=True)
        times.append(dt)
        
        #apply 7x7 median filter
        from scipy import signal
        img=signal.medfilt2d(currentdata, kernel_size=7).astype(float)
        
        # Apply horizon
        ind = np.where(el < el_cutoff)
        
        img[ind] = np.nan
        if coord_nan:
            lat[ind] = np.nan
            lon[ind] = np.nan
        IM3D[:,:,i] = img

    print("Processing finished.")
    return IM3D, lat, lon, times, ht, observation_point, success

#prepare_airglow_bu('mho', 2018, 336, height = 250, el_cutoff = 15)

EAST_REGULAR  = None
NORTH_REGULAR = None
UP_REGULAR    = None
east_nonans   = None
north_nonans  = None
IM3Dfilt      = None
IM3Dfull      = None
nonans        = None
use_full      = None

def project_data(frameNo):

#    print(".", end = '')
#    print(frameNo)

    global IM3Dfilt
    global IM3Dfull
    global east_nonans
    global north_nonans
    global EAST_REGULAR
    global NORTH_REGULAR
    global nonans
    global use_full

    if use_full is True:
        data_mf         = IM3Dfull[:,:,frameNo]
    else:
        data_mf         = IM3Dfilt[:,:,frameNo]
    
#     print('data_mf size', data_mf.shape)

    data_mf_nonans          = np.zeros(data_mf.shape)
    data_mf_nonans[nonans]  = data_mf[nonans]

    DATA_REGULAR_MF = interpolate.griddata((np.ravel(east_nonans), np.ravel(north_nonans)), np.ravel(data_mf_nonans),\
                           (np.ravel(EAST_REGULAR), np.ravel(NORTH_REGULAR)), method='linear')

    DATA_REGULAR_MF = np.reshape(DATA_REGULAR_MF, EAST_REGULAR.shape)

    return DATA_REGULAR_MF


def remap_airglow_uniform(IM3Dfiltarg, lat, lon, emissionHeight, observation_point, returnReg = False, IM3DFullarg = None):
#     print('IM3Dfiltarg shape', IM3Dfiltarg.shape)
    global east_nonans
    global north_nonans
    global EAST_REGULAR
    global NORTH_REGULAR
    global UP_REGULAR
    global IM3Dfilt
    global IM3Dfull
    global nonans
    global use_full

    use_full = False
    IM3Dfilt = IM3Dfiltarg
    IM3Dfull = IM3DFullarg

    x,y,z = latlon2xyz(lat,lon, emissionHeight)
    
#     print('IM3Dfilt shape', IM3Dfilt.shape)

    # Take the sampling rate as the average 'dr' spacing between pixels.

    # Using adjacent column
    x2 = np.zeros(x.shape); y2 = np.zeros(y.shape); z2 = np.zeros(z.shape)
    x2[:,0:-1] = x[:,1:]; x2[:,-1] = x2[:,-2]
    y2[:,0:-1] = y[:,1:]; y2[:,-1] = y2[:,-2]
    z2[:,0:-1] = z[:,1:]; z2[:,-1] = z2[:,-2]

    dx = np.abs(x2 - x)
    dy = np.abs(y2 - y)
    dz = np.abs(z2 - z)

    dr = np.sqrt(dx**2 + dy**2 + dz**2)

    # Using adjacent row
    x3 = np.zeros(x.shape); y3 = np.zeros(y.shape); z3 = np.zeros(z.shape)
    x3[0:-1,:] = x[1::,:]; x3[-1,:] = x3[-2,:]
    y3[0:-1,:] = y[1::,:]; y3[-1,:] = y3[-2,:]
    z3[0:-1,:] = z[1::,:]; z3[-1,:] = z3[-2,:]

    dx2 = np.abs(x3 - x)
    dy2 = np.abs(y3 - y)
    dz2 = np.abs(z3 - z)

    dr2 = np.sqrt(dx2**2 + dy2**2 + dz2**2)

    a1 = np.nanmean(dr)
    a2 = np.nanmean(dr2)

    drmin = np.nanmin(dr)
    drmax = np.nanmax(dr)

    dr = np.nanmean([a1, a2])

    print('drmin: ', drmin, 'm')
    print('drmax: ', drmax, 'm')

    n_lat             = lon.shape[0]; n_lon = lon.shape[1]
    RANGE_VECTORS_ENU = np.zeros([n_lat, n_lon, 3])

    # ECEF, VEN matrix (rotation matrix) of observation point
    localOriginX, localOriginY, localOriginZ = latlon2xyz(observation_point[0],observation_point[1], 0)
    VEN                                      = GPS_VEN(observation_point[0], observation_point[1], 0, False)

    localOriginVertical, localOriginEast, localOriginNorth = \
    np.dot(VEN, np.array([localOriginX, localOriginY, localOriginZ]))/1000. # km

    print("Remapping grid...")
    pool = mp.Pool(processes = NUM_PROCESSORS)
    arglist = []
    for i in range(0, n_lat):
        for j in range(0, n_lon):
            arglist.append((lat[i,j], lon[i,j], emissionHeight, localOriginX, \
                localOriginY, localOriginZ, VEN))

    RANGE_VECTORS_ENU_LIST = pool.map(projection_to_regular, arglist)
    pool.close()
    pool.join()

    k = 0
    for i in range(0, n_lat):
        for j in range(0, n_lon):
            RANGE_VECTORS_ENU[i,j,:] = RANGE_VECTORS_ENU_LIST[k]
            k = k + 1
    
    
    east  = RANGE_VECTORS_ENU[:,:,0]/1000. # km
    north = RANGE_VECTORS_ENU[:,:,1]/1000. # km
    up    = RANGE_VECTORS_ENU[:,:,2]/1000. # km
    
#     print('east shape', east.shape)

    # remove nans
    nonans = np.where(~np.isnan(east))

    north_nonans = np.zeros(north.shape)
    east_nonans  = np.zeros(east.shape)
    up_nonans    = np.zeros(up.shape)

    north_nonans[nonans] = north[nonans]
    east_nonans[nonans]  = east[nonans]
    up_nonans[nonans]    = up[nonans]

    from scipy.interpolate import griddata
    east_regular  = np.arange(np.nanmin(east),  np.nanmax(east),  dr/1000.)
    north_regular = np.arange(np.nanmin(north), np.nanmax(north), dr/1000.)

    EAST_REGULAR, NORTH_REGULAR = np.meshgrid(east_regular, north_regular)
    UP_REGULAR                  = interpolate.griddata((np.ravel(east_nonans), np.ravel(north_nonans)), np.ravel(up_nonans),\
                           (np.ravel(EAST_REGULAR), np.ravel(NORTH_REGULAR)), method='cubic')

    UP_REGULAR = np.reshape(UP_REGULAR, EAST_REGULAR.shape)

#   matplotlib.rcParams['savefig.dpi']    = 150
#   matplotlib.rcParams['figure.figsize'] = (7,3)
#   matplotlib.rcParams['font.size']      = 6

#   a = plt.figure()
#   plt.subplot(121)
#   ax2 = plt.gca()
#   warped = ax2.scatter(east.flatten(), north.flatten(), c = up.flatten(), linewidths = 0, color = 'k', s = 0.1, vmin = 220, vmax = 250)
#   ax2.scatter(0,0, c='r', linewidths = 0, label = 'Origin (%0.2f, %0.2f)' % (observation_point[0], observation_point[1]), s = 40)
#   ax2.scatter(east.flatten()[::40], north.flatten()[::40], c = 'k', linewidths = 0, s = 1, alpha = 0.5)
#   plt.xlabel('East [km]'); plt.ylabel('North [km]'); a.colorbar(warped)
#   ax2.legend(loc = 'upper left', prop = {'size': 8}); plt.xlim(-600, 600); plt.ylim(-600, 600);
#   plt.title('Original Sampling')

#   plt.subplot(122)
#   ax3 = plt.gca()
#   euclid = ax3.scatter(EAST_REGULAR.flatten(), NORTH_REGULAR.flatten(), c = UP_REGULAR.flatten(), linewidths = 0, color = 'k', s = 0.1, vmin = 220, vmax = 250)
#   ax3.scatter(0,0, c='r', linewidths = 0, label = 'Origin (%0.2f, %0.2f)' % (observation_point[0], observation_point[1]), s = 40)
#   ax3.scatter(EAST_REGULAR.flatten()[::40], NORTH_REGULAR.flatten()[::40], c = 'k', linewidths = 0, s = 1, alpha = 0.5)
#   plt.xlabel('East [km]'); plt.ylabel('North [km]'); cb = plt.colorbar(euclid)
#   cb.set_label('Up [km]')
#   ax3.legend(loc = 'upper left', prop = {'size': 8}); plt.xlim(-600, 600); plt.ylim(-600, 600);
#   plt.title('Resampling')
#   plt.savefig('sampling_grid.png')
#   plt.close()

    K = np.shape(IM3Dfilt)[2]
    frames = range(0, K-1)

    # interpolate images onto uniform grid and load the data into a container
    data_frames  = []
    nans         = np.where(np.isnan(east))

## THESE LINES SERIALIZE AND CAN BE COMMENTED IF USING MULTICORE
##    for frame in frames:
##       print (frame)
##       data_mf         = IM3Dfilt[:,:,frame]
##       data_mf_nonans  = data_mf[nonans]
##
##       DATA_REGULAR_MF = interpolate.griddata((np.ravel(east_nonans), np.ravel(north_nonans)), np.ravel(data_mf_nonans),\
##                          (np.ravel(EAST_REGULAR), np.ravel(NORTH_REGULAR)), method='linear')
##       DATA_REGULAR_MF = np.reshape(DATA_REGULAR_MF, EAST_REGULAR.shape)
##
##       data_frames.append(DATA_REGULAR_MF)

    pool = mp.Pool(processes = NUM_PROCESSORS)
    data_frames = pool.map(project_data, frames)
    pool.close()
    pool.join()

    if returnReg is True and IM3DFullarg is not None:
        print("\nRemapping unfiltered grid...")

        use_full = True

        pool            = mp.Pool(processes = NUM_PROCESSORS)
        data_frames_reg = pool.map(project_data, frames)
        pool.close()
        pool.join()

        print("Finished.")
        print('dr: ' + str(dr) + ' m')

        return data_frames, EAST_REGULAR, NORTH_REGULAR, UP_REGULAR, dr, data_frames_reg

    print("Finished.")
    print('dr: ' + str(dr) + ' m')

    return data_frames, EAST_REGULAR, NORTH_REGULAR, UP_REGULAR, dr

def initialize_airglow_filter(ntaps, Tlo, Thi, times, t_irregular = 3, raise_irregular = False):

    sample_periods = np.array([(times[i+1]-times[i]).total_seconds()/60.0 for i in range(1, len(times)-1)]) # minutes

    Tst            = np.median(sample_periods)

    print("median sampling rate: " + str(Tst) + " min")

    if raise_irregular is True:
        if any((sample_periods - Tst) > t_irregular):
            print(sample_periods - Tst)
            raise Exception('Sampling rate too irregular.')

    fs        = 1./(Tst*60.)
    cutofflo  = (1./(Thi*60))/fs
    cutoffhi  = (1./(Tlo*60))/fs

    print("lo %.2f; hi %.2f" % (cutofflo, cutoffhi))

    #b           = signal.firwin(ntaps, [cutofflo, cutoffhi], pass_zero=False, nyq = 0.5)
    b           = signal.firwin(ntaps, [cutofflo], pass_zero=False, nyq = 0.5)
    w,h         = signal.freqz(b,[1])

    # Frequency Response
    h_dB    = 20*np.log10(np.abs(h)); frequencies = w/(2*np.pi)*fs*1000
    h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))

#    plt.subplot(211)
#    plt.plot(frequencies, 20*np.log10(abs(h)))
#    locs, labels = plt.xticks()
#    ticks_min = 1/(locs[1:]/1000)/60.
#
#    for i,f in enumerate(ticks_min):
#         ticks_min[i] = '%.2f' % ticks_min[i]
#
#    labels[1:] = ticks_min; labels[0] = '$\infty$'; plt.xticks(locs, labels)
#    plt.xlim(0, frequencies[-1])
#    plt.ylabel('Magnitude (db)')
#    plt.title(r'Frequency response')
#    plt.ylim(-20, 20)
#    plt.grid()
#
#    plt.subplot(212)
#    plt.plot(frequencies, h_Phase)
#    plt.xticks(locs, labels)
#    plt.xlim(0, frequencies[-1])
#    plt.ylabel('Phase (radians)')
#    plt.xlabel(r'Period (min)')
#    plt.title(r'Phase response')
#    plt.grid(); plt.tight_layout()
#    plt.savefig('AGFilterResponse.png')
#    plt.close()

    return b

def conv_singleframe(args):

    global IM3D
    b, k, K, ntaps = args
    imtmp = np.zeros(np.shape(IM3D)[0:2])

    args = [(b, i, K, ntaps) for i in range(K)]

    for bi in range(len(b)):

        kk = k + bi - (ntaps-1)//2

        if kk < 0:
            kk = 0
        elif kk >= K:
            kk = K-1

        imtmp += b[bi] * IM3D[:,:,kk]

    return imtmp

IM3D = None

def filter_airglow(IM3Darg, b, ntaps):

    global IM3D
    IM3D  = IM3Darg

    IM3Dfilt  = np.zeros(np.shape(IM3D))
    K = np.shape(IM3D)[2]

    print("Temporally convolving...")
    pool = mp.Pool(processes = NUM_PROCESSORS)

    args = [(b, i, K, ntaps) for i in range(K)]

    conv_singleframe(args[0])

    IM3Dfilt_LIST = pool.map(conv_singleframe, args)
    pool.close()
    pool.join()

    for k in range(K):
        IM3Dfilt[:,:,k] = IM3Dfilt_LIST[k]

    print("finished.")

    return IM3Dfilt
