from matplotlib import ticker,gridspec,dates
from cartopy import crs,feature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import numpy as np
from datetime import datetime,timedelta
#from apexpy import Apex
#from apexpy.apex import ApexHeightError
import os,shutil
import sys
from glob import glob
from pandas import date_range
from scipy import interpolate
import pytz,cv2
import FPI,ASI
import gc
import matplotlib
import fpiinfo

def changeHistogram(data1,nbins=12):
    ind=~np.isnan(data1)
    nbins=nbins
    imhist,bins = np.histogram(data1[ind],bins=2**nbins,normed=False)
    x=0.5*(bins[1:]+bins[:-1])
    y=imhist
    rangen=0
    rangei=0
    paso1=False
    paso2=False
    h=list(range(len(y)))
    for i in range(len(y)-2):
        h[i]=(y[i+1]-y[i])*1.0/(x[i+1]-x[i])
    h[len(h)-1]=0
    h[len(h)-2]=0
    for i in range(len(h)-1,0,-1):
        if abs(h[i])>=2. and paso1==False:
            rangen=int(x[i])
            paso1=True
        if y[i]>=1000 and paso1:
            paso2=True
        if y[i]<200 and paso2:
            rangei=int(x[i])
            break
    im=(data1.astype(float)-rangei)/(rangen-rangei)*(2**8-1)
    ind=im>2**8-1
    im[ind]=2**8-1
    ind=im<0
    im[ind]=0
    return im

def changeHistogram2(gdata,nbins=12):
    numberBins = 10000  # A good balance between time and space complexity, and well as precision
    contrast = 95.
    gflattenedImageData = np.array(gdata).flatten()
    imageHistogram, bins = np.histogram(gflattenedImageData, numberBins)
    imageHistogram = imageHistogram[1:]
    bins = bins[1:]
    cdf = np.cumsum(imageHistogram)
    # spliced to cut off non-image area
    cdf = cdf[:9996]
    max_cdf = max(cdf)
    maxIndex = np.argmin(abs(cdf - contrast/100 * max_cdf))
    minIndex = np.argmin(abs(cdf - (100 - contrast)/100 * max_cdf))
    vmax = float(bins[maxIndex])
    vmin = float(bins[minIndex])
    lowValueIndices = gflattenedImageData < vmin
    gflattenedImageData[lowValueIndices] = vmin
    highValueIndices = gflattenedImageData > vmax
    gflattenedImageData[highValueIndices] = vmax
    imageHistogram, bins = np.histogram(gflattenedImageData, numberBins)
    gdata = gflattenedImageData.reshape(gdata.shape)
    return gdata

def magnetic_gridlines(ax,xlocs=None,ylocs=None):
    if not xlocs:
        xlocs = np.arange(0,360,30)
    if not ylocs:
        ylocs = np.arange(-90,90,15)
    A = Apex(2017)
    gridlines = []
    gridlabels = []
    for mlat in ylocs:
        try:
            gdlat, gdlon = A.convert(mlat,np.linspace(0,360,100),'apex','geo',height=0)
            ax.plot(gdlon,gdlat,transform=crs.Geodetic(),linewidth=2, color='gray', alpha=0.5, linestyle='--')
            gridlines.append(np.array([gdlon,gdlat]))
            if mlat<0:
                gridlabels.append('{}magS'.format(abs(mlat)))
            else:
                gridlabels.append('{}magN'.format(mlat))
        except ApexHeightError as _:
            continue
    for mlon in xlocs:
        try:
            gdlat, gdlon = A.convert(np.linspace(-90,90,100),mlon,'apex','geo',height=0)
            ax.plot(gdlon,gdlat,transform=crs.Geodetic(),linewidth=2, color='gray', alpha=0.5, linestyle='--')
            gridlines.append(np.array([gdlon,gdlat]))
            if mlon<0:
                gridlabels.append('{}magW'.format(abs(mlon)))
            else:
                gridlabels.append('{}magE'.format(mlon))
        except ApexHeightError as _:
            continue

    return gridlines, gridlabels
def map_ticks(ax, lines, labels):
    xticks = []
    yticks = []
    xticklabels = []
    yticklabels = []

    for line, label in zip(lines,labels):
        tick_locations = tick_location(ax,line)
        if label.endswith('E') or label.endswith('W'):
            xticks.extend(tick_locations['bottom'])
            xticklabels.extend([label]*len(tick_locations['bottom']))

        if label.endswith('N') or label.endswith('S'):
            yticks.extend(tick_locations['left'])
            yticklabels.extend([label]*len(tick_locations['left']))

    return xticks, yticks, xticklabels, yticklabels
def tick_location(ax,line,edges_with_ticks=['left','bottom']):

    # convert line from geodetic coordinates to map projection coordinates
    line_map = ax.projection.transform_points(crs.Geodetic(),line[0],line[1]).T

    # parameters specific for finding ticks on each edge of plot
    edge_params = {'left':{'axis_i':0,'axis_d':1,'edge':ax.viewLim.x0,'bounds':[ax.viewLim.y0,ax.viewLim.y1]},
                   'right':{'axis_i':0,'axis_d':1,'edge':ax.viewLim.x1,'bounds':[ax.viewLim.y0,ax.viewLim.y1]},
                   'bottom':{'axis_i':1,'axis_d':0,'edge':ax.viewLim.y0,'bounds':[ax.viewLim.x0,ax.viewLim.x1]},
                   'top':{'axis_i':1,'axis_d':0,'edge':ax.viewLim.y1,'bounds':[ax.viewLim.x0,ax.viewLim.x1]}}

    # initialize empty dictionary to be returned wiht tick locations
    tick_locations = {}

    for e in edges_with_ticks:
        axis = edge_params[e]['axis_i']
        axisd = edge_params[e]['axis_d']
        edge = edge_params[e]['edge']
        bounds = edge_params[e]['bounds']
        tick_locations[e] = []

        # find indicies where the line crosses the edge
        line_map_shift = line_map[axis]-edge
        args = np.argwhere(line_map_shift[:-1]*line_map_shift[1:]<0).flatten()

        # for each crossing, interpolate to edge to find the tick location
        for a in args:
            l = line_map[0:2,a:a+2]
            l = l[:,l[axis].argsort()]

            tick = np.interp(edge,l[axis],l[axisd])
            # if tick is located within the plot, add it to list of tick locations
            if tick>bounds[0] and tick<bounds[1]:
                tick_locations[e].append(tick)

    return tick_locations

def readfile(ipath,data=False):
    import h5py
    fg = h5py.File(ipath,'r')
    gdata = fg['image']
    latitude=fg['image'].attrs['latitude']
    longitude=fg['image'].attrs['longitude']
    start_time=datetime(1970,1,1)+timedelta(seconds=int(fg['image'].attrs['start_time']))
    raw0=None
    if data:
        raw0=gdata[:,:]
    fg.close()
    return start_time,latitude,longitude,raw0

def getHorizontalDataFromNPZfile(f, sky_line_tag='X', chosenReference=None):
    '''
    Input:
        sky_line_tag: X or XG
        chosenReference: 'laser' or 'zenith'. See FPI.DopplerReference
    Output:
        Col  0: UT datetime object
        Col  1: Horizontal wind
        Col  2: error
        Col  3: Temperature
        Col  4: error
        Col  5: Intensity sky
        Col  6: error
        Col  7: Background sky
        Col  8: error
        Col  9: azimuthal angle
        Col 10: zenith angle
        Col 11: latitude for 250km
        Col 12: longitude for 250km
        Col 13: wind quality flag, 0(good), 1(caution), 2(bad)
        Col 14: temp quality flag, 0(good), 1(caution), 2(bad)
        Col 15: flag for corrected i.e. None, North or East
        Col 16: direction label i.e. minime91.nzk.cv_nzk_a3o_1
        Col 17: clouds
    '''
    def decodekeys(old):
        new={}
        for k,v in old.items():

            if isinstance(v,dict):
                val = decodekeys(v)
            elif isinstance(v,bytes):
                val = v.decode('utf-8')
            elif isinstance(v,np.ndarray):
                val = np.copy(v)
            elif isinstance(v,list):
                newlist=[]
                for iv in v:
                    if isinstance(iv,bytes):
                        newlist.append(iv.decode('utf-8'))
                    else:
                        newlist.append(iv)
                val = newlist
            else:
                val = v

            if isinstance(k,bytes):
                new[k.decode('utf-8')]=val
            else:
                new[k]=val

        return new

    # Read in the file
#    print( 'Accesing npz file %s...'%f, )
    try:
        npzfile = np.load(f,allow_pickle=True,encoding='bytes')
        FPI_Results = decodekeys( npzfile['FPI_Results'].item() )
        FPI_instrument = decodekeys( npzfile['instrument'].item() )
        FPI_site = decodekeys( npzfile['site'].item() )
        npzfile.close()

    except Exception as e:
        print(e)
        return None,None,None
#    print( 'OK')

    if 'Abbreviation' in FPI_instrument.keys():
        FPI_instrument_abbreviation=FPI_instrument['Abbreviation']
    else:
        FPI_instrument_abbreviation=f.split("_")[0]

    if chosenReference is None:
        (ref_Dop, e_ref_Dop) = FPI.DopplerReference(FPI_Results,reference=FPI_Results['reference'])
    else:
        (ref_Dop, e_ref_Dop) = FPI.DopplerReference(FPI_Results,reference=chosenReference)

    # Calculate the vertical wind and interpolate it
    ind = FPI.all_indices('Zenith',FPI_Results['direction'])
    w = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]) # LOS is away from instrument
    sigma_w = FPI_Results['sigma_fit_LOSwind'][ind]
    dt = []
    for x in FPI_Results['sky_times'][ind]:
        diff = (x - FPI_Results['sky_times'][0])
        dt.append(diff.seconds+diff.days*86400.)
    dt = np.array(dt)

    # Remove outliers
    ind = abs(w) < 200.

    if np.count_nonzero(ind) <= 1:
    # No good data, just use all ind
        ind = np.ones(len(w),dtype=np.bool)

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

    forbidden=['Laser','Unknown']
    l = [x for x in np.unique(FPI_Results['direction']) if x not in forbidden ]

    #choose height of emission
    if sky_line_tag=='X':
        height_emission=250
    elif sky_line_tag=='XG':
        height_emission=95

    DATA=None

    for x in l:

        ind = FPI.all_indices(x,FPI_Results['direction'])

        if x == 'Zenith':
            Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind])
            Doppler_Error = np.sqrt(FPI_Results['sigma_fit_LOSwind'][ind]**2)
        else:
            Doppler_Wind = (FPI_Results['LOSwind'][ind]-ref_Dop[ind]-w[ind]*np.cos(FPI_Results['ze'][ind]*np.pi/180.))/np.sin(FPI_Results['ze'][ind]*np.pi/180.)
            Doppler_Error = np.sqrt(FPI_Results['sigma_fit_LOSwind'][ind]**2+sigma_w[ind]**2)

        if ('Clouds' in FPI_Results.keys()) is False:
            FPI_Results['Clouds'] = None

        name="%s.%s.%s"%(FPI_instrument_abbreviation,FPI_site['Abbreviation'],x)

        az,ze=FPI_site['Directions'][x]['az'],FPI_site['Directions'][x]['ze']
        if ze < 0:
            az = az + 180
            ze = -ze

        if x == 'South':
            Doppler_Wind = -Doppler_Wind
        elif x == 'West':
            Doppler_Wind = -Doppler_Wind

        lat, lon = ASI.ConvertAzEl2LatLon(az, 90-ze, height_emission, FPI_site['Location'][0], FPI_site['Location'][1])

        lt2utobj=lambda lt:lt.astimezone(pytz.utc).replace(tzinfo=None)
        uts=list(map(lt2utobj,FPI_Results['sky_times'][ind]))
        data=np.array(list(zip(uts,
                        Doppler_Wind,Doppler_Error,
                        FPI_Results['T'][ind],FPI_Results['sigma_T'][ind],
                        FPI_Results['skyI'][ind],FPI_Results['sigma_skyI'][ind],
                        FPI_Results['skyB'][ind],FPI_Results['sigma_skyB'][ind],
                        [az,]*len(ind),
                        [ze,]*len(ind),
                        [lat,]*len(ind),
                        [lon,]*len(ind),
                        FPI_Results['wind_quality_flag'][ind],
                        FPI_Results['temp_quality_flag'][ind],
                        [name,]*len(ind),
                        FPI_Results['Clouds']['mean'][ind],
                        )))

        _i_=data[:,-1]<-900
        data[_i_,-1]=np.nan

        if DATA is None:
            DATA=data
        else:
            DATA=np.concatenate((DATA,data),axis=0)

    return DATA

def readHorizontalNPZdata(d0,dn,pathID,sky_line_tag='X',**fpikwargs):
    DATA=None
    path=__REPOSITORY_FPI__+"/%s*%s.npz" % (pathID, d0.strftime('%Y%m%d'))
#    print path
    if (sky_line_tag is not None) and (sky_line_tag !='X'):
        path=path[:-4]+"_%s"%(sky_line_tag.lower())+".npz"
    for dt in date_range(d0,dn).to_pydatetime():
        fmt=dt.strftime(path)
        paths=glob(fmt)
        for f in paths:

            try:
                data=getHorizontalDataFromNPZfile(f, sky_line_tag=sky_line_tag, **fpikwargs)
            except:
                continue

            if data is None:
                continue

            if DATA is None:
                DATA=data
            else:
                DATA=np.concatenate((DATA,data),axis=0)
    return DATA

def readHorizontalSiteNPZdata(d0,dn,site,sky_line_tag='X',**fpikwargs):
    DATA=None
    path=__REPOSITORY_FPI__+"/*%s*%s.npz" % (site, d0.strftime('%Y%m%d'))
#    print path
    if (sky_line_tag is not None) and (sky_line_tag !='X'):
        path=path[:-4]+"_%s"%(sky_line_tag.lower())+".npz"
    for dt in date_range(d0,dn).to_pydatetime():
        fmt=dt.strftime(path)
        paths=glob(fmt)
        for f in paths:

            try:
                data=getHorizontalDataFromNPZfile(f, sky_line_tag=sky_line_tag, **fpikwargs)
            except:
                continue

            if data is None:
                continue

            if DATA is None:
                DATA=data
            else:
                DATA=np.concatenate((DATA,data),axis=0)
    return DATA


#def getDASI02data(d0,dn,**fpikwargs):
#    DATA=readHorizontalNPZdata(d0,dn,pathID="minime11",**fpikwargs)
#    return DATA
#def getDASI01data(d0,dn,**fpikwargs):
#    DATA=readHorizontalNPZdata(d0,dn,pathID="minime12",**fpikwargs)
#    return DATA
#def getDASI03data(d0,dn,**fpikwargs):
#    DATA=readHorizontalNPZdata(d0,dn,pathID="minime10",**fpikwargs)
#    return DATA

def splitInCardinal(data):
    '''
    Returns dictionary for east,west,north,south,zenith data
    '''
    DATA={}
    for label in ['Zenith','East','West','North','South']:
        idx=np.array(list(map(lambda row:label.lower()in row[15].lower(),data)),dtype=np.bool)
        if np.count_nonzero(idx)==0:
            continue
        DATA[label]=np.copy(data[idx,:])
    return DATA


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



def generate_plot4(d0,dn,extent=[-130,-90,20,55],sky_line_tag="X",site_tags=['blo','cfs','cvo','low'],):

    from matplotlib import pyplot as plt

    data_crs = crs.PlateCarree()

    if sky_line_tag=='XG':
        fig = plt.figure(figsize=(8,4))
        spec=gridspec.GridSpec(ncols=2,nrows=2,figure=fig,
                                left=0.04,right=0.94,bottom=0.06,top=0.85,
                                wspace=0.05,hspace=0.25,
                                width_ratios=[1.5,2,])
    elif sky_line_tag=='X':
        fig = plt.figure(figsize=(8,4))
        spec=gridspec.GridSpec(ncols=2,nrows=2,figure=fig,
                                left=0.04,right=0.94,bottom=0.07,top=0.85,
                                wspace=0.05,hspace=0.25,
                                width_ratios=[1.5,2,])

    axes00 = fig.add_subplot(spec[:,0], projection=crs.Orthographic(np.nanmean(extent[:2]),np.nanmean(extent[2:])))
    axes01 = fig.add_subplot(spec[0,1])
    axes11 = fig.add_subplot(spec[1,1])

    #===============================================================================
    # setting up map
    #===============================================================================
    axes00.add_feature(feature.COASTLINE)
    axes00.add_feature(feature.STATES,alpha=0.2)
    axes00.set_extent(extent,crs=crs.PlateCarree())
    axes00.set_aspect('auto')

    xticks = []
    yticks = []
    xticklabels = []
    yticklabels = []
    xtickcolor = []
    ytickcolor = []
    if sky_line_tag=='X':
        mlatlines = range(-90,90,5)
        mlonlines = range(-180,180,10)
    else:
        mlatlines = range(-90,90,5)
        mlonlines = range(-180,180,10)
#    gridlines, gridlabels = magnetic_gridlines(axes00,xlocs=mlonlines,ylocs=mlatlines)
#    xt, yt, xtl, ytl = map_ticks(axes00,gridlines,gridlabels)
#    print( ytl)
#    print( yt)
#    if sky_line_tag=='X':
#        pass
#        #ytl=[ytl[2],ytl[4],ytl[6],ytl[7]]
#        #yt=[yt[2],yt[4],yt[6],yt[7]]
#    else:
#        ytl=[ytl[-3],ytl[-2],ytl[-1]]
#        yt=[yt[-3],yt[-2],yt[-1]]
#    xticks.extend(xt)
#    yticks.extend(yt)
#    xticklabels.extend(xtl)
#    yticklabels.extend(ytl)
#    xtickcolor.extend(['black']*len(xt))
#    ytickcolor.extend(['black']*len(yt))
#    # put gridline markers on plot
#    axes00.set_xticks(xticks)
#    axes00.set_xticklabels(xticklabels)
#    for t, c in zip(axes00.xaxis.get_ticklabels(),xtickcolor):
#        t.set_color(c)
#    axes00.set_yticks(yticks)
#    axes00.set_yticklabels(yticklabels,rotation='vertical',va='center')
#    for t, c in zip(axes00.yaxis.get_ticklabels(),ytickcolor):
#        t.set_color(c)





    #===============================================================================
    # adding stuff on top
    #===============================================================================
    if sky_line_tag=='XG':
        axes01.set_title("GREEN LINE\nEASTWARD WIND (M/S)",fontsize=10,pad=2)
    elif sky_line_tag=='X':
        axes01.set_title("RED LINE\nEASTWARD WIND (M/S)",fontsize=10,pad=2)
    axes11.set_title("NORTHWARD WIND (M/S)",fontsize=10,pad=2)

    err1=axes00.errorbar([],[],yerr=[],marker='o',color='black',markerfacecolor='white',markersize=4.5,linewidth=0,capsize=0,elinewidth=0.,markeredgewidth=0.5)
    err2=axes00.errorbar([],[],yerr=[],marker='^',color='black',markerfacecolor='white',markersize=6,linewidth=0,capsize=0,elinewidth=0.,markeredgewidth=0.5)
    err3=axes00.errorbar([],[],yerr=[],marker='s',color='black',markerfacecolor='white',markersize=6,linewidth=0,capsize=0,elinewidth=0.,markeredgewidth=0.5)
    handlesA=[err1,err2,err3]
    labelsA=["LOW, AZ","BLO, UT","CVO, OR"]
    axes00.legend(handlesA,labelsA,ncol=3,frameon=False, bbox_to_anchor=(0.5,0.97), loc='lower center',prop={'size':11},columnspacing=1,handletextpad=0.)


    dt2hr=lambda _dt:_dt.hour+_dt.minute/60.+_dt.second/3600.
    row2hr2=lambda row:dt2hr(row[0])-24 if row[0].hour>15 else dt2hr(row[0])

    #===============================================================================
    # adding FPI red
    #===============================================================================

    colors={'East':'tab:blue','West':'tab:orange','North':'tab:blue','South':'tab:orange'}
    marker={'minime11':'o','minime12':'^','minime10':'s'}

    red_winds={}

#    for getfn,instr in [(getDASI01data,'minime11'),(getDASI02data,'minime12'),(getDASI03data,'minime10')]:
    for s in site_tags:
        instr = fpiinfo.get_instr_at(s,dn)
        if instr == []:
            continue
        else:
            instr = instr[0]

        fpikwargs={'sky_line_tag':sky_line_tag}
#        data=getfn(d0,dn,**fpikwargs)
        data=readHorizontalSiteNPZdata(d0,dn,s,**fpikwargs)

        if data is None:
            continue

        data_los=splitInCardinal(data)
        for key,item in iter(data_los.items()):
            item=np.insert(item, 1, list(map(row2hr2,item)), axis=1)
            data_los[key]=item

        red_winds[instr]={k:np.copy(i) for k,i in iter(data_los.items())}

        for key,item in iter(data_los.items()):
            if key=='Zenith':
                continue
            elif key in ['East','West']:
                ax=axes01
            elif key in ['North','South']:
                ax=axes11
            color=colors[key]
            idx=item[:,14]==2
            ax.plot(item[idx,0],item[idx,2],marker=marker[instr],color=color,markersize=4,alpha=0.8,linewidth=0,markeredgewidth=0.5,markerfacecolor='None',markeredgecolor=color)
            idx=item[:,14]==1
            ax.plot(item[idx,0],item[idx,2],marker=marker[instr],color=color,markersize=4,alpha=0.8,linewidth=0,markeredgewidth=0.5,markerfacecolor='None',markeredgecolor=color)
            idx=item[:,14]==0
            ax.plot(item[idx,0],item[idx,2],marker=marker[instr],color=color,markersize=4,linewidth=0,markeredgewidth=0.5,markerfacecolor='None',markeredgecolor=color)
            ax.plot(item[:,0],item[:,2],linewidth=0.5,color=color,marker='None')

            lat,lon=np.nanmedian(item[:,12].astype(np.float)),np.nanmedian(item[:,13].astype(np.float))
            markersize=4.5 if instr=='minime11' else 5.5
            axes00.plot([lon,],[lat,],marker=marker[instr],color=color,markersize=markersize,markerfacecolor='None',markeredgewidth=1.5,zorder=100,transform=data_crs)

    for ax in [axes01,axes11]:
        ax.set_ylim(-200,200)
        ax.xaxis.set_major_formatter(dates.DateFormatter("%H"))
        ax.xaxis.set_minor_locator(dates.HourLocator(byhour=np.arange(0,24,1)))
        ax.xaxis.set_major_locator(dates.HourLocator(byhour=np.arange(0,24,2)))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(50))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
        ax.axhline(y=0,color='black',linewidth=0.5)
        ax.tick_params(axis='both',which='major',size=6,direction='in',right=True,top=True,labelright=True,labelleft=False)
        ax.tick_params(axis='both',which='minor',size=3,direction='in',right=True,top=True,labelright=True,labelleft=False)


    xlims=np.array([ax.get_xlim() for ax in [axes01,axes11,]])
    xlims=[np.nanmin(xlims[:,0]),np.nanmax(xlims[:,1])]
    for ax in [axes01,axes11,]:
        ax.set_xlim(*xlims)

    return fig,red_winds





def add_ASI2(fig,chunk,green_winds=None,sky_line_tag='XG',site_tags_ASI=['cfs','cvo','blo','low'],site_tags_FPI=['cvo','blo','low']):

    fpitags=site_tags_FPI

    green_fns={}
    if green_winds is not None:
        for tag in fpitags:
            if tag=='blo':
                inst='minime12'
            elif tag=='low':
                inst='minime11'
            elif tag=='cvo':
                inst='minime10'
            if inst not in green_winds.keys():
                continue
            pp={}
            for key,item in iter(green_winds[inst].items()):
                idx=item[:,14]<2
                if np.count_nonzero(idx)==0:
                    continue
                xx,yy=dates.date2num(item[idx,0]),item[idx,2].astype(np.float)
                pp[key]=lambda x,_yy=yy,_xx=xx:np.interp(x, _xx, _yy, left=np.nan, right=np.nan, period=None)
            green_fns[tag]=pp

    #===========================================================================
    # show ASI data
    #===========================================================================
    if sky_line_tag=='X':
        tags=site_tags_ASI
        tagline="redline"
    elif sky_line_tag=='XG':
        tags=site_tags_ASI
        tagline="greenline"

    DTS=None
    gbg=[]
    for tag in tags:
        if tag not in chunk.keys():
            continue
        path=chunk[tag]

        print( "Accesing hdf5 file %s..."%path,)
        dt,lat0,lon0,currentdata=readfile(path,data=True)
        print( 'OK')
        DTS=dt

#        daily_dir=os.path.dirname(os.path.dirname(path))

        #Calibration path
        calibration_path=__DAILY_DIR__+'/calibration_%s_%s_*.npz'%(tag.lower(),tagline,)
        nbins=12#for my way of histogramming!
        #print calibration_path

        calibration_path=glob(calibration_path)

        if len(calibration_path)>0:
            calibration_path=calibration_path[0]
            print( 'Loading calibration file ',calibration_path)
            obj=np.load(calibration_path,allow_pickle=True,encoding='bytes')
            calibration_params = {key:obj[key] for key in obj.keys()}
            obj.close()

        if tagline=='redline':
            cobermax=1000
            ht=250
        else:
            cobermax=300
            ht=95

        #getting nan mask
        cleared=changeHistogram2(currentdata,nbins=nbins)
        obj=ASIimage(path,calibration_params=calibration_params)
        obj.data=cleared
        obj.geometricalCalibration(cobermax=cobermax,H=ht,latlon0=[lat0,lon0])
        F=(np.nanmax(obj.data)-np.nanmin(obj.data))/256
        obj.data=(obj.data-np.nanmin(obj.data))/F
        #
        if tagline=='redline':
            _,mask=cv2.threshold(obj.data.astype('uint8'),0,256,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
            mask=mask.astype(np.float)
            idx=mask==0
            mask[idx]=np.nan
            idx=obj.geoelev<13
            mask[idx]=np.nan
        else:
            _,mask=cv2.threshold(obj.data.astype('uint8'),0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
            mask=mask.astype(np.float)
            idx=mask==0
            mask[idx]=np.nan
            idx=obj.geoelev<15.
            mask[idx]=np.nan
        #doing geometric calibration
        obj=ASIimage(path,calibration_params=calibration_params)
        obj.data=cleared
        obj.geometricalCalibration(cobermax=cobermax,H=ht,latlon0=[lat0,lon0])
        obj.data=obj.data*mask

        #draw ASI map
        ax=fig.axes[0]
        alpha=0.5 if len(chunk.keys())>1 else 1
        ax.pcolormesh(obj.geolon, obj.geolat, obj.data, transform=crs.PlateCarree(),cmap='gray',alpha=alpha)

        gbg.append([np.copy(obj.geolon), np.copy(obj.geolat), np.copy(obj.data)])

    #=======================================================================
    # show FPI arrows
    #=======================================================================

    quiverobj1=None

    for tag in fpitags:
        ax=fig.axes[0]

        if sky_line_tag=='XG':
            qcolor = 'yellowgreen'
            scale=25.
            xoff,yoff=8e4, -1.3e5
            width=0.015
        elif sky_line_tag=='X':
            qcolor='red'
            scale=50.
            xoff,yoff=1.5e5, -1.3e5
            width=0.01

        headwidth = 4
        headlength = 4
        headaxislength= headlength-1
        minshaft = 2
        sc = ax.bbox.width/fig.dpi/(10.*scale) # 1/10th of the plot width <==> 250 m/s, for quiver

        #draw FPI quiver
        if len(green_fns.keys())>0 and \
            tag in green_fns.keys() and \
            'East' in green_fns[tag].keys() and \
            'West' in green_fns[tag].keys() and \
            'North' in green_fns[tag].keys() and \
            'South' in green_fns[tag].keys():


            dtnum=dates.date2num(DTS)


            u=np.nanmean([green_fns[tag]['East'](dtnum),green_fns[tag]['West'](dtnum)])
            v=np.nanmean([green_fns[tag]['North'](dtnum),green_fns[tag]['South'](dtnum)])

            u=sc*u
            v=sc*v

            if tag=='blo':
                inst='minime12'
            elif tag=='low':
                inst='minime11'
            elif tag=='cvo':
                inst='minime10'

            lat0=np.nanmedian(green_winds[inst]['Zenith'][:,12].astype(np.float))
            lon0=np.nanmedian(green_winds[inst]['Zenith'][:,13].astype(np.float))

            obj = ax.quiver(np.array([lon0]),np.array([lat0]),np.array([u]),np.array([v]),
                            angles='uv', scale_units='inches', scale=1,width=width,
                            pivot='tail', headwidth=headwidth, headlength=headlength,
                            minshaft=minshaft, headaxislength=headaxislength, color=qcolor, transform=crs.PlateCarree(),zorder=1000)

            if quiverobj1 is None:
                x0,x1,y0,y1=ax.get_extent()
                quiverobj1=ax.quiverkey(obj, x0+xoff, y1+yoff, sc*scale, r'$%i\,\frac{m}{s}$'%scale,labelsep=0,color="black", coordinates='data',labelpos='N', transform=crs.PlateCarree(),)


    xlims_old=[dates.num2date(xl).replace(tzinfo=None) for xl in fig.axes[2].get_xlim()]

    if not (DTS>=xlims_old[0] and DTS<=xlims_old[1]):
        xlims=sorted(xlims_old+[DTS])
        for ax in fig.axes[1:]:
            ax.set_xlim(xlims[0],xlims[-1])

    lines=[]
    for ax in fig.axes[1:]:
        l=ax.axvline(x=DTS,color='red',linewidth=0.5)
        lines.append(l)

    title="%s"%( DTS.strftime("%d %b %Y, %H:%M UT") )
    fig.suptitle(title)


def add_DT_ASI4(fig,chunk,green_winds=None,sky_line_tag='XG',daily_dir=None,site_tags_ASI=['cfs','cvo','blo','low'],site_tags_FPI=['cvo','blo','low']):
    '''
    Function to add DT images
    '''
    fpitags=site_tags_FPI

    green_fns={}
    if green_winds is not None:
        for tag in fpitags:
            if tag=='blo':
                inst='minime12'
            elif tag=='low':
                inst='minime11'
            elif tag=='cvo':
                inst='minime10'
            if inst not in green_winds.keys():
                continue
            pp={}
            for key,item in iter(green_winds[inst].items()):
                idx=item[:,14]<2
                if np.count_nonzero(idx)==0:
                    continue
                xx,yy=dates.date2num(item[idx,0]),item[idx,2].astype(np.float)
                pp[key]=lambda x,_yy=yy,_xx=xx:np.interp(x, _xx, _yy, left=np.nan, right=np.nan, period=None)
            green_fns[tag]=pp

    #===========================================================================
    # show ASI data
    #===========================================================================
    if sky_line_tag=='X':
        tags=site_tags_ASI
        tagline="redline"
    elif sky_line_tag=='XG':
        tags=site_tags_ASI
        tagline="greenline"

    DTS=None
    gbg=[]
    for tag in tags:
        if tag not in chunk.keys():
            continue
        path0,path1=chunk[tag]

#        print( "Accesing hdf5 file %s..."%path0,)
        prevdt,lat0,lon0,prevdata=readfile(path0,data=True)
#        print( 'OK')

#        print( "Accesing hdf5 file %s..."%path1,)
        postdt,lat0,lon0,postdata=readfile(path1,data=True)
#        print( 'OK')

        if prevdata.shape!=postdata.shape:
            continue

        DTS=prevdt+timedelta(seconds=(postdt-prevdt).total_seconds()/2.)

#        daily_dir=os.path.dirname(os.path.dirname(path0))

        #Calibration path
        calibration_path=__DAILY_DIR__+'/calibration_%s_%s_*.npz'%(tag, tagline,)
#        print calibration_path
        nbins=12#for my histogram!

        calibration_path=glob(calibration_path)

        if len(calibration_path)>0:
            calibration_path=calibration_path[0]
#            print( 'Loading calibration file ',calibration_path)
            obj=np.load(calibration_path,allow_pickle=True,encoding='bytes')
            calibration_params = {key:obj[key] for key in obj.keys()}
            obj.close()

        if tagline=='redline':
            cobermax=1000
            ht=250
        else:
            cobermax=300
            ht=95

        #getting nan mask
        cleared=changeHistogram(prevdata,nbins=nbins)
        obj=ASIimage(path0,calibration_params=calibration_params)
        obj.data=cleared
        obj.geometricalCalibration(cobermax=cobermax,H=ht,latlon0=[lat0,lon0])
        #
        if tagline == 'redline':
            _,mask=cv2.threshold(obj.data.astype('uint8'),0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
            mask=mask.astype(np.float)
            idx=mask==0
            mask[idx]=np.nan
        else:
            _,mask=cv2.threshold(obj.data.astype('uint8'),0,255,cv2.THRESH_BINARY)
            mask=mask.astype(np.float)
            idx=mask==0
            mask[idx]=np.nan
            idx=obj.geoelev<13.
            mask[idx]=np.nan
        #doing geometric calibration
        obj=ASIimage(path0,calibration_params=calibration_params)
        obj.data=postdata-prevdata
        obj.geometricalCalibration(cobermax=cobermax,H=ht,latlon0=[lat0,lon0])
        obj.data=obj.data*mask

        #apply 7x7 median filter
        from scipy import signal
        obj.data=signal.medfilt2d(obj.data, kernel_size=7)

        #draw ASI map
        ax=fig.axes[0]
        alpha=0.5 if len(chunk.keys())>1 else 1
        ax.pcolormesh(obj.geolon, obj.geolat, obj.data, transform=crs.PlateCarree(),cmap='gray',alpha=alpha)

        gbg.append([np.copy(obj.geolon), np.copy(obj.geolat), np.copy(obj.data)])

    #=======================================================================
    # show FPI arrows
    #=======================================================================

    quiverobj1=None

    for tag in fpitags:
        ax=fig.axes[0]

        if sky_line_tag=='XG':
            qcolor = 'yellowgreen'
            scale=25.
            xoff,yoff=1e5, 5e4
            width=0.015
        elif sky_line_tag=='X':
            qcolor='red'
            scale=50.
            xoff,yoff=1.5e5, -1.3e5
            width=0.01

        headwidth = 4
        headlength = 4
        headaxislength= headlength-1
        minshaft = 2
        sc = ax.bbox.width/fig.dpi/(10.*scale) # 1/10th of the plot width <==> 250 m/s, for quiver

        #draw FPI quiver
        if len(green_fns.keys())>0 and \
            tag in green_fns.keys() and \
            'East' in green_fns[tag].keys() and \
            'West' in green_fns[tag].keys() and \
            'North' in green_fns[tag].keys() and \
            'South' in green_fns[tag].keys():


            dtnum=dates.date2num(DTS)


            u=np.nanmean([green_fns[tag]['East'](dtnum),green_fns[tag]['West'](dtnum)])
            v=np.nanmean([green_fns[tag]['North'](dtnum),green_fns[tag]['South'](dtnum)])

            u=sc*u
            v=sc*v

            if tag=='blo':
                inst='minime12'
            elif tag=='low':
                inst='minime11'
            elif tag=='cvo':
                inst='minime10'

            lat0=np.nanmedian(green_winds[inst]['Zenith'][:,12].astype(np.float))
            lon0=np.nanmedian(green_winds[inst]['Zenith'][:,13].astype(np.float))

            obj = ax.quiver(np.array([lon0]),np.array([lat0]),np.array([u]),np.array([v]),
                            angles='uv', scale_units='inches', scale=1,width=width,
                            pivot='tail', headwidth=headwidth, headlength=headlength,
                            minshaft=minshaft, headaxislength=headaxislength, color=qcolor, transform=crs.PlateCarree(),zorder=1000)

            if quiverobj1 is None:
                x0,x1,y0,y1=ax.get_extent()
                quiverobj1=ax.quiverkey(obj, x0+xoff, y0+yoff, sc*scale, r'$%i\,\frac{m}{s}$'%scale,labelsep=0,color="black", coordinates='data',labelpos='N', transform=crs.PlateCarree(),)

    if len(fpitags) > 0:
        xlims_old=[dates.num2date(xl).replace(tzinfo=None) for xl in fig.axes[2].get_xlim()]

        if not (DTS>=xlims_old[0] and DTS<=xlims_old[1]):
            xlims=sorted(xlims_old+[DTS])
            for ax in fig.axes[1:]:
                ax.set_xlim(xlims[0],xlims[-1])

        lines=[]
        for ax in fig.axes[1:]:
            l=ax.axvline(x=DTS,color='red',linewidth=0.5)
            lines.append(l)

    title="%s"%( DTS.strftime("%d %b %Y, %H:%M UT") )
    fig.suptitle(title)

def set_FPI_REPOSITORY(fpi_repo):
    __REPOSITORY_FPI__ = fpi_repo

def set_ASI_REPOSITORY(asi_repo):
    __REPOSITORY_ASI__ = asi_repo

def set_OUTFOLDER(outf):
    __OUTFOLDER__ = outf

def set_DAILY_DIR(daily_d):
    __DAILY_DIR__ = daily_d
