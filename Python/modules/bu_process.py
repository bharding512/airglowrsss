import matplotlib
#matplotlib.use('Agg')
import numpy as np
from scipy.io import loadmat, savemat
import matplotlib.pyplot as plt
import TifImagePlugin
import datetime
import os
import array
import glob
import ASI

BASE_IMAGING_DIRECTORY = '/rdata/airglow/zama/asi'

def extract_timestamp(fn, date):
    """
    Returns a datetime corresponding to the BU image filename (the filename format is
    described in the documentation for get_datetime_paths)
    """

    N = len(fn)
    n1 = N - fn[::-1].find('.')
    timestamp = fn[n1-8:n1-2]
    
    hour = int(timestamp[0:2])
    minute = int(timestamp[2:4])
    second = int(timestamp[4:6])

    dt = datetime.datetime(date.year, date.month, date.day, hour, minute, second)

    return dt

def sort_key(e):
    return e[1]

sao_month_codes_2018 = {}
sao_month_codes_2018['Jan'] = 'January'
sao_month_codes_2018['Feb'] = 'February'
sao_month_codes_2018['Mar'] = 'March'
sao_month_codes_2018['Apr'] = 'April'
sao_month_codes_2018['May'] = 'May'
sao_month_codes_2018['Jun'] = 'June'
sao_month_codes_2018['Jul'] = 'July'
sao_month_codes_2018['Aug'] = 'August'
sao_month_codes_2018['Sep'] = 'September'
sao_month_codes_2018['Oct'] = 'October'
sao_month_codes_2018['Nov'] = 'Nov'
sao_month_codes_2018['Dec'] = 'Dec'

sao_month_codes_2019 = {}
sao_month_codes_2019['Jan'] = 'Jan'
sao_month_codes_2019['Feb'] = 'Feb'
sao_month_codes_2019['Mar'] = 'Mar'
sao_month_codes_2019['Apr'] = 'April'
sao_month_codes_2019['May'] = 'May'
sao_month_codes_2019['Jun'] = 'Jun'
sao_month_codes_2019['Jul'] = 'Jul'
sao_month_codes_2019['Aug'] = 'Aug'
sao_month_codes_2019['Sep'] = 'Sept'
sao_month_codes_2019['Oct'] = 'Oct'
sao_month_codes_2019['Nov'] = 'Nov'
sao_month_codes_2019['Dec'] = 'Dec'

def get_datetime_paths(date, station = 'mho'):
    """
    get_date_path will return the directory of bu imaging data stored within BASE_IMAGING_DIRECTORY
    for the specified datetime and for the specified station, if data exists for the specified date.
    If data is not available for the specified date, False is returned. The data should be stored
    within BASE_IMAGING_DIRECTORY with the following filestructure:

    BASE_IMAGING_DIRECTORY/sss/yyyy/mmmddyy/XHHMMSSF.DOY

    where sss is the 3-letter site abbreviation (sao, mho, etc), yyyy is the data year, mmmddyy
    is the month, day, year code (e.g., December 02, 2018 <--> Dec0218), X is the single-letter
    station code (S = SAO, M = MHO), HH is the hour, MM is the minute and SS is the second
    (e.g., 21:15:00 <--> 211500), and F is the image type (D, B, C).
    """

    # check and see if the directory corresponding to the passed date exists

    if station == 'mho':
        dir_path = os.path.join(BASE_IMAGING_DIRECTORY, '%s' % station, str(date.year),
        '%s%02d%s' % (date.strftime('%b'), date.day, str(date.year)[2:]))
    elif station == 'sao':
        mo = date.strftime('%b')

        if date.year == 2018:
            month = sao_month_codes_2018[mo]
        else:
            month = sao_month_codes_2019[mo]

        dir_path = os.path.join(BASE_IMAGING_DIRECTORY, '%s' % station, str(date.year), month,
        '%s%02d%s' % (mo, date.day, str(date.year)[2:]))

    if not os.path.isdir(dir_path):
        return [], [], [], [], [], []
    else:
        doy = date.timetuple().tm_yday
        Bfiles = glob.glob(os.path.join(dir_path, '*B.%03d' % doy))
        Cfiles = glob.glob(os.path.join(dir_path, '*C.%03d' % doy))
        Dfiles = glob.glob(os.path.join(dir_path, '*D.%03d' % doy))

        # extract the timestamps from each of the groups
        times_B = []
        times_C = []
        times_D = []

        for fn in Bfiles:
            times_B.append(extract_timestamp(fn, date))

        for fn in Cfiles:
            times_C.append(extract_timestamp(fn, date))

        for fn in Dfiles:
            times_D.append(extract_timestamp(fn, date))

        # sort the files (ascending order, wrt time)
        tB = sorted(list(tuple(zip(Bfiles, times_B))), key = sort_key)
        tC = sorted(list(tuple(zip(Cfiles, times_C))), key = sort_key)
        tD = sorted(list(tuple(zip(Dfiles, times_D))), key = sort_key)

        try:
                Bfiles_sorted, times_B_sorted = list(zip(*tB))
        except ValueError as e:
                print("No images of type 'B' found for %s" % date)
                Bfiles_sorted = []
                times_B_sorted = []

        try:
                Cfiles_sorted, times_C_sorted = list(zip(*tC))
        except ValueError as e:
                print("No images of type 'C' found for %s" % date)
                Cfiles_sorted = []
                times_C_sorted = []

        try:
                Dfiles_sorted, times_D_sorted = list(zip(*tD))
        except ValueError as e:
                print("No images of type 'D' found for %s" % date)
                Dfiles_sorted = []
                times_D_sorted = []

        return Bfiles_sorted, Cfiles_sorted, Dfiles_sorted, times_B_sorted, times_C_sorted, times_D_sorted

def load_bu(fn):
    """
    Script that loads in imaging data from the Boston University Imaging Science data repository:

    http://sirius.bu.edu/dataview/

    and returns the image as a 2D numpy array.
    """

    arr = array.array("h")

    with open(fn, "rb") as f:

        # skip header which is 128 bytes 
        f.seek(128,0)

        # append rest of file to arr structure defined; size of file is 512 x 512
        arr.fromfile(f, 512*512)

        # conversion to unsigned?
        tmp = np.array(arr)
        ind = np.where(tmp < 0)[0]
        tmp[ind] = tmp[ind] + 2**16

        # put in correct format
        img = np.reshape(tmp, (512,512))
        img = np.fliplr(img)

    return img

def bu_latlon_cal(fn_az_cal, fn_el_cal, latlon_imager, projection_height = 250, horizon = 15, make_plot = False):

    star_az = loadmat(fn_az_cal)['ccd_az']
    star_el = loadmat(fn_el_cal)['ccd_el']

    LATS, LONS = ASI.ConvertAzEl2LatLon(star_az.flatten(), star_el.flatten(), projection_height, latlon_imager[0],\
    latlon_imager[1], horizon = horizon)

    LATS = np.reshape(LATS, star_az.shape)
    LONS = np.reshape(LONS, star_el.shape)

    return LATS, LONS
