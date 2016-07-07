# This module is an interface for obtaining site and instrument information
# for all of our ASIs. It also houses all this information.
#
# When information needs to be added or updated, the "sites" and
# "dates" global variables defined below should be changed. Also, the function
# "get_instr_info" may need to be changed if a parameter of an instrument
# depends on the date (e.g., if a ccd is changed or a new calibration is used)

# All users should access this information through the functions in this module.

import ephem
import subprocess
import os
import datetime

# Define the start and stop dates of the instruments at each site.
# (Dates are inclusive, meaning that the instrument is at the site
#  on both its start and stop dates)
# HISTORY:
#   23 Jul 2013 jjm: Initialized dates by looking at data on remote2
_dates = {}

_dates['casi01'] = {}
_dates['casi01']['hka'] = { 'start': datetime.datetime(2001,10,19),
                            'stop': None,}

_dates['cnfi01'] = {}
_dates['cnfi01']['hka'] = { 'start': datetime.datetime(2001,10,23),
                            'stop': None,}

_dates['picasso01'] = {}
_dates['picasso01']['nso'] = { 'start': datetime.datetime(2010,5,14),
                            'stop': datetime.datetime(2012,6,25),}

_dates['picasso02'] = {}
_dates['picasso02']['cto'] = { 'start': datetime.datetime(2010,8,10),
                            'stop': None,}

_dates['picasso03'] = {}
_dates['picasso03']['bon'] = { 'start': datetime.datetime(2009,3,22),
                            'stop': datetime.datetime(2011,7,29),}
_dates['picasso03']['nso'] = { 'start': datetime.datetime(2016,6,1),
                            'stop': None,}

_dates['picasso04'] = {}
_dates['picasso04']['sgt'] = { 'start': datetime.datetime(2008,8,8),
                            'stop': datetime.datetime(2012,7,12),}
_dates['picasso04']['mor'] = { 'start': datetime.datetime(2013,11,4),
                            'stop': None,}

_dates['picasso05'] = {}
_dates['picasso05']['caj'] = { 'start': datetime.datetime(2009,7,31),
                            'stop': None,}

_dates['picasso06'] = {}
_dates['picasso06']['tht'] = {'start': datetime.datetime(2014,1,29),
                            'stop': None,}

_dates['swenson01'] = {}
_dates['swenson01']['soc'] = { 'start': datetime.datetime(2010,5,12),
                            'stop': datetime.datetime(2012,6,25),}

# Define the site information.
# HISTORY:
#   23 Jul 2013 jjm - Initialized by using info from ASI processing code

_sites = {}

_sites['hka'] = {
        'Location':     (20.71, 203.74, 3000.),
        'Name':         'Haleakala, Hawaii',
        'Abbreviation': 'HKA',
        'Timezone':     'US/Hawaii',
        'sql_id':       3,
        'share':        False,
        'borders':      True,
    }

_sites['nso'] = {
        'Location':     (32.79,254.19, 2804.),
        'Name':         'National Solar Observatory, New Mexico',
        'Abbreviation': 'NSO',
        'Timezone':     'US/Mountain',
        'sql_id':       12,
        'share':        False,
        'borders':      True,
    }

_sites['cto'] = {
        'Location':     (-30.71, 289.19, 2207.),
        'Name':         'Cerro Tololo Observatory, Chile',
        'Abbreviation': 'CTO',
        'Timezone':     'America/Santiago',
        'sql_id':       5,
        'share':        False,
        'borders':      True,
    }

_sites['bon'] = {
        'Location':     (12.19,291.76,0.),
        'Name':         'Sunspots, Bonaire',
        'Abbreviation': 'BON',
        'Timezone':     'America/Kralendijk',
        'sql_id':       9,
        'borders':      True,
    }

_sites['sgt'] = {
        'Location':     (10.57,298.83,0.),
        'Name':         'Sangre Grande, Trinidad',
        'Abbreviation': 'SGT',
        'Timezone':     'America/Port_of_Spain',
        'sql_id':       8,
        'share':        False,
        'borders':      True,
    }

_sites['caj'] = {
        'Location':     (-6.87,321.44,0.),
        'Name':         'Cajazeiras, Brazil',
        'Abbreviation': 'CAJ',
        'Timezone':     'America/Recife',
        'sql_id':       11,
        'share':        False,
        'borders':      True,
    }

_sites['soc'] = {
        'Location':     (34.05,253.08,1396.),
        'Name':         'Socorro, New Mexico',
        'Abbreviation': 'SOC',
        'Timezone':     'US/Mountain',
        'sql_id':       7,
        'share':        False,
        'borders':      True
    }
    
_sites['mor'] = {
        'Location':     (31.206,352.134,2700.),
        'Name':         'Oukaimeden, Morocco',
        'Abbreviation': 'MOR',
        'Timezone':     'UTC',
        'sql_id':       24,
        'share':        True,
        'borders':      False,
    }

_sites['tht'] = {
        'Location':     (-17.5667,210.426,350),
        'Name':         'Tahiti',
        'Abbreviation': 'THT',
        'Timezone':     'Pacific/Tahiti',
        'sql_id':       25,
        'share':        False,
        'borders':      True,
    }

# Instrument definitions
_instruments = {}
#TODO: lat, lon and cal file are site dependant ANDt insturment dependent (filter type)  ??
# SOLUTION: NEVER TOUCH THESE AGAIN UNLESS TO ADD NEW INSTRUMENT.  CHANGES ARE MADE IN GET_INSTR_INFO BELOW
_instruments['casi01'] = {
        'horizon':          20.,
        'filters':          ['6300','5577','7774','BGND','5890','Na'],
        'filter_names':     ['630.0 nm','557.7 nm', '777.4 nm', 'Background', '589.0 nm', '589.0 nm'],
        'unwarp_ht':        [250., 97., 350., 250., 92., 92.],
        'sql_inst_id':      [21,20,22,24,23,23],
        't_lat':            [20.71, 20.71, 20.71, 20.71, 20.71,20.71],
        't_lon':            [203.74, 203.74, 203.74, 203.74, 203.74,203.74],
        'cal_file':         '/rdata/airglow/imaging/calibrations/CASIelaz_0482002.npz',
        'kernel_size':      3,
        'ignore_dark':      False,
        'ccd_temp_set':     -35,
    }

_instruments['cnfi01'] = {
        'horizon':          0.,
        'filters':          ['6300','5577','7774','BGND','NONE'],
        'filter_names':     ['630.0 nm','557.7 nm', '777.4 nm', 'Background', 'none'],
        'unwarp_ht':        [250., 97., 350., 250., 250.],
        'sql_inst_id':      [9,10,11,12,19],
        't_lat':            [8.0, 12.0, 8.0, 8.0, 8.0],
        't_lon':            [203.74, 203.74, 203.74, 203.74, 203.74],
        'cal_file':         '/rdata/airglow/imaging/calibrations/CNFIelaz_0022002.npz',
        'kernel_size':      3,
        'ignore_dark':      False,
        'ccd_temp_set':     -40,
    }

_instruments['picasso01'] = {
        'horizon':          15.,
        'filters':          ['5511','5577','6300','8655'],
        'filter_names':     ['Background','557.7 nm','630.0 nm','OH'],
        'unwarp_ht':        [97., 97., 250., 95.],
        'sql_inst_id':      [54,4,3,55],
        't_lat':            [32.788,32.788,32.788,32.788],
        't_lon':            [254.194,254.194,254.194,254.194],
        'cal_file':         '/rdata/airglow/imaging/calibrations/PICASSO1elaz_1342010.npz',
        'kernel_size':      3,
        'ignore_dark':      False,
        'ccd_temp_set':     -30,
    }

_instruments['picasso02'] = {
        'horizon':          0.,
        'filters':          ['6300','7774','6343'],
        'filter_names':     ['630.0 nm','777.4 nm','Background'],
        'unwarp_ht':        [250., 350., 250.],
        'sql_inst_id':      [27,28,29],
        't_lat':            [-20.,-19.,-20.],
        't_lon':            [289.194,289.194,289.194],
        'cal_file':         '/rdata/airglow/imaging/calibrations/PICASSO2elaz_2342006.npz',
        'kernel_size':      3,
        'ignore_dark':      False,
        'ccd_temp_set':     -20,
    }

_instruments['picasso03'] = {
        'horizon':          15.,
        'filters':          ['6300','7774','BGND'],
        'filter_names':     ['630.0 nm','777.4 nm','Background'],
        'unwarp_ht':        [250., 350., 250.],
        'sql_inst_id':      [34,33,35],
        't_lat':            [12.19,12.19,12.19],
        't_lon':            [291.76,291.76,291.76],
        'cal_file':         '/rdata/airglow/imaging/calibrations/PICASSO3_0852009_approximate.npz',
        'kernel_size':      3,
        'ignore_dark':      False,
        'ccd_temp_set':     -60,
    }

_instruments['picasso04'] = {
        'horizon':          15.,
        'filters':          ['6300','7774','6343'],
        'filter_names':     ['630.0 nm','777.4 nm','Background'],
        'unwarp_ht':        [250., 350., 250.],
        'sql_inst_id':      [36,38,37],
        't_lat':            [10.57,10.57,10.57],
        't_lon':            [298.83,298.83,298.83],
        'cal_file':         '/rdata/airglow/imaging/calibrations/PICASSO4elaz_2202008.npz',
        'kernel_size':      3,
        'ignore_dark':      False,
        'ccd_temp_set':     -30,  # UNSURE
    }

_instruments['picasso05'] = {
        'horizon':          20.,
        'filters':          ['6300','BGND','7774','OH'],
        'filter_names':     ['630.0 nm','Background','777.4 nm','OH'],
        'unwarp_ht':        [250., 250., 350., 95.],
        'sql_inst_id':      [30,31,32,40],
        't_lat':            [-6.8726,-6.8726,-6.8726,-6.8726],
        't_lon':            [321.4434,321.4434,321.4434,321.4434],
        'cal_file':         '/rdata/airglow/imaging/calibrations/PICASSO5elaz_2122009.npz',
        'kernel_size':      3,
        'ignore_dark':      False,
        'ccd_temp_set':     -30,
    }

_instruments['picasso06'] = { # TODO: check all of these
        'horizon':          20.,
        'filters':          ['6300','6343','7774'],
        'filter_names':     ['630.0 nm','Background','777.4 nm'],
        'unwarp_ht':        [250., 250., 350., 95.],
        'sql_inst_id':      [101, 102, 103],
        't_lat':            [-17.56, -17.56, -17.56],
        't_lon':            [210.43, 210.43, 210.43],
#        'cal_file':         '/rdata/airglow/imaging/calibrations/PICASSO06elaz_0302014.npz',
	'cal_file':	    '/rdata/airglow/imaging/calibrations/picasso06_tht_20140130.npz',
        'kernel_size':      5,
        'ignore_dark':      False,
        'ccd_temp_set':     -69,
    }

_instruments['swenson01'] = {
        'horizon':          15.,
        'filters':          ['5510','5577','6300','OH'],
        'filter_names':     ['Background','557.7 nm','630.0 nm','OH'],
        'unwarp_ht':        [97., 97., 250., 95.],
        'sql_inst_id':      [59,60,61,62],
        't_lat':            [34.05,34.05,34.05,34.05],
        't_lon':            [253.08,253.08,253.08,253.08],
        'cal_file':         '/rdata/airglow/imaging/calibrations/SWENSON01_1322010.npz',
        'kernel_size':      3,
        'ignore_dark':      False,
        'ccd_temp_set':     -30,
    }

################## Site information functions ####################

def get_site_info(site_name):
    '''
    Return a dictionary with information about a site.

    INPUT:
        site_name - abbreviated site name (e.g., 'uao')
    OUTPUT:
        site_info - dictionary with information about the site.
        (keys are 'Location', 'TimeZone', 'Name', etc.)
    '''
    
    try:
        site_info = _sites[site_name].copy()
    except KeyError:
        raise Exception('Site name ("%s") not recognized. Try one of %s.' % \
                         (site_name, str(_sites.keys())))         
    
    return _sites[site_name].copy()

def get_all_sites_info():
    '''
    Return a dictionary of all site dictionaries, keyed by the site name.

    OUTPUT:
        all_site_dicts - dictionary whose keys are site names (e.g., 'uao') and values
        are dictionaries with information about the site (see get_site_info).
        For example, all_site_dicts['uao'] is the same as get_site_info('uao')
    '''
    
    return _sites.copy()




################### Instrument information functions #######################

def get_site_of(instr_name, dn):
    '''
    Return the name of the site where the instrument with name
    instr_name was on the day with date number dn. Return None
    if this instrument was nowhere on this date.
    INPUTS:
        instr_name - str, e.g., 'minime05'
        dn - datetime.datetime
    OUTPUTS:
        site_name - str, site name on this day, or None
    ''' 
      
    try:
        startstop = _dates[instr_name]
    except KeyError:
        raise Exception('Instrument name "%s" not recognized. Try one of: %s' \
                        % ( instr_name, str(sorted(_instruments.keys())) ))
    site_name = None
    for sn in startstop:
        start = startstop[sn]['start']
        stop  = startstop[sn]['stop' ]
        if start <= dn and (stop is None or dn <= stop):
            # We have a match. Make sure that we haven't found > one matches.
            if site_name is not None:
                raise Exception('Two sites ("%s", "%s") found for instrument "%s" on "%s")' \
                                % (site_name, sn, instr_name, str(dn)))
            site_name = sn

    return site_name
        

def get_instr_at(site_name, dn):
    '''
    Return the names of the instruments at site site_name
    on the day with date number dn.
    INPUTS:
        site_name - str, e.g., 'uao'
        dn - datetime.datetime
    OUTPUTS:
        instr_names - list of strs, comprising instrument names
    '''

    instr_names = []
    for n in _dates: # loop over possible instrument names
        if get_site_of(n, dn) == site_name:
            # add it to the list
            instr_names.append(n)
    return instr_names
        

def get_instr_info(instr_name, dn):
    '''
    Return an instrument dictionary for instrument with name instr_name.
    This may depend on the date number dn.
    INPUTS:
        instr_name - str, e.g., 'minime01'
        dn - datetime.datetime
    OUTPUTS:
        instrument - dictionary specifying instrument parameters
    '''
    try:
        instrument = _instruments[instr_name].copy()
    except KeyError:
        raise Exception('Instrument name "%s" not recognized. Try one of: %s' \
                        % ( instr_name, str(sorted(_instruments.keys())) ))


    if instr_name == 'casi01':
        if dn < datetime.datetime(2005,1,1) + datetime.timedelta(days = 253):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CASIelaz_0482002.npz'
        elif dn < datetime.datetime(2006,1,1) + datetime.timedelta(days = 261):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CASIelaz_02542005.npz'
        elif dn < datetime.datetime(2007,1,1) + datetime.timedelta(days = 56):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CASI01elaz_2622006.npz'
        elif dn < datetime.datetime(2009,1,1) + datetime.timedelta(days = 267):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CASIelaz_2622006.npz'
        elif dn < datetime.datetime(2011,1,1) + datetime.timedelta(days = 69):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CASIelaz_2682009.npz'
        else:
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CASIelaz_0702011.npz'
    elif instr_name == 'cnfi01':
        if dn < datetime.datetime(2002,1,1) + datetime.timedelta(days = 65):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CNFIelaz_0022002.npz'
        elif dn < datetime.datetime(2002,1,1) + datetime.timedelta(days = 77):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CNFIelaz_0662002.npz'
        elif dn < datetime.datetime(2005,1,1) + datetime.timedelta(days = 253):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CNFIelaz_0782002.npz'
        elif dn < datetime.datetime(2015,10,15):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CNFIelaz_2542005.npz'
        else:
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/CNFIelaz_3342015.npz'
    elif instr_name == 'picasso01':
        if dn < datetime.datetime(2010,5,14):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO1elaz_1342010.npz'
        else:
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO1elaz_1352010.npz'
    elif instr_name == 'picasso02':
        if dn < datetime.datetime(2007,1,20):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO2elaz_2342006.npz'
        elif dn < datetime.datetime(2015,4,26):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO2elaz_0202007.npz'
        elif dn < datetime.datetime(2015,7,27):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO02elaz_1182015.npz'
        else:
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO02elaz_2092015.npz'
    elif instr_name == 'picasso03': #BON
        if dn < datetime.datetime(2016,6,1):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO3_0852009_approximate.npz'
        else: #NSO
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO03elaz_1552016.npz'
            instrument['t_lat'] = [33.,33.,33.]
            instrument['t_lon'] = [254.,254.,254.]
    elif instr_name == 'picasso04': #SGT
        if dn < datetime.datetime(2008,1,1) + datetime.timedelta(days = 269):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO4elaz_2202008.npz'
        elif dn < datetime.datetime(2009,1,1) + datetime.timedelta(days = 61):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO4elaz_2702008.npz'
        elif dn < datetime.datetime(2009,1,1) + datetime.timedelta(days = 268):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO4elaz_0622009.npz'
        elif dn < datetime.datetime(2010,1,1) + datetime.timedelta(days = 17):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO4elaz_2692009.npz'
        elif dn < datetime.datetime(2011,1,1) + datetime.timedelta(days = 5):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO4elaz_0182010.npz'
        elif dn < datetime.datetime(2013,11,1):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO4elaz_0062011.npz'
        elif dn < datetime.datetime(2014,3,3): #MOR
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO4elaz_3122013.npz'
            instrument['t_lat'] = [31.2,31.2,31.2]
            instrument['t_lon'] = [352.1,352.1,352.1]
            instrument['kernel_size'] = 5
            instrument['ignore_dark'] = True
        elif dn < datetime.datetime(2014,7,1):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO04elaz_0622014.npz'
            instrument['t_lat'] = [31.2,31.2,31.2]
            instrument['t_lon'] = [352.1,352.1,352.1]
            instrument['kernel_size'] = 5
            instrument['ignore_dark'] = True
        elif dn < datetime.datetime(2015,3,27):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO04elaz_1822014.npz'
            instrument['t_lat'] = [31.2,31.2,31.2]
            instrument['t_lon'] = [352.1,352.1,352.1]
            instrument['kernel_size'] = 5
            instrument['ignore_dark'] = True
        else:
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO04elaz_0862015.npz'
            instrument['t_lat'] = [31.2,31.2,31.2]
            instrument['t_lon'] = [352.1,352.1,352.1]
            instrument['kernel_size'] = 3  
            instrument['ignore_dark'] = True
    elif instr_name == 'picasso05':
        if dn < datetime.datetime(2009,1,1) + datetime.timedelta(days = 264):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO5elaz_2122009.npz'
        elif dn < datetime.datetime(2012,1,1) + datetime.timedelta(days = 207):
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO5elaz_2652009.npz'
        else:
            instrument['cal_file'] = '/rdata/airglow/imaging/calibrations/PICASSO5elaz_2082012.npz'
    
    return instrument
    
    
def get_all_instr_names():
    '''
    Return a list of all of the instrument names (strings).
    '''
    return _instruments.keys()
