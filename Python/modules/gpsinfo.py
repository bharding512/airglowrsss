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
#   05 Feb 2014 djf: Reformated asiinfo for gpsinfo
_dates = {}

_dates['casi01'] = {}
_dates['casi01']['hka'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}

_dates['scinda01'] = {}
_dates['scinda01']['sgt'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}

_dates['scinda02'] = {}
_dates['scinda02']['caj'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}
                            
_dates['scinda03'] = {}
_dates['scinda03']['bon'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}
                            
_dates['scintmonK'] = {}
_dates['scintmonK']['hka'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}
                            
_dates['scintmonL'] = {}
_dates['scintmonL']['cto'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}

_dates['scintmonM'] = {}
_dates['scintmonM']['cto'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}

_dates['scintmonO'] = {}
_dates['scintmonO']['bog'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}           
                            
_dates['scintmonS'] = {}
_dates['scintmonS']['car'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}
                            
_dates['scintmonT'] = {}
_dates['scintmonT']['car'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}

_dates['scintmonU'] = {}
_dates['scintmonU']['caj'] = { 'start': datetime.datetime(2001,2,16),
                            'stop': None,}


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
    }

_sites['cto'] = {
        'Location':     (-30.71, 289.19, 2207.),
        'Name':         'Cerro Tololo Observatory, Chile',
        'Abbreviation': 'CTO',
        'Timezone':     'America/Santiago',
        'sql_id':       5,
        'share':        False,
    }

_sites['bog'] = {
        'Location':     (4.64, 285.92, 2560.),
        'Name':         'Bogota, Colombia',
        'Abbreviation': 'BOG',
        'Timezone':     'America/Bogota',
        'sql_id':       18,
        'share':        False,
    }

_sites['bon'] = {
        'Location':     (12.19,291.76,0.),
        'Name':         'Sunspots, Bonaire',
        'Abbreviation': 'BON',
        'Timezone':     'America/Kralendijk',
        'sql_id':       9,
    }

_sites['sgt'] = {
        'Location':     (10.57,298.83,0.),
        'Name':         'Sangre Grande, Trinidad',
        'Abbreviation': 'SGT',
        'Timezone':     'America/Port_of_Spain',
        'sql_id':       8,
        'share':        False,
    }

_sites['caj'] = {
        'Location':     (-6.87,321.44,0.),
        'Name':         'Cajazeiras, Brazil',
        'Abbreviation': 'CAJ',
        'Timezone':     'America/Recife',
        'sql_id':       11,
        'share':        False,
    }

_sites['car'] = {
        'Location':     (-7.38, -36.52, 460.),
        'Name':         'Cariri, Brazil',
        'Abbreviation': 'CAR',
        'Timezone':     'America/Recife',
        'sql_id':       10,
        'share':        False,
    }

# Instrument definitions
_instruments = {}
#TODO: sql_inst_id is [s4,tec] should it be separated?

_instruments['cases01'] = {
        'horizon':          15.,
        'sql_inst_id':      63,
    }

_instruments['scinda01'] = {
        'horizon':          25.,
        'sql_inst_id':      [49,48],
    }
 
_instruments['scinda02'] = {
        'horizon':          25.,
        'sql_inst_id':      [47,46],
    }
 
_instruments['scinda03'] = {
        'horizon':          25.,
        'sql_inst_id':      [50,51],
    }
   
_instruments['scintmonK'] = {
        'horizon':          25.,
        'sql_inst_id':      45,
    }
    
_instruments['scintmonL'] = {
        'horizon':          25.,
        'sql_inst_id':      43,
    }

_instruments['scintmonM'] = {
        'horizon':          25.,
        'sql_inst_id':      44,
    }
 
_instruments['scintmonS'] = {
        'horizon':          25.,
        'sql_inst_id':      39,
    }
    
_instruments['scintmonT'] = {
        'horizon':          25.,
        'sql_inst_id':      41,
    }

_instruments['scintmonU'] = {
        'horizon':          25.,
        'sql_inst_id':      42,
    }
    
_instruments['scintmonO'] = {
        'horizon':          25.,
        'sql_inst_id':      71,
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

    
    return instrument
