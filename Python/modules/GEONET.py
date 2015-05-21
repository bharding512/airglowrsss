import pandas as pd
import MySQLdb as mdb
import pandas.io.sql as psql
import calendar
from datetime import datetime, timedelta
from matplotlib import dates

hostname = 'airglow.ece.illinois.edu'

def get_rx_data(rxID, tstart, tstop):
    '''
    Query the GEONET database for data from a receiver in a certain time span.
    INPUTS:
        rxID - integer, receiver ID
        tstart - datetime, start of time window
        tstop  - datetime, end of time window
    OUTPUTS:
        frame - an object containing the desired GPS data:
                rxID, time [POSIX], latitude [deg], longitude [deg], and altitude [m]
    '''
    
    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')
    
    # This is the UT POSIX time for the start/stop times of interest
    startUT = calendar.timegm(tstart.timetuple())
    stopUT =  calendar.timegm(tstop.timetuple())
    
    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, UT, lat, lon, alt ' + \
              'FROM rxposition WHERE rxID = %s AND UT >= %d AND UT <= %d' % (rxID, startUT, stopUT)
    frame = psql.frame_query(sql_cmd, con=con)
    con.close()
    
    return frame

def get_rx_tec_data(rxID, tstart, tstop):
    '''
    Query the GEONET database for tec data from a receiver in a certain time span.
    INPUTS:
        rxID - integer, receiver ID
        tstart - datetime, start of time window
        tstop  - datetime, end of time window
    OUTPUTS:
        frame - an object containing the desired GPS data:
                rxID, time [POSIX], svPRN, slant TEC [tecu], el [deg], and az [deg]
    '''

    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')

    # This is the UT POSIX time for the start/stop times of interest
    startUT = calendar.timegm(tstart.timetuple())
    stopUT =  calendar.timegm(tstop.timetuple())

    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, UT, svPRN, sTEC, el, az ' + \
              'FROM tecdata WHERE rxID = %s AND UT >= %d AND UT <= %d' % (rxID, startUT, stopUT)
    frame = psql.frame_query(sql_cmd, con=con)
    con.close()

    return frame

def get_rx_ftec_data(rxID, tstart, tstop):
    '''
    Query the GEONET database for filtered tec data from a receiver in a certain time span.
    INPUTS:
        rxID - integer, receiver ID
        tstart - datetime, start of time window
        tstop  - datetime, end of time window
    OUTPUTS:
        frame - an object containing the desired GPS data:
                rxID, time [POSIX], svPRN, filtered TEC [tecu], el [deg], and az [deg]
    '''

    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')

    # This is the UT POSIX time for the start/stop times of interest
    startUT = calendar.timegm(tstart.timetuple())
    stopUT =  calendar.timegm(tstop.timetuple())

    # Grab the slant TEC data (need the el/az from this)
    sql_cmd = 'SELECT rxID, UT, svPRN, sTEC, el, az ' + \
              'FROM tecdata WHERE rxID = %s AND UT >= %d AND UT <= %d' % (rxID, startUT, stopUT)
    frame1 = psql.frame_query(sql_cmd, con=con)
    
    # Grab the filtered TEC
    sql_cmd = 'SELECT rxID, UT, svPRN, fTEC ' +\
              'FROM ftecdata where rxID = %s and UT >= %d AND UT <= %d' % (rxID, startUT, stopUT)
    frame2 = psql.frame_query(sql_cmd, con=con)

    con.close()

    # Join the two frames
    frame = pd.merge(frame1, frame2, how='inner', on=['UT', 'rxID', 'svPRN']);

    return frame

def get_geonet_data(tstart, tstop):
    '''
    Query the GEONET database for data for all receivers in a certain time span.
    INPUTS:
        tstart - datetime, start of time window
        tstop  - datetime, end of time window
    OUTPUTS:
        frame - an object containing the desired GPS data:
                rxID, time [POSIX], latitude [deg], longitude [deg], and altitude [m]
    '''
    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')
    
    # This is the UT POSIX time for the start/stop times of interest
    startUT = calendar.timegm(tstart.timetuple())
    stopUT =  calendar.timegm(tstop.timetuple())
    
    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, UT, lat, lon, alt FROM rxposition WHERE UT >= %d AND UT <= %d' \
              % (startUT, stopUT)
    frame = psql.frame_query(sql_cmd, con=con)
    con.close()
    
    return frame

def get_geonet_tec_data(tstart, tstop):
    '''
    Query the GEONET TEC database for data for all receivers in a certain time span.
    INPUTS:
        tstart - datetime, start of time window
        tstop  - datetime, end of time window
    OUTPUTS:
        frame - an object containing the desired GPS data:
		rxID, time [POSIX], svPRN, slant TEC [tecu], el [deg], and az [deg]
    '''
    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')

    # This is the UT POSIX time for the start/stop times of interest
    startUT = calendar.timegm(tstart.timetuple())
    stopUT =  calendar.timegm(tstop.timetuple())

    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, UT, svPRN, sTEC, el, az FROM tecdata WHERE UT >= %d AND UT <= %d' \
              % (startUT, stopUT)
    frame = psql.frame_query(sql_cmd, con=con)
    con.close()

    return frame

def get_geonet_data_at_time(t):
    '''
    Query the GEONET database for data for all receivers at a certain time.
    INPUTS:
        t - datetime, desired sample time
    OUTPUTS:
        frame - an object containing the desired GPS data:
                rxID, time [POSIX], latitude [deg], longitude [deg], and altitude [m]
    '''
    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')
    
    # This is the UT POSIX time for the time of interest
    UT = calendar.timegm(t.timetuple())
    
    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, UT, lat, lon, alt FROM rxposition WHERE UT = %d' \
              % (UT)
    frame = psql.frame_query(sql_cmd, con=con)
    con.close()
    
    return frame

def get_geonet_tec_data_at_time(t):
    '''
    Query the GEONET TEC database for data for all receivers at a certain time.
    INPUTS:
        t - datetime, desired sample time
    OUTPUTS:
        frame - an object containing the desired GPS data:
		rxID, time [POSIX], svPRN, slant TEC [tecu], el [deg], and az [deg]
    '''
    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')

    # This is the UT POSIX time for the time of interest
    UT = calendar.timegm(t.timetuple())

    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, UT, svPRN, sTEC, el, az FROM tecdata WHERE UT = %d' \
              % (UT)
    frame = psql.frame_query(sql_cmd, con=con)
    con.close()

    return frame

def get_geonet_ftec_data_at_time(t):
    '''
    Query the GEONET filtered TEC database for data for all receivers at a certain time.
    INPUTS:
        t - datetime, desired sample time
    OUTPUTS:
        frame - an object containing the desired GPS data:
                rxID, time [POSIX], svPRN, slat TEC [tecu], filtered TEC [tecu], el [deg], and az [deg]
    '''
    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')

    # This is the UT POSIX time for the time of interest
    UT = calendar.timegm(t.timetuple())

    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, UT, svPRN, sTEC, el, az FROM tecdata WHERE UT = %d' \
              % (UT)
    frame1 = psql.frame_query(sql_cmd, con=con)

    sql_cmd = 'SELECT rxID, UT, svPRN, fTEC FROM ftecdata WHERE UT = %d' % (UT)
    frame2 = psql.frame_query(sql_cmd, con=con)

    # Join
    frame = pd.merge(frame1,frame2,how='inner',on=['UT','rxID','svPRN'])

    con.close()

    return frame

def get_geonet_ftec_data_sv(sv, tstart, tstop):
    '''
    Query the GEONET database for filtered tec data from all receivers in a certain time span
    for a given satellite.
    INPUTS:
        sv - integer, satellite PRN id
        tstart - datetime, start of time window
        tstop  - datetime, end of time window
    OUTPUTS:
        frame - an object containing the desired GPS data:
                rxID, time [POSIX], svPRN, slant TEC [tecu], ,filtered TEC [tecu], el [deg], and az [deg]
    '''

    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')

    # This is the UT POSIX time for the start/stop times of interest
    startUT = calendar.timegm(tstart.timetuple())
    stopUT =  calendar.timegm(tstop.timetuple())

    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, UT, svPRN, sTEC, el, az ' + \
              'FROM tecdata WHERE svPRN = %d AND UT >= %d AND UT <= %d' % (sv, startUT, stopUT)
    frame1 = psql.frame_query(sql_cmd, con=con)
    
    # Get the filtered TEC data
    sql_cmd = 'SELECT rxID, UT, svPRN, fTEC FROM ftecdata WHERE svPRN = %d AND UT >= %d and UT <= %d' %(sv, startUT, stopUT)
    frame2 = psql.frame_query(sql_cmd, con=con)

    # Join
    frame = pd.merge(frame1,frame2,how='inner',on=['UT','rxID','svPRN'])

    con.close()

    return frame

def get_rx_info(rxID):
    '''
    Query the GEONET info database for data for all receivers at a certain time.
    INPUTS:
	rxID - integer, receiver ID 
    OUTPUTS:
        frame - an object containing the desired GPS data:
		rxID, lat [deg], lon [deg], alt [m], ECEFx [m], ECEFy [m], ECEFz [m]
    '''
    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')

    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, lat, lon, alt, ECEFx, ECEFy, ECEFz from rxinfo WHERE rxID = %d' % (rxID)
    frame = psql.frame_query(sql_cmd, con=con)
    con.close()

    return frame

def get_all_rx_info():
    '''
    Query the GEONET info database for data for all receivers at a certain time.
    INPUTS:
    OUTPUTS:
        frame - an object containing the desired GPS data:
                rxID, lat [deg], lon [deg], alt [m], ECEFx [m], ECEFy [m], ECEFz [m]
    '''
    # Open a connection to the database
    con = mdb.connect(host=hostname,user='bduser',passwd='bdpass',db='gpsdatabase')

    # First find out if the entry is in there (i.e., we are just updating the png and avi file)
    sql_cmd = 'SELECT rxID, lat, lon, alt, ECEFx, ECEFy, ECEFz from rxinfo'
    frame = psql.frame_query(sql_cmd, con=con)
    con.close()

    return frame
