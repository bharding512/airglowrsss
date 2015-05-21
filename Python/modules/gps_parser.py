


def parse_tec(my_file):
    '''
      data = parse_tec(my_file)

    Parses a TEC data file
    
    Input
    -----
      my_file - a string for the *.11_TEC file

    Output
    ------
      data - a dictionary with the following keys:
               'header' - contains header information
                          for XYZ and lat, lon, alt
                          of reciever
                '[1-32]' - a string representing
                           the satellite PRN.
                           These contain a dictionary:
                             'dns' - a list of times
                             'tec' - a list of tec values
    History
    -------
      1/10/14 - initially written by Timothy Duly 
                (duly2@illinois.edu)


    '''
    import numpy as np
    from datetime import datetime, timedelta 
    import pytz

    fid = open(my_file,'r')

    # initialize output data structure
    data = {}
    for kk in range(32):
        #data[str(kk+1)] = {
        #'dns':[],
        #'tec':[], 
        #}
        data[str(kk+1)] = {}
    data['header'] = {
            'xyz':[],
            'lla':[],
            }
    
    # we'll change this variable once we get past the main header
    past_header = False

    # main loop through the file
    for line in fid:
        if "APPROX POSITION XYZ" in line:
            # parse xyz:
            data['header']['xyz'] = \
                    np.array([float(aa) for aa in line.split()[:3]])
        if "POSITION LAT LON ALT" in line:
            # parse lat, lon, alt:
            data['header']['lla'] = \
                    np.array([float(aa) for aa in line.split()[:3]])
        if past_header:
            # now we're looking for a sub header describing
            # the next few lines of TEC data
            if "G" in line:
                line_info = line.split('G')
                year              = int(line_info[0].split()[0]) + 2000
                month             = int(line_info[0].split()[1])
                day               = int(line_info[0].split()[2])
                hour              = int(line_info[0].split()[3])
                minute            = int(line_info[0].split()[4])
                sec               = round(float(line_info[0].split()[5]))  ## TO HANDLE 29.999 entries
                number_of_entries = int(line_info[0].split()[7])
    
                dn = datetime(year, month, day, hour, minute,\
                        tzinfo=pytz.utc)\
                        + timedelta(seconds=sec)
    
                prns = [int(prn) for prn in line.split('G')[1:]]
    
                # loop through the PRNs in this set...
                for kk in range(number_of_entries):
                    tec_info = fid.next()
                    tec      = float(tec_info.split()[0])
                    tec_flag = int(tec_info.split()[1])
                    # ... and add the TEC & datenum (dn)
                    #     to the data structure
                    if tec_flag==0: # healthy measurement
                        #data[str(prns[kk])]['dns'].append(dn)
                        #data[str(prns[kk])]['tec'].append(tec)
                        data[str(prns[kk])][dn] = {'tec':tec}

        # we found the end of the main header:
        if "END OF HEADER" in line: past_header=True
    fid.close()
    return data

def tec_at_time(dn_search, data):
    '''
      tec, minutes_off = tec_at_time(dn_search, data)

     search through "data" for time "dn_search"

    Input
    -----
      dn_search - a datetime to search
                  * must be timezone aware *
      data - the dictionary created by 
             'parse_rinex' function

    Output
    ------
        tec - a list of TEC values that are
              close to dn_search.
              The index refers to PRN
        minutes_off - a list of minutes that
                      describe how far off
                      the TEC time values are 
                      from dn_search

    History
    -------
      1/10/14 - initially written by Timothy Duly 
                (duly2@illinois.edu)


    '''
    minutes_off = []
    tec = []

    for kk in range(32):
        prn = str(kk+1)
        diffs = []
        for dn in data[prn]['dns']:
            # calculate the number of seconds from dn_search
            # and add it to an array
            diff = np.abs((dn_search-dn).total_seconds()) / 60.
            diffs.append(diff)
        # make sure we have values in the diff array:
        if len(diffs)==0:
            tec.append(np.nan)
            minutes_off.append(np.nan)
        else:
            # from the list of time differences, find the
            # one that has the shortest time
            index_match = np.argmin(diffs)
            # add this to our output array
            tec.append(
                    data[prn]['tec'][index_match]
                    )
            minutes_off.append(
                    np.abs((dn_search-data[prn]['dns'][index_match])
                        .total_seconds())/60.
                    )
    return tec, minutes_off


def add_pos_stat(my_file, data):
    import numpy as np
    from datetime import datetime, timedelta 
    import pytz
    
    
    fid = open(my_file,'r')
    
    gps_epoch = datetime(1980,1,6, tzinfo=pytz.utc)
    
    
    #data = {}
    #for kk in range(32):
    #    data[str(kk+1)] = {}
    
    for line in fid:
        info = line.split(',')
        if info[0] == "$SAT":
            weeks = int(info[1])
            secs = float(info[2])
            satno = int(info[3])
            az = float(info[5])
            el = float(info[6])
    
            dn = gps_epoch + timedelta(weeks=weeks) +\
                    timedelta(seconds=secs) +\
                    timedelta(seconds=0)
            try:
                data[str(satno)][dn]['el'] = el
                data[str(satno)][dn]['az'] = az
            except KeyError:
                data[str(satno)][dn] = {}
                data[str(satno)][dn]['tec'] = np.nan
                data[str(satno)][dn]['el'] = el
                data[str(satno)][dn]['az'] = az

    fid.close()
    #return data
    return

def parse_rx_pos(my_file):
    import numpy as np
    from datetime import datetime, timedelta
    import pytz

    # input: a *.pos file

    data = {}

    fid = open(my_file,'r')

    for line in fid:
        if line[0] != '%':
            info = line.split()
            date = info[0]
            time = info[1]
            year  = int(date.split('/')[0])
            month = int(date.split('/')[1])
            day   = int(date.split('/')[2])
            hour  = int(time.split(':')[0])
            minute= int(time.split(':')[1])
            second= float(time.split(':')[2])

            dn = datetime(year,month,day,hour,minute,\
                    tzinfo=pytz.utc)+\
                    timedelta(seconds=second)
            lat = float(info[2])
            lon = float(info[3])
            alt = float(info[4])

            data[dn] = {'lla': \
                    np.array([lat,lon,alt])}

    fid.close()
    return data


def parse_sv_pos(my_file, data):
    import numpy as np
    from datetime import datetime, timedelta
    import pytz

    # input: a *.sp3 file

    fid = open(my_file,'r')

    for line in fid:
        if line[0] == "*":
            year  = int(line.split()[1])
            month = int(line.split()[2])
            day   = int(line.split()[3])
            hour  = int(line.split()[4])
            minute= int(line.split()[5])
            second= float(line.split()[6])
            dn = datetime(year,month,day,hour,minute,\
                    tzinfo=pytz.utc) + \
                            timedelta(seconds=second)
            for kk in range(32):
                prn = kk + 1
                subline = fid.next()
                x = float(subline.split()[1])
                y = float(subline.split()[2])
                z = float(subline.split()[3])
                data[str(prn)][dn] = {'xyz' : np.array([x,y,z])}

    fid.close()
    return


def create_interp_fs(sv_pos):
    fxs, fys, fzs = [], [], []
    for kk in range(32):
        prn = kk + 1
        dns = sv_pos[str(prn)].keys()
        dns.sort()
        dts = [] # coarse
        x, y, z = [], [], [] # coarse
        epoch = dns[0]
        for dn in dns:
            dts.append((dn-epoch).total_seconds())
            x.append(sv_pos[str(prn)][dn]['xyz'][0])
            y.append(sv_pos[str(prn)][dn]['xyz'][1])
            z.append(sv_pos[str(prn)][dn]['xyz'][2])
        fx = interp1d(dts, x, kind='cubic', bounds_error=False)
        fy = interp1d(dts, y, kind='cubic', bounds_error=False)
        fz = interp1d(dts, z, kind='cubic', bounds_error=False)
        fxs.append(fx)
        fys.append(fy)
        fzs.append(fz)
    return (fxs,fys,fzs), epoch

def sv_pos_at(dn,prn, fs, epoch):
    fxs, fys, fzs = fs
    dt = (dn-epoch).total_seconds()
    return np.array([fxs[prn-1](dt),fys[prn-1](dt),fzs[prn-1](dt)])

if __name__ == "__main__":
    from matplotlib.pylab import *
    import numpy as np
    from datetime import datetime, timedelta
    import glob
    import pytz
    from scipy.interpolate import UnivariateSpline
    from scipy.interpolate import interp1d

    #entry = "30990700"
    entry = "02130700"
    path = "/Users/duly/data/GEONET/test/"

    print "looking at ",path+entry

    tec = parse_tec(path+entry+".11_TEC")
    rx_pos = parse_rx_pos(path+entry+".pos")

    add_pos_stat(path+entry+".pos.stat", tec)

    dns = tec['2'].keys()
    dns.sort()
    dn = dns[500]
    print "tec['2'][dn=%s] = " % (dn)
    print tec['2'][dn]


    sv_pos = {}
    for kk in range(32):
        sv_pos[str(kk+1)] = {}
    parse_sv_pos(path+"igr16264.sp3", sv_pos)
    parse_sv_pos(path+"igr16265.sp3", sv_pos)

    # create the interpolation functions 
    fs, epoch = create_interp_fs(sv_pos)


    # parse a receiver position file:
    rx_pos = parse_rx_pos(path+entry+".pos")

    dns = rx_pos.keys()
    dns.sort()
    figure(1); clf()
    plot(dns,
            [rx_pos[dn]['lla'][2] for dn in dns],
            )
    draw(); show()
    grid()
            



    # test out interpolation:
    dns = sv_pos['2'].keys()
    dns.sort()
    dn_coarse = []
    x_coarse, y_coarse, z_coarse = [], [], []
    for dn in dns:
        dn_coarse.append(dn)
        x_coarse.append(sv_pos['2'][dn]['xyz'][0])
        y_coarse.append(sv_pos['2'][dn]['xyz'][1])
        z_coarse.append(sv_pos['2'][dn]['xyz'][2])
    N = len(dns)
    dts = ((dns[1]-dns[0]).total_seconds() ) / 10.
    dn = dns[0]
    dn_fine = []
    x_fine, y_fine, z_fine = [], [], []
    while dn < (dns[-1]+timedelta(hours=12)):
        dn = dn + timedelta(seconds = dts)
        dn_fine.append(dn)
        xyz = sv_pos_at(dn,2,fs, epoch)
        x_fine.append(xyz[0])
        y_fine.append(xyz[1])
        z_fine.append(xyz[2])

    figure(2); clf()
    plot(dn_coarse,y_coarse,'o')
    plot(dn_fine,  y_fine,'.g')
    draw(); show()
    

    # out of bounds should get nan's:
    dn = datetime(2011,3,11,1,3,4,34,tzinfo=pytz.utc)
    print sv_pos_at(dn,2, fs, epoch)
    dn = dn + timedelta(days=20)
    print sv_pos_at(dn,2, fs, epoch)
    


    


    




'''
    #print "data has the following keys:"
    #print data.keys()
    ##for key in data.keys(): print "  " + key

    #print "data['header'] = "
    #print data['header']
    #print "for example, data['4'] has %3i dn and TEC values" \
    #        % (len(data['4']['dns']))
    #print "first time of data['4'] = ", data['4']['dns'][0]
    #print "last  time of data['4'] = ", data['4']['dns'][-1]



    #print "\n\nDemo of searching for a specific time"
    #print "---------------------------------------"


    ## the time to find:
    #dn_search = datetime(2011,3,11,21,23,15,tzinfo=pytz.utc)

    #print " lets search for ",dn_search," in 'data'"

    #tec, minutes_off = tec_at_time(dn_search, data)

    #print "for ",dn_search
    #print " the tec for each satellite is "
    #print tec
    #print "However, these data are within "
    #print minutes_off
    #print " minutes from ", dn_search



    #for dn, tec in zip(data['2']['dns'],data['2']['tec']):
    #    try:
    #        s = pos['2'][dn]
    #    except KeyError:
    #        s = 'DN NOT FOUND'
    #    print dn, tec, s
'''
