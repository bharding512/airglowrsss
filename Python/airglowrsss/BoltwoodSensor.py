
import numpy as np
from datetime import datetime, timedelta
#import matplotlib.pyplot as plt
import pytz
from pytz import timezone

def ReadTempLog_newformat(file, tz):
    """
    Function dns, sky_temp, amb_temp = ReadTempLog_newformat(file)
    ----------------------------------------------------
    Parses Boltwood Cloud Sensors log files

    INPUT:
        file : the name of the file to be parsed

    OUTPUT:
        dns : array of datetime object of the parsed file
        sky_temp : array of sky temperatures
        amb_temp : array of ambient temperatures
        
    History:
        8/20/12 : Written by Timothy Duly (duly2@illinois.edu)
        11/17/12 : Added timezone localization (jmakela@illinois.edu)
    """

    # Get the local timezone
    local = pytz.timezone(tz)

    try:
        fid = open(file,'r')
    except IOError as e:
        print file, 'does not exist'
        return [],[],[]
    data = []
    for line in fid:
        single_line = line.split()
        data.append(single_line)
    fid.close()

    N = len(data) # number of lines in file

    dns=[]
    sky_temp=[]
    amb_temp=[]

    for k in range(N):
        if len(data[k]) > 2 and data[k][2] is 'M' and data[k][3] == "~D": # healthy line
            ''' parse time '''
            date = data[k][0].split('-')
            year = int(date[0])
            month = int(date[1])
            day = int(date[2])

            # Depending on locale on the remote machine, decimals may be commas
            time = data[k][1].replace(',','.').split(':')
            hour = int(time[0])
            minute = int(time[1])
            second = int(float(time[2]))

            total_seconds = hour*3600. + minute*60. + second

            dn = local.localize(datetime(year,month,day) + timedelta(seconds=total_seconds))
            dns.append(dn)

            ''' parse temp '''
            try:
                sky_temp.append(float(data[k][10]))
            except:
                sky_temp.append(float('nan'))

            try:
                amb_temp.append(float(data[k][11]))
            except:
                amb_temp.append(float('nan'))

    # Finally, make numpy arrays:
    dns = np.array(dns)
    sky_temp = np.array(sky_temp)
    amb_temp = np.array(amb_temp)

    sky_temp[sky_temp == -998] = 0 # replace wet values with 0 (cloudy)
    sky_temp[sky_temp == -999] = float('nan') # replace bad values with nan

    # TODO: get rid of other bogus values... but which ones are bogus?

    return dns, sky_temp, amb_temp


def ReadRawTempLog(file, tz):
    """
        Function dns, sky_temp, amb_temp = ReadRawTempLog(file)
        ----------------------------------------------------
        Parses Boltwood Cloud Sensors Rawlog files

        INPUT:
        file : the name of the file to be parsed

        OUTPUT:
        dns : array of datetime object of the parsed file
        sky_temp : array of sky temperatures
        amb_temp : array of ambient temperatures

        History:
        8/28/12 : Written by Timothy Duly (duly2@illinois.edu)
                11/26/12: Added tz support
        12/13/12: updated to include Old format
        """

    # Get the local timezone
    local = pytz.timezone(tz)

    fid = open(file,'r')
    data = []
    for line in fid:
        single_line = line.split()
        data.append(single_line)
    fid.close()

    N = len(data) # number of lines in file

    dns=[]
    sky_temp=[]
    amb_temp=[]

    #print "N=",N
    for k in range(N):
        if len(data[k]) > 2: # healthy line
            #print data[k]
            if len(data[k]) == 11:
                # lines that have length=11 correspond to the old format"
                i_sky_temp = 3
                i_amb_temp = 4
            else:
                # new format
                i_sky_temp = 4
                i_amb_temp = 5
            ''' parse time '''
            date = data[k][0].split('-')
            year = int(date[0])
            month = int(date[1])
            day = int(date[2])

            time = data[k][1].split(':')
            hour = int(time[0])
            minute = int(time[1])
            second = int(float(time[2]))

            total_seconds = hour*3600. + minute*60. + second

            dn = local.localize(datetime(year,month,day) + timedelta(seconds=total_seconds))
            print 'here'
            dns.append(dn)

            ''' parse temp '''
            try:
                sky_temp.append(float(data[k][i_sky_temp]))
            except:
                sky_temp.append(float('nan'))

            try:
                amb_temp.append(float(data[k][i_amb_temp]))
            except:
                amb_temp.append(float('nan'))

    # Finally, make numpy arrays:
    dns = np.array(dns)
    sky_temp = np.array(sky_temp)
    amb_temp = np.array(amb_temp)

    sky_temp[sky_temp == -998] = 0 # replace wet values with 0 (cloudy)
    sky_temp[sky_temp == -999] = float('nan') # Replace bad with nan (unknown)

    return dns, sky_temp, amb_temp
    
def ReadTempLog_oldformat(file,tz):
    """
    Function dns, sky_temp, amb_temp = ReadTempLog_oldformat(file)
    ----------------------------------------------------
    Parses Boltwood Cloud Sensors log files

    INPUT:
        file : the name of the file to be parsed

    OUTPUT:
        dns : array of datetime object of the parsed file
        sky_temp : array of sky temperatures
        amb_temp : array of ambient temperatures
        
    History:
        11/29/12 : Written by Timothy Duly (duly2@illinois.edu)
    """

    # Get the local timezone
    local = pytz.timezone(tz)
    
    try:
        fid = open(file,'r')
    except IOError as e:
        print file, 'does not exist'
        return [],[],[]

    data = []
    for line in fid:
        single_line = line.split()
        data.append(single_line)
    fid.close()

    N = len(data) # number of lines in file

    dns=[]
    sky_temp=[]
    amb_temp=[]

    for k in range(N):
        line = data[k]
        if line[0] != 'Date' and 'illegal' not in line: # healthy line
            year, month, day = line[0].split('-')
            year = int(year)
            month = int(month)
            day = int(day)

            hour, minute, second = line[1].split(':')
            hour = int(hour)
            minute = int(minute)
            second = int(second)

            dns.append( local.localize(datetime(year,month,day,hour,minute,second)) )
            
            mess = line[10]
            sky = float(line[3])
            amb = float(line[4])
            sky_temp.append( sky )
            amb_temp.append( amb )
            
    # Finally, make numpy arrays:
    dns = np.array(dns)
    sky_temp = np.array(sky_temp)
    amb_temp = np.array(amb_temp)
    
    sky_temp[sky_temp == -998] = float('nan') # replace bad values with NaNs

    return dns, sky_temp, amb_temp

def DetermineIfOldFormat(file):
    fid = open(file,'r')
    old = 0
    for line in fid:
        # this messy string is only in the old format:
        if "FvrSkyTemAmbBlkHeaFlgC1sC1tC1aC1bC2sC2tC2aC2bC3sC3tC3aC3bTetAmtRqtSmtDitVctSeqChk" in line:
            old = 1
            return old
    return old

def ReadTempLog(file,tz):
    try:
	with open(file) as f: pass
    except:
    	return ([],[],[])
	
    if DetermineIfOldFormat(file):
        #print "old format"
        dns, sky_temp, amb_temp = ReadTempLog_oldformat(file,tz)
    else:
        #print "new format"
        dns, sky_temp, amb_temp = ReadTempLog_newformat(file,tz)
    return dns, sky_temp, amb_temp


def BoltwoodReduce(file,dn):
    """
    writes a day temperature file from a big temperature file
    BoltwoodReduce(file,dn)

    currently only old format is supported (maybe?)
    
    12/20/12 -- Timothy Duly (duly2@illinois.edu)
    """
    fid = open(file,'r')
    out_file_name = "%s%02d%02d_dailytemp.txt" % (dn.year, dn.month, dn.day)
    fid_out = open(out_file_name,'w')
        
    for line in fid:
        data = line.split()
        if data[0] not in 'Date':
            date = data[0].split('-')
            year = int(date[0])
            month = int(date[1])
            day = int(date[2])
            dn_line = datetime(year, month, day)
            if dn == dn_line:
                fid_out.write(line)
    print "created daily temp log: %s" % (out_file_name)
    fid.close()
    fid_out.close()
    
    
if __name__ == '__main__':
    #dns, sky_temp, amb_temp = ReadRawTempLog("/Users/duly/data/FPIData/temps/Cloud_raw.txt")
    
    #file1 = "/Users/duly/data/FPIData/temps/clarity_log.txt"
    #dns, sky_temp, amb_temp = ReadTempLog(file1,'US/Eastern')

    #file2 = "/Users/duly/data/FPIData/temps/Cloud_UAO_20120811.txt"
    #dns, sky_temp, amb_temp = ReadTempLog(file2,'US/Eastern')

    file1 = "/data/FPIData/temps/PARCloud.txt"
    file2 = "/data/FPIData/temps/EKUCloud.txt"

    dns1, sky_temp1, amb_temp1 = ReadRawTempLog(file1,'US/Eastern')
    dns2, sky_temp2, amb_temp2 = ReadRawTempLog(file2,'US/Eastern')

    #data = ReadRawTempLog(file1,'US/Eastern')






