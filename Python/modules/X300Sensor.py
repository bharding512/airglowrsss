
#import numpy as np 
from datetime import datetime, timedelta
import pytz
from pytz import timezone
#import matplotlib
#matplotlib.use('AGG')
#import matplotlib.pyplot as plt 

def ReadTempLog(file, tz):
	"""
	Function dns, sky_temp, amb_temp = ReadTempLog(file)
	----------------------------------------------------
	Parses X300 log files

	INPUT:
		file : the name of the file to be parsed

	OUTPUT:
		dns : array of datetime object of the parsed file
		temp1 : temp1
		temp2 : temp2

	History:
		8/23/12 : Written by Timothy Duly (duly2@illinois.edu)
		02/26/13: Added timezone localization (jmakela@illinois.edu)
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

	if len(data) == 0:
	    # No data
       	    return [],[],[]

	# Check if this is a thermoHID
        print "[ReadTempLog]: file, tz = ",file,tz
	if data[0][0] == '#Output':
                print file, 'is a thermoHID file'
                return [], [], []

	N = len(data) # number of lines in file

	dns = []
	temp1 = []
	temp2 = []

	for k in range(N):
		if 'x' not in data[k][0]: # check to make sure we don't have a bogus line

			''' parse time ''' 
			date = data[k][0].split('/')
			month = int(date[0])
			day = int(date[1])
			year = int(date[2])

			time = data[k][1].split(',')[0].split(':')
			hour = int(time[0])
			minute = int(time[1])
			second = int(time[2])

			total_seconds = hour*3600. + minute*60. + second

			dn = local.localize(datetime(year,month,day) + timedelta(seconds=total_seconds))
			dns.append(dn)

			''' parse temp '''
			temp1.append(float(data[k][1].split(',')[1]))
			temp2.append(float(data[k][1].split(',')[2]))


	# let's sort the data by dns. 
	# First, create a tuple of data
	N = len(dns)

	my_list = []
	for k in range(N):
		my_list.append( (dns[k], temp1[k], temp2[k]) )
	
	# Next, sort
	my_list_sorted = sorted(my_list, key=lambda datum: datum[0])

	# Finally, place data back in list
	dns = []; temp1 = []; temp2=[]
	for k in range(N):
		dns.append(my_list_sorted[k][0])
		temp1.append(my_list_sorted[k][1])
		temp2.append(my_list_sorted[k][2])

	# make numpy arrays:
	import numpy as np
	dns = np.array(dns)
	temp1 = np.array(temp1)
	temp2 = np.array(temp2)

	return np.array(dns), np.array(temp1), np.array(temp2)

def WriteLogFromRaw(dn,read_path,write_path,site):
	"""
	Function WriteLogFromRaw(dn,read_path,write_path)
	----------------------------------------------------
	Creates a TempL_UAO_YYYYMMDD.txt file from the raw 
	file

	INPUT:
		dn : a datetime object for the day created

	OUTPUT:
		<none>, 
		but a file is created in write_path

	History:
		8/24/12 : Written by Timothy Duly (duly2@illinois.edu)
	"""
	# open new file:
	fname_new = "%s/TempL_%s_%04d%02d%02d.txt" % (write_path,site,dn.year,dn.month,dn.day)
	fid_new = open(fname_new,'w')

	# Load in raw data:
	fname = "%s/X300_temp_log.txt" % (read_path)
	fid = open(fname,'r')
	data = []
	for line in fid:
		(month, day, year)= line.split()[0].split("/")
		if 'x' not in month: # healthy line
			dn_line = datetime(int(year),int(month),int(day))
			if dn == dn_line: # the day we want
				fid_new.write(line)

	fid.close()
	fid_new.close()
	

if __name__ == '__main__':
	# ''' test out ReadTempLog '''
	# test_file = '/Users/duly/data/FPIData/temps/TempL_UAO_20120808.txt'
	# dns, temp1, temp2 = ReadTempLog(test_file)

	# fig = plt.figure(1); plt.clf()
	# plt.plot(dns,temp1,'-',dns,temp2,'-')
	# plt.title('Temperature Data')
	# plt.xlabel('time')
	# plt.ylabel('temperature [C]')
	# plt.legend(('temp1','temp2'))
	# plt.grid('on')
	# fig.autofmt_xdate()
	# plt.draw(); plt.show()


	''' test out WriteLogFromRaw '''
	read_path = '/Users/duly/data/FPIData/temps/raw/'
	write_path = '/Users/duly/data/FPIData/temps/'
	dn = datetime(2012,8,12)
	out = WriteLogFromRaw(dn,read_path,write_path)


