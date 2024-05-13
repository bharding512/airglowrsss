
#import numpy as np 
from datetime import datetime, timedelta
import pytz
from pytz import timezone
#import matplotlib
#matplotlib.use('AGG')
#import matplotlib.pyplot as plt 

def getTempVals(d0,dn,tz,site,num=None):
    import numpy as np
    dts,temp1,temp2=None,None,None
    idt=d0

    while idt<=dn:
        path='/rdata/airglow/templogs/x300/%s/TempL%02s_%s_%s.txt'%(site,num,site,idt.strftime("%Y%m%d"))
        idts,itemp1,itemp2=ReadTempLog(path,tz)
        if len(idts)==0:
            idt=idt+timedelta(days=1)
            continue
        if dts is None:
            dts=idts
            temp1=itemp1
            temp2=itemp2
        else:
            dts=np.concatenate((dts,idts))
            temp1=np.concatenate((temp1,itemp1))
            if temp2 is not None:
                temp2=np.concatenate((temp2,itemp2))
        idt=idt+timedelta(days=1)

    return dts,temp1,temp2

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
		temp2 : temp2 (None for rpi and hass)

	History:
		8/23/12 : Written by Timothy Duly (duly2@illinois.edu)
		02/26/13: Added timezone localization (jmakela@illinois.edu)
		08/26/21: Added formats rpi and hass (lnav@illinois.edu)
	"""

	# Get the local timezone
	local = pytz.timezone(tz)

        try:
            fid = open(file,'r')
            data = fid.readlines()
            fid.close()
        except IOError as e:
            print file, 'does not exist'
            return [],[],[]

	if len(data) == 0:
	    # No data
	    return [],[],[]

	# Check if this is a thermoHID
        print "[ReadTempLog]: file, tz = ",file,tz
	if '#Output' in data[0]:
                print file, 'is a thermoHID file'
                return [], [], []

        import numpy as np
        #Parsing each line
        def row2list(row,localtz=local):
            words=row.split(',')
            try:#hass/rpi format: dt is UT/LT
                dn=datetime.strptime(words[0],"%Y/%m/%d %H:%M:%S")
                dn=dn.replace(tzinfo=pytz.utc)
            except:
                try:#old format: dt is in local time
                    dn=datetime.strptime(words[0],"%m/%d/%Y %H:%M:%S")
                    dn=local.localize(dn)
                except:
                    dn=None
            return [dn,]+[np.float(w) for w in words[1:]]
        my_list=np.array(map(row2list,data))

	# Next, sort
	my_list_sorted = np.array(sorted(my_list.tolist()))

        #Choosing format: dn,dome,inside
        if my_list_sorted.shape[1]==3:
            #old format
            return my_list_sorted[:,0],my_list_sorted[:,1],my_list_sorted[:,2]
        elif my_list_sorted.shape[1]==2:
            #hass format
            return my_list_sorted[:,0],my_list_sorted[:,1],None
        elif my_list_sorted.shape[1]==5:
            #rpi format
            return my_list_sorted[:,0],my_list_sorted[:,3]-14.43,None


def WriteLogFromRaw(dn,read_path,write_path,site,num=None,dtfmt="%m/%d/%Y"):
	"""
	Function WriteLogFromRaw(dn,read_path,write_path)
	----------------------------------------------------
        Creates a log file from the raw file (located in read_path/X300_temp_log.txt)
        By default, the name of this file is write_path/TempL_uao_YYYYMMDD.txt
        if num=="01", then name of this is write_path/TempL01_uao_YYYYMMDD.txt

        INPUT:
            dn : a datetime object for the day created
            read_path : a string indicating path folder where X300_temp_log.txt is located
            write_path : a string indicating path folder where output file will be written
            num (optional): a string indicating the number of sensing device. For example, 00 for rpi info 01 for zigbee's homeassistant
            dtfmt (optional): a string date format used to read dates within the raw files.

	History:
		8/24/12 : Written by Timothy Duly (duly2@illinois.edu)
		05/21 : Added number of sending device (lnav@illinois.edu)
	"""
	# open new file:
        if num is None:
            fname_new = "%s/TempL_%s_%04d%02d%02d.txt" % (write_path,site,dn.year,dn.month,dn.day)
        else:
            fname_new = "%s/TempL%s_%s_%04d%02d%02d.txt" % (write_path,num,site,dn.year,dn.month,dn.day)
        fid_new = open(fname_new,'w')

	# Load in raw data:
        fname = "%s/X300_temp_log.txt" % (read_path)
        fid = open(fname,'r')
        data = []
        for line in fid:
            try:
                dn_line=datetime.strptime(line.split()[0],dtfmt)
            except:
                continue
            if dn == dn_line: # the day we want
                fid_new.write(line)

        fid.close()
        fid_new.close()
	

if __name__ == '__main__':
	# ''' test out ReadTempLog '''
	#test_file = '/rdata/airglow/templogs/x300/uao/TempL_uao_20210523.txt'
	#test_file = '/rdata/airglow/templogs/x300/low/TempL00_low_20210823.txt'
        test_file = '/rdata/airglow/templogs/x300/low/TempL01_low_20210823.txt'
	dns, temp1, temp2 = ReadTempLog(test_file,tz='US/Mountain')
        print temp1
        print temp2
        print dns.shape
	# fig = plt.figure(1); plt.clf()
	# plt.plot(dns,temp1,'-',dns,temp2,'-')
	# plt.title('Temperature Data')
	# plt.xlabel('time')
	# plt.ylabel('temperature [C]')
	# plt.legend(('temp1','temp2'))
	# plt.grid('on')
	# fig.autofmt_xdate()
	# plt.draw(); plt.show()


	#''' test out WriteLogFromRaw '''
	#read_path = '/Users/duly/data/FPIData/temps/raw/'
	#read_path = '/Users/duly/data/FPIData/temps/raw/'
	#dn = datetime(2012,8,12)
	#out = WriteLogFromRaw(dn,read_path,write_path)


