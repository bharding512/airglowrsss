
import sys
import datetime

from optparse import OptionParser


'''
Parse command line:
'''

usage = "usage: cleanup_X300 -y YEAR -d DOY -s SITE"
parser = OptionParser(usage=usage)

# Get yesterday's date information as a default
dn = datetime.date.today()-datetime.timedelta(days = 1)
yyyy = dn.year
ddd = dn.timetuple().tm_yday

# Parse the command line
parser.add_option("-y", "--year", dest="year", help="year to be processed",
		metavar="YEAR", type="int",default=yyyy)
parser.add_option("-d", "--doy", dest="doy", help="day of year to be processed",
		metavar="DOY", type="int",default=ddd)
parser.add_option("-s", "--site", dest="site", help="site to be processed",
		metavar="SITE", type="str",default=ddd)

(options, args) = parser.parse_args()
year = options.year
doy = options.doy
site = options.site

# Grab the datetime of our input
dn = datetime.datetime(year,1,1) + datetime.timedelta(days=doy-1)


'''
Done parsing command line
'''

read = {}
write = {}

read['DAN'] = "/home/fisher/airglowrsss/Python/"
write['DAN'] = "/home/fisher/Video/"

read['UAO'] = "C:/Scripts/Python/"
write['UAO'] = "C:/Sending/"

read['EKU'] = "C:/Scripts/Python/"
write['EKU'] = "C:/Sending/"

read['PAR'] = "C:/Scripts/Python/"
write['PAR'] = "C:/Sending/"

read['ANN'] = "C:/Scripts/Python/"
write['ANN'] = "C:/Sending/"

read['CAJ'] = "C:/Scripts/Python/"
write['CAJ'] = "C:/Sending/"
read['CAR'] = "C:/Scripts/Python/"
write['CAR'] = "F:/Sending/"

print "\n" + " system path is:"
sys.path = sys.path + [read[site]+'modules/']
print sys.path
print "\n"
import X300Sensor

read_path = read[site]+'programs/'
write_path = write[site]

print "\n"
import sys
print "the system path is:"
print sys.path
print "\n"
print datetime.datetime.now()
print dn
print read_path
print write_path
print site
X300Sensor.WriteLogFromRaw(dn,read_path,write_path,site)
print "Done with X300Sensor.WriteLogFromRaw()"

