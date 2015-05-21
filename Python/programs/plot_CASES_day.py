#!/usr/bin/python
# Filename: plot_CASES_day.py
#
# Run as plot_CASES_day -y YEAR -d DOY
#
# If no options given, run for yesterday's data

import matplotlib as mpl
mpl.use('Agg')
import CASES as cases
from optparse import OptionParser
import datetime
import os
import tarfile

usage = "usage: plot_CASES_day -y YEAR -d DOY"

parser = OptionParser(usage=usage)

# Get yesterday's date information as a default
dn = datetime.date.today()-datetime.timedelta(days = 1)
yyyy = dn.year
ddd = dn.timetuple().tm_yday

# Paths to use
datapath = '/mnt/CASESData/Hawaii/streaming/'
pngpath = '/data/Hawaii/s4pngs/'
archivepath = '/mnt/ScintMonData/CASES/Hawaii/'
log_fname = '/home/jmakela/Database/Hawaii_CASES_log.txt'
log_cmd = '/usr/bin/perl /home/jmakela/Database/load_HawaiiCASES_logfile.pl'

# Database info
site_id = 3;
s4_id = 63;

# Parse the command line
parser.add_option("-y", "--year", dest="year", help="year to be processed",
		metavar="YEAR", type="int",default=yyyy)
parser.add_option("-d", "--doy", dest="doy", help="day of year to be processed",
		metavar="DOY", type="int",default=ddd)

(options, args) = parser.parse_args()

# Create the filename of the data to be parsed
fname = '{:s}dataout_{:04d}_{:03d}.bin'.format(datapath,options.year,options.doy) 

# Run binflate
os.popen('/usr/local/sbin/binflate -i ' + fname)

# Load txinfo.log file
txfname = 'txinfo.log'
txprn,txdn, el, az, txsystem = cases.load_txinfo(txfname)

# Load scint.log file
s4fname = 'scint.log'
s4prn, s4dn, s4, s4system = cases.load_scint(s4fname)

# Create plots
dn = datetime.date(options.year,1,1)+datetime.timedelta(days=options.doy-1)
s4fname = '{:s}{:s}H.png'.format(pngpath,dn.strftime('%y%m%d'))
cases.plot_s4summary(txprn,txdn,el,az,txsystem,s4prn,s4dn,s4,s4system,s4fname)

# Write the logfile
fid = open(log_fname,'w')
fid.writelines('Site,Instrument,StartUTTime,StopUTTime,SummaryImage,MovieFile\n')
fid.writelines('{:d},{:d},{:s},{:s},{:s}H.png'.format(site_id,s4_id,
		s4dn[0].strftime('%Y-%m-%d %H:%M:%S'),
		s4dn[-1].strftime('%Y-%m-%d %H:%M:%S'),
		dn.strftime('%y%m%d')))
fid.close()

# Load the log file into the database
os.popen(log_cmd)

# Move the data
os.popen('mv {:s}dataout*_{:04d}_{:03d}.bin .'.format(datapath,options.year,options.doy
) )
tar = tarfile.open('{:s}.tgz'.format(dn.strftime('%y%m%d')), 'w:gz')

tar.add('dataout_{:04d}_{:03d}.bin'.format(options.year,options.doy))

if os.path.exists('dataoutiq_{:04d}_{:03d}.bin'.format(options.year,options.doy)):
   tar.add('dataoutiq_{:04d}_{:03d}.bin'.format(options.year,options.doy))

tar.close()

# Clean up files
#os.popen('tar czvf {:s}.tgz dataout*_{:04d}_{:03d}.bin'.format(dn.strftime('%y%m%d'), options.year, options.doy))
os.popen('rm dataout*_{:04d}_{:03d}.bin'.format(options.year,options.doy))
os.popen('mv {:s}.tgz {:s}{:04d}'.format(dn.strftime('%y%m%d'),archivepath,options.year))
os.popen('rm channel.log')
os.popen('rm iono.log')
os.popen('rm navsol.log')
os.popen('rm scint.log')
os.popen('rm txinfo.log')

