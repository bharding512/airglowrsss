#!/usr/bin/python

'''
This code is run at 20:05 CT to allow Cases01_HKA to process with Sorting scripts
It makes a copy of the data to use (which enables accurate tracking).
Note it relys on current server time.  This assumption is okay because the data should be realtime transmitted.

Creator: Dan Fisher (dfisher2@illinois.edu)
'''

import datetime as dt
from glob import glob
import os

# Create "True" filename
dn = dt.datetime.now()
doy = (dn-dt.datetime(dn.year,1,1)).days+1
truefn = '/rdata/airglow/rx/cas01_hka_{:04d}{:02d}{:02d}.tar.gz'.format(dn.year,dn.month,dn.day)

# Grab cases file from yesterday (in UT)
datafn = '/rdata/airglow/gps/cases01/hka/streaming/dataout_{:04d}_{:03d}.bin'.format(dn.year,doy)
files = glob(datafn)

# if found, copy file to RX folder
if files:
    os.system('cp '+datafn+' '+truefn)

