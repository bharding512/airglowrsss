#!/usr/bin/python
# Filename: plot_CASES_RT.py

import matplotlib as mpl
mpl.use('Agg')
import CASES as cases

# Load txinfo.log file
txfname = '/mnt/CASESData/Hawaii/streaming/txinfo.log'
txprn,txdn, el, az, txsystem = cases.load_txinfo(txfname)

# Load scint.log file
s4fname = '/mnt/CASESData/Hawaii/streaming/scint.log'
s4prn, s4dn, s4, s4system = cases.load_scint(s4fname)

# Create plots
elazfname = '/data/Hawaii/s4pngs/elaz.png'
cases.plot_elaz(txprn,txdn,el,az,txsystem,s4prn,s4dn,s4,s4system,elazfname,dt = 15, el_mask = 15)

s4fname = '/data/Hawaii/s4pngs/realtime.png'
cases.plot_s4summary(txprn,txdn,el,az,txsystem,s4prn,s4dn,s4,s4system,s4fname)
