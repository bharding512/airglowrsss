#!/usr/bin/python
'''
Script to send files back to Airglow
  -s to add SITE

History: 27 Sep 2012 - initial script written; Daniel J. Fisher (dfisher2@illinois.edu)
         12 Feb 2014 - Updated to v3.0 - txtcheck; Daniel J. Fisher (dfisher2@illinois.edu)

'''

# Import required modules
import os
import sys
from glob import glob

# Set up correct folder options
folder = {}
folder[0] = '/cygdrive/c/Sending/'
folder[1] = '/cygdrive/d/Sending/'
folder[2] = '/cygdrive/f/Sending/'
folder[3] = '/home/gps/Sending/'
folder[4] = '/home/airglow/Sending/'
folder[5] = '/home/scintmon/Sending/'
folder[6] = '/data/Sending/'
folder[7] = 'C:/Sending/'
folder[8] = 'D:/Sending/'
folder[9] = 'F:/Sending/'
destiny = 'tx@remote2.ece.illinois.edu:/rdata/airglow/rx/.'
mfs = 70 #minimum file size to send

# Go through all possible folders for data
files = []
for f in range(len(folder)):
    files = files + glob(folder[f]+'*.tar.gz*')
    files = files + glob(folder[f]+'*.txt')

# Send all files one by one; remove if sent successfully
for f in [fx for fx in files if os.stat(fx).st_size > mfs]:
    os.system('chmod 774 '+f)
    flag = os.system('scp ' + f + ' ' + destiny)
    if flag == 0:
        os.remove(f)
        print 'Completed...'

print 'All Sending is now Complete!'

