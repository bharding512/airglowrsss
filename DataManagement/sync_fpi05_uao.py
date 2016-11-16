#!/usr/bin/python2.7
'''
Script to send files back to remote2 and remove already copied data to free space
NOTE: this is for nfi only

History: 05 Aug 2014 - initial script written based on Tom Gehrels sh script; Daniel J. Fisher (dfisher2@illinois.edu)
'''

# Import required modules
import os
import sys
import datetime as dt
from glob import glob

print "You will need to enter the scp password..."

# Filestuff
#TODO: puth this in an array so it works for all FPI sites
DATA = '/cygdrive/c/FPI_Data/'
AIRGLOW = 'airglow@remote2.ece.illinois.edu:/rdata/airglow/fpi/minime05/uao/'
FILENAME = DATA + 'foldersizes.txt'
r = {}
l = {}

# SCP txt file (from list*.sh on remote2)
zelda = os.system('scp ' + AIRGLOW + '/foldersizes.txt ' + FILENAME)

# Parse txt file
if zelda == 0:
    with open(FILENAME, 'r') as f:
        for line in f:
            a,b,c = line.strip().split()
            r[dt.datetime(int(a),int(b[:2]),int(b[2:]))] = [a+b,c]
    print "remote file parsed"
    
    # Grab local files
    days = glob(DATA + '*/*')
    for d in days:
        try:
            t = dt.datetime.strptime(d[21:25]+' '+d[25:27]+' '+d[27:29],"%Y %b %d")
            size = 0
            for f in glob(d+'/*.img'):
                size += os.stat(f).st_size
            l[t] = [d,size]
        except:
            print 'not data folder'
    print 'local files parsed'
        
    # Do comparison
    for k in l.keys():
        # if file exists on remote...
        if k in r.keys():
            # if local is larger than remote
            if int(l[k][1]) > int(r[k][1]):
                print 'Not completely sent'
                # send file
                flag = os.system('scp ' + l[k][0]+'/*.tif ' + AIRGLOW + '%04i/%03s/.'%(k.year,r[k][0]))
                '''
                #Should I remove it once sent???
                if flag == 0:
                    #os.remove(l[k][0])
                    os.system('rm -rf ' + l[k][0])
                '''
                
            # if local is smaller or equal to remote
            else:
                print 'File okay to remove'
                # delete file
                #os.remove(l[k][0])
                os.system('rm -rf ' + l[k][0])

        # if file not on remote
        else:
            print 'Send new data'
            # send file
            flag = os.system('scp ' + l[k][0]+'/*.tif ' + AIRGLOW + '%04i/%03i/.'%(k.year,k.timetuple().tm_yday))
            '''        
            #Should I remove it once sent???
            if flag == 0:
                #os.remove(l[k][0])
                os.system('rm -rf ' + l[k][0])
            '''

        
else:
    print 'oh dear, this is no good!'
    
