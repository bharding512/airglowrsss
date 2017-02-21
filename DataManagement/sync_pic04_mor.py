#!/usr/local/bin/python
'''
Script to send files back to remote2 and remove already copied data to free space
NOTE: this is for pic04 only
TODO: MAKE one script for imagers and one for FPIs

History: 05 Aug 2014 - initial script written based on Tom Gehrels sh script; Daniel J. Fisher (dfisher2@illinois.edu)
'''

# Import required modules
import os
import sys
import datetime as dt
from glob import glob

# Why add this?  You need to type in the scp for every day to send data back
# This is because we don't want site computers to have passwordless login to our server
send = input("Do you want to send back files [1] or just delete duplicate data [0]?")
if send not in [0,1]:
    print 'Bad input'
    exit()
            
print "You will need to enter the scp password..."

# Filestuff
DATA = '/data/'
AIRGLOW = 'airglow@remote2.ece.illinois.edu:/rdata/airglow/imaging/picasso04/mor/'
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
            r[dt.datetime(int(a),1,1)+dt.timedelta(days=int(b)-1)] = [b,c]
    print "remote file parsed"
    
    # Grab local files
    days = glob(DATA + '*/*')
    for d in days:
        try:
            t = dt.datetime(int(d[-8:-4]),1,1)+dt.timedelta(days=int(d[-3:])-1)
            size = 0
            for f in glob(d+'/*.tif'):
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
            if int(l[k][1]) > int(r[k][1]) and send:
                print 'Not completely sent'
                # send file
                flag = os.system('scp ' + l[k][0]+'/*.tif ' + AIRGLOW + '%04i/%03s/.'%(k.year,r[k][0]))
                '''
                #Should I remove it once sent???
                if flag == 0:
                    #os.remove(l[k][0])
                    os.system('rm -rf ' + l[k][0])
                '''
            
            elif int(l[k][1]) > int(r[k][1]) and not(send):
                print "File needs to be sent later"

            # if local is smaller or equal to remote
            else:
                print 'File okay to remove'
                # delete file
                #os.remove(l[k][0])
                os.system('rm -rf ' + l[k][0])

        # if file not on remote
        elif send:
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
            print "Send new data later - %s"%k
        
else:
    print 'oh dear, this is no good!'
    
