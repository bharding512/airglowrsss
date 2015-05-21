#!/usr/local/bin/python2.7
'''
Script to send files back to remote2 and remove already copied data to free space
NOTE: this is for scn only, pic use sync*.sh

History: 05 Aug 2014 - initial script written based on Tom Gehrels sh script; Daniel J. Fisher (dfisher2@illinois.edu)
'''

# Import required modules
import os
import sys
from glob import glob

print "You will need to enter the scp password..."

# Filestuff
DATA = '/home/scintmon/cascade-1.62/'
AIRGLOW = 'airglow@remote2.ece.illinois.edu:/rdata/airglow/gps/scintmonU/caj/'
FILENAME = DATA + 'foldersizes.txt'
letter = 'U'
century = '20'
r = {}
l = {}

# SCP txt file (from list*.sh on remote2)
zelda = os.system('scp ' + AIRGLOW + '/foldersizes.txt ' + FILENAME)

# Parse txt file
if zelda == 0:
    with open(FILENAME, 'r') as f:
        for line in f:
            key,val = line.strip().split()
            r[key] = val
    print "text file parsed"
    
    # Grab local files
    files = []
    files = files + glob(DATA + '*'+letter+'*.fsl')
    files = files + glob(DATA + '*'+letter+'*.nav')
    files = files + glob(DATA + '*'+letter+'*.obs')
    
    # Grab local sizes
    for f in files:
        s = os.stat(f).st_size
        l[f[len(DATA):]] = s,f
        
    # Do comparison
    for k in l.keys():
        # if file exists on remote...
        if k in r.keys():
            # if local is larger than remote
            if int(l[k][0]) > int(r[k]):
                print 'Not completely sent'
                # send file
                flag = os.system('scp ' + l[k][1] + ' ' + AIRGLOW + century+k[:2] + '/raw_data/.')
                
                #Should I remove it once sent???
                if flag == 0:
                    os.remove(l[k][1])
                
            # if local is smaller or equal to remote
            else:
                print 'File okay to remove'
                # delete file
                os.system('rm -f ' + l[k][1])

        # if file not on remote
        else:
            print 'Send new data'
            # send file
            flag = os.system('scp ' + l[k][1] + ' ' + AIRGLOW + century+k[:2] + '/raw_data/.')
            
            #Should I remove it once sent???
            if flag == 0:
                os.remove(l[k][1])
            
        
        
else:
    print 'oh dear, this is no good!'
