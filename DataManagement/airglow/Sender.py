#!/usr/bin/python
'''
Script to send files back to Airglow
  -s to add SITE

History: 27 Sep 2012 - initial script written; Daniel J. Fisher (dfisher2@illinois.edu)
         12 Feb 2014 - Updated to v3.0 - txtcheck; Daniel J. Fisher (dfisher2@illinois.edu)

'''

# Import required modules
import os
from glob import glob


def main():
    print("Hello World!")

    # Set up correct folder options
    folders = {
        '/cygdrive/c/Sending/',
        '/cygdrive/d/Sending/',
        '/cygdrive/f/Sending/',
        '/home/gps/Sending/',
        '/home/airglow/Sending/',
        '/home/scintmon/Sending/',
        '/data/Sending/',
        'C:/Sending/',
        'D:/Sending/',
        'F:/Sending/'
    }

    # destiny = 'tx@remote2.ece.illinois.edu:/rdata/airglow/rx/.'
    destiny = 'tx@remote2.ece.illinois.edu:/home/tx/rx/.'

    mfs = 70  # minimum file size to send

    send(folders, destiny, mfs)


def send(folders: set[str], destiny: str, mfs: int):
    # Go through all possible folders for data
    files = []
    for f in folders:
        f = f+"/" if f[-1] != "/" else f
        files = files + glob(f+'*.tar.gz*')
        files = files + glob(f+'*.txt')

    # Send all files one by one; remove if sent successfully
    for f in [fx for fx in files if os.stat(fx).st_size > mfs]:
        os.system('chmod 774 '+f)
        flag = os.system('scp ' + f + ' ' + destiny)
        if flag == 0:
            os.remove(f)
            print('Completed...')

    print('All Sending is now Complete!')
