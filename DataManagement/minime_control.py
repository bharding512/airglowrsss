#!/usr/bin/python
'''
Script to run minime instruments:
    -r --run = runs instrument from start to finish
    -s --stop = stops and turns off instrument
    #-n --run-now = runs instrument NOW (power outage option)
    -f --focus = focus TODO
    -c --calibrate = sun calibration TODO
    -t --test_SS = test cardinal sky scanner positions TODO
    -h --help = help

History: 15 Jun 2017 - initial script written

Written by Daniel J Fisher (dfisher2@illinois.edu)
        by Brian J Harding (bhardin2@illinois.edu)
'''

# Import modules
import os
import sys
import glob as glob
from collections import defaultdict
import argparse
#import Andor


# Scheduling is done with at_night, which has run, do nothing, and stop options.


##### OPERATIONS #####

def run(now=False):
    '''
    Run typical night observations
    
    Input:
        now [False=default] - flag to force immediately taking data

    '''
    print 'in'
    #Schedule start/stop time
    scheduler(now=now)
    #run system
    master()
    return()



def focus():

    print 'TODO'

    return()

# To do a sun calibration
def cal():

    print "WARNING: cover CCD (need paper)"

    print 'TODO'

    return()


def test():

    print 'TODO'
    return()


##### OPERATION MODULES #####

# Turns on peripherals
def periperals_on(LASER=True,SS=True,CCD=True):

    if LASER:
        # Turn on Laser
    
        # Wait 30 minutes

    if SS:
        # Turn on Sky Scanner

    if CCD:
        # Turn on CCD

        # Begin cooling CCD

    return()


# Main data-taking script (dynamic_NATION)
def master():

    print 'todo master'

    return()


# Turns off peripherals
def peripherals_off(LASER=True,SS=True,CCD=True):

    if SS:
        # Turn off SkyScanner

    if LASER:
        # Turn off laser

    if CCD:
        # Cool CCD

        # Turn off CCD

    return()


# Zips and sends previous nights data
def tx_data()

    # Zip FPI

    # Zip Cloud

    # Zip Temps (if applicable)

    # Send over

    return()


if __name__=="__main__":
    '''
    Main module run on command line
    '''
    
    # Close program if already running
    pid = str(os.getpid())
    pidfile = "/tmp/minime_control.pid"
    if os.path.isfile(pidfile):
        print "Already running, will exit"
        sys.exit()
    else:
        file(pidfile, 'w').write(pid)
    
    # Parse command
    parser = argparse.ArgumentParser()
    #parser.add_option("-h","--help")
    parser.add_argument("-r","--run",action="store_true",\
            help="Automatically runs instrument")
    parser.add_argument("-s","--stop",action="store_true",\
            help="Automatically stops instrument")
    #parser.add_argument("-n","--run-now",action="store_true",\
    #        help="Immediately runs instrument (recover from powerloss)")
    parser.add_argument("-c","--calibrate",action="store_true",\
            help="Enters sun calibration mode")
    parser.add_argument("-f","--focus",action="store_true",\
            help="Does autofocus")
    parser.add_argument("-t","--test",action="store_true",\
            help="Tests SkyScanner positions")
    args = parser.parse_args()

    # Branch to function to run
    if args.run:
        print args.run
        print 'running'
        peripherals_on()
        master()

    if args.stop:
        print 'stopping'
        peripherals_off()
        tx_data()

    if args.run-now:
        print 'run NOW'

    elif args.calibrate:
        print 'todo'
        # Turn on SS
        peripherals_on(CCD=False,LASER=False)  
        # SUN CALIBRATE FUNCTION
        #sunseeker()
        # Turn off SS
        peripherals_off(CCD=False,LASER=False)

    elif args.focus:
        print 'todo'
        peripherals_on()
        # focus()
        peripherals_off()

    elif args.test:
        print 'todo'
        # Turn on SS
        peripherals_on(CCD=False,LASER=False)  
        # SUN CALIBRATE FUNCTION
        #ss_cardinal_test()
        # Turn off SS
        peripherals_off(CCD=False,LASER=False)

    print "Function Off"
