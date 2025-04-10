#!/usr/bin/python
# Filename: CASES.py

import matplotlib
import numpy as np
import datetime
from matplotlib.pyplot import figure, show, rc, grid, setp
from math import pi
from pylab import *

def load_txinfo(txfname, N = None):
    if N is None:
        # Load all data from txinfo
        wk, gpsSec, gpsFracSec, az, el, health, txsystem, txprn = genfromtxt(txfname,skip_footer=1, unpack=True)
    else:
        # Load last N entries (N = 15 min * 60 s/min * 12 channels = 10800)
        import os
        lines = os.popen('tail -' + str(N) + ' ' + txfname).readlines()

        wk = []
        gpsSec = []
        az = []
        el = []
        health = []
        txsystem = []
        txprn = []
        for line in lines[0:-1]:
            temp = line.split()
            wk.append(int(temp[0]))
            gpsSec.append(int(temp[1]))
            az.append(float(temp[3]))
            el.append(float(temp[4]))
            health.append(int(temp[5]))
            txsystem.append(int(temp[6]))
            txprn.append(int(temp[7]))
            
        # Cast data
        txsystem = np.array(txsystem)
        txprn = np.array(txprn)
        el = np.array(el)
        az = np.array(az)

    # Convert in to GPS Time datetime structure
    dn = []
    for i in range(0,len(wk)):
        dn.append(datetime.datetime(1980,1,6,0,0,0)+datetime.timedelta(wk[i]*7,gpsSec[i]))

    txdn = np.array(dn)
    
    return txprn, txdn, el, az, txsystem
    
def load_scint(s4fname):
    # Load data from scintinfo
    wk, gpsSec, gpsFracSec, s4, s4system, s4prn = np.genfromtxt(s4fname,skip_footer=1, unpack=True, usecols=[0,1,2,4,13,14])

    # Convert in to GPS Time datetime structure
    dn = []
    for i in range(0,wk.size):
        dn.append(datetime.datetime(1980,1,6,0,0,0)+datetime.timedelta(wk[i]*7,gpsSec[i]))

    s4dn = np.array(dn)
    
    return s4prn, s4dn, s4, s4system
    
def plot_elaz(txprn,txdn,el,az,txsystem,s4prn,s4dn,s4,s4system,fname, el_mask = None, dt = None):
    # Elevation mask to use
    if el_mask is None:
        el_mask = 25
        
    if dt is None:
        dt = 15
    
    # The target date and timedelta to plot
    target_dn = s4dn[-1];
    target_dt = datetime.timedelta(0,dt*60);
    
    # Create a polar plot with north to the top and east to the right with zenith
    # angles running from 0 to 90
    fig = figure(figsize=[8,6])
    cm = matplotlib.cm.get_cmap('jet')
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(pi/2)
    ax.set_theta_direction(-1)
    ax.set_rticks([0,30,60,90])
    ax.set_rgrids([30,60,90], labels=('', '', ''))
    ax.set_rlim(0,90)
    ax.set_thetagrids((0,45,90,135,180,225,270,315),labels=('N','','E','','S','','W',''))
    
    # Grab only data for L1 within the time frame of interest
    txidx = (txsystem == 0) & (abs(txdn-target_dn) < target_dt)
    my_txdn = txdn[txidx]
    my_txprn = txprn[txidx]
    my_el = el[txidx]
    my_az = az[txidx]
    
    s4idx = (s4system == 0) & (abs(s4dn-target_dn) < target_dt)
    my_s4prn = s4prn[s4idx]
    my_s4dn = s4dn[s4idx]
    my_s4 = s4[s4idx]
    
    c = None
    
    # Make the plot
    for i in range(1,33):
        # Grab indexes for the current satellite
        txidx = (my_txprn.T == i) & (my_el.T > el_mask)
        s4idx = (my_s4prn.T == i)
        
        if any(txidx) and any(s4idx):
            
            # Just save this satellite's s4 index values
            mys4 = my_s4[s4idx]
        
            # Find the intersection of the time indexes 
            idx_int = set(my_s4dn[s4idx]) & set(my_txdn[txidx])
    
            # Sort the list
            idx_int = np.sort(list(idx_int))
            
            for j in range(0,idx_int.size):
                # Find matching times for this satellite
                idx = (my_txdn == idx_int[j]) & (my_txprn == i)
        
                if any(idx):
                    # Plot the data on a scatter plot with the colorscale going from 0 to 1
                    c = scatter(my_az[idx]*pi/180,90-my_el[idx],c=mys4[j],s=20,edgecolor='none',vmin=0.0,vmax=1.0)
        
            # Label the last point of each satellite with the satellite PRN
            if any(idx_int):
                # Plot the last point a little larger
                c = scatter(my_az[idx]*pi/180,90-my_el[idx],c=mys4[j],s=40,edgecolor='none',vmin=0.0,vmax=1.0)
                c.set_alpha(0.75)
                    
                temp = my_el[txidx]
                if temp[-1] > temp[-2]:
                    # Satellite is setting
                    text(my_az[idx]*pi/180,90-my_el[idx]-7.5,i,
                         horizontalalignment='center',
                         verticalalignment='center')
                else:
                    # Satellite is rise
                    text(my_az[idx]*pi/180,90-my_el[idx]+7.5,i,
                         horizontalalignment='center',
                         verticalalignment='center')
            
    # Set the title with the target time
    title('CASES at Hawaii: ' + target_dn.strftime('%d %b %Y %H:%M') + ' UT')
    
    # The colorbar
    if c is not None:
        cb = colorbar(c,ticks=[0,0.2,0.4,0.6,0.8,1.0])
        cb.set_label('S$_4$')

    fig.savefig(fname)
    
def plot_s4summary(txprn,txdn,el,az,txsystem,s4prn,s4dn,s4,s4system,fname):
    # Create a polar plot with north to the top and east to the right with zenith
    # angles running from 0 to 90
    fig = figure(figsize=[6,8])
    subplots_adjust(hspace=0,wspace=0)

    # Grab only data for L1 within the time frame of interest
    txidx = (txsystem == 0)
    my_txdn = txdn[txidx]
    my_txprn = txprn[txidx]
    my_el = el[txidx]
    my_az = az[txidx]

    s4idx = (s4system == 0)
    my_s4prn = s4prn[s4idx]
    my_s4dn = s4dn[s4idx]
    my_s4 = s4[s4idx]

    d0 = s4dn[len(s4dn)/2].replace(hour = 0, minute = 0, second = 0)
    d1 = d0+datetime.timedelta(1)

    #d0 = s4dn[0].replace(hour = 0, minute = 0, second = 0)
    #d1 = s4dn[0]+datetime.timedelta(days = 1)

    # Make the plot
    for i in range(1,33):
        # Set up the axis for this subplot
        ax = subplot(19,2,i)
        ax.set_ylim(0,1.2)
        ax.yaxis.set_major_locator( FixedLocator([0.0,0.6]) )
        ax.yaxis.set_minor_locator( FixedLocator([0.3,0.9]) )
        ax.yaxis.set_major_formatter( FormatStrFormatter('%.1f') )
        setp( ax.get_yticklabels(), fontsize=10)
        ax.set_xlim(d0,d1)
        ax.xaxis.set_major_locator( HourLocator(byhour=range(0,24,6)) )
        ax.xaxis.set_minor_locator( HourLocator(byhour=range(0,24,3)) )
        setp( ax.get_xticklabels(), visible=False)
        
        # turn off first y-axis on even plots
        if mod(i,2) == 0:
            setp( ax.get_yticklabels(), visible=False)
        
        # Grab indexes for the current satellite
        txidx = (my_txprn.T == i)
        s4idx = (my_s4prn.T == i)
        
        # Plot the s4 data
        if any(txidx) and any(s4idx):
            scatter(my_s4dn[s4idx],my_s4[s4idx],color='blue',edgecolor='none',s=4)
            
        # Put the PRN label
        text(d0+datetime.timedelta(minutes = 30),1.0,'PRN ' + str(i),
            horizontalalignment='left',
            verticalalignment='center',
	    fontsize=8)
            
        # Set up the axis for the el axis
        ax2 = twinx()
        ax2.set_ylim(0,90)
        ax2.yaxis.set_major_locator( FixedLocator([0,45]) )
        ax2.yaxis.set_minor_locator( FixedLocator([15,30,60,75]) )
        ax2.yaxis.set_major_formatter( FormatStrFormatter('%.1f') )
        setp( ax2.get_yticklabels(), fontsize=10)
        ax2.set_xlim(d0,d1)
        ax2.xaxis.set_major_locator( HourLocator(byhour=range(0,24,6)) )
        ax2.xaxis.set_minor_locator( HourLocator(byhour=range(0,24,3)) )
        setp( ax2.get_xticklabels(), visible=False)
        
        # turn off second y-axis on odd plots
        if mod(i,2) == 1:
            setp( ax2.get_yticklabels(), visible=False)
        
        # Plot elevation angle
        if any(txidx) and any(s4idx):
            scatter(my_txdn[txidx],my_el[txidx],color='red',edgecolor='none',s=4)

        # Special subplot labels
        if i == 1:
            title('CASES at Hawaii', fontsize=10)
            ax.yaxis.set_major_locator( FixedLocator([0.0,0.6,1.2]) )

        if i == 2:
            title(s4dn[0].strftime('%d %b %Y %H:%M') + '-' + s4dn[-1].strftime('%H:%M') + ' UT', fontsize=10)
            ax2.yaxis.set_major_locator( FixedLocator([0,45,90]) )
        
        if (i == 31) or (i == 32):
            ax.xaxis.set_major_formatter( DateFormatter('%H') )
            setp( ax.get_xticklabels(), visible=True, fontsize=10)
            ax.set_xlabel('UT (hrs)')
            
        if i == 15:
            ax.set_ylabel('S$_4$')
    
        if i == 16:
            ax2.set_ylabel('Elevation')

        # Plot the summary s4 plot at the bottom
        if any(txidx) and any(s4idx):
            ax = subplot(12,1,12)
            scatter(my_s4dn[s4idx],my_s4[s4idx],color=cm.jet(1.*i/32),edgecolor='none',s=4)

    # Labels for the summary s4 plot    
    ax = subplot(12,1,12)
    ax.set_ylim(0,1.2)
    ax.yaxis.set_major_locator( FixedLocator([0.0,0.5,1.0]) )
    ax.yaxis.set_minor_locator( FixedLocator([0.25,0.75]) )
    ax.yaxis.set_major_formatter( FormatStrFormatter('%.1f') )
    ax.set_ylabel('S$_4$',fontsize=10)
    ax.set_xlim(d0,d1)
    ax.xaxis.set_major_locator( HourLocator(byhour=range(0,24,6)) )
    ax.xaxis.set_minor_locator( HourLocator(byhour=range(0,24,3)) )
    ax.xaxis.set_major_formatter( DateFormatter('%H') )
    setp( ax.get_xticklabels(), visible=True)
    ax.set_xlabel('UT (hrs)')
    setp( ax.get_xticklabels(), visible=True)

    fig.savefig(fname)
