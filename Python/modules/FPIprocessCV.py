
import sys
import numpy as np
import FPI
from datetime import datetime
import time
from scipy import interpolate
#from matplotlib.pyplot import *


class Data:
    def __init__(self, dn, project):
        self.dn = dn
        self.project = project
        self.key = ""
        self.f = ""
        self.u  = np.array([]) 
        self.ue = np.array([])
        self.v  = np.array([])
        self.ve = np.array([])
        self.w  = np.array([])
        self.we = np.array([])

        self.cloud = np.array([])

        # interpolated stuff:
        self.it  = np.array([])
        self.iw  = np.array([])
        self.iwe = np.array([])

        self.T  = np.array([])
        self.Te = np.array([])
        self.t1 = np.array([])
        self.t2 = np.array([])
        self.log = ""
        self.notes = ""
        self.length = 0
        self.rev = "??"
        self.error = False
    def __str__(self):
        string = ""
        string += "%11s" % "dn = "     + self.dn.strftime("%Y-%m-%d") + "\n"
        string += "%11s" % "project = " + self.project + "\n"
        string += "%11s" % "f = " + self.f + "\n"
        string += "%11s" % "key = " + self.key + "\n"
        string += "%11s" % "log = " + self.log + "\n"
        string += "%11s" % "notes = " + self.notes + "\n"
        string += "%11s" % "length = " + "%3i" % self.length + "\n"
        return string
    def plot(self, switch_onefig=True):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        from datetime import timedelta
        import matplotlib as mpl
        mpl.rcParams['font.family'] = 'monospace'

        switch_plot_u = True
        switch_plot_v = True
        switch_plot_w = True
        switch_plot_iw = True

        if len(self.u) <3:
            switch_plot_u = False
        if len(self.v) <3:
            switch_plot_v = False
        if len(self.w) <3:
            switch_plot_w = False
        if len(self.iw) <3:
            switch_plot_iw = False

        if (not (switch_plot_u or switch_plot_v or switch_plot_w)):
            return None
        if self.error:
            return None

        if switch_onefig:
            fig = plt.figure(1); 
        else:
            fig = plt.figure();
        plt.clf()

        ax = fig.add_axes((.1,.2,.8,.7)) # left, bottom, width, height
        
        if switch_plot_u:
            plt.errorbar(self.t1, self.u, yerr=self.ue, \
                    color='b', marker='o', label='u')
        if switch_plot_v:
            plt.errorbar(self.t1, self.v, yerr=self.ve, \
                    color='g', marker='+', label='v')
        if switch_plot_w:
            plt.errorbar(self.t1, self.w, yerr=self.we, \
                    color='r', marker='*', label='w')

        if switch_plot_iw:
            plt.errorbar(self.it, self.iw, yerr=self.iwe, \
                color='k', label='iw')

        dnp1 = self.dn + timedelta(days=1)
        plt.xlim( [datetime(self.dn.year, self.dn.month, self.dn.day, 20), datetime(dnp1.year, dnp1.month, dnp1.day, 12)] )
        plt.ylim([-200.,200.]) 
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        plt.legend()
        plt.grid()
        fig.text(.1,.05,self.notes)
        datestr = self.dn.strftime("%Y-%m-%d")
        fig.text(.1,.92,"%10s, %12s, %10s" % (self.project, self.key, datestr))
        
        if len(self.cloud) > 0:
            ax1 = fig.gca()
            ax2 = ax1.twinx()
            ax2.plot(self.it, self.cloud, 'c--', linewidth=.5)
            ax2.set_ylim([-50., -15.])
        

        fig.text(.7,.030, self.log)
        plt.draw();
        #plt.show()
        return 0


def MasterDictionary():
    Sites = { \
        'UAO': {'Name':'Urbana Atmospheric Observatory',
            'Project': 'NATION',
            #'Location': (40.133, -88.2, 200),
            'Combos': ('ANN','EKU')},
        'ANN': {'Name': 'Peach Mountain',
            'Project': 'NATION',
            'Combos': ('UAO','EKU')},
        'EKU': {'Name': 'Eastern Kentucky University',
            'Project': 'NATION',
            'Combos': ('UAO','ANN')},
        'PAR': {'Name': 'Pisgah Aeronomical Research Institute',
            'Project': 'NATION',
            'Combos': ('EKU',)},
        'CAJ': {'Name': 'Cajazeiras',
            'Project': 'RENOIR',
            'Combos': ('CAR',)},
        'CAR': {'Name': 'Cariri',
            'Project': 'RENOIR',
            'Combos': ('CAJ',)},
        'MRH': {'Name': 'Merihill',
            'Project': 'PERU',
            'Combos': ('NZK',)},
        'NZK': {'Name': 'Nazca',
            'Project': 'PERU',
            'Combos': ('MRH',)}, 
        'KAF': {'Name': 'Kirtland Airforce Base',
            'Project': 'NATION',
            'Combos': ()},
        'NATION': ['UAO','ANN','EKU', 'PAR', 'KAF'], # CLEAN THIS UP LATER
        'RENOIR': ['CAJ','CAR'] # CLEAN THIS UP LATER
        }
    return Sites

def CardFinder(dn,site1,project, reference):
    '''
    Summary
    -------

    data = CardFinder(dn, site1, project)
    
    Returns data for cardinal mode points,

    Inputs
    ------
        dn = datetime day
        site1 = site, e.g. 'EKU' or 'PAR'

    Outputs
    -------
        data =xxxx dictionary of data whose keys are the look directions of the site
              xxxxx The entries of these keys are the
              xxxx    time, u/ue, v/ve, w/we, T/Te, clouds?
                 xxx time, winds, Temps, clouds, etc...
              xxx Run data['EKU_N'].keys() to see full list of entries

    History
    -------
    3/26/13 -- Written by DJF (dfisher2@illionis.edu),
                        & TMD (duly2@illinois.edu)

    '''
    import os
    import copy
    #print "CardFinder: dn=",dn,"site1=",site1

    # Ouput variable is an instances of class Data()
    d = Data(dn, project)

    d.key = site1
    
    # get Master Dictionary
    Sites = MasterDictionary()

    stub = '/mnt/FPIData/Results/'

    # load in data:
    f = '%s/%s_%s_%4i%02d%02d.npz' % (stub,project, site1, dn.year, dn.month, dn.day)
    if not os.path.isfile(f):
        d.log += "%s is not found \n" % f
        d.error = True
        return [d]
    d.f += "\n%s" % f
    npzfile = np.load(f,allow_pickle=True)
    r1 = npzfile['FPI_Results'].ravel()[0]
    #r1['el'] = 90 - r1['ze']  # Elevation angles
    npzfile.close()
    
    # get reference winds:
    zind1 = FPI.all_indices('Zenith',r1['direction'])

    # removes erroneous values from std:
    zind1 = np.delete(zind1, np.where(np.abs(r1['sigma_LOSwind'][zind1]) > 50000.)) 
    
    # check to see if we have some vertical wind measurements
    # our CV winds depend on these
    if len(zind1)==0:
        d.log += "no vertical wind measurements found\n"
        d.error = True
        return [d]

    # reference winds:
    ref1 = np.array(FPI.DopplerReference(r1, reference=reference))
    d.notes += "Reference is " + reference + "\n"
    
    # get vertical winds:
    w1 = (r1['LOSwind'][zind1]-ref1[0][zind1]) / -np.cos(r1['ze'][zind1]*np.pi/180.)
    zind1 = np.delete(zind1, np.where(np.abs(w1) > 300.)) # removes erroneous values actual wind values
    
    # for the interpolation, get rid of indice
    # values that correspond to cloudy observations
    if 'Clouds' in r1.keys():
        if r1['Clouds'] is not None:
            zind1 = np.delete(zind1, np.where( r1['Clouds']['mean'][zind1] > -25. ))
            d.cloud = r1['Clouds']['mean']

    # now we have to re-run the vertical wind with the updated indices:
    w1 = (r1['LOSwind'][zind1]-ref1[0][zind1]) / -np.cos(r1['ze'][zind1]*np.pi/180.)
    we1 = r1['sigma_LOSwind'][zind1]

    # if we only have less than 4 data points, we can't do a cubic spline,
    # so just exit with none:
    if len(zind1) < 4:
        d.log += "less than 4 data points, can't do cubic spline \n"
        d.error = True
        return [d]
        
    # interpolate vertical winds/errors for all times:
    tck1 =    interpolate.splrep(np.array([time.mktime(dn.timetuple()) for dn in  r1['sky_times'][zind1]]),w1)
    w1_all =  interpolate.splev([time.mktime(dn.timetuple()) for dn in r1['sky_times']], tck1)

    # bad spline interpolation:
    #tck1 =    interpolate.splrep(np.array([time.mktime(dn.timetuple()) for dn in  r1['sky_times'][zind1]]),we1)
    #we1_all = interpolate.splev([time.mktime(dn.timetuple()) for dn in r1['sky_times']], tck1)

    f1 = interpolate.interp1d( np.array( [time.mktime(dn.timetuple()) for dn in  r1['sky_times'][zind1]] ), we1,\
            bounds_error=False, fill_value=0.0)
    we1_all = f1([time.mktime(dn.timetuple()) for dn in r1['sky_times']])

    d.it = r1['sky_times']
    d.iw = w1_all
    d.iwe = we1_all
    
    # keep cardinal looks only and make it unique:
    looks = list(set([val for val in r1['direction'] if val in ['Zenith','North','East','South','West']]))

    # ------------------------------------------------
    # loop thru for different cardinal directions
    # ------------------------------------------------
    ds = []
    for pt in looks:
        
        # copy the data instance with 
        # information we have so far:
        d_loop = copy.deepcopy(d)

        # reset output:
        u = np.array([]); ue = np.array([])
        v = np.array([]); ve = np.array([])
        w = np.array([]); we = np.array([])

        #print "pt=",pt,
        ind1 = FPI.all_indices(pt,r1['direction'])
        
        # Record look times
        t1 = r1['sky_times'][ind1]
    
        if 'Z' in pt[0]:
            #print "zenith measurement"
            # ------------------
            # Zenith measurement
            # ------------------
    
            t1 = r1['sky_times'][zind1]
            w = w1
            we = we1
            d_loop.notes += 'Vertical wind is measurement\n'           
            
        elif 'E' in pt[0]:
            # ------------------
            # Eastern Zonal measurement
            # ------------------

	    # Calculated Horizontal Winds
            u = (r1['LOSwind'][ind1]-ref1[0][ind1]+w1_all[ind1]*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    -np.sin(r1['ze'][ind1]*np.pi/180.)
            ue = np.sqrt( r1['sigma_LOSwind'][ind1]**2+(we1_all[ind1])**2*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    np.sin(r1['ze'][ind1]*np.pi/180.)
            d_loop.notes += 'Vertical wind is interpolated\n'
            
        elif 'W' in pt[0]:
            # ------------------
            # Western Zonal measurement
            # ------------------

	    # Calculated Horizontal Winds
            u = (r1['LOSwind'][ind1]-ref1[0][ind1]+w1_all[ind1]*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    np.sin(r1['ze'][ind1]*np.pi/180.)
            ue = np.sqrt( r1['sigma_LOSwind'][ind1]**2+(we1_all[ind1])**2*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    np.sin(r1['ze'][ind1]*np.pi/180.)
            d_loop.notes += 'Vertical wind is interpolated\n'
            
        elif 'N' in pt[0]:
            # ----------------------
            # Northern Meridional measurement
            # ----------------------

	    # Calculated Horizontal Winds
            v = (r1['LOSwind'][ind1]-ref1[0][ind1]+w1_all[ind1]*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    -np.sin(r1['ze'][ind1]*np.pi/180.)
            ve = np.sqrt( r1['sigma_LOSwind'][ind1]**2+(we1_all[ind1])**2*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    np.sin(r1['ze'][ind1]*np.pi/180.)
            d_loop.notes += 'Vertical wind is interpolated\n'
            
        else:
            # ----------------------
            # Sourthern Meridional measurement
            # ----------------------

	    # Calculated Horizontal Winds
            v = (r1['LOSwind'][ind1]-ref1[0][ind1]+w1_all[ind1]*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    np.sin(r1['ze'][ind1]*np.pi/180.)
            ve = np.sqrt( r1['sigma_LOSwind'][ind1]**2+(we1_all[ind1])**2*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    np.sin(r1['ze'][ind1]*np.pi/180.)
            d_loop.notes += 'Vertical wind is interpolated\n'


        # ------------
        #  end pt loop
        # ------------


        # Save information
        d_loop.key = "%s_%s" % (site1, pt)

        d_loop.u = u
        d_loop.ue = ue
        d_loop.v = v
        d_loop.ve = ve
        d_loop.w = w
        d_loop.we = we
        d_loop.T = r1['T'][ind1]
        d_loop.Te = r1['sigma_T'][ind1]
        d_loop.t1 = t1
        d_loop.length = len(t1)

        ds.append(d_loop)

    return ds

def CVFinder(dn,site1,site2,project, reference):
    '''
    Summary
    -------

    data = CVFinder(dn, site1, site2)
    
    Returns data for common value points,

    Inputs
    ------
        dn = datetime day
        site1, site2 = two sites of comparision, e.g. 'EKU' and 'PAR'

    Outputs
    -------
        data = xxxxdictionary of data whose keys are the look directions between 
               xxxxxxthe two sites, e.g. 'CV_EKU_PAR_2',
               xxxxxThe entries of these keys are the
                  meridional wind/error
                  zonal wind/error
                  times of measurements
               Run data['CV_EKU_PAR_2'].keys() to see full list of entries

    History
    -------
    3/11/13 -- Written by DJF (dfisher2@illionis.edu),
                        & TMD (duly2@illinois.edu)

    '''
    #print "dn=",dn,'site1=',site1,'site2=',site2
    import os 
    import copy


    stub = '/mnt/FPIData/Results/'

    # Ouput variable is an instances of class Data()
    d = Data(dn, project)

    d.log = "\ncreated with CVFinder\n"
    d.log += "site1 = "+site1+", site2 = "+site2 + "\n"

    d.key = "%s_%s" % (site1, site2)

    # get Master Dictionary
    Sites = MasterDictionary()

    # Check to see that we have a valid pair:
    if (site1 in Sites.keys() and site2 in Sites[site1]['Combos']):
        project = Sites[site1]['Project']
    else:
        sys.exit("Bad Site Pair")

    # load in data:

    # for site1:
    f = '%s%s_%s_%4i%02d%02d.npz' % (stub,project, site1, dn.year, dn.month, dn.day)
    # if file is not found, log and return
    if not os.path.isfile(f):
        d.log += "%s is not found \n" % f
        d.error = True
        return [d]
    d.f += "\n%s" % f
    npzfile = np.load(f,allow_pickle=True)
    r1 = npzfile['FPI_Results'].ravel()[0]
    # r1['az'] = 90 - r1['ze']
    npzfile.close()
    
    # for site2:
    f = '%s%s_%s_%4i%02d%02d.npz' % (stub,project, site2, dn.year, dn.month, dn.day)
    # if file is not found, log and return
    if not os.path.isfile(f):
        d.log += "%s is not found \n" % f
        d.error = True
        return [d]
    d.f += "\n%s" % f
    npzfile = np.load(f,allow_pickle=True)
    r2 = npzfile['FPI_Results'].ravel()[0]
    # r2['az'] = 90 - r2['ze']
    npzfile.close()
    
    # get ----??----reference winds:
    zind1 = FPI.all_indices('Zenith',r1['direction'])
    zind1 = np.delete(zind1, np.where(np.abs(r1['sigma_LOSwind'][zind1]) > 50000.)) # removes erroneous values from std
    
    zind2 = FPI.all_indices('Zenith',r2['direction'])
    zind2 = np.delete(zind2, np.where(r2['sigma_LOSwind'][zind2] > 50000.)) # removes erroneous values

    # check to see if we have some vertical wind measurements
    # our CV winds depend on these
    if ((len(zind1)==0) or (len(zind2)==0)):
        d.log += "no vertical wind measurements found\n"
        d.error = True
        return [d]
    
    # reference winds:
    ref1 = np.array(FPI.DopplerReference(r1, reference=reference))
    ref2 = np.array(FPI.DopplerReference(r2, reference=reference))
    d.notes += "Reference is " + reference + "\n"
    
    # get vertical winds:
    w1 = (r1['LOSwind'][zind1]-ref1[0][zind1])/-np.cos(r1['ze'][zind1]*np.pi/180.)
    w2 = (r2['LOSwind'][zind2]-ref2[0][zind2])/-np.cos(r2['ze'][zind2]*np.pi/180.)

    # removes erroneous values actual wind values
    zind1 = np.delete(zind1, np.where(np.abs(w1) > 300.)) 
    zind2 = np.delete(zind2, np.where(np.abs(w2) > 300.))

    # for the interpolation, get rid of indice
    # values that correspond to cloudy observations
    if (('Clouds' in r1.keys()) and ('Clouds' in r2.keys())):
        if ((r1['Clouds'] is not None) and (r2['Clouds'] is not None)):
            zind1 = np.delete(zind1, np.where( r1['Clouds']['mean'][zind1] > -25. ))
            zind2 = np.delete(zind2, np.where( r2['Clouds']['mean'][zind2] > -25. ))
            d.cloud = r1['Clouds']['mean']
    
    # now we have to re-run the vertical wind with the updated indices:
    w1 = (r1['LOSwind'][zind1]-ref1[0][zind1])/-np.cos(r1['ze'][zind1]*np.pi/180.)
    w2 = (r2['LOSwind'][zind2]-ref2[0][zind2])/-np.cos(r2['ze'][zind2]*np.pi/180.)
    we1 = r1['sigma_LOSwind'][zind1]
    we2 = r2['sigma_LOSwind'][zind2]
    
    # if we only have less than 4 data points, we can't do a cubic spline,
    # so just exit with none:
    if (len(zind1) < 4) or (len(zind2) < 4):
        d.log += "less than 4 data points, can't do cubic spline \n"
        d.error = True
        return [d]

    # interpolate vertical winds/erros for all times:
    tck1 = interpolate.splrep(np.array([time.mktime(dn.timetuple()) for dn in  r1['sky_times'][zind1]]),w1)
    w1_all = interpolate.splev([time.mktime(dn.timetuple()) for dn in r1['sky_times']], tck1)
    #tck1 = interpolate.splrep(np.array([time.mktime(dn.timetuple()) for dn in  r1['sky_times'][zind1]]),we1)
    #we1_all = interpolate.splev([time.mktime(dn.timetuple()) for dn in r1['sky_times']], tck1)
    
    tck2 = interpolate.splrep(np.array([time.mktime(dn.timetuple()) for dn in  r2['sky_times'][zind2]]),w2)
    w2_all = interpolate.splev([time.mktime(dn.timetuple()) for dn in r2['sky_times']], tck2)

    #tck2 = interpolate.splrep(np.array([time.mktime(dn.timetuple()) for dn in  r2['sky_times'][zind2]]),we2)
    #we2_all = interpolate.splev([time.mktime(dn.timetuple()) for dn in r2['sky_times']], tck2)
    
    f1 = interpolate.interp1d( np.array( [time.mktime(dn.timetuple()) for dn in  r1['sky_times'][zind1]] ), we1,\
            bounds_error=False, fill_value=0.0)
    we1_all = f1([time.mktime(dn.timetuple()) for dn in r1['sky_times']])

    f2 = interpolate.interp1d( np.array( [time.mktime(dn.timetuple()) for dn in  r2['sky_times'][zind2]] ), we2,\
            bounds_error=False, fill_value=0.0)
    we2_all = f2([time.mktime(dn.timetuple()) for dn in r2['sky_times']])

    d.it = r1['sky_times']
    d.iw = w1_all
    d.iwe = we1_all
    
    # get common value locations between the 2 sites:
    common_pair = [val for val in r1['direction'] if val in r2['direction']]
    
    # get rid of cardinal modes in our list, and make it unique:
    common_pair = list(set([val for val in common_pair if val not in ['Zenith','North','East','South','West','None']]))
    
    # check to make sure we have some common pairs. 
    # if not, exit
    if len(common_pair) == 0:
        d.log += "no common pairs found \n"
        d.error = True
        return [d]

    # ------------------------------------------------
    # loop thru for different common volume directions
    # ------------------------------------------------
    ds = []
    for cv in common_pair:
        
        # copy the data instance with 
        # information we have so far:
        d_loop = copy.deepcopy(d)
        
        # reset output:
        u = np.array([]); ue = np.array([]) 
        v = np.array([]); ve = np.array([]) 
        w = np.array([]); we = np.array([]) 

        ind1 = FPI.all_indices(cv,r1['direction'])
        ind2 = FPI.all_indices(cv,r2['direction'])
    
        # -----------------------------------------------------
        # this section of code gets rid of times that do 
        # not match their partner CV times
        t1_list = r1['sky_times'][ind1]
        t2_list = r2['sky_times'][ind2]
    
        bad_t = []
        for kk, t1 in enumerate(t1_list):
            flag = False
            for t2 in t2_list:
                dt = (t1-t2)
                if np.abs(dt.seconds+dt.days*86400.) < 3*60.:
                    flag = True
            if not flag:
                bad_t.append(kk)
        ind1 = np.delete(ind1, bad_t)
    
    
        bad_t = []
        for kk, t2 in enumerate(t2_list):
            flag = False
            for t1 in t1_list:
                dt = (t1-t2)
                if np.abs(dt.seconds+dt.days*86400.) < 3*60. :
                    flag = True
            if not flag:
                bad_t.append(kk)
        ind2 = np.delete(ind2, bad_t)
        # -----------------------------------------------------
        
        # Keep coincident look times
        t1 = r1['sky_times'][ind1]
        t2 = r2['sky_times'][ind2]
    
        # Get Average Temperature/Error at pair pt
        T = np.average(np.vstack((r1['T'][ind1],r2['T'][ind2])), \
                axis=0, \
                weights=1./np.vstack((r1['sigma_T'][ind1],r2['sigma_T'][ind2])))
        Te = np.average(np.vstack((r1['sigma_T'][ind1],r2['sigma_T'][ind2])), \
                axis=0, \
                weights=1./np.vstack((r1['sigma_T'][ind1],r2['sigma_T'][ind2])))
        
        if 'IN' in cv[:2]:
            # ------------------
            # INline measurement
            # ------------------
    
            # FIX cos[ze] is from look direction, not at the point in the sky where they are equal...
            w = (r1['LOSwind'][ind1]-ref1[0][ind1]+r2['LOSwind'][ind2]-ref2[0][ind2])/ \
                    (-np.cos(r1['ze'][ind1]*np.pi/180.) - np.cos(r2['ze'][ind2]*np.pi/180.))
            we = np.sqrt(r1['sigma_LOSwind'][ind1]**2+r2['sigma_LOSwind'][ind2]**2)/ \
                    ( np.cos(r1['ze'][ind1]*np.pi/180.) + np.cos(r2['ze'][ind2]*np.pi/180.))
            d_loop.notes += 'Vertical wind is inline measurement\n'
            
        else:
            # ------------------
            # CV measurement
            # ------------------
            
            # Calculated Average vertical wind at cv points
            w = np.average(np.vstack((w1_all[ind1],w2_all[ind2])), \
                    axis=0, \
                    weights=1./np.vstack((we1_all[ind1],we2_all[ind2])))
            we = np.average(np.vstack((we1_all[ind1],we2_all[ind2])), \
                    axis=0, \
                    weights=1./np.vstack((we1_all[ind1],we2_all[ind2])))
            d_loop.notes += 'Vertical wind is interpolated\n'
            
            # FIX cos[ze] is from look direction, not at the point in the sky where they are equal...
            # vh = velocity horizontal
            vh1 = (r1['LOSwind'][ind1]-ref1[0][ind1]+w*np.cos(r1['ze'][ind1]*np.pi/180.))/ \
                    -np.sin(r1['ze'][ind1]*np.pi/180.)
            vh2 = (r2['LOSwind'][ind2]-ref2[0][ind2]+w*np.cos(r2['ze'][ind2]*np.pi/180.))/\
                    -np.sin(r2['ze'][ind2]*np.pi/180.)

            vh1e = np.sqrt(r1['sigma_LOSwind'][ind1]**2+(we*np.cos(r1['ze'][ind1]*np.pi/180.))**2)/ \
                    np.sin(r1['ze'][ind1]*np.pi/180.)
            vh2e = np.sqrt(r2['sigma_LOSwind'][ind2]**2+(we*np.cos(r2['ze'][ind2]*np.pi/180.))**2)/ \
                    np.sin(r2['ze'][ind2]*np.pi/180.)

            az1 = np.mean(r1['az'][ind1])
            az2 = np.mean(r2['az'][ind2])
            
            # Calculate winds
            M = np.array([[np.sin(az1*np.pi/180.0),np.cos(az1*np.pi/180.0)],[np.sin(az2*np.pi/180.0),np.cos(az2*np.pi/180.0)]])
            u = [] ; ue = []
            v = [] ; ve = []
            for (myt,a) in zip(t1,vh1):
                i = min(range(len(t2)), key=lambda i: abs(t2[i]-myt))
                b = vh2[i]
                temp = np.linalg.solve(M,np.array([[a],[b]]))
                u.append(temp[0][0])
                v.append(temp[1][0])
                temp = np.sqrt(np.linalg.solve(M**2,np.array([[vh1e[i]**2],[vh2e[i]**2]])))
                ue.append(temp[0][0])
                ve.append(temp[1][0])

        # Save information
        d_loop.key = cv

        d_loop.u = np.array(u)
        d_loop.ue = np.array(ue)
        d_loop.v = np.array(v)
        d_loop.ve = np.array(ve)
        d_loop.w = w
        d_loop.we = we
        d_loop.T = T
        d_loop.Te = Te
        d_loop.t1 = t1
        d_loop.t2 = t2
        d_loop.length = len(t1)

        ds.append(d_loop)

    return ds

def GetDataForADay(dn,project):
    import glob

    datas = []
    Sites = MasterDictionary()
    #stub = "/mnt/FPIData/Results/"
    #path = "%s/%s*_%s%02d%02d.npz" % (stub,project, dn.year, dn.month, dn.day)
    #fs = glob.glob(path)
    #
    #Sites = MasterDictionary()

    #datas = []
    #print fs
    #for f in fs:
    #    site1 = f[-16:-13]
    #    #print f, site1

    #    #print site1, Sites[site1]['Combos']
    #    #print 'site1=',site1, 'combos=',Sites[site1]['Combos']


    #    '''
    #    OK SO YOU JUST ASSIGNED SITE1, BUT THEN YOU COMPLETELY
    #    DISREGARDED THAT BELOW BECAUSE YOU LOOP THRU 
    #    SITES[PROJECT].....
    #    '''

    reference = 'Laser'
    #reference = 'Zenith'

    for site1 in Sites[project]:
        temps= CardFinder(dn, site1, project, reference)
        for temp in temps:
            datas.append(temp)
                

        for combo in Sites[site1]['Combos']:
            temps = CVFinder(dn, site1, combo, project, reference)
            # flatten it to append to out:
            for temp in temps:
                datas.append(temp)

    return datas
    #out = {}
    #for data in datas:
    #    for a_key in data.keys():
    #        out[a_key] = data[a_key]
    #return out


def simple_test():
    dn = datetime(2013,3,14)
    #site1 = "ANN"
    #site2 = "EKU"
    return CardFinder(dn, "PAR", "NATION", 'Zenith')

    #return GetDataForADay(dn, "NATION")

def complicated_test():
    import os
    print "removing pngs..."; os.system("rm -rf /home/duly/cv/*.png")
    #dn_start = datetime(2013,3,14)
    #dn_end = datetime(2013,3,15)
    dn_start = datetime(2013,3,2)
    dn_end = datetime(2013,3,3)
    dns = [dn_start + timedelta(days=k) for k in range((dn_end-dn_start).days)]

    ds = []
    for dn in dns:
        out = GetDataForADay(dn, "RENOIR" )

        for o in out:
            ds.append(o)

        #plt.close('all')
        for kk, d in enumerate(out):
            datestr = d.dn.strftime("%Y-%m-%d")
            #print "-------------------------",datestr
            p = d.plot()
            if p is not None:
                file_name = "%s_%s_%s" % (d.project, datestr, d.key)
                print "saving:", file_name
                plt.savefig("/home/duly/cv/%s.png" % file_name)
            #print d
            #print "len(d.t1)=",len(d.t1)
            #print "len(d.t2)=",len(d.t2)
            #print "len(d.u)=",len(d.u)
            #print "len(d.v)=",len(d.v)
            #print "len(d.w)=",len(d.w)
        #out[-1].plot()

    #count = 0
    #for a_d in ds:
    #    if a_d.length > 10: count+=1
    #print count
    return ds
    #return None

if __name__=="__main__":
    import matplotlib.pyplot as plt
    from datetime import timedelta

    #out = GetDataForADay( datetime(2013,1,2), "NATION" )
    #out = GetDataForADay( datetime(2013,3,1), "NATION" )
    #out = GetDataForADay( datetime(2013,2,5), "NATION" )
    #out = simple_test()
    out = complicated_test()

    #count = {}
    #for o in out:
    #    if o.length > 10:
    #        if o.key in count.keys():
    #            count[o.key] += 1
    #        else:
    #            count[o.key] = 1

    
    #ds = CVFinder( datetime( 2013, 3, 3), "PAR", "EKU", "NATION")
    #for d in ds:
    #    d.plot()


