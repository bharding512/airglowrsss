#!/usr/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.mlab import prctile
import numpy as np
import Image
import TifImagePlugin
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from sys import argv, exit
import time

# RESIZE FIGURES HERE...
mpl.rcParams['savefig.dpi'] = 400
mpl.rcParams['figure.figsize'] = (8,5)


class img:
    def __init__(self,path,filename):
        new = Image.open(path+'/'+filename)
        self.size = new.size
        self.data = np.reshape(new.getdata(), self.size)
        self.pmax = prctile(new.getdata(),p=95)
        self.pmin = prctile(new.getdata(),p=5)
        self.title= 'New range = %i-%i'%(self.pmin,self.pmax)
        self.label= 'UT: %s\nExp. Times: %.1f sec\nTemp: %.1fC' % (new.info['UniversalTime'].strftime('%m-%d-%y %H:%M:%S'),new.info['ExposureTime'],new.info['CCDTemperature'])

class FileModifiedHandler(FileSystemEventHandler):

    def __init__(self,path,filename):
        self.filename = filename
        self.path = path
        print 'Ready...'

    def on_modified(self,event):
        global flag
        if not event.is_directory and event.src_path.endswith(self.filename):
            print "File Modified"
            flag = True


# Call this function as ./view_asi_focus.py PATH IMAGENAME
# Function will automatically plot new image and save the previous one.
# Resizing figures currently doesn't work
# Written by Dan Fisher - dfisher2@illinois.edu
if __name__ == "__main__":
    if not len(argv) == 3:
        print("no path and file specified")
        exit(1)
    path = argv[1]
    filename = argv[2]

    # Setup Plots
    fig, (ax1, ax2) = plt.subplots(1, 2,sharey=False)
    fig.show()
    plt.gcf().subplots_adjust(bottom=0.2)
    old = []

    print path,filename
    event_handler = FileModifiedHandler(path,filename)
    observer = Observer()
    observer.schedule(event_handler,path,recursive=False)
    observer.start()
    flag = False 

    try:
        while True:
            time.sleep(1)
            if flag:
                # move figures if not first image
                if old:
                    print 'Moving Old Plot'
                    ax1.imshow(old.data,vmin=old.pmin,vmax=old.pmax,cmap = mpl.cm.gray)
                    ax1.set_xlim([0,old.size[1]])
                    ax1.set_ylim([0,old.size[0]])
                    ax1.set_title(old.title)
                    ax1.set_xlabel(old.label)
                    plt.draw()
                    #time.sleep(2)

                # add new image
                print 'Adding New Plot'
                new = img(path,filename)
                ax2.imshow(new.data,vmin=new.pmin,vmax=new.pmax,cmap = mpl.cm.gray)
                ax2.set_xlim([0,new.size[1]])
                ax2.set_ylim([0,new.size[0]])
                ax2.set_title(new.title)
                ax2.set_xlabel(new.label)
                
                #plt.ioff()
                plt.draw()
                flag = False
                old = new

    except KeyboardInterrupt:
        observer.stop()
