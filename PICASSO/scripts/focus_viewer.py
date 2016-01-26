#!/usr/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.mlab import prctile
import numpy as np
import Image
import TifImagePlugin
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler
from watchdog.events import FileSystemEventHandler
from sys import argv, exit
import time

# RESIZE FIGURES HERE...
mpl.rcParams['savefig.dpi'] = 400
mpl.rcParams['figure.figsize'] = (8,5)


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
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    fig.show()
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
                    pmin = prctile(old.getdata(),p=5)
                    pmax = prctile(old.getdata(),p=95)
                    ax1.imshow(np.reshape(old.getdata(), old.size),vmin=pmin,vmax=pmax,cmap = mpl.cm.gray)
                    ax1.set_title('Old range = %i-%i'%(pmin,pmax))
                    ax1.set_xlabel('UT: %s\nExp. Time: %.1f sec\nTemp: %.1fC' % (old.info['UniversalTime'].strftime('%m-%d-%y %H:%M:%S'),old.info['ExposureTime'],old.info['CCDTemperature']))

                # add new image
                print 'Adding New Plot'
                new = Image.open(path+'/'+filename)
                pmin = prctile(new.getdata(),p=5)
                pmax = prctile(new.getdata(),p=95)   
                ax2.imshow(np.reshape(new.getdata(), new.size),vmin=pmin,vmax=pmax,cmap = mpl.cm.gray)
                ax2.set_title('New range = %i-%i'%(pmin,pmax))
                ax2.set_xlabel('UT: %s\nExp. Times: %.1f sec\nTemp: %.1fC' % (new.info['UniversalTime'].strftime('%m-%d-%y %H:%M:%S'),new.info['ExposureTime'],new.info['CCDTemperature']))

                # Plot stuff
                old = new
                plt.gcf().subplots_adjust(bottom=0.2)
                plt.ion()
                plt.draw()
                flag = False

    except KeyboardInterrupt:
        observer.stop()
