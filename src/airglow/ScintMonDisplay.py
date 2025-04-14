import ScintMon
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as dates

def PlotDay(sum_file, el_mask = 20.):
# Plots a day's s4 index with the requested el_mask

    t, s4, rx_pos, sx, sy, sz, site = ScintMon.parse_sum(sum_file)
    
    (M,N) = np.shape(s4)
    el = np.ones_like(s4) * np.nan
    
    for prn in range(N):
        sats = np.array([sx[:,prn],sy[:,prn],sz[:,prn]])
        (e,az) = ScintMon.elaz(rx_pos,sats)
        el[:,prn] = e
        
    data = np.ma.array(s4,mask=(el<el_mask) | (np.isnan(s4)))
    print(el_mask)

    s4_fig = plt.figure()
    s4_graph = s4_fig.add_subplot(111)
    s4_graph.hold(True)
    
    s4_graph.plot(t,data,linewidth=0.15)
    s4_graph.plot(t,data.max(axis=1),'k',linewidth=1)
    s4_graph.set_ylim([0,1.2])
    s4_graph.set_ylabel('s4')
    s4_graph.set_xlabel('UT')
    if t[0].day != t[-1].day:
        s4_graph.set_title(site + ': ' + t[0].strftime('%d') + t[-1].strftime('-%d %b %Y'))
    else:
        s4_graph.set_title(site + ': ' + t[0].strftime('%d %b %Y'))   
    s4_graph.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))

    return min(t), max(t), s4_fig
