import psutil

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import strpdate2num
import pandas as pd

def plotLog(filename, watch=False, **kargs):
    fig, ax = plt.subplots(1,1)
    fig.subplots_adjust(hspace=0)
    plt.show(block=False)
    while True:
        with open(filename, 'r') as f:
            names = f.readline().rstrip().split('\t')
        df = pd.read_csv(filename, skiprows=1, delimiter='\t', parse_dates=[0], na_values="<undefined>", names=names)
        #df = df.dropna()
        ax2 = df.plot("Time", subplots=True, ax=ax, sharex=True)
        for a in ax2:
            a.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
            a.grid()
        plt.minorticks_off()

        mypause(3)
        if not watch:
            break
                
def mypause(interval):
    backend = plt.rcParams['backend']
    if backend in mpl.rcsetup.interactive_bk:
        figManager = mpl._pylab_helpers.Gcf.get_active()
        if figManager is not None:
            canvas = figManager.canvas
            if canvas.figure.stale:
                canvas.draw()
            canvas.start_event_loop(interval)
            return

if __name__ == '__main__':
    if len(sys.args)>0:
        filename = sys.argv[1]
        plotLog(filename, watch=False, debug=options.debug)
    else:
        F = Fpanel()
        F.plotLog(watch=True)
