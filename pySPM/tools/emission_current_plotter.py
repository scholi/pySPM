import psutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import strpdate2num
import pandas as pd
import sys
from pySPM.tools.fpanel import Fpanel

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
        dt = df.iloc[-1,0]-df.iloc[0,0]
        for a in ax2:
            if dt.seconds < 15*60:
                a.xaxis.set_major_locator(mpl.dates.MinuteLocator(interval=1))
            elif dt.seconds < 3*60*60:
                a.xaxis.set_major_locator(mpl.dates.MinuteLocator(interval=5))
            else:
                a.xaxis.set_major_locator(mpl.dates.MinuteLocator(interval=15))
            a.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
            a.grid()
        plt.minorticks_off()

        if watch:
            mypause(3)
        else:
            plt.show()
                
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

def main():
    if len(sys.argv)>1:
        filename = sys.argv[1]
        print("Plot file \"{}\"".format(filename))
        plotLog(filename, watch=False)
    else:
        F = Fpanel()
        logfile = F.getLogFile()
        plotLog(logfile, watch=True)

if __name__ == '__main__':
    main()
