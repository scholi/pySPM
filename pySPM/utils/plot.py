# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Various helper function for plotting data with matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt

def plotMask(ax, mask, color, **kargs):
    """
    Create an overlay of a given color where mask is True.
    The transparency and other parameters can be adjusted with **kargs. see help(plt.imshow)
    """
    import copy
    m = np.ma.masked_array(mask, ~mask)
    palette = copy.copy(plt.cm.gray)
    palette.set_over(color, 1.0)
    ax.imshow(m, cmap=palette, vmin=0, vmax=0.5, **kargs)
    
def Xdist(ax,left, right, y, color='r', linestyle=':', fmt='.2f', xtransf=lambda x: x, **kargs):
    ax.axvline(left,color=color, linestyle=linestyle)
    ax.axvline(right,color=color, linestyle=linestyle)
    s = "{:"+fmt+"}"+kargs.get('unit','')
    ax.annotate(s.format(xtransf(right-left)),(.5*(left+right),y),(0,2),textcoords='offset pixels',va='bottom',ha='center')
    ax.annotate("",(left,y),(right,y),arrowprops=dict(arrowstyle=kargs.get('arrowstyle','<->')))
    
def DualPlot(ax, col1='C0',col2='C1'):
    axb = ax.twinx()
    axb.spines['left'].set_color(col1)
    axb.spines['right'].set_color(col2)
    ax.yaxis.label.set_color(col1)
    axb.yaxis.label.set_color(col2)
    ax.tick_params(axis='y', colors=col1)
    axb.tick_params(axis='y', colors=col2)
    return axb

def sp(M,N=1,W=21):
    """
    Shortcut for creating subplots with max width = W (default 21,
        which seems to correspond to 100% of the width in jupyter).
    Height is calculated in order to have square subplots
    """
    fig, ax = plt.subplots(N,M,figsize=(W,N*W/M))
    return ax