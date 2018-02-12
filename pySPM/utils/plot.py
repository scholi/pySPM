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
    
def Xdist(ax,left, right, y, color='r', linestyle=':', fmt="{dist:.1f}{unit}", xtransf=lambda x: x, va='bottom', ha='center', offset=(0,2), **kargs):
    ax.axvline(left, color=color, linestyle=linestyle)
    ax.axvline(right, color=color, linestyle=linestyle)
    ann = dict(va=va, ha=ha, color=color)
    ann.update({k[3:]:kargs[k] for k in kargs if k.startswith('an_')})
    ax.annotate(fmt.format(dist=xtransf(right-left),unit=kargs.get('unit','')), ({'center':.5*(left+right),'left':right,'right':left}[ann['ha']], y), offset, textcoords='offset pixels', **ann)
    arr = dict(arrowstyle='<->',color=color)
    arr.update({k[4:]:kargs[k] for k in kargs if k.startswith('arr_')})
    ax.annotate("", (left, y), (right, y), arrowprops=arr)
    
def DualPlot(ax, col1='C0',col2='C1'):
    axb = ax.twinx()
    axb.spines['left'].set_color(col1)
    axb.spines['right'].set_color(col2)
    ax.yaxis.label.set_color(col1)
    axb.yaxis.label.set_color(col2)
    ax.tick_params(axis='y', colors=col1)
    axb.tick_params(axis='y', colors=col2)
    return axb

def sp(M, N=1, W=21, ravel=True, fig=False):
    """
    Shortcut for creating subplots with max width = W (default 21,
        which seems to correspond to 100% of the width in jupyter).
    Height is calculated in order to have square subplots
    """
    f, ax = plt.subplots(N, M, figsize=(W, N*W/M))
    if ravel:
        if fig:
            return np.ravel(ax), f
        return np.ravel(ax)
    if fig:
        return ax, f
    return ax
    
def get_rect(img, bottom, top, left, right, ax, color='r'):
    ax.axhline(bottom, color=color)
    ax.axhline(top, color=color)
    ax.axvline(left, color=color)
    ax.axvline(right, color=color)
    return img[bottom:top,left:right]
    
def sublegend(*ax, labels=None, color='white', margin=9, titles=None, fontsize=14):
    props = dict(boxstyle='round', facecolor=color, alpha=1)  
    if not hasattr(margin, "__getitem__") and not hasattr(margin, "__iter__"):
        margin = (margin, margin)

    if labels is None:
        labels = [chr(ord('a')+i) for i,_ in enumerate(np.ravel(ax))]
    for i,a in enumerate(np.ravel(ax)):
        if titles is False:
            a.set_title("")
        a.annotate(labels[i],(0, 1),xytext=(margin[0],-margin[1]),fontsize=fontsize,verticalalignment='top', bbox=props, xycoords='axes fraction',textcoords='offset pixels');

    