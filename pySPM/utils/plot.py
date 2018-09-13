# -- coding: utf-8 --
# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Various helper function for plotting data with matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from .misc import dec_debug, do_debug

def plotMask(ax, mask, color, **kargs):
    """
    Create an overlay of a given color where mask is True.
    The transparency and other parameters can be adjusted with **kargs. see help(plt.imshow)
    
    ax: matplotlib axis
    mask: An array with value True where the mask should be plotted
    color: color of the mask
    **kargs: dictionnary of arguments which can be passed to imshow. Useful is mainly: alpha
    """
    import copy
    m = np.ma.masked_array(mask, ~mask)
    palette = copy.copy(plt.cm.gray)
    palette.set_over(color, 1.0)
    ax.imshow(m, cmap=palette, vmin=0, vmax=0.5, **kargs)

def offset_coord(xy, offset=(0,0), ax=None, fig=None, unit='px'):
    if ax is None:
        ax = plt.gca()

    if fig is None:
        fig = plt.gcf()

    tr = ax.transData
    tri = tr.inverted()

    if unit is 'px':
        offset = np.array(offset)
    elif unit is 'ax':
        offset = ax.transAxes.transform(offset)
    elif unit is 'fig':
        offset = fig.transFigure.transform(offset)
    return tri.transform(tr.transform(xy)+offset)

def Xdist(ax,left, right, y, color='r', linestyle=':', fmt="{dist:.1f}{unit}", xtransf=lambda x: x, va='bottom', ha='center', offset=(0,2), **kargs):
    """
    Show the distance between to x-values
    
    ax: matplotlib axis
    left: x-value of the left mark
    right: x-value of the right mark
    y: height of the label
    color: color of the lines / labels
    linestyle: linestyle of the vertical lines
    fmt: Format string to parse the distance
    xtransf: a function can be passed to display the distance in an other unit that the x axis.
        example: xtransf=lambda x: x*1e3 to display the distance in nm on an axis that uses micro-meters.
    va: vertical alignment of the distance label
    ha: horizontal alignment of the distance label
    offset: (x,y) offset of the label in pixels
    **kargs: additional arguments.
        arguments starting with amn_ will be sent to the annotation
        arguments starting with arr_ will be sent to the arrow
        for example arr_color='r', ann_color='b' will display a red arrow with the distance in blue
    """
    
    ax.axvline(left, color=color, linestyle=linestyle)
    ax.axvline(right, color=color, linestyle=linestyle)
    ann = dict(va=va, ha=ha, color=color)
    ann.update({k[3:]:kargs[k] for k in kargs if k.startswith('an_')})
    ax.annotate(fmt.format(dist=xtransf(right-left),unit=kargs.get('unit','')), ({'center':.5*(left+right),'left':right,'right':left}[ann['ha']], y), offset, textcoords='offset pixels', **ann)
    arr = dict(arrowstyle='<->',color=color)
    arr.update({k[4:]:kargs[k] for k in kargs if k.startswith('arr_')})
    ax.annotate("", (left, y), (right, y), arrowprops=arr)

def Ydist(ax, down, up, x, color='r', linestyle=':', fmt="{dist:.2f}{unit}", ytransf=lambda y: y, rotation=90, va='center', ha='right', offset=(-2,0), **kargs):
    """
    Show the distance between to x-values
    
    ax: matplotlib axis
    down: y-value of the bottom mark
    up: y-value of the top mark
    x: gorizontal position of the label
    color: color of the lines / labels
    linestyle: linestyle of the vertical lines
    fmt: Format string to parse the distance
    xtransf: a function can be passed to display the distance in an other unit that the x axis.
        example: xtransf=lambda x: x*1e3 to display the distance in nm on an axis that uses micro-meters.
    va: vertical alignment of the distance label
    ha: horizontal alignment of the distance label
    offset: (x,y) offset of the label in pixels
    **kargs: additional arguments.
        arguments starting with amn_ will be sent to the annotation
        arguments starting with arr_ will be sent to the arrow
        for example arr_color='r', ann_color='b' will display a red arrow with the distance in blue
    """
    
    ax.axhline(down, color=color, linestyle=linestyle)
    ax.axhline(up, color=color, linestyle=linestyle)
    ann = dict(va=va, ha=ha, color=color)
    ann.update({k[3:]:kargs[k] for k in kargs if k.startswith('an_')})
    ax.annotate(fmt.format(dist=ytransf(up-down),unit=kargs.get('unit','')), (x,{'center':.5*(down+up),'top':up,'bottom':down}[ann['va']]), offset, textcoords='offset pixels', rotation=rotation, **ann)
    arr = dict(arrowstyle='<->',color=color)
    arr.update({k[4:]:kargs[k] for k in kargs if k.startswith('arr_')})
    ax.annotate("", (x, down), (x, up), arrowprops=arr)
    
def DualPlot(ax, col1='C0',col2='C1'):
    """
    Create a dual axis from ax and tune the color of the left/right axes to col1, col2 resp.
    """
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
    Height is calculated in order to have square subplots.
    The layout is given by M columns and N lines.
    If N is negative, abs(N) is interpreted as the number of elements and the number of line will be calculated automatically in function of M.
    """
    if N < 0:
        tot = -N
        N = 1+(tot-1)//M
    f, ax = plt.subplots(N, M, figsize=(W, N*W/M))
    if ravel:
        if fig:
            return np.ravel(ax), f
        return np.ravel(ax)
    if fig:
        return ax, f
    return ax
    
def get_rect(img, bottom, top, left, right, ax, color='r'):
    """
    return a sub-area on a 2D-array img and display the boundaries on a matplotlib axis ax
    """
    ax.axhline(bottom, color=color)
    ax.axhline(top, color=color)
    ax.axvline(left, color=color)
    ax.axvline(right, color=color)
    return img[bottom:top, left:right]
    
def sublegend(*ax, **kargs):
    labels   = kargs.get('labels', None)
    color    = kargs.get('color', 'white')
    margin   = kargs.get('margin', 9)
    titles   = kargs.get('titles', None)
    fontsize = kargs.get('fontsize',14)
    """
    ax: list of axes
    labels: If None, the will be a, b, c, d, ...
    color: background color of the rectangle
    margin: margin in pixel from the axis border
    titles: set to False to remove all titles. (Useful to keep the set_title info in the code to remember what is plotted)
    fontsize: The font size of the labels
    """
    props = dict(boxstyle='round', facecolor=color, alpha=1)  
    if not hasattr(margin, "__getitem__") and not hasattr(margin, "__iter__"):
        margin = (margin, margin)

    if labels is None:
        labels = [chr(ord('a')+i) for i,_ in enumerate(np.ravel(ax))]
    for i,a in enumerate(np.ravel(ax)):
        if titles is False:
            a.set_title("")
        a.annotate(labels[i],(0, 1),xytext=(margin[0],-margin[1]),fontsize=fontsize,verticalalignment='top', bbox=props, xycoords='axes fraction',textcoords='offset pixels');

def formula(x):
    import re
    x = re.sub('_([0-9]+)',r'$_{\1}$',x)
    x = re.sub(r'\^([0-9\+\-]+)',r'$^{\1}$',x)
    return x

def __points_in_bbox(x,y,bbox):
    x_in = np.logical_and(x>=bbox.xmin, x<=bbox.xmax)
    y_in = np.logical_and(y>=bbox.ymin, y<=bbox.ymax)
    mask = x_in & y_in
    return np.any(mask), mask

def __bbox_overlap(bbox1, bbox2):
    x_overlap = (bbox1.xmax >= bbox2.xmin) and (bbox1.xmin <= bbox2.xmax)
    y_overlap = (bbox1.ymax >= bbox2.ymin) and (bbox1.ymin <= bbox2.ymax)
    return x_overlap and y_overlap

__cached_points = {}
__cached_pointsd = {}

def __overlap(ax, obj, objs, debug=False):
    global __cached_points
    fig = ax.get_figure()
    renderer = fig.canvas.get_renderer()
    tr = ax.transData
#    tri = tr.inverted()
    r = obj.get_window_extent(renderer) # got object bbox
    for o in objs:
        if type(o) is mpl.lines.Line2D:
            if o in __cached_points:
                xy = __cached_points[o]
                xyd = __cached_pointsd[o]
            else:
                xyd = np.vstack(o.get_data()).T # retriev data in data coordinates
                xy = tr.transform(xyd) # retrieve data in pixel coordinates
                __cached_points[o] = xy
                __cached_pointsd[o] = xyd
            isin, ins = __points_in_bbox(xy[:,0],xy[:,1],r)
            if isin:
                if debug: ax.plot(xyd[:,0][ins],xyd[:,1][ins],'rx')
                return o
        else:
            ro = o.get_window_extent(renderer)
            if __bbox_overlap(r,ro):
                return o
    return False

def __get_yextend_of_curve(ax, bbox, curve):
    tri = ax.transData.inverted()
    rd = bbox.transformed(tri) # bbox in data coordinates
    x, y = curve.get_data()
    mask = (x>=rd.xmin)*(x<=rd.xmax)
    if len(y[mask])==0:
        return (0,0)
    return (np.min(y[mask]), np.max(y[mask]))


def put_Xlabels(ax, pos, labels, colors='rgb', debug=False, save=False, bbox=False, **kargs):
    fig = ax.get_figure()
    renderer = fig.canvas.get_renderer()
    P = list(zip(pos, labels))
    P.sort(key=lambda x: x[0]) # sort labels by ascending x
    labs = []
    ax.draw(renderer) # make sure to draw the new object
    ylim = ax.get_ylim()
    objs = [o for o in ax.get_children() if type(o) not in [mpl.patches.Rectangle]]
    coli = [None, None]
    for i, (x, lab) in enumerate(P):
        col = colors[i%len(colors)]
        pos = offset_coord((x, ylim[0]), (3, 0), ax=ax)
        labs.append(ax.text(pos[0], pos[1], lab, va='bottom', ha='left', rotation=90, color=col))
        ax.axvline(x, color=col, linestyle=kargs.get('linestyle','-'), alpha=kargs.get('line_alpha',.5))
        ax.draw(renderer) # make sure to draw the new object
        r = labs[-1].get_window_extent(renderer) # get the Bbox of the last label
        ov = __overlap(ax, labs[-1], objs+labs[:-1])
        if save:
            if ov:
                r1 = r.transformed(ax.transData.inverted())
                r2 = ov.get_window_extent(renderer).transformed(ax.transData.inverted())
                print(repr(ov), r2.ymin, r2.ymax)
                coli[0] = mpl.patches.Rectangle((r1.xmin, r1.ymin), r1.xmax-r1.xmin, r1.ymax-r1.ymin, ec='b', fill=False)                
                coli[1] = mpl.patches.Rectangle((r2.xmin, r2.ymin), r2.xmax-r2.xmin, r2.ymax-r2.ymin, ec='r', fill=False)

                ax.add_patch(coli[0])
                ax.add_patch(coli[1])
            fig.savefig("put_Xlabels_{:03d}_000.png".format(i))
            if ov:
                coli[0].remove()
                coli[1].remove()
        last_ov = ov
        y_offset = 0
        it = 0
        while ov:
            it += 1
            if debug: print("Label \"{}\" overlap with {}".format(labs[-1].get_text(), repr(ov)))
            if type(ov) is mpl.lines.Line2D:
                new_y = __get_yextend_of_curve(ax, r, ov)[1]
            else:
                rov = ov.get_window_extent(renderer).transformed(ax.transData.inverted())
                new_y = rov.ymax
            if ov == last_ov:
                y_offset += kargs.get('offset_step', 5)
                if y_offset > kargs.get('offset_step', 5)*kargs.get('max_iteration', 6):
                    break
            else:
                y_offset = kargs.get('offset_step', 5)
            last_ov = ov
            new_xy_px = ax.transData.transform((x,new_y))+np.array([3,y_offset])
            new_xy = ax.transData.inverted().transform(new_xy_px) # Offset the new_y of 5 pixels
            labs[-1].set_position(new_xy)
            ax.draw(renderer)
            ov = __overlap(ax, labs[-1], objs+labs[:-1], debug=debug)
            if save:
                ax.draw(renderer)
                fig.savefig("put_Xlabels_{:03d}_{:03d}.png".format(i,it))
                if ov:
                    r2 = ov.get_window_extent(renderer) # get the Bbox of the last label
                    mpl.patches.draw_bbox(r, renderer, color='r')
                    mpl.patches.draw_bbox(r2, renderer, color='r')
    if bbox:
        ax.draw(renderer)
        for o in objs:
            if type(o) is not mpl.lines.Line2D:
                r = o.get_window_extent(renderer).transformed(ax.transData.inverted())
                ax.add_patch(
                        mpl.patches.Rectangle((r.xmin,r.ymin), r.xmax-r.xmin, r.ymax-r.ymin, ec='r', fill=False)
                        )

        for l in labs:
            r = l.get_window_extent(renderer).transformed(ax.transData.inverted())
            ax.add_patch(
                    mpl.patches.Rectangle((r.xmin,r.ymin), r.xmax-r.xmin, r.ymax-r.ymin, ec='b', fill=False)
                    )
