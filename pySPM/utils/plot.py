import numpy as np
import matplotlib.pyplot as plt

def plotMask(ax, mask, color, **kargs):
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
