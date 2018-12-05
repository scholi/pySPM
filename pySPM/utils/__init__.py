# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

import numpy as np
from scipy import stats, optimize as opt
import re
import sqlite3
from .math import *
from .elts import *
from .constants import *
from .spectra import *
from .plot import *
from . import fit, misc, colors
from .save import *
from .restoration import *
from .misc import alias, PB

def funit(value, unit=None):
    """
    Convert a value and unit to proper value
    e.g. 0.01m will be converted to 10mm
    e.g. 1000A will be converted to 1kA
    etc.

    Parameters
    ----------
    value : int, float or dictionary with 'value' and 'unit' key
        numerical value to work with
    unit : string or None
        name of the unit

    Returns
    -------
    dictionary with formated 'value' and 'unit'.
    The value will be: ≥1 and <1000

    Examples
    --------
    >>> funit(0.01,'m')
    {'value': 10.0, 'unit': 'mm'}
    >>> funit(1, 'cm')
    {'value': 1.0, 'unit': 'cm'}
    >>> funit({'value':2340, 'unit': 'um'})
    {'value': 2.34, 'unit': 'mm'}
    """
    if unit == None:
        unit = value['unit']
        value = value['value']
    import math
    if value == 0:
        return {'value': value, 'unit': unit}
    shift = int(math.floor(math.log10(value)/3.0))  # base10 exponent
    mag = u'afpnum1kMGTPE'
    index_mag = mag.index('1')
    # Test if unit has a scale factor (n,m,k,M,G, etc.)
    if len(unit) > 1 and unit[0] in mag and unit[1:] and index_mag:
        index_mag = mag.index(unit[0])
        unit = unit[1:]
    value /= 10**(3*shift)
    index_mag += shift
    if index_mag < 0:
        value *= 10**(index_mag*3)
        index_mag = 0
    elif index_mag >= len(mag):
        value *= 10**((index_mag-len(mag)+1)*3)
        index_mag = len(mag)-1
    unit_prefix = mag[index_mag]
    if unit_prefix == '1':
        unit_prefix = ''
    return {'value': value, 'unit': u'{mag}{unit}'.format(mag=unit_prefix, unit=unit)}
    
def get_shifts_bbox(shifts, shape):
    """
    Return the bounding-box of the overlapping area of scans which have a thermal drift given by shifts.
    Return: dictionary with the top, left, right, bottom position
    
    Parameters
    ----------
    shifts: list of tuple [(dx0, dy0), (dx1, dy1), ..., (dxN, dyN)]
        list of tuple given the shift for each scan in pixel in the x- and y-direction.
        The number of element should equal the number of scans (so N=Nscan-1)
    shape: tuple (Height, Width)
        shape of the image.
    """
    from .geometry import Bbox
    right = shape[1]+min(0,np.min([x[0] for x in shifts]))
    left = np.max([x[0] for x in shifts])
    top = shape[0]+min(0,np.min([x[1] for x in shifts]))
    bottom = np.max([x[1] for x in shifts])
    return Bbox(right=right, left=left, top=top, bottom=bottom)
    
def s2hms(s):
    """Convert seconds to hour, minutes or seconds"""
    M = np.max(s)
    if M>120*60:
        return s/(60*60), 'h'
    if M>300:
        return s/60, 'min'
    return s, 's'

def time2hms(s, string=True):
    h = int(s//3600)
    s -= h*3600
    m = int(s//60)
    s -= m*60
    if string:
        return "{:02d}:{:02d}:{:02.2f}".format(h,m,s)
    return (h,m,s)

@alias("fitSpectrum")
def fit_spectrum(t, m, error=False):
    """
    fit and return the sf and k0 parameters for given known times t and masses m
    """
    M = np.sqrt(np.array(m))
    T = np.hstack([np.array(t)[:,None],np.ones((len(t),1))])
    x = np.linalg.lstsq(T,M,rcond=-1)[0]
    sf = 1/x[0]
    k0 = -sf*x[1]
    res = [sf, k0]
    if error:
        mfit = np.dot(T, x)
        mresid = M - mfit
        X = np.dot(T.T, T)
        epsvar = np.var(mresid, ddof=2)
        xvar = np.linalg.inv(X) * epsvar # correlation matrix
        D = np.sqrt(np.diag(xvar)) # error vector
        Dsf = np.abs(D[0]/x[0]**2)
        Dk0 = np.sqrt((D[1]/x[0])**2+(x[1]*Dsf)**2)
        res += [Dsf, Dk0]
    return res

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def mass2time(m, sf, k0):
    """
    Here the time is actually the number of channels.
    To convert it to real time, you should multiply the answer by the Time Resolution which can be obtained by pySPM.ITM.get_value("Registration.TimeResolution")
    """
    r = m*0+k0
    if not hasattr(r, '__setitem__'):
        if m<0:
            return r
        return sf*np.sqrt(m)+k0
    r[m>=0] = sf*np.sqrt(m[m>=0])+k0
    return r

def time2mass(t, sf, k0):
    """
    Here the time t is actually the number of channels.
    To convert a real time to the number of channel t, you should divide the real time by the Time Resolution which can be obtained by pySPM.ITM.get_value("Registration.TimeResolution")
    """
    return ((t-k0)/sf)**2

def show_table(t):
    from IPython.core.display import display, HTML
    display(HTML(html_table(t)))

def html_table(t, header=False):
    S = "<table>"
    if header:
        S += "<tr>"+"".join(["<th>{0}</th>".format(x) for x in t[0]])+"</tr>"
        t = t[1:]
    for x in t:
        S += "<tr>"+"".join(["<td>{0}</td>".format(y) for y in x])+"</tr>"
    S += "</table>"
    return S

def aa_table(t, header=False):
    """
    print a list of list in a nice ascii-art
    """
    Ncols = len(t[0])
    Lcol = [0]*Ncols
    for x in t:
        for i in range(Ncols):
            Lcol[i] = max(Lcol[i], len(repr(x[i])))
    if header:
        print(
            "  ".join([u"{: <"+str(Lcol[i]+4)+"}" for i in range(Ncols)]).format(*t[0]))
        print("="*sum(Lcol))
        t = t[1:]
    for j, x in enumerate(t):
        print("  ".join([u"{:"+['.', '_'][j % 2]+"<" +
                         str(Lcol[i]+4)+"}" for i in range(Ncols)]).format(*x))

def dict_update(d, u):
    import collections
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = dict_update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d

def htmlTable(t, show=True, header=False):
    import html
    s = "<table><tr>"
    if header:
        s += "<th>"+"</th><th>".join([html.escape(str(y))
                                      for y in t[0]])+"</th></tr><tr>"
        t = t[1:]
    s += "".join([
        "<tr><td>"+"</td><td>".join([
            html.escape(str(y)) for y in x]) + "</td></tr>" for x in t])+"</table>"
    if show:
        from IPython.display import display, HTML
        display(HTML(s))
    else:
        return s

def zfit_level(A, processes=4):
    import multiprocessing as mp
    level = np.array(mp.Pool(processes).map(fitCDF1line, (A[:,i,:] for i in range(A.shape[1]))))
    return level

def getToFimg(I, N=100, prog=False):
    """
    Simulate obtained counts for N scans from an image
    """
    R = np.zeros(I.shape)
    L = range(N)
    if prog:
        L = PB(L)
    for i in L:
        R += (np.random.rand(*I.shape)<(1-np.exp(-I)))
    return R

def getToFsimg(I, N=[10, 50, 100, 300, 500], prog=False):
    Ns = np.insert(np.diff(N),0,N[0])
    T = [getToFimg(I, n) for n in Ns]
    return dict(zip(N,np.cumsum(T, axis=0)))

def centered_meshgrid(A):
    w = A.shape[1]
    h = A.shape[0]
    Y, X = np.mgrid[-h//2:h//2,-w//2:w//2]
    return Y,X
