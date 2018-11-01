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
from . import fit
from .save import *
from .restoration import *

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
    
def fitSpectrum(t, m, error=False):
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
        Dsf = D[0]/x[0]**2
        Dk0 = np.sqrt((D[1]/x[0])**2+(x[1]*D[0]/x[0]**2)**2)
        res.append(Dsf, Dk0)
    return res

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def mass2time(m, sf, k0):
    return sf*np.sqrt(m)+k0

def time2mass(t, sf, k0):
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

def in_ipynb():
    try:
        cfg = get_ipython().config 
        if 'IPKernelApp' in cfg:
            return True
        else:
            return False
    except NameError:
        return False

def getToFimg(I, N=100, prog=False):
    """
    Simulate obtained counts for N scans from an image
    """
    R = np.zeros(I.shape)
    L = range(N)
    if prog:
        from tqdm import tqdm_notebook as tqdm
        L = tqdm(L)
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
