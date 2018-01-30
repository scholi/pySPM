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

def fitSpectrum(t, m, error=False):
    """
    fit and return the sf and k0 parameters for given known times t and masses m
    """
    M = np.sqrt(np.array(m))
    T = np.hstack([np.array(t)[:,None],np.ones((len(t),1))])
    x = np.linalg.lstsq(T,M)[0]
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
    return k0+sf*np.sqrt(m)
    
def time2mass(t, sf, k0):
    return ((t-k0)/sf)**2
    
def getMass(elt):
    from warnings import warn
    warn("Function getMass is deprecated. Please use \"get_mass\" instead")
    return get_mass(elt)
    
def get_mass(elt):
    import os
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    m = 0
    for A,x,n in re.findall('(?:\\^([0-9]+))?([A-Z][a-z]?)_?([0-9]*)',elt):
        if n == '':
            n = 1
        else:
            n = int(n)
        if A == '':
            c.execute("SELECT mass from elements where symbol='{sym}' and abund=(SELECT max(abund) from elements where symbol='{sym}')".format(sym=x))
        else:
            c.execute("SELECT mass from elements where symbol='{sym}' and A={A}".format(sym=x, A=A))
        res = c.fetchone()
        if res is None:
            raise Exception("Cannot fetch mass of {}".format(x))
        m += n*res[0]
    m -= me*elt.count('+')
    return m
    
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
