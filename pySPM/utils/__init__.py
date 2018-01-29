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

def fitSpectrum(t, m, error=False, dev=False):
    """
    fit and return the sf and k0 parameters fôr given know times t and masses m
    """
    M = np.sqrt(np.array(m))
    T = np.hstack([np.array(t)[:,None],np.ones((len(t),1))])
    x = np.linalg.lstsq(T,M)
    sf = 1/x[0][0]
    k0 = -sf*x[0][1]
    res = [sf, k0]
    if error:
        res.append(x[1][0])
    if dev:
        res.append(np.matmul(T,x[0])-M)
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


def RL(x, image, psf):
    from scipy.signal import convolve
    I = convolve(x, psf, mode='same') # reconvoluted estimation.
    relative_blur = np.copy(image)
    relative_blur[I>0] = relative_blur[I>0]/ I[I>0]
    return x * convolve(relative_blur, psf[::-1,::-1], mode='same') # Correlation is the convolution of mirrored psf
            
def dampedRL(x, image, psf, T=1, N=3):
    from scipy.signal import convolve
    I = convolve(x, psf, mode='same') # reconvoluted estimation.
    ratio = np.zeros(image.shape)
    ratio[image>0] = I[image>0] / image[image>0]
    logarithm = np.zeros(image.shape)
    logarithm[ratio>0] = np.log(ratio[ratio>0])
    U = -2*(image*logarithm-I+image)/T**2
    Ut = np.minimum(U, 1)
    relative_blur = np.ones(image.shape)
    relative_blur[I>0] += (Ut[I>0]**(N-1)) * (N-(N-1)*Ut[I>0]) * (image[I>0]-I[I>0])/I[I>0]
    return x * convolve(relative_blur, psf, mode='same')
    
def richardson_lucy(image, psf, iterations, damped=0, N=10, **kargs):
    import scipy
    from scipy.signal import convolve
    
    image = image.astype(np.float)
    psf = psf.astype(np.float)
    psf /= np.sum(psf) # Normalize the psf ⇒ ∫∫ psf(x,y) dx dy = 1
    
    # im_deconv = 0.5 * np.ones(image.shape)
    im_deconv = kargs.get('x0', np.mean(image) * np.ones(image.shape))
            
    # As the image and the psf are both ≥ 0, if the image value is 0, then the result should also be 0 at this position
    #im_deconv[image == 0] = 0

    # Is iterations a number of a list of number?
    dict_output = True
    if type(iterations) is int:
        dict_output = False
        iterations = [iterations]
  
    N = max(iterations)
        
    results = {}

    for i in range(N):
        if damped == 0:            
            im_deconv = RL(im_deconv, image, psf)
        else:
            im_deconv = dampedRL(im_deconv, image, psf, T=damped, N=N)
        if i+1 in iterations:
            results[i+1] = np.copy(im_deconv)
    if dict_output:
        return results
    return results[N]

def accelerated_richardson_lucy(image, psf, iterations, **kargs):    
    def _RL(x):
        T = kargs.get('T', kargs.get('damped', 0))
        if T==0:   
            return RL(x, image, psf)
        return dampedRL(x, image, psf, N=kargs.get('N',3), T=T)
    
    image = image.astype(np.float)
    psf   = psf.astype(np.float)
    psf  /= np.sum(psf)
    
    x0 = kargs.get('x0', np.mean(image) * np.ones(image.shape))    
    x1 = kargs.get('x1', _RL(x0))

    h1 = x1 - x0
    y0 = x0
    y1 = x1
    g0 = _RL(y0) - y0
    g1 = _RL(y1) - y1
        
    dict_output = True
    if type(iterations) is int:
        dict_output = False
        iterations = [iterations]
  
    N = max(iterations)
        
    results = {}

    for i in range(2, N):
        x2 = y1 + g1 # x_{k+1} = y_k + g_k
        h2 = x2 - x1
        a2 = np.sum(g1*g0)/np.sum(g0**2)
        y2 = x2 + a2*h2
        g2 = _RL(y2)-y2
        
        # Shift data
        x0, x1 = x1, x2
        y0, y1 = y1, y2
        g0, g1 = g1, g2
        h0, h1 = h1, h2
        if i+1 in iterations:
            results[i+1] = np.copy(x1)
    if dict_output:
        return results
    return results[N]