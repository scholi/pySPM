import numpy as np
from scipy import stats, optimize as opt
import re
import sqlite3

# electron mass
me = 0.00054858 # u

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
    import os
    this_dir, this_filename = os.path.split(__file__)
    DB_PATH = os.path.join(this_dir, "elements.db")
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
            c.execute("SELECT mass from elements where symbol='{sym}' and A={A}".format(sym=x,A=A))
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


def fact(x):
    if x < 2:
        return x
    import math
    f = []
    i = 2
    while True:
        while x % i == 0:
            f.append(i)
            x /= i
        i += 1
        if x == 1:
            return f


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

def moving_average(x, N):
    import numpy as np
    assert len(x) > N
    c = np.cumsum(x)
    return (c[N:]-c[:-N])/N


def butter_lowpass(cutOff, fs, order=5):
    import scipy.signal
    nyq = 0.5 * fs
    normalCutoff = cutOff / nyq
    b, a = scipy.signal.butter(order, normalCutoff, btype='low', analog=True)
    return b, a


def butter_lowpass_filter(data, cutOff, fs, order=4):
    import scipy.signal
    b, a = butter_lowpass(cutOff, fs, order=order)
    y = scipy.signal.lfilter(b, a, data)
    return y

def Gauss(x, x0, s, A=None):
    R = np.exp(-(x-x0)**2/(2*s**2))
    if A is None:
        return R /(s*np.sqrt(2*np.pi))
    return A*R
    
def Lorentz(x, x0, gamma, A=None):
    R = 1/((x-x0)**2+(.5*gamma)**2)
    if A is None:
        return .5*gamma*R/np.pi
    return A*R*(.5*gamma)**2

def CDF(x,mu,sig, lg=0):
    from scipy.special import erf
    g = sig*np.sqrt(2*np.log(2))
    return lg*(.5+np.arctan2(x-mu,g)/np.pi)+(1-lg)*.5*(1+erf((x-mu)/(sig*np.sqrt(2))))
    
def LG(x, x0, sig=None, Amp=1, lg=.5, FWHM=None):
    assert sig is not None or FWHM is not None
    if FWHM is None:
        FWHM = 2*np.sqrt(2*np.log(2))*sig
    if sig is None:
        sig = FWHM/(2*np.sqrt(2*np.log(2)))
    return Amp*((1-lg)*Gauss(x,x0,sig,A=1)+lg*Lorentz(x,x0,FWHM,A=1))

def logistic(x, lower=0, upper=1, growth=1, x0=0, nu=1, C=1):
    return lower+(upper-lower)/(C+np.exp(-growth*(x-x0)))**(1/nu)    
    
def fitCDF1line(A):
    line = np.zeros(A.shape[1])
    for x in range(A.shape[1]):
        popt, pcov = opt.curve_fit(CDF, np.arange(A.shape[0]), A[:,x], (10, A.shape[0]/2, 1))
        line[x] = popt[1]
    return line
    
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

def binning(data, N=2, axis=0, ufunc=np.sum):
    w = int(np.floor(data.shape[axis]/N))
    r = np.copy(data)
    size = list(data.shape)
    size[axis] = w
    size = size[:axis+1]+[N]+size[axis+1:]
    r.resize(size)
    return ufunc(r, axis=axis+1)
    
class ProgressBar:
    def __init__(self, seq=None, min=0, max=100, value=0):        
        if seq is not None:
            self.min = 0
            self.max = len(seq)-1
            self.value = 0
        else:
            self.min = min
            self.max = max
            self.value = value
        if in_ipynb():
            from ipywidgets import IntProgress
            from IPython.display import display
            self.pb = IntProgress(min=self.min,max=self.max,value=self.value)
            display(self.pb)
        else:
            self.pb = None
            
    def set_value(self, value):
        self.value = value
        if self.pb:
            self.pb.value = value