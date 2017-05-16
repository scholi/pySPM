import numpy as np
from scipy import stats, optimize as opt

Elts = {
    'B': [10, 11],
    'C': 12,
    'N': 14,
    'O': 16,
    'OH': 17,
    'F': 19,
    'Na': 23,
    'Mg': 24,
    'CN': 26,
    'Al': 27,
    'Si': 28,
    'P': 31,
    'S': 32,
    'SH': 33,
    'Cl': [35, 37],
    'TiO': [64, 65],
    'TiO2': [80, 81],
    'Cs': 133,
    'Au': 197
}


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

def CDF(x,a,x0,sig):
    return a*stats.norm.cdf((x-x0),scale=sig)

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