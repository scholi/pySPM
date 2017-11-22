import numpy as np

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
    
def LG(x, x0, sig=None, Amp=1, lg=.5, FWHM=None, asym=1):
    assert sig is not None or FWHM is not None
      
    if FWHM is None:
        FWHM = 2*np.sqrt(2*np.log(2))*sig
    if sig is None:
        sig = FWHM/(2*np.sqrt(2*np.log(2)))
    Y = Amp*((1-lg)*Gauss(x,x0,sig,A=1)+lg*Lorentz(x,x0,FWHM,A=1))
    if asym!=1:
        Yr = Amp*((1-lg)*Gauss(x,x0,sig*asym,A=1)+lg*Lorentz(x,x0,FWHM*asym,A=1))
        Y[x>x0] = Yr[x>x0]
    return Y

def logistic(x, lower=0, upper=1, growth=1, x0=0, nu=1, C=1):
    return lower+(upper-lower)/(C+np.exp(-growth*(x-x0)))**(1/nu)    
    
def fitCDF1line(A):
    line = np.zeros(A.shape[1])
    for x in range(A.shape[1]):
        popt, pcov = opt.curve_fit(CDF, np.arange(A.shape[0]), A[:,x], (10, A.shape[0]/2, 1))
        line[x] = popt[1]
    return line
    
def binning(data, N=2, axis=0, ufunc=np.sum):
    w = int(np.floor(data.shape[axis]/N))
    r = np.copy(data)
    size = list(data.shape)
    size[axis] = w
    size = size[:axis+1]+[N]+size[axis+1:]
    r.resize(size)
    return ufunc(r, axis=axis+1)
    