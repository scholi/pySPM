# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Provides some useful mathematical functions which are not present in numpy.
"""

import numpy as np

def fact(x):
    """
    Return the factors of an integer as a list.
    Warning: This function is not a factorial!
    """
    if x < 0 or type(x) is not int:
        raise ValueError("input must be a positive integer")
    if x < 2:
        return x    
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

def FT(x, ufunc=np.real, real=False):
    """
    shortcut for 1D/2D Fourier Transform (real and centered)
    """
    assert isinstance(x, np.ndarray)
    if len(x.shape) == 1:
        if real:
            F = np.fft.rfft(x)
        else:
            F = np.fft.fft(x)
    elif len(x.shape) == 2:
        if real:
            F = np.fft.rfft2(x)
        else:
            F = np.fft.fft2(x)
    else:
        raise TypeError("The array should be 1D or 2D")
    return ufunc(np.fft.fftshift(F))
    
def binning(data, N=2, axis=0, ufunc=np.sum):
    w = int(np.floor(data.shape[axis]/N))
    r = np.copy(data)
    size = list(data.shape)
    size[axis] = w
    size = size[:axis+1]+[N]+size[axis+1:]
    r.resize(size)
    return ufunc(r, axis=axis+1)
    
def stat_info(data):
    import matplotlib.pyplot as plt
    D = np.ravel(data)
    U = np.unique(D)
    if len(U)>1:
        sep = np.min(U[1:]-U[:-1])
        N = min(100, int(np.ceil((np.max(D)-np.min(D))/sep)))
    else:
        N = 1
    
    mean = np.mean(D)
    std = np.std(D)
    
    fig, ax = plt.subplots(2,1,figsize=(21,4))
    ax[0].boxplot(D, 0, 'ro', 0);
    ax[1].hist(D, N, density=True);
    ax[1].axvline(mean, color='r', label='mean')
    ax[1].axvline(mean+std, color='r', linestyle='--', label='1$\\sigma$')
    ax[1].axvline(mean-std, color='r', linestyle='--', label='1$\\sigma$')
    if mean-2*std >= U[0]:
        ax[1].axvline(mean-2*std, color='r', linestyle=':', label='2$\\sigma$')
    if mean+2*std <= U[-1]:
        ax[1].axvline(mean+2*std, color='r', linestyle=':', label='2$\\sigma$')
    ax[1].legend();
    print("Stats")
    print("\tAverage:", mean)
    print("\tStandard-deviation:", std)
    print("\tMinimum:", np.min(D))
    print("\tQ1:", np.percentile(D, 25))
    print("\tMedian:", np.percentile(D, 50))
    print("\tQ3:", np.percentile(D, 75))
    print("\tMaximum:", np.max(D))
    
def LG2D(XY, A, a, sx, sy, x0, y0, lgx, lgy):
    """
    Return a 2D Lorentz-Gauss.
    XY: (X,Y) tuple
    A: amplitude
    a: angle
    sx: sigma for x-axis
    sy: sigma fdor y-axis
    x0,y0 : center coordinates of the peak
    lgx, lgy: Lorentz-Gauss proportion (for x,y axis)
    """
    X1 = (XY[0]-x0)*np.cos(a) - (XY[1]-y0)*np.sin(a)
    Y1 = (XY[0]-x0)*np.sin(a) + (XY[1]-y0)*np.cos(a)
    return A*LG(X1, 0, sx, 1, lgx)*LG(Y1,0,sy,1, lgy)