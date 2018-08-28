# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

import scipy.optimize as opt
"""
Provides some useful mathematical functions which are not present in numpy.
"""

import numpy as np
    
def strictly_positify(x):
    """
    Make the result strictly positive by setting the minimum value to
    the lower allowed float value
    """
    return np.fmax(x, np.finfo(x.dtype).eps)
    
def positify(x):
    """
    Set to zero all negative values
    """
    return np.fmax(0, x)

def clip01(x):
    """
    clip data x between 0 and 1
    """
    return np.fmax(np.fmin(x, 1), 0)
    
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
    old = np.geterr()
    np.seterr(all='ignore')
    R = np.exp(-(x-x0)**2/(2*s**2))
    if A is None:
        R /= (s*np.sqrt(2*np.pi))
    else:
        R *= A
    R[s==0] = (x[s==0]==x0)*1.0
    np.seterr(**old)
    return R

def Lorentz(x, x0, gamma, A=None):
    R = 1/((x-x0)**2+(.5*gamma)**2)
    if A is None:
        return .5*gamma*R/np.pi
    return A*R*(.5*gamma)**2

def CDF(x,mu,sig, Amp=1, lg=0):
    from scipy.special import erf
    g = sig*np.sqrt(2*np.log(2))
    return Amp*lg*(.5+np.arctan2(x-mu,g)/np.pi)+(1-lg)*Amp*.5*(1+erf((x-mu)/(sig*np.sqrt(2))))
    
def LG(x, x0, sig=None, Amp=1, lg=.5, asym=1, FWHM=None):
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
        popt, pcov = opt.curve_fit(CDF,
            np.arange(A.shape[0]),
            A[:,x],
            (A.shape[0]/2, 1, np.max(A[:,x])),
            bounds=(0,(A.shape[1],np.inf,np.inf))
            )
        line[x] = popt[0]
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

def ellipse(a,b,phi):
    """
    return the x,y coordinates of an ellipse with major axis=a and minor axis=b at angle phi(radians)
    """
    r = np.zeros(phi.shape) # if a & b is zero, result is zero
    m = np.logical_and(a!=0,b!=0)
    r[m] = a*b/np.sqrt((b*np.cos(phi[m]))**2+(a*np.sin(phi[m]))**2)
    return r
    
def asymm_ellipse(left, right, upper, lower, phi):
    phi = np.divmod(phi+2*np.pi,2*np.pi)[1] # Be sure phi ∈ [0,2π]
    b = np.where(phi<=np.pi, upper, lower)
    a = np.where(np.logical_or(phi<=np.pi/2, phi>=3*np.pi/2), right, left)
    m = np.logical_and(a!=0,b!=0)
    r = np.zeros(phi.shape) # if a & b is zero, result is zero
    r[m] = a[m]*b[m]/np.sqrt((b[m]*np.cos(phi[m]))**2+(a[m]*np.sin(phi[m]))**2)
    return r

def LG2D(XY, amplitude=1, angle=0, sig_x=10, sig_y=10, x0=None, y0=None, LG_x=0, LG_y=0, assym_x=1, assym_y=1, bg=0):
    """
    Return a 2D Lorentz-Gauss.
    XY: (X,Y) tuple
    A: amplitude
    a: angle
    sx: sigma for x-axis
    sy: sigma fdor y-axis
    x0,y0 : center coordinates of the peak
    lgx, lgy: Lorentz-Gauss proportion (for x,y axis)
    assym_x, assym_y: The assymetry in sig_x/sig_y for the left/right or upper/lower part of the curve
    """
    if x0 is None:
        x0 = XY[0][0,XY[0].shape[1]//2]
    if y0 is None:
        y0 = XY[1][XY[1].shape[0]//2,0]
    X1 = (XY[0]-x0)*np.cos(angle) - (XY[1]-y0)*np.sin(angle)
    Y1 = (XY[0]-x0)*np.sin(angle) + (XY[1]-y0)*np.cos(angle)
    
    R1 = np.sqrt(X1**2+Y1**2)
    angle = np.arctan2(Y1,X1)
    sig = asymm_ellipse(sig_x, sig_x*assym_x, sig_y, sig_y*assym_y, angle)
    gamma = np.sqrt(2*np.log(2))*sig
    LG = ellipse(LG_x+1, LG_y+1, angle)-1
    
    Gxy = Gauss(R1, 0, sig, 1)
    Lxy = 1/((R1/gamma)**2+1)
    
    f = (1-LG)*Gxy+LG*Lxy
    out = bg+amplitude*f
    return out
    
def LG2Da(XY, amplitude=1, angle=0, sigN=10, sigS=None, sigE=10, sigW=None, x0=None, y0=None, LGN=0, LGS=None, LGE=0, LGW=None, bg=0):
    if x0 is None:
        x0 = XY[0][0,XY[0].shape[1]//2]
    if y0 is None:
        y0 = XY[1][XY[1].shape[0]//2,0]
    if sigS is None:
        sigS = sigN
    if sigW is None:
        sigW = sigE
    if LGS is None:
        LGS = LGN
    if LGW is None:
        LGW = LGE
        
    X1 = (XY[0]-x0)*np.cos(angle) - (XY[1]-y0)*np.sin(angle)
    Y1 = (XY[0]-x0)*np.sin(angle) + (XY[1]-y0)*np.cos(angle)
    
    R1 = np.sqrt(X1**2+Y1**2)
    angle = np.arctan2(Y1, X1)
    sig = asymm_ellipse(sigW, sigE, sigN, sigS, angle)
    gamma = np.sqrt(2*np.log(2))*sig # HFHM
    LG = asymm_ellipse(LGW+1, LGE+1, LGN+1, LGS+1, angle)-1
    
    Gxy = Gauss(R1, 0, sig, 1)
    Lxy = 1/((R1/gamma)**2+1)
    
    f = (1-LG)*Gxy+LG*Lxy
    out = bg+amplitude*f
    return out

def MaxwellBoltzmann(E,T):
    from . import constants as const
    return 2*const.qe*np.sqrt(E/np.pi)*np.exp(-E/(const.kb*T))/(const.kb*T)**1.5
    
def Voigt(x, x0, sig, gamma, A=1):
    import scipy
    L = Lorentz(x, x0, gamma, A=1)
    G = Gauss(x, x0, sig)
    out = scipy.signal.convolve(L, G, 'same')
    out /= np.max(out)
    return A*out
