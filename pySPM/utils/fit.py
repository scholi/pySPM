# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Provides useful functions for fitting with scipy.optimize.curve_fit
"""

from . import math
import numpy as np

def CDF(x, bg, *args, **kargs):
    """
    Return the sum of several CDFs.
    x: values used to evaluate the function
    bg: background (value on the very left)
    args: suits of one or several Ai, xi, si describing each a CDF
        Ai: Amplitude
        xi: position
        si: sigma of CDF
                     _  
    ex: step edge: _| |_ between 0.2 and 0.7 in y and -1 and 1 in x with sigma=5
        CDF(x, 0.2, 0.5, -1, 5, -0.5, 1, 5)
    """
    if 'A' in kargs and 'x0' in kargs and 'sig' in kargs:
        args = [x for y in zip(kargs['A'], kargs['x0'], kargs['sig']) for x in y]
        
    r = bg * np.ones(x.shape)
    if len(args)%3 != 0:
        raise Exception("Invalid number of arguments. see help.")
        return
        
    for i in range(len(args)//3):
        r += args[3*i]*math.CDF(x, args[3*i+1], args[3*i+2])
    return r

def lgCDF(x, bg, lg=0, *args, **kargs):
    """
    Return the sum of several CDFs.
    x: values used to evaluate the function
    bg: background (value on the very left)
    lg: Lorentz-Gauss factor (0=Gauss, 1=Lorentz, in-between = mix)
    args: suits of one or several Ai, xi, si describing each a CDF
        Ai: Amplitude
        xi: position
        si: sigma of CDF
                     _  
    ex: step edge: _| |_ between 0.2 and 0.7 in y and -1 and 1 in x with sigma=5
        lgCDF(x, lg, 0.2, 0.5, -1, 5, -0.5, 1, 5)
    """
    
    if 'A' in kargs and 'x0' in kargs and 'sig' in kargs:
        args = [x for y in zip(kargs['A'], kargs['x0'], kargs['sig']) for x in y]
    
    r = bg * np.ones(x.shape)
    if len(args)%3 != 0:
        raise Exception("Invalid number of arguments. see help.")
        return
        
    for i in range(len(args)//3):
        r += args[3*i]*math.CDF(x, args[3*i+1], args[3*i+2], lg)
    return r

def lgCDF_fit(x,y, p0, dic=False):
    import scipy.optimize as opt
    n = (len(p0)-2)//3
    bounds=[
        [0,0]+[-np.inf,-np.inf,0]*n,
        [np.inf,1]+[np.inf,np.inf,np.inf]*n]
    popt, pcov = opt.curve_fit(lgCDF, x, y, p0, bounds=bounds)
    if dic:
        d = np.diag(pcov)
        return dict(bg=popt[0],lg=popt[1],A=popt[2::3],x0=popt[3::3],sig=popt[4::3]), dict(bg=d[0],lg=d[1],A=d[2::3],x0=d[3::3],sig=d[4::3])
    return popt, pcov

def CDF_fit(x,y, p0, dic=False):
    import scipy.optimize as opt
    n = (len(p0)-1)//3
    bounds=[
        [0]+[-np.inf,-np.inf,0]*n,
        [np.inf]+[np.inf,np.inf,np.inf]*n]
    popt, pcov = opt.curve_fit(CDF, x, y, p0, bounds=bounds)
    if dic:
        d = np.diag(pcov)
        return dict(bg=popt[0], A=popt[1::3],x0=popt[2::3],sig=popt[3::3]), dict(bg=d[0],A=d[2::3],x0=d[3::3],sig=d[4::3])
    return popt, pcov
    
def LG2D(A, Rweight=None, sigma=None, dic=False, **kargs):
    """
    Perform a 2D Lorentz-Gauss fitting on a 2D array A
    dic: if True return the solution as a dictionary
    
    output: list of popt and pcov array.
    The parameter order is:
    Amplitude, Angle, Sigma_x, Sigma_y, Center_x, Center_y, LGx, LGy
    """
    from scipy.optimize import curve_fit
    params = ['amplitude', 'angle', 'sig_x', 'sig_y', 'x0', 'y0', 'LG_x', 'LG_y','bg']
    def fit(XY, *args):
        # Add default parameters to the argument list at the right position
        j = 0
        a = []
        for x in params:
            if x in kargs:
                a.append(kargs[x])
            else:
                a.append(args[j])
                j += 1
        # Calculate the 2D Lorentz-Gauss
        return a[-1]+math.LG2D(XY, *a[:-1]).ravel()
        
    # First guess for parameters
    Amplitude = np.max(A)
    Center = np.unravel_index(np.argmax(A), A.shape)
    p0_def = [Amplitude, 0, 10, 10, Center[1], Center[0], 0, 0, 0]
    bounds_def = [
        [0     , np.radians(-45),       0,       0, -np.inf, -np.inf, 0, 0, 0],
        [np.inf, np.radians(45) , +np.inf, +np.inf,  np.inf,  np.inf, 1, 1, np.inf]
        ]
    p0 = [p0_def[i] for i,x in enumerate(params) if x not in kargs]
    bounds = [
        [bounds_def[0][i] for i,x in enumerate(params) if x not in kargs],
        [bounds_def[1][i] for i,x in enumerate(params) if x not in kargs]
        ]
    # Kick out arguments for the given default parameters

    x = np.arange(A.shape[1])
    y = np.arange(A.shape[0])
    X, Y = np.meshgrid(x, y)
    R = np.sqrt((X-Center[1])**2+(Y-Center[0])**2)
    if Rweight is not None:
        if Rweight is 'R':
            sigma = 0.001+R
        else:
            sigma = Rweight(R)
    if sigma is not None:
        sigma = np.ravel(sigma)
    pfit = curve_fit(fit, (X, Y), np.ravel(A), p0=p0, sigma = sigma, bounds=bounds)        
    
    # Re-add default parameters to the output
    pout = []
    dpout = []
    j = 0
    for i,x in enumerate(params):
        if x in kargs:
            pout.append(kargs[x])
            dpout.append(0)
        else:
            pout.append(pfit[0][j])
            dpout.append(np.sqrt(pfit[1][j,j]))
            j += 1
            
    if dic:
       return dict(zip(params, pout)), dict(zip(params, dpout))
    return pfit
