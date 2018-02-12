# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Provides useful functions for fitting with scipy.optimize.curve_fit
"""

from . import math
import numpy as np

def CDF(x, bg, *args):
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
    
    r = bg * np.ones(x.shape)
    if len(args)%3 != 0:
        raise Exception("Invalid number of arguments. see help.")
        return
        
    for i in range(len(args)//3):
        r += args[3*i]*math.CDF(x, args[3*i+1], args[3*i+2])
    return r
    
def lgCDF(x, bg, lg=0, *args):
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
    
    r = bg * np.ones(x.shape)
    if len(args)%3 != 0:
        raise Exception("Invalid number of arguments. see help.")
        return
        
    for i in range(len(args)//3):
        r += args[3*i]*math.CDF(x, args[3*i+1], args[3*i+2], lg)
    return r
    
def LG2D(A, Rweight=None, sigma=None, dic=False):
    """
    Perform a 2D Lorentz-Gauss fitting on a 2D array A
    dic: if True return the solution as a dictionary
    
    output: list of popt and pcov array.
    The parameter order is:
    Amplitude, Angle, Sigma_x, Sigma_y, Center_x, Center_y, LGx, LGy
    """
    from scipy.optimize import curve_fit
    def fit(XY, *args):
        return math.LG2D(XY, *args).ravel()
        
    # First guess for parameters
    Amplitude = np.max(A)
    Center = np.unravel_index(np.argmax(A), A.shape)
    p0 = [Amplitude, 0, 10, 10, Center[1], Center[0], 0, 0]
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
    pfit = curve_fit(fit, (X, Y), np.ravel(A), p0=p0, sigma = sigma,
        bounds=[
        (0     , np.radians(-45),       0,       0, -np.inf, -np.inf, 0, 0),
        (np.inf, np.radians(45) , +np.inf, +np.inf,  np.inf,  np.inf, 1, 1)
        ])
    if dic:
       params = ['amplitude','angle','sig_x','sig_y','x0','y0','LG_x','LG_y']
       return dict(zip(params, pfit[0])), {x:np.sqrt(pfit[1][i,i]) for i,x in enumerate(params)}
    return pfit
