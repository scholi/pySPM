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
        r += args[3*i]*math.CDF(x,args[3*i+1],args[3*i+2])
    return r