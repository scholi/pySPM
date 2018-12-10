import numpy as np
from warnings import warn

def dec_debug(debug):
    if debug>0:
        return debug-1
    if debug<0:
        return debug+1
    return False

def do_debug(debug):
    return debug==-1 or debug>0

def adaptive(img):
    import skimage
    mi = np.min(img)
    img = np.asarray(256**2*(img-mi)/(np.max(img)-mi), dtype=np.uint16)
    img = skimage.exposure.equalize_adapthist(img, clip_limit=0.03)
    return img

def smiley(width, height=None, ratio=.9, eye=.15, thick=.1, eye_sep=.35, eye_height=.3, mouth_rad=.5, mouth_thick=.1):
    if height is None:
        height = width
    size = min(width, height)
    Y, X = np.mgrid[-height//2:height//2, -width//2:width//2]
    R = np.sqrt(X**2+Y**2)/(size/2)
    smiley = (R<ratio)*(R>((ratio-thick)))
    smiley += np.sqrt((X-eye_sep*ratio*size/2)**2+(Y+eye_height*ratio*size/2)**2)<ratio*eye*size/2
    smiley += np.sqrt((X+eye_sep*ratio*size/2)**2+(Y+eye_height*ratio*size/2)**2)<ratio*eye*size/2
    smiley += (R<ratio*mouth_rad)*(R>(ratio*mouth_rad-ratio*mouth_thick))*(Y>0)
    return smiley*1.0
import functools

def aliased(cls):
    original_methods = cls.__dict__.copy()
    for name in original_methods:
        method = original_methods[name]
        if hasattr(method, '_aliases'):
            for alias in method._aliases:
                setattr(cls, *alias)
    return cls

class alias(object):
    def __init__(self, *aliases):
       self.aliases = aliases
       
    def __call__(self, func):
        if not hasattr(func, '_aliases'):
            func._aliases = []
        func._aliases += [(x, func) for x in self.aliases]
        import sys
        module = sys.modules[func.__module__]
        for alias in self.aliases:
            setattr(module, alias, func)
        return func

class deprecated(object):
    def __init__(self, alias):
        self.alias = alias
        
    def __call__(self, func):
        import sys
        import inspect
        msg = "the function {} is DEPRECATED. Please use {} instead".format(self.alias, func.__name__)
        module = sys.modules[func.__module__]
        @functools.wraps(func)
        def wrapper(*args, **kargs):
            from warnings import warn
            warn(msg)
            return func(*args, **kargs)
        setattr(module, self.alias, wrapper)
        if not hasattr(func, '_aliases'):
            func._aliases = []
        func._aliases += [(self.alias, wrapper)]
        return func

@alias("get_bam", "get_BAM")
def getBAM(x, x0, N=10, least_one=False):
    """
    get the profile of the GaAlAs layer from a BAM-L200 sample.
    
    Parameters
    ----------
    x: numpy.ndarray
        x-axis in nm
    x0: float
        offset of the x-axis in nm (ie. where the first edge of the BAM is located)
    N: int
        The number of headings grating (seems to be 13 on the SEM image in the doc, but 10 in reality?)
    least_one: bool
        if True each GaAlAs strip are at least 1px wide.
    """
    P = x*0
    l = 0
    
    # Here are the certified known distances
    d = dict(W1=691, W2=691, W3=293, W4=294, W5=19.5, W6=195, W7=195, W8=38, W9=3.6, W10=14.2, W11=3.5, W12=96, W13=5, W14=1, P1=587,P2=389,P3=273,P4=193,P5=136,P6=97,P7=67.5,P8=48.5,P9=76.5,P10=57,P11=42,P12=31,P13=23,P14=17.5,P15=13.3,P16=9.4,P17=6.9,P18=4.6,P19=3,P20=2)
    # The prelines
    lines = [(-147*i, 80) for i in range(N,0,-1)]
        
    # Estimated coordinates of the starting of lines
    pos = dict(W1=0, P1=d['W1']+d['W2'], W5=2683, P9=8760, W8=7300, W9=7800, W10=8250, W11=9500, W12=9840, W14=8100)
    pos_rel = dict(W5=dict(P2=397, P3=1380, P4=2140, P5=2730, P6=3200, P7=3535, P8=3730),P9=dict(P10=231, P11=402, P12=539, P13=642, P14=710),W11=dict(P15=60, P16=104, P17=141, P18=162, P19=182, P20=191))
    for rel in pos_rel:
        for key in pos_rel[rel]:
            pos[key] = pos[rel]+pos_rel[rel][key]

    for key in sorted(pos):
        width = [d[key],d[key]//2][key[0] == 'P']
        lines.append((pos[key], width))
        if key[0] == 'P': lines.append((pos[key]+d[key], width))        
            
    for pos, width in lines:
        P[((x-x0)>=pos)*((x-x0)<=pos+width)] = 1
        if least_one: P[np.argmin(abs(x-x0-pos))] = 1 # At least 1 pixel set ?
    return P

def in_ipynb():
    try:
        cfg = get_ipython().config 
        if 'IPKernelApp' in cfg:
            return True
        else:
            return False
    except NameError:
        return False
        
if in_ipynb():
    try:
        from tqdm import tqdm_notebook as tqdm
    except:
        try:
            from tqdm import tqdm
        except:
            warn("the library tqdm cannot be found. All progressbar will be disabled.")
            tqdm = lambda x: x
    PB = tqdm
else:
    try:
        from tqdm import tqdm
    except:
        warn("the library tqdm cannot be found. All progressbar will be disabled.")
        tqdm = lambda x: x
    PB = tqdm