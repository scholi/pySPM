import numpy as np

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

def getBAM(x, x0, N=10, edges=False):
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
    edges: boolean
        If False (default) return the binary profile (0 for GaAs and 1 for GaAlAs)
        If True returns the location of the edges as a binary vector (0 everywhere and 1 at the position of each edge). Note that at least a 1 will be set for each edge, but this might lead to wrong distances if the input x is given at low resolution.
    """
    P = x*0
    l = 0
    Edges = []
    
    # Here are the position and width of each Al stripe
    for xx in [(68,80)]*N+[(67,691),
            (695,293), (294,293),
            (380,19.5),
            (420,195), (195,195),
            (415,135), (135,135),
            (360,96), (96,96),
            (300,68), (68,68),
            (260,48), (49,48),
            (210,38), (39,38),
            (105,24), (24,24),
            (800,38),
            (494,3.6),
            (495,14.2)]:
        l += xx[0]
        Edges.append(x0+l)
        P = P + ((x-x0)>=l)*((x-x0)<=l+xx[1])
        if x0+l>=x[0] and x0+l<=x[-1]:
            P[np.argmin(abs(x-x0-l))] = 1 # At least 1 pixel set
        l += xx[1]
        Edges.append(x0+l)
    if edges:
        return P, Edges
    return P
        