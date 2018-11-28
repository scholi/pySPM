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
    