import numpy as np
import scipy.signal.signaltools as sig
from scipy.signal import convolve
import numpy as np
import scipy.signal.signaltools as sig
                           
def psf(img, sx, sy=None, angle=0):
    """
    Return a Gaussian PSF of the same size as img.

    img: image (reference for the output size)
    sx: sigma value for the long axis
    sy: sigma value for the short axis. If None take the same value as sx [default]
    angle: geometric angle (in radian) of the long axis. [default: 0]
    """
    from .math import Gauss
    if sy is None:
        sy = sx
    x = np.arange(img.shape[1])
    y = np.arange(img.shape[0])
    X, Y = np.meshgrid(x,y)
    X -= img.shape[1]//2
    Y -= img.shape[0]//2
    if angle != 0:
        Xp = X*np.cos(angle) - Y*np.sin(angle)
        Yp = X*np.sin(angle) + Y*np.cos(angle)
    else:
        Xp = X
        Yp = Y
    return Gauss(Xp, 0, sx)*Gauss(Yp, 0, sy)

def strictly_positify(x):
    """
    Make the result strictly positive ba setting the minimum value to
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
        
def _rl(x, image, psf):
    """
    Richardson-Lucy core algorithm
    Reference: L. B. Lucy / The Astronomical Journal / vol. 79 / No. 6 / June 1974 / pp. 745-754

    By giving an estimate x_k this function returns the next estimate x_{k+1}.
    
    x: x_k estimate
    image: input image to enhance
    psf: point spread functional
    """
    from scipy.signal import convolve
    I = strictly_positify(convolve(x, psf, mode='same')) # reconvoluted estimation.
    relative_blur = image / I
    return x * convolve(relative_blur, psf[::-1,::-1], mode='same') # Correlation is the convolution of mirrored psf

def _rl_fft(x, image, psf):
    from scipy.signal import fftconvolve
    I = strictly_positify(fftconvolve(x, psf, mode='same')) # reconvoluted estimation.
    relative_blur = image / I
    return x * fftconvolve(relative_blur, psf[::-1,::-1], mode='same') # Correlation is the convolution of mirrored psf

def _rl_accelerate(x, x1, x2, g1=None, g2=None, order=1):
    """
    Accelerated Richardson-Lucy algorithm.
    Reference: David S. C. Biggs and Mark Andrews, Appl. Opt./ Vol. 36 / No. 8 / 10 March 1997 / pp. 1766-1775
    
    Notation in reference to paper:
    x  = x_k
    x1 = x_{k-1}
    x2 = x_{k_2}
    g1 = g_{k-1}
    g2 = g_{k-2}
    y  = y_k
    """
    
    if g2 is None:
        alpha = 0 # Initialization
    else:
        alpha = np.sum(g1*g2)/strictly_positify(np.sum(g2**2)) # Eq. 10
        alpha = clip01(alpha) # be sure α∈[0,1]
    if alpha == 0:
        return x # the prediction is the same as x (initialization)
    h1 = x - x1 # Eq. 7
    y = x + alpha * h1 # Eq. 6
    if order>1:
        h2 = x - 2*x1 + x2 # Eq. 17
        y += h2 * alpha**2 / 2 # Eq. 14
    return y

def _get_fshape_slice(image, psf):
    """This is necessary for the fast Richardson-Lucy Algorithm"""
    s1 = np.array(image.shape)
    s2 = np.array(psf.shape)
    assert (s1 >= s2).all()
    shape = s1 + s2 - 1
    # Speed up FFT by padding to optimal size for FFTPACK
    fshape = [sig.fftpack.helper.next_fast_len(int(d)) for d in shape]
    fslice = tuple([slice(0, int(sz)) for sz in shape])
    return fshape, fslice
    
def _rl_fast(x, image, otf, iotf, fshape, fslice, **kargs):
    """
    fast algorithm for R-L.
    """
    import scipy.signal.signaltools as sig
    I = sig._centered(strictly_positify(np.fft.irfftn(np.fft.rfftn(x, fshape) * otf, fshape)[fslice]), image.shape) # reconvoluted estimation.
    im_ratio = image / I
    return x *  sig._centered(np.fft.irfftn(np.fft.rfftn(im_ratio, fshape) * iotf, fshape)[fslice], image.shape)
    
def _rl_damped(x, image, psf, T=1, N=3):
    """
    Core for damped Richardson-Lucy
    Only first draft here. It's not yet functional.
    """
    from scipy.signal import convolve
    I = convolve(x, psf, mode='same') # reconvoluted estimation.
    ratio = I / strictly_positify(image)
    logarithm = np.log(strictly_positify(ratio))
    U = -2*(image*logarithm-I+image)/T**2
    Ut = clip01(U)
    relative_blur = 1 + (Ut**(N-1)) * (N-(N-1)*Ut) * (image-I)/strictly_positify(I)
    return x * convolve(relative_blur, psf, mode='same')
    
def richardson_lucy(image, psf, iterations, T=0, N=10,
        core='default', acceleration=1, init='mean'):
    """
    Richardson-Lucy algorithm
    image: the image to enhance (numpy 2d array)
    psf: the Point Spread Function (numpy 2d array)
    iterations:
        The number of iterations to perform. It can be either an integer or a list of integer.
        For the later case, the returned solution is a dictionary with keys K and value being the enhancement after K interations.
    T: Damping factor ( to be used with core='damped' )
    N: N factor used with core='damped'
    core:
        default: default R-L algorithm using convolve from scipy.signal
        fft: performs a fftconvolution
        damped: perform a damped R-L (beta)
    acceleration:
        0: (No acceleration. standard R-L)
        1: First order acceleration     
        2: Second order acceleration
        higher orders are not yet implemented
    init:
        'mean': the default. The start value for x is the mean of the image
        'image': the start value x is the image itself
        numpy array: if init is a 2d numpy array, its value will be used as init value for x
        
    TODO: implement the damped algorithm
    """

    rl_args = dict(image=image, psf=psf)
    if core is 'default':
        core = _rl
    elif core is 'fft':
        core = _rl_fft
    elif core is 'damped':
        core = _rl_damped
        rl_args['T'] = T
        rl_args['N'] = N
    elif core is 'fast':
        core = _rl_fast
        fshape, fslice = _get_fshape_slice(image, psf)
        otf = np.fft.rfftn(psf, fshape)
        iotf = np.fft.rfftn(psf[::-1,::-1], fshape)
        rl_args.update(dict(otf=otf, iotf=iotf, fshape=fshape, fslice=fslice))
    else:
        raise ValueError("Invalid core value. See help for detail.")
        
    image = image.astype(np.float)
    psf = psf.astype(np.float)
    psf /= np.sum(psf) # Normalize the psf ⇒ ∫∫ psf(x,y) dx dy = 1
    
    if init is 'mean':
        x = 0.5 * np.ones(image.shape)
    elif init is 'image':
        x = image
    else:
        x = init
            
    # As the image and the psf are both ≥ 0, if the image value is 0, then the result should also be 0 at this position
    #im_deconv[image == 0] = 0

    # Is iterations a number of a list of number?
    dict_output = True
    if type(iterations) is int:
        dict_output = False
        iterations = [iterations]
  
    N = max(iterations)
        
    results = {}
    x1 = x2 = None
    g1 = g2 = None
    
    for i in range(N):
        if acceleration:
            y = _rl_accelerate(x, x1, x2, g1, g2, order=acceleration)
        else:
            y = x
        x_new = core(positify(y), **rl_args)
        g2 = g1
        g1 = x_new - y
        x, x1, x2 = x_new, x, x1 # rotate elements for nect iteration
        if i+1 in iterations:
            results[i+1] = np.copy(x)
    if dict_output:
        return results
    return results[N]