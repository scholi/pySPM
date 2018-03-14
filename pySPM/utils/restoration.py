import numpy as np
import scipy.signal.signaltools as sig
from scipy.signal import convolve
import numpy as np
import scipy.signal.signaltools as sig
import copy

def sign(abs_var, sign_var):
    return abs(abs_var) * (1 - np.where(sign_var < 0, 2*sign_var, sign_var))
    
def hfilter(diff_image, var_image, threshold=1, ndamp=10):
    him = htrans(diff_image).astype(np.float)
    dvarim = htrans(var_image, var_image=True).astype(np.float)
    
    sqhim = ((him/threshold)**2)/dvarim
    index = np.where(sqhim < 1)
    
    if len(index[0]) == 0:
        return diff_image
    
    sqhim = sqhim[index] * (ndamp * sqhim[index]**(ndamp-1) - (ndamp-1)*sqhim[index]**ndamp)
    
    him[index] = sign(threshold*np.sqrt(dvarim[index] * sqhim), him[index])
    
    return hinv(him)
    
def htrans(image, var_image=False):    
    out_image = copy.deepcopy(image).astype(np.float)
    
    y_size, x_size = image.shape
        
    h0, hx, hy, hc = hreduce(out_image)

    xstop = x_size
    ystop = y_size
    
    xstart = (xstop+1) // 2
    ystart = (ystop+1) // 2
    
    if not var_image:
        #-- bottom right quadrant of image
        out_image[0:ystart, xstart:xstop] = hx
        #-- top left quadrant
        out_image[ystart:ystop, 0:xstart] = hy
        #-- top right quadrant
        out_image[ystart:ystop, xstart:xstop] = hc
    else:
        #-- Increase variance at edge if section size is odd
        #-- edge pixels are used in the sum twice
        if y_size%2:
            out_image[-1] *= 2
        if x_size%2:
            out_image[:, -1] *= 2
            
        #-- bottom right quadrant of image
        out_image[0:ystart, xstart:xstop] = h0[:, 0:xstop-xstart]
        #-- top left quadrant
        out_image[ystart:ystop, 0:xstart] = h0[0:ystop-ystart, :]
        #-- top right quadrant
        out_image[ystart:ystop, xstart:xstop] = h0[0:ystop-ystart, 0:xstop-xstart]
    
    while xstart > 1 or ystart > 1: 
        h0, hx, hy, hc = hreduce(h0)
        
        xstop = xstart 
        xstart = (xstart+1) // 2
        
        ystop = ystart 
        ystart = (ystart+1) // 2
        
        
        if not var_image:
            if xstop >= xstart:
                out_image[0:ystart, xstart:xstop] = hx/2
            if ystop >= ystart:
                out_image[ystart:ystop, 0:xstart] = hy/2
            if xstop >= xstart and ystop >= ystart:
                out_image[ystart:ystop, xstart:xstop] = hc/2
            h0 = (h0/2).reshape(ystart, xstart)
        else:
            #-- Increase variance at edge if section size is odd
            #-- edge pixels are used in the sum twice
            if ystart*2 != ystop:
                h0[-1] *= 2
            if xstart*2 != xstop:
                h0[:, -1] *= 2     
                
            if xstop >= xstart:
                out_image[0:ystart, xstart:xstop] = h0[:, 0:xstop-xstart]/4
            if ystop >= ystart:
                out_image[ystart:ystop, 0:xstart] = h0[0:ystop-ystart, :]/4
            if xstop >= xstart and ystop >= ystart:
                out_image[ystart:ystop, xstart:xstop] = h0[0:ystop-ystart, 0:xstop-xstart]/4
            h0 = h0/4  
            
    out_image[0, 0] = h0
    
    return out_image
    
def hinv(image):
    y_size, x_size = image.shape
    
    xstop = [x_size]
    ystop = [y_size]
    
    xstart = [(xstop[0]+1) // 2]
    ystart = [(ystop[0]+1) // 2]
    
    while xstart[0] > 1 or ystart[0] > 1:
        xstop = [xstart[0]] + xstop
        ystop = [ystart[0]] + ystop
        
        xstart = [(xstart[0] + 1)//2] + xstart
        ystart = [(ystart[0] + 1)//2] + ystart
        
    h = copy.deepcopy(image).astype(np.float)
    
    for i in range(len(ystart)-1):
        h0 = (h[0:ystart[i], 0:xstart[i]]).reshape(ystart[i], xstart[i])
        
        if xstop >= xstart:
            hx = (h[0:ystart[i], xstart[i]:xstop[i]]).reshape(ystart[i], xstop[i]-xstart[i])
        if ystop >= ystart:
            hy = (h[ystart[i]:ystop[i], 0:xstart[i]]).reshape(ystop[i] - ystart[i], xstart[i])
        if xstop >= xstart and ystop >= ystart:
            hc = (h[ystart[i]:ystop[i], xstart[i]:xstop[i]]).reshape(ystop[i]-ystart[i], xstop[i]-xstart[i])

        h[0:ystop[i], 0:xstop[i]] = hexpand(h0, hx, hy, hc) / 2
        
    #-- continue from last
    i += 1
    h0 = (h[0:ystart[i], 0:xstart[i]]).reshape(ystart[i], xstart[i])
    if xstop >= xstart:
        hx = (h[0:ystart[i], xstart[i]:xstop[i]]).reshape(ystart[i], xstop[i]-xstart[i])
    if ystop >= ystart:
        hy = (h[ystart[i]:ystop[i], 0:xstart[i]]).reshape(ystop[i] - ystart[i], xstart[i])
    if xstop >= xstart and ystop >= ystart:
        hc = (h[ystart[i]:ystop[i], xstart[i]:xstop[i]]).reshape(ystop[i]-ystart[i], xstop[i]-xstart[i])

    h[0:ystop[i], 0:xstop[i]] = hexpand(h0, hx, hy, hc) / 4
        
    return h
    
def hexpand(h0, hy, hx, hc):
    y_length, x_length = h0.shape
    
    #assert h0.shape == hy.shape == hx.shape == hc.shape, "{} {} {} {}".format(h0.shape, hy.shape, hx.shape, hc.shape)
    
    padded_x = False
    padded_y = False
    
    out_image = np.zeros((y_length*2, x_length*2)).astype(np.float)
    
    #if y_length%2:
    #    padded_y = True
    #    out_image = np.pad(out_image, ((0, 1), (0, 0)), mode='constant', constant_values=0)
    #if x_length%2:
    #    padded_x = True
    #    out_image = np.pad(out_image, ((0, 0), (0, 1)), mode='constant', constant_values=0)   
    
    if h0.shape != hy.shape:
        hy = np.pad(hy, ((0, 0),(0, 1)), mode='constant', constant_values=0)
        padded_y = True
    if h0.shape != hx.shape:
        hx = np.pad(hx, ((0, 1),(0, 0)), mode='constant', constant_values=0)
        padded_x = True
    if h0.shape != hc.shape:
        hc = np.pad(hc, ((0, 1),(0, 1)), mode='constant', constant_values=0)
        paddex_x = padded_y = True
    
    out_image[1::2, 1::2] = h0 + hx + hy + hc
    out_image[1::2, 0::2] = h0 + hx - hy - hc
    out_image[0::2, 1::2] = h0 - hx + hy - hc
    out_image[0::2, 0::2] = h0 - hx - hy + hc
    
    if padded_y and y_length > 1:
        out_image = out_image[:-1]
    if padded_x and x_length > 1:
        out_image = out_image[:, :-1]
        
    return out_image #out_image[:y_length*2, :x_length*2]
    
def hreduce(image):
    
    y_length, x_length = image.shape
    
    padded_x = False
    padded_y = False
    
    if y_length%2:
        padded_y = True
        #image = np.pad(image, ((0, 1), (0, 0)), mode='constant', constant_values=0)
        image = np.pad(image, ((0, 1), (0, 0)), mode='edge')
    if x_length%2:
        padded_x = True
        image = np.pad(image, ((0, 0), (0, 1)), mode='edge')      
        
    a11 = image[1::2, 1::2].astype(np.float)
    a10 = image[1::2, 0::2].astype(np.float)
    a01 = image[0::2, 1::2].astype(np.float)
    a00 = image[0::2, 0::2].astype(np.float)
    
    ### Each factor may be missing division by two according to internet
    ### switch hy and hx for aggreement with IDL version
    h0 = (a11 + a10 + a01 + a00)
    hy = (a11 + a10 - a01 - a00)
    hx = (a11 - a10 + a01 - a00)
    hc = (a11 - a10 - a01 + a00)
    
    if padded_y and y_length > 2:
        hy = hy[:-1]
        hc = hc[:-1]
    if padded_x and x_length > 2:
        hx = hx[:, :-1]
        hc = hc[:, :-1]
        
    return h0, hx, hy, hc
    
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
    
def _rl_accurate(x, image, psf):
    """
    Richardson-Lucy core algorithm
    Reference: L. B. Lucy / The Astronomical Journal / vol. 79 / No. 6 / June 1974 / pp. 745-754

    By giving an estimate x_k this function returns the next estimate x_{k+1}.
    
    x: x_k estimate
    image: input image to enhance
    psf: point spread functional
    """
    from scipy.signal import convolve2d as convolve
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
    
def _rl_damped_old(x, image, psf, T=1, N=3):
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
    
def _rl_damped(x, image, psf, gain=1, con_var=1, T=1, N=10, conv=None):
    """ Perform damped lucy richardson algorithm
    
    Parameters
    ----------
    gain: float, int
        CCD gain (relic?)
    con_var: float, int, np.ndarray
        Noise value or image
    threshold: float, int
        noise sigma threshold for dampening
    ndamp: float, int
        order of the dampening
    """
    
    if conv is None:
        from scipy.signal import convolve
        conv = convolve
    I = conv(x, psf, mode='same') # reconvoluted estimation.
    psf_mirror = psf[::-1, ::-1]
    rrr = image - I
    rrr = hfilter(rrr, (I+con_var)/gain, T, ndamp=N)
    rrr[np.isnan(rrr)] = 0
    ratio = gain*(1 + rrr / (I+con_var))
    return x * conv(ratio, psf_mirror, 'same')
    
def richardson_lucy(image, psf, iterations, T=0, N=10,
        core='default', acceleration=1, init='mean', clip=False):
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
    elif core is 'accurate':
        core = _rl_accurate
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
        x, x1, x2 = x_new, x, x1 # rotate elements for next iteration
        if clip:
            x[x<0] = 0
            x[x>clip] = clip
        if i+1 in iterations:
            results[i+1] = np.copy(x)
    if dict_output:
        return results
    return results[N]