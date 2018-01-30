def RL(x, image, psf):
    from scipy.signal import convolve
    I = convolve(x, psf, mode='same') # reconvoluted estimation.
    relative_blur = np.copy(image)
    relative_blur[I>0] = relative_blur[I>0]/ I[I>0]
    return x * convolve(relative_blur, psf[::-1,::-1], mode='same') # Correlation is the convolution of mirrored psf
            
def dampedRL(x, image, psf, T=1, N=3):
    from scipy.signal import convolve
    I = convolve(x, psf, mode='same') # reconvoluted estimation.
    ratio = np.zeros(image.shape)
    ratio[image>0] = I[image>0] / image[image>0]
    logarithm = np.zeros(image.shape)
    logarithm[ratio>0] = np.log(ratio[ratio>0])
    U = -2*(image*logarithm-I+image)/T**2
    Ut = np.minimum(U, 1)
    relative_blur = np.ones(image.shape)
    relative_blur[I>0] += (Ut[I>0]**(N-1)) * (N-(N-1)*Ut[I>0]) * (image[I>0]-I[I>0])/I[I>0]
    return x * convolve(relative_blur, psf, mode='same')
    
def richardson_lucy(image, psf, iterations, damped=0, N=10, **kargs):
    import scipy
    from scipy.signal import convolve
    
    image = image.astype(np.float)
    psf = psf.astype(np.float)
    psf /= np.sum(psf) # Normalize the psf ⇒ ∫∫ psf(x,y) dx dy = 1
    
    # im_deconv = 0.5 * np.ones(image.shape)
    im_deconv = kargs.get('x0', np.mean(image) * np.ones(image.shape))
            
    # As the image and the psf are both ≥ 0, if the image value is 0, then the result should also be 0 at this position
    #im_deconv[image == 0] = 0

    # Is iterations a number of a list of number?
    dict_output = True
    if type(iterations) is int:
        dict_output = False
        iterations = [iterations]
  
    N = max(iterations)
        
    results = {}

    for i in range(N):
        if damped == 0:            
            im_deconv = RL(im_deconv, image, psf)
        else:
            im_deconv = dampedRL(im_deconv, image, psf, T=damped, N=N)
        if i+1 in iterations:
            results[i+1] = np.copy(im_deconv)
    if dict_output:
        return results
    return results[N]

def accelerated_richardson_lucy(image, psf, iterations, **kargs):    
    def _RL(x):
        T = kargs.get('T', kargs.get('damped', 0))
        if T==0:   
            return RL(x, image, psf)
        return dampedRL(x, image, psf, N=kargs.get('N',3), T=T)
    
    image = image.astype(np.float)
    psf   = psf.astype(np.float)
    psf  /= np.sum(psf)
    
    x0 = kargs.get('x0', np.mean(image) * np.ones(image.shape))    
    x1 = kargs.get('x1', _RL(x0))

    h1 = x1 - x0
    y0 = x0
    y1 = x1
    g0 = _RL(y0) - y0
    g1 = _RL(y1) - y1
        
    dict_output = True
    if type(iterations) is int:
        dict_output = False
        iterations = [iterations]
  
    N = max(iterations)
        
    results = {}

    for i in range(2, N):
        x2 = y1 + g1 # x_{k+1} = y_k + g_k
        h2 = x2 - x1
        a2 = np.sum(g1*g0)/np.sum(g0**2)
        y2 = x2 + a2*h2
        g2 = _RL(y2)-y2
        
        # Shift data
        x0, x1 = x1, x2
        y0, y1 = y1, y2
        g0, g1 = g1, g2
        h0, h1 = h1, h2
        if i+1 in iterations:
            results[i+1] = np.copy(x1)
    if dict_output:
        return results
    return results[N]