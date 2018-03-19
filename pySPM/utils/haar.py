import numpy as np
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