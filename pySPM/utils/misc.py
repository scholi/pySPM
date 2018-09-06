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
