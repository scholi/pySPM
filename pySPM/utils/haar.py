import numpy as np
import copy
import pywt


def sign(abs_var, sign_var):
    return abs(abs_var) * (1 - np.where(sign_var < 0, 2*sign_var, sign_var))
    
def hfilter(diff_image, var_image, threshold=1, ndamp=10):
    """
    This code was inspired from: https://github.com/spacetelescope/sprint_notebooks/blob/master/lucy_damped_haar.ipynb
    I believe it was initially written by Justin Ely: https://github.com/justincely
    It was buggy and not working properly with every image sizes.
    I have thus exchanged it by using pyWavelet (pywt) and a custom function htrans
    to calculate the matrix for the var_image.
    """
    him, coeff_slices = pywt.coeffs_to_array(pywt.wavedec2(diff_image.astype(np.float), 'haar'), padding=0)
    dvarim = htrans(var_image.astype(np.float))
    
    sqhim = ((him/threshold)**2)/dvarim
    index = np.where(sqhim < 1)
    
    if len(index[0]) == 0:
        return diff_image
    
    # Eq. 8 of White is derived leading to N*x^(N-1)-(N-1)*x^N  :DOI: 10.1117/12.176819
    sqhim = sqhim[index] * (ndamp * sqhim[index]**(ndamp-1) - (ndamp-1)*sqhim[index]**ndamp)
    him[index] = sign(threshold*np.sqrt(dvarim[index] * sqhim), him[index])
    
    return pywt.waverec2(pywt.array_to_coeffs(him, coeff_slices, output_format='wavedec2'), 'haar')[:diff_image.shape[0],:diff_image.shape[1]]

def htrans(A):
    h0 = A
    res = []
    while h0.shape[0]>1 and h0.shape[1]>1:
        h0, (hx, hy, hc) = pywt.dwt2(h0, 'haar')
        res = [(h0, h0, h0)]+res
    out, _ = pywt.coeffs_to_array([h0]+res, padding=1)
    return out
