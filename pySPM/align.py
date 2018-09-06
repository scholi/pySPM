# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
This module is used in order to align two different images.
Usually one from an SPM and the other from the ToF-SIMS.

This module also gives the ability to perform shift correction on images
which is used in order to align the different scans from ToF-SIMS images.
"""

from __future__ import print_function
import numpy as np
from skimage import transform as tf
from scipy.ndimage.filters import gaussian_filter


class Aligner:
    def __init__(self, fixed, other, prog=True, FFT=True):
        self.fixed = np.copy(fixed)
        self.other = np.copy(other)
        self.size = fixed.shape
        self. FFT = FFT
        self.trans = [0, 0]
        self.scale = [1, 1]
        self.rotation = 0

        self.initIDX = self.getMatchingIndex()
    
    def compute(self, prog=False):
        # Compute the correlation to find the best shift as first guess
        if prog:
            print("Progress [1/4] Improve Shift...", end='\r')
        self.ImproveShift()
        if prog:
            print("Progress [2/6] Improve Scale X...", end='\r')
        self.ImproveScaleX()
        if prog:
            print("Progress [3/6] Improve Scale Y...", end='\r')
        self.ImproveScaleY()
        if prog:
            print("Progress [4/6] Improve Rotation...", end='\r')
        self.ImproveRotation()
        if prog:
            print("Progress [5/6] Improve Scale X...", end='\r')
        self.ImproveScaleX()
        if prog:
            print("Progress [6/6] Improve Scale Y...", end='\r')
        self.ImproveScaleY()
        if prog:
            print("", end='\r')

    def ImproveShift(self, verbose=False, **kargs):
        tform = tf.AffineTransform(
            scale=self.scale, rotation=self.rotation, translation=(0,0))
        
        O = tf.warp(self.other, tform, output_shape=self.other.shape, preserve_range=True)
        if self.FFT:
            Corr = np.real(np.fft.fftshift(np.fft.ifft2(
                np.conj(np.fft.fft2(self.fixed))*np.fft.fft2(O))))
            cord = np.unravel_index(np.argmax(Corr), self.fixed.shape)
            self.trans = [cord[1]-self.size[0]/2, cord[0]-self.size[1]/2]
        else:
            shift, D = AutoShift(self.fixed, O, shift=[-self.trans[0], self.trans[1]],**kargs)
            if verbose:
                print("ImproveShift",shift,D)
            self.trans = [-shift[0], shift[1]]

    def getTf(self, verbose=False):
        """
        Get the Afdfine transform.
        You can apply it to a pySPM Image (img) with: img.align(this.getTf())
        """
        if verbose:
            print("Transpose: {0[0]}, {0[1]}".format(self.trans))
        return tf.AffineTransform(scale=self.scale, rotation=self.rotation, translation=self.trans)

    def getMatchingIndex(self, power=1):
        img = tf.warp(self.other, self.getTf(), output_shape=self.other.shape, preserve_range=True)[:self.fixed.shape[0],:self.fixed.shape[1]]
        return np.sum(np.abs(self.fixed-img)**power)

    def ImproveScaleX(self, fact=.1, count=0, verbose=False):
        old = self.scale
        IDX1 = self.getMatchingIndex()
        self.scale[0] *= 1+fact
        self.ImproveShift()
        IDX2 = self.getMatchingIndex()
        print("ImproveX",count, IDX1, IDX2)
        if IDX2 > IDX1:
            self.scale[0] = old[0]*(1-fact)
            self.ImproveShift()
            IDX2 = self.getMatchingIndex()
            if IDX2 > IDX1:
                self.scale[0] = old[0]
                if fact > 1e-3:
                    self.ImproveScaleY(fact)
            elif count <= 11:
                self.ImproveScaleX(fact, count+1)
        elif count <= 11:
            self.ImproveScaleX(fact, count+1)

    def ImproveRotation(self, delta=.1, count=0, prog=False):
        IDX1 = self.getMatchingIndex()
        IDX2 = IDX1
        init = self.rotation
        if prog:
            print("Progress [{i}/4] Improve Rotation. Passes: {count} (max 11)".format(
                i=prog, count=count+1), end='\r')
        while IDX2 <= IDX1:
            IDX1 = IDX2
            count += 1
            self.rotation += delta
            self.ImproveShift()
            IDX2 = self.getMatchingIndex()
        self.rotation -= delta
        if delta > 0:
            self.ImproveRotation(-delta, count=count+1)
        elif abs(delta) > 1e-5:
            self.ImproveRotation(abs(delta)/10., count=count+1)

    def ImproveScaleY(self, fact=.1, count=0):
        old = self.scale
        IDX1 = self.getMatchingIndex()
        self.scale[1] *= 1+fact
        self.ImproveShift()
        IDX2 = self.getMatchingIndex()
        if IDX2 > IDX1:
            self.scale[1] = old[1]*(1-fact)
            self.ImproveShift()
            IDX2 = self.getMatchingIndex()
            if IDX2 > IDX1:
                self.scale[1] = old[1]
                if fact > 1e-3:
                    self.ImproveScaleX(fact/10)
            elif count <= 11:
                self.ImproveScaleY(fact, count+1)
        elif count <= 11:
            self.ImproveScaleY(fact, count+1)

    def __repr__(self):
        return u"Scale: ({scale[0]},{scale[1]})\nRotation: {rot:.6f} deg.\nTranslation: ({trans[0]},{trans[1]})".format(rot=self.rotation, scale=self.scale, trans=self.trans)

def ApplyShift(Img, shift):
    dx, dy = [-int(x) for x in shift]
    return np.pad(Img,((max(0,dy),max(0,-dy)),(max(0,-dx),max(0,dx))),
                mode='constant', constant_values=0)[max(0,-dy):max(0,-dy)+Img.shape[0],
                                                    max(0,dx):max(0,dx)+Img.shape[1]]
def ShiftScore(Ref, Img,  shift, gauss = 5, mean=True, norm=False, debug=False, normData=False):
    assert Ref.shape[0] <= Img.shape[0]
    assert Ref.shape[1] <= Img.shape[1]
    if mean:
        Ref = Ref - np.mean(Ref)
        Img = Img - np.mean(Img)
    if normData:
       assert np.std(Ref)>0
       assert np.std(Img)>0
       Ref /= np.std(Ref)
       Img /= np.std(Img)
    # blur the images for score calculation
    if gauss in [0,None,False]:
        im1 = Ref
        im2 = Img
    else:
        im1 = gaussian_filter( Ref, gauss)
        im2 = gaussian_filter( Img, gauss)
        
    corr2 = ApplyShift(im2, shift)
    dx, dy = shift
    # Create a copy of the reference and erase parts which are not overlaping with img
    DSX = Img.shape[1] - Ref.shape[1]
    DSY = Img.shape[0] - Ref.shape[0]
    Or = np.copy(im1)
    if dy < 0:
        Or[:-dy, :] = 0
    elif DSY-dy < 0:
        Or[DSY-dy:, :] = 0
    if dx > 0:
        Or[:, :dx] = 0
    elif DSX+dx < 0:
        Or[:, dx+DSX:] = 0
        
    corr2 = corr2[:Ref.shape[0],:Ref.shape[1]]
    # calculate the score: absolute of the differendces normed by the overlaping area
    D = np.sum(np.abs( Or - corr2 ))
    if norm:
        D /= ((Ref.shape[0]-2*dy)*(Ref.shape[1]-2*dx))
    if debug:
        return D, Or, corr2
    return D

def AutoShift(Ref, Img, Delta = 50, shift=(0,0), step=5, gauss=5, mean=True, test=False, norm=False, normData=False):
    """Function to find the best shift between two images by using brute force
    It shift the two iumages and calculate a difference score between the two.
    The function will return the shift which gives the lowerst score (least difference)
    The score is the norm of the difference between the two images where all non-overlaping parts of the images
    due to the shifts are set to 0. The score is then normes by the effective area.
    In order to avoid the errors due to shot-noise, the images are gaussian blured.
  
    Delta: shifts between shift[0/1]-Delta and shift[0/1]+Delta will be tested
    step: The step between tested delta values
    gauss: For noisy image it is better to use a gaussian filter in order to improve the score accuracy.
           The value is the gaussian size.
    mean: If True, the mean value of each image are subtracted. This is prefered when the intensities of the two images does not match perfectly.
          Set it to False if you know that the intensities of your two images are identical
    Note: This function was developed as the maximum in FFT correlation does not seems to give very acurate
          result for images with low counts. If the shifts is expected to be big, the first guess shift can be calculated
          by FFT correlation. ie.:
          s = np.fft.fftshift( np.abs( np.fft.ifft2( np.fft.fft2(Reference) * np.conj(np.fft.fft2(Image)))))
          shift = [x-s.shape[i]/2 for i,x in enumerate(np.unravel_index(np.argmax(s), s.shape))]
    """
    assert Ref.shape[0] <= Img.shape[0]
    assert Ref.shape[1] <= Img.shape[1]
    if mean:
        Ref = Ref - np.mean(Ref)
        Img = Img - np.mean(Img)
    if normData:
       assert np.std(Ref)>0
       assert np.std(Img)>0
       Ref /= np.std(Ref)
       Img /= np.std(Img)
    # blur the images for score calculation
    if gauss in [0,None,False]:
        im1 = Ref
        im2 = Img
    else:
        im1 = gaussian_filter( Ref, gauss)
        im2 = gaussian_filter( Img, gauss)
        
    # the following two variables save the best score
    best = (0,0)
    Dbest = Ref.shape[0]*Ref.shape[1]*max(np.max(im2),np.max(im1))
    
    tested = np.zeros((int(2*Delta/step)+1,int(2*Delta/step)+1))
    # Sweep through all possible shifts (brute force)
    for iy,Dy in enumerate(np.arange(shift[1]-Delta, shift[1]+Delta+1, step)):
        dy = int(Dy)
        for ix,Dx in enumerate(np.arange(shift[0]-Delta, shift[0]+Delta+1, step)):
            dx = int(Dx)
            corr2 = ApplyShift(im2, (dx,dy))
            DSX = Img.shape[1] - Ref.shape[1]
            DSY = Img.shape[0] - Ref.shape[0]
            # Create a copy of the reference and erase parts which are not overlaping with img
            Or = np.copy(im1)
            if dy < 0:
                Or[:-dy, :] = 0
            elif DSY-dy < 0:
                Or[DSY-dy:, :] = 0
            if dx > 0:
                Or[:, :dx] = 0
            elif DSX+dx < 0:
                Or[:, dx+DSX:] = 0
            corr2 = corr2[:Ref.shape[0],:Ref.shape[1]]
            # calculate the score: absolute of the differendces normed by the overlaping area
            D = np.sum(np.abs( Or - corr2 ))
            if norm:
                D /= ((Ref.shape[0]-2*dy)*(Ref.shape[1]-2*dx))
            if test:
                tested[iy,ix] = D
            if D < Dbest:
                Dbest = D
                best = (dx,dy)
    if test:
        return best, Dbest, tested
    return best, Dbest
