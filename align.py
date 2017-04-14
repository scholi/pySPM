# -- coding: utf-8 --

import numpy as np
from skimage import transform as tf


class Aligner:

    def __init__(self, fixed, other, prog=True):
        self.fixed = fixed
        self.other = other
        self.size = fixed.shape

        self.trans = [0, 0]
        self.scale = [1, 1]
        self.rotation = 0

        self.initIDX = self.getMatchingIndex()

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

    def ImproveShift(self):
        tform = tf.AffineTransform(
            scale=self.scale, rotation=self.rotation, translation=(0, 0))
        O = tf.warp(self.other, tform, output_shape=self.fixed.shape)
        Corr = np.real(np.fft.fftshift(np.fft.ifft2(
            np.conj(np.fft.fft2(self.fixed))*np.fft.fft2(O))))
        cord = np.unravel_index(np.argmax(Corr), self.fixed.shape)
        self.trans = (cord[1]-self.size[0]/2, cord[0]-self.size[1]/2)

    def getTf(self):
        return tf.AffineTransform(scale=self.scale, rotation=self.rotation, translation=self.trans)

    def getMatchingIndex(self):
        return np.sum((self.fixed-tf.warp(self.other, self.getTf(), output_shape=self.fixed.shape))**2)

    def ImproveScaleX(self, fact=.1, count=0):
        old = self.scale
        IDX1 = self.getMatchingIndex()
        self.scale[0] *= 1+fact
        self.ImproveShift()
        IDX2 = self.getMatchingIndex()
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
