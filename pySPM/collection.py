# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
collection is a module to handle collection of images.
This is specially useful for SPM data which store several channels for the same measurement
"""

import copy
import re
import numpy as np
import matplotlib.pyplot as plt
from .SPM import SPM_image
from . import SPM


def atoi(text):
    """
    Convert string to integer
    """
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    """
    return [atoi(c) for c in re.split('(\\d+)', text)]


class Collection:
    """
    Class to handle a collection of SPM images
    """

    def __init__(self, sx=None, sy=None, unit='px', name='RawData', cls=None):
        """
        Create a new collection.
        You should provide a size.
        sx for the x-axis
        sy for the y-axis
        and unit for the name of the dimensional units
        """

        if isinstance(cls, Collection):
            self.size = cls.size
            self.name = cls.name

        if sy is None and sx is not None:
            sy = sx

        self.name = name
        self.channels = {}
        self.size = {'x': sx, 'y': sy, 'unit': unit}

    def add(self, Img, name, force=False):
        """
        Add a new Img to the collection
        """

        if len(self.channels) == 0 and self.size['unit'] == 'px':
            if isinstance(Img, SPM_image):
                self.size = {'x': Img.pixels.shape[
                    1], 'y': Img.pixels.shape[0]}
            else:
                self.size['x'] = len(Img[0])
                self.size['y'] = len(Img)
        if name in self.channels:
            if force:
                del self.channels[name]
            else:
                raise KeyError('The channel {} is already present in '
                           'the collection. Please delete it before')
                return
        self.channels[name] = Img
        
    def create_image(self, img, key=None):
        return SPM_image(img, _type=self.name, real=self.size, channel=key)

    def __len__(self):
        return len(self.channels)
        
    def __getitem__(self, key):
        if key not in self.channels:
            return None
        if isinstance(self.channels[key], SPM_image):
            return self.channels[key]
        return self.create_image(self.channels[key], key=key)
        
    def __delitem__(self, key):
        del self.channels[key]

    def __iter__(self):
        self.iterator = iter(self.channels)
        return self

    def __next__(self):
        it = self.iterator.__next__()
        return self.__getitem__(it)

    def __setitem__(self, key, value):
        self.add(value, key, force=True)

    def show(self, ax=None, channels=None, cmap='hot', ncols=3, width=20, **kargs):
        """
        Display the images of the collection in a matplotlib plot.

        ax: is the axis of the matplotib plot.
            If None is provided, the current axis will be retrieved (or new one)
        channels: The channels to plot.
            If None, all will be plotted
        cmap: The colormap to use (hot)
        **kargs: further arguments will be passed to the show function of the generated SPM_Image
        """
        if channels is None:
            channels = list(self.channels.keys())
        channels_number = len(channels)
        assert channels_number > 0
        channels.sort(key=natural_keys)
        if ax is None:
            if channels_number == 4 and ncols==3:
                fig, ax = plt.subplots(2, 2, figsize=(width,
                                                      self[channels[0]].pixels.shape[
                                                          0]*width
                                                      / self[channels[0]].pixels.shape[1]))
            else:
                Ny = (channels_number-1)//ncols+1
                Nx = min(ncols, channels_number)
                fig, ax = plt.subplots(Ny, Nx,
                                       figsize=(width, ((channels_number-1)//ncols+1)*width/Nx))
        if type(ax) is not list:
            ax = np.array(ax).ravel()
        for i, x in enumerate(channels):
            self[x].show(ax=ax[i], cmap=cmap, **kargs)
        plt.tight_layout()

    def get_multivariate(self, channels=None):
        """
        Create and return a (pandas) DataFrame of the collection
        channels: List of the channels to use (default: None => all channels)

        Note: The images will be unraveled (flattened)
        """
        import pandas as pd
        if channels is None:
            channels = self.channels.keys()
        if isinstance(list(self.channels.values())[0], SPM_image):
            data = pd.DataFrame(
                {k: self.channels[k].pixels.ravel() for k in channels})
        else:
            data = pd.DataFrame(
                {k: self.channels[k].ravel() for k in channels})
        return data

    def overlay(self, channel_names, colors=[[1, 0, 0], [0,1,0],[0,0,1]],
                sig=None, vmin=None, vmax=None, **kargs):
        """
        Create an overlay (in color) of several channels.
        channel_names: a list of the channels to select for the overlay
        colors: List of [Red,Green,Blue] color to use.
            (Its length should be identical to channel_names)
        """
        assert len(channel_names) >= 2
        assert len(colors) >= len(channel_names)
        data = [SPM.normalize(
            self[ch].pixels, sig=sig, vmin=vmin, vmax=vmax)
            for ch in channel_names]
        layers = [np.stack([data[i]*colors[i][j] for j in range(3)], axis=2)
                  for i in range(len(channel_names))]
        overlay = np.clip(np.sum(layers, axis=0), 0, 1)
        o = self.create_image(overlay, key="overlay")
        ch = [self.create_image(x, key=self[channel_names[i]].channel) for i,x in enumerate(layers)]
        if 'ax' in kargs:
            o.show(**kargs)
        return o, ch

    def stitch_correction(self, channel, stitches):
        """
        Function to correct for anomalies seen in stitched image.
        The function will calculate an average distribution of the stitched field
        and average the image with it and return a new collection with the result

        channel: name of the channel used as a reference
            (take one with homogeneous intensities in the whole image)
        stitches: a tuple/list containing the number of stitches in the image (x,y)
        """
        result = copy.deepcopy(self)
        del result.channels
        result.channels = {}
        size = list(self.channels.values())[0]
        S = np.zeros((size[0]/stitches[0], size[1]/stitches[1]))
        for i in range(stitches[0]):
            for j in range(stitches[1]):
                S += self.channels[channel][sy*i:sy*(i+1), sx*j:sx*(j+1)]
        S[S == 0] = 1
        sy, sx = S.shape
        for x in self.channels:
            F = np.zeros(size)
            for i in range(stitches[0]):
                for j in range(stitches[1]):
                    F[sy*i:sy*(i+1), sx*j:sx*(j+1)] \
                        = self.channels[x][sy*i:sy*(i+1), sx*j:sx*(j+1)]/S
            result.add(F, x)
        return result


def __Tsign(p1, p2, p3):
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])


def PointInTriangle(pt, v1, v2, v3):
    """
    Is a point pt in the triangle formed by the vertices v1, v2 and v3?
    pt,v1,v2,v3: tuple/list contai9ng the (x,y) coordinate of each vertices/point
    """
    b1 = __Tsign(pt, v1, v2) < 0
    b2 = __Tsign(pt, v2, v3) < 0
    b3 = __Tsign(pt, v3, v1) < 0
    return (b1 == b2) * (b2 == b3)


def overlay_triangle(channel_names, colors=[[1,0,0],[0,1,0],[0,0,1]], radius=.8,proportion = .8,ax=None, size=512, bgcolor=[0,0,0], textcolor='white',fontsize=20):
    """
    Create the image of a triangle, where the color is a gradient between three colors
    (one for each vertex).

    channel_names: a list of 3 channels
    colors: a list of 3 [R,G,B] list color
    ax: the matplotlib axis to plot in (if none use plt.gca())
    size: size of the rastered image generated.
        As there are in RGB no more than 256 values,
        the default value of 512 should be good in most of the cases.
    """
    assert len(channel_names) == 3
    assert len(colors) == 3
    if ax is None:
        ax = plt.gca()
    RGB = [bgcolor[i]*np.ones((size, size)) for i in range(3)]
    distance = 2*radius*proportion*np.sin(np.radians(120))

    x = np.linspace(-.7, 1.1, size)
    y = np.linspace(-1, 1, size)
    X, Y = np.meshgrid(x, y)
    centers = [(radius*proportion*np.cos(np.radians(120*i)), radius *
                proportion*np.sin(np.radians(120*i))) for i in range(3)]
    dist = [np.sqrt((X-centers[i][0])**2+(Y-centers[i][1])**2)
            for i in range(3)]

    for j in range(3):
        RGB[j] = sum(
            [colors[i][j]*np.maximum((distance-dist[i])/distance, 0) for i in range(3)])
        ax.annotate(channel_names[j], (radius*np.cos(np.radians(120*j)), radius*np.sin(np.radians(120*j))),
                    color=textcolor,
                    fontsize=fontsize,
                    va="center",
                    ha="center")
        RGB[j][PointInTriangle([X, Y], *centers) == 0] = bgcolor[j]
    image = np.stack(RGB, axis=2)
    ax.imshow(image, extent=[x[0], x[-1], y[-1], y[0]])
    ax.set_xticks([])
    ax.set_yticks([])
