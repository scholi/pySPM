# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Library to handle SPM data.
This is the core module of all images retrieved by SPM and ToF-SIMS.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage
import scipy.optimize
import skimage
import skimage.exposure
import skimage.filters
import scipy.interpolate
from skimage import transform as tf
import copy
from .utils import CDF
import sys

try:
    from tqdm import tqdm_notebook as tqdm
except:
    try:
        from tqdm import tqdm
    except:
        tqdm = lambda x: x
        
import matplotlib as mpl
import warnings

try:
    from skimage.filters import threshold_local
except:
    # For compatibility with old versions of skimage
    from skimage.filters import threshold_adaptive as threshold_local


def funit(value, unit=None):
    """
    Convert a value and unit to proper value
    e.g. 0.01m will be converted to 10mm
    e.g. 1000A will be converted to 1kA
    etc.

    Parameters
    ----------
    value : int, float or dictionary with 'value' and 'unit' key
        numerical value to work with
    unit : string or None
        name of the unit

    Returns
    -------
    dictionary with formated 'value' and 'unit'.
    The value will be: ≥1 and <1000

    Examples
    --------
    >>> funit(0.01,'m')
    {'value': 10.0, 'unit': 'mm'}
    >>> funit(1, 'cm')
    {'value': 1.0, 'unit': 'cm'}
    >>> funit({'value':2340, 'unit': 'um'})
    {'value': 2.34, 'unit': 'mm'}
    """
    if unit == None:
        unit = value['unit']
        value = value['value']
    import math
    shift = int(math.floor(math.log10(value)/3.0))  # base10 exponent
    mag = u'afpnum1kMGTPE'
    index_mag = mag.index('1')
    # Test if unit has a scale factor (n,m,k,M,G, etc.)
    if len(unit) > 1 and unit[0] in mag and unit[1:] and index_mag:
        index_mag = mag.index(unit[0])
        unit = unit[1:]
    value /= 10**(3*shift)
    index_mag += shift
    if index_mag < 0:
        value *= 10**(index_mag*3)
        index_mag = 0
    elif index_mag >= len(mag):
        value *= 10**((index_mag-len(mag)+1)*3)
        index_mag = len(mag)-1
    unit_prefix = mag[index_mag]
    if unit_prefix == '1':
        unit_prefix = ''
    return {'value': value, 'unit': u'{mag}{unit}'.format(mag=unit_prefix, unit=unit)}
   
class SPM_image:
    """
    Main class to handle SPM images.
    This class contains the pixels data of the images as well as it's real size.
    It also provides a lot of tools to correct and perform various analysis and tasks on the image.
    """
    
    def __init__(self, BIN, channel='Topography',
                 corr=None, real=None, zscale='?', _type='Unknown'):
        """
        Create a new SPM_image

        Parameters
        ----------
        BIN : 2D numpy array
            The pixel values of the image as a 2D numpy array
        channel : string
            The name of the channel. What does the image represents?
        corr : string or None
            'slope' : correct the SPM image for its slope (see pySPM.SPM.SPM_image.correct_slope)
            'lines' : correct the SPM image for its lines (see pySPM.SPM.SPM_image.correct_lines)
            'plane' : correct the SPM image by plane fitting (see pySPM.SPM.SPM_image.correct_plane)
        real : None or dictionary
            Information about the real size of the image {'x':width,'y':height,'unit':unit_name}
        zscale : string
            Unit used to describe the z-scale. (units of the data of BIN)
        _type : string
            represent the type of measurement
        """
        self.channel = channel
        self.direction = 'Unknown'
        self.size = {'pixels': {'x': BIN.shape[1], 'y': BIN.shape[0]}}
        if not real is None:
            self.size['real'] = real
        else:
            self.size['real'] = {'unit': 'pixels',
                                 'x': BIN.shape[1], 'y': BIN.shape[0]}
        if not 'unit' in self.size['real']:
            self.size['real']['unit'] = 'px'
        self.pixels = BIN
        self.type = _type
        self.zscale = zscale
        if corr is not None:
            if corr.lower() == 'slope':
                self.correct_slope()
            elif corr.lower() == 'lines':
                self.correct_lines()
            elif corr.lower() == 'plane':
                self.correct_plane()
            
    def __add__(self, b):
        """
        Add up two images. This is a low level function and no check is performed to proof that both images have the same size.
        """
        New = copy.deepcopy(self)
        New.pixels += b.pixels
        New.channel += " + "+b.channel
        return New
    
    def pxs(self):
        """
        Return the pixel size
        """
        fxy = {xy: funit(self.size['real'][xy], self.size['real']['unit']) for xy in 'xy'}
        return [(fxy[xy]['value']/self.size['pixels'][xy], fxy[xy]['unit']) for xy in 'xy']
        
    def add_scale(self, length, ax=None, height=20, color='w', loc=4, text=True, pixels=True, fontsize=20):
        """
        Display a scale marker on an existing image

        Parameters
        ----------

        length : float
            The length of the scale in real units
        ax : matplotlib axis
            if None the current axis will be taken (plt.gca())
        height : int
            The height of the scale bar in pixels
        color : string
            The color used to display the scale bar
        loc : int
            The location of the scale bar.
            1 : top right
            2 : top left
            3 : bottom left
            4 : bottom right
        text : bool
            display the size of the scale on top of it?
        pixels : bool
            Is the image plotted in ax with a x/y scale in pixels?
        fontsize : float
            The fontsize used to display the text

        Example
        -------
        >>> img = pySPM.SPM_image()
        >>> img.show()
        >>> img.add_scale(50e-6, pixels=False);
        Add a scale of 50 μm on an image displayed with real units

        >>> img = pySPM.SPM_image()
        >>> img.show(pixels=True)
        >>> img.add_scale(50e-6);
        Add a scale of 50 μm on an image displayed in pixels
        """
        import matplotlib.patches
        L = length*self.size['pixels']['x']/self.size['real']['x']
        ref = [height, height]
        if loc == 1 or loc == 4:
            ref[0] = self.size['pixels']['x'] - ref[0] - L
        if loc == 3 or loc == 4:
            ref[1] = self.size['pixels']['y'] - ref[1] - height
        if ax is None:
            ax = plt.gca()
        if text and loc in [1,2]:
            ref[1] += fontsize + height
        if not pixels:
            x,y = self.px2real(ref[0]+L,ref[1]+height)
            ref = self.px2real(*ref)
            L = x-ref[0]
            height = y-ref[1]
        ax.add_patch(matplotlib.patches.Rectangle(
            ref, width=L, height=height, color=color))
        if text:
            r = funit(length, self.size['real']['unit'])
            if r['unit'][0] == 'u':
                r['unit'] = '$\\mu$' + r['unit'][1:]
            ax.annotate("{value:.01f} {unit}".format(**r),
                        (ref[0]+L/2, ref[1]), color=color,
                        fontsize=fontsize, va="bottom", ha="center")

    def offset(self, profiles, width=1, ax=None, col='w', inline=True, **kargs):
        """
        Correct an image by offsetting each row individually in order that the lines passed as argument in "profiles" becomes flat.

        Parameters
        ----------
        profiles: list of list
            each sublist represent a line as [x1, y1, x2, y2] in pixels known to be flat
        width : int, float
            the line width in pixels used for better statistics
        ax : matplotlib axis or None
            If not None, axis in which the profiles will be plotted in
        inline : bool
            If True perform the correction on the current object, otherwise return a new image
        col : string
            matrplotlib color used to plot the profiles (if ax is not None)
        labels : bool
            display a label number with each profile
        **kargs: arguments passed further to get_row_profile.
            axPixels: set to True if you axis "ax" have the data plotted in pixel instead of real distance

        Example
        -------
        Exampel if the data are plotted in pixels:
        >>> topo = pySPM.SPM_image(...)
        >>> fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        >>> topoC = topo.offset([[150, 0, 220, 255]], inline=False,axPixels=True)
        >>> topo.show(pixels=True, ax=ax[0])
        >>> topoC.show(ax=ax[1]);

        Example if the data are plotted with real units
        >>> topo = pySPM.SPM_image(...)
        >>> fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        >>> topoC = topo.offset([[150, 0, 220, 255]], inline=False)
        >>> topo.show(ax=ax[0])
        >>> topoC.show(ax=ax[1]);
        """
        offset = np.zeros(self.pixels.shape[0])
        counts = np.zeros(self.pixels.shape[0])
        for i, p in enumerate(profiles):
            if kargs.get('labels', False):            
                y, D = self.get_row_profile(*p, width=width, ax=ax, col=col, label=str(i), **kargs)
            else:  
                y, D = self.get_row_profile(*p, width=width, ax=ax, col=col, **kargs)
            counts[y] += 1
            offset[y[1:]] += np.diff(D)
        counts[counts == 0] = 1
        offset = offset/counts
        offset = np.cumsum(offset)
        offset = offset.reshape((self.pixels.shape[0], 1))
        if inline:
            self.pixels = self.pixels - \
                np.flipud(np.repeat(offset, self.pixels.shape[1], axis=1))
            return self
        else:
            C = copy.deepcopy(self)
            C.pixels = self.pixels - \
                np.flipud(np.repeat(offset, self.pixels.shape[1], axis=1))
            return C
            
    def pxRect2Real(self, xy, width, height):
        """
        Transform a xy, width, height data in pixels to an equivalentz one with real units
        """
        ll = self.px2real(xy[0],xy[1])
        ur = self.px2real(xy[0]+width,xy[1]+height)
        return ll,ur[0]-ll[0],ur[1]-ll[1]
    
    def get_row_profile(self, x1, y1, x2, y2, width=1, col='C1', ax=None, alpha=0, **kargs):
        """
        Get a profile per row along a given line. This function is mainly useful for the function offset.
        
        x1, y1, x2, y2: int
            coordinates of the line.
        width : int
            the width of the line used for statistics (in pixels)
        col: string
            color used to plot the line position
        ax : matplotlib axis
            axis in which the lines position will plotted
        alpha : float
            The alpha channel of the line color (≥0 and ≤1)
        **kargs:
            line style arguments: linewidth, color and linestyle
            axis units: axPixels set to True if ax has the image plotted in pixels.
            
        Returns
        -------
        Y coordinates : 1D numpy array
            distance along the profile starting at 0
        Z coordinates : 1D numpy array
            profile
        """
        plotargs = { key: kargs[key] for key in ['linewidth', 'color', 'linestyle'] if key in kargs }
        if y2 < y1:
            x1, y1, x2, y2 = x2, y2, x1, y1
        if ax is not None:
            d = np.sqrt((x2-x1)**2+(y2-y1)**2)
            dx = -width/2*(y2-y1)/d
            dy = width/2*(x2-x1)/d
            if kargs.get('axPixels', False):
                ax.plot([x1-dx, x1+dx], [y1-dy, y1+dy], col)
                ax.plot([x2-dx, x2+dx], [y2-dy, y2+dy], col)
                ax.plot((x1, x2), (y1, y2), col, **plotargs)
                if kargs.get('label', False):
                    ax.annotate(kargs.get('label'), (.5*(x1+x2),.5*(y1+y2)), color=col)
                if alpha>0:
                    import matplotlib.patches
                    ax.add_patch(matplotlib.patches.Rectangle((x1+dx,y1+dy),width, d, -np.degrees(np.arctan2(x2-x1,y2-y1)), color=col, alpha=alpha))
            else:
                h = self.pixels.shape[0]
                pxs = self.size['real']['x'] / self.pixels.shape[1]
                pys = self.size['real']['y'] / h
                ax.plot([(x1-dx)*pxs, (x1+dx)*pxs], [(h-(y1-dy))*pys, (h-(y1+dy))*pys], col)
                ax.plot([(x2-dx)*pxs, (x2+dx)*pxs], [(h-(y2-dy))*pys, (h-(y2+dy))*pys], col)                
                ax.plot((x1*pxs, x2*pxs), ((h-y1)*pys, (h-y2)*pys), col, **plotargs)
                if kargs.get('label', False):
                    ax.annotate(kargs.get('label'), (.5*(x1+x2)*pxs,.5*(2*h-y1-y2)*pys), color=col)
                if alpha>0:
                    import matplotlib.patches
                    W = np.sqrt((2*dx*pxs)**2+(2*dy*pys)**2)
                    L = np.sqrt(((x2-x1)*pxs)**2+((y2-y1)*pys)**2)
                    ax.add_patch(matplotlib.patches.Rectangle(((x1+dx)*pxs,(y1+dy)*pys), W, L, -np.degrees(np.arctan2((x2-x1)*pxs,(y2-y1)*pys)), color=col, alpha=alpha))

        x = np.arange(self.pixels.shape[1])
        y = np.arange(self.pixels.shape[0])
        I = scipy.interpolate.interp2d(x, y, np.flipud(self.pixels))

        Y = np.arange(y1, y2+1)
        V = np.zeros(len(Y))
        for w in np.arange(width):
            xl = np.linspace(x1-(width-1)/2.+w, x2-(width-1)/2.+w, len(Y))
            for i in range(len(Y)):
                Z = I(xl[i], Y[i])
                V[i] += Z
        return Y, V/width

    def correct_median_diff(self, inline=True):
        """
        Correct the image with the median difference
        """
        N = self.pixels
        # Difference of the pixel between two consecutive row
        N2 = np.vstack([N[1:, :], N[-1:, :]])-N
        # Take the median of the difference and cumsum them
        C = np.cumsum(np.median(N2, axis=1))
        # Extend the vector to a matrix (row copy)
        D = np.tile(C, (N.shape[0], 1)).T
        if inline:
            self.pixels = N-D
        else:
            New = copy.deepcopy(self)
            New.pixels = N-D
            return New

    def correct_slope(self, inline=True):
        """
        Correct the image by subtracting a fitted slope along the y-axis
        """
        s = np.mean(self.pixels, axis=1)
        i = np.arange(len(s))
        fit = np.polyfit(i, s, 1)
        if inline:
            self.pixels -= np.tile(np.polyval(fit, i).reshape(len(i), 1), len(i))
            return self
        else:
            New = copy.deepcopy(self)
            New.pixels -= np.tile(np.polyval(fit, i).reshape(len(i), 1), len(i))
            return New

    def correct_plane(self, inline=True, mask=None):
        """
        Correct the image by subtracting a fitted 2D-plane on the data

        Parameters
        ----------
        inline : bool
            If True the data of the current image will be updated otherwise a new image is created
        mask : None or 2D numpy array
            If not None define on which pixels the data should be taken.
        """
        x = np.arange(self.pixels.shape[1])
        y = np.arange(self.pixels.shape[0])
        X0, Y0 = np.meshgrid(x, y)
        Z0 = self.pixels
        if mask is not None:
            X = X0[mask]
            Y = Y0[mask]
            Z = Z0[mask]
        else:
            X = X0
            Y = Y0
            Z = Z0
        A = np.column_stack((np.ones(Z.ravel().size), X.ravel(), Y.ravel()))
        c, resid, rank, sigma = np.linalg.lstsq(A, Z.ravel(), rcond=-1)
        if inline:
            self.pixels -= c[0] * \
                np.ones(self.pixels.shape) + c[1] * X0 + c[2] * Y0
            return self
        else:
            New = copy.deepcopy(self)
            New.pixels -= c[0]*np.ones(self.pixels.shape) + c[1] * X0 + c[2] * Y0
            return New

    def correct_lines(self, inline=True):
        """
        Subtract the average of each line for the image.

        if inline is True the current data are updated otherwise a new image with the corrected data is returned
        """
        if inline:
            self.pixels -= np.tile(np.mean(self.pixels, axis=1).T, (self.pixels.shape[0], 1)).T
            return self
        else:
            New = copy.deepcopy(self)
            New.pixels -= np.tile(np.mean(self.pixels, axis=1).T, (self.pixels.shape[0], 1)).T
            return New

    def dist_v2(self, pixel=False):
        """
        Return a 2D array with the distance between each pixel and the closest border.
        Might be usefull for FFT filtering
        """
        if pixel:
            dx = 1
            dy = 1
        else:
            dx = self.size['real']['x']/self.size['pixels']['x']
            dy = self.size['real']['y']/self.size['pixels']['y']
        x2 = np.arange(self.size['pixels']['x'])
        x2 = (np.minimum(x2, self.size['pixels']['x']-x2) * dx)**2
        y2 = np.arange(self.size['pixels']['y'])
        y2 = (np.minimum(y2, self.size['pixels']['y'] - y2) * dy)**2
        X, Y = np.meshgrid(x2, y2)
        return np.sqrt(X+Y)
    
    def inv_calc_flat(self, d, l=0.1):
        """
        Function used for inverse MFM calculation (inspired from http://qmfm.empa.ch/qmfm/)
        The function is in its early devlopment stage as not used by the developed.

        Parameters
        ----------
        d : float
            Height distance in the input data
        l : float
            Tikhonov parameter for the deconvolution
        """
        work_image = self.pixels
        ny, nx = self.pixels.shape
        dx = self.size['real']['x']/self.size['pixels']['x']
        dy = self.size['real']['y']/self.size['pixels']['y']

        k = self.dist_v2()
        k[0, 0] = 1e-10

        tf = np.exp(-d*k)
        tf[0, 0] = np.mean(tf)
        tf /= 2
        tf *= 1-np.exp(-d * k)

        recon_tf = np.ones(tf.shape) / (tf+l*np.ones(tf.shape) / np.conj(tf))
        tf *= recon_tf
        return np.real(np.fft.ifft2(np.fft.fft2(work_image)*recon_tf))

    def get_extent(self):
        """
        Get the image extent in real data
        """
        W = self.size['recorded']['real']['x']
        H = self.size['recorded']['real']['y']
        return (0, W, 0, H)

    def show(self, ax=None, sig=None, cmap=None, title=None,
             adaptive=False, dmin=0, dmax=0, pixels=False, flip=False, wrap=None, mul=1, symmetric=False, **kargs):
        """
        Function to display the image with a lot of parametrization

        Parameters
        ----------
        ax : matplotlib axis or None
            matplotlib axis if given otherwise current axis will be used (plt.gca())
        sig : float
            sigma values to adjust the contrast range around the mean ±sig times the standard-deviation
        cmap : string
            colormap name used. By default a gray map is used. If the zscale of the data are in 'meter' (i.e. topography data) the 'hot' colormap is used
        title : string
            The title of the plot. By default is the channel name
        adaptive : bool
            The color scale used is linear. If adaptive is True a non linear color scale is used in order that each color is used with the same amount.
        dmin : float
            minimum value adjustment used for the colorscale
        dmax: float
            maximum value adjustment used for the colorscale
        pixels : bool
            Display the image with x/y-labels with real unit. If pixels is True, the axes are in pixels
        flip : bool
            Flip the image upside-down
        wrap : Nont or int
            wrap the title to a width of wrap chars
        symmetric : bool
            If True will place the middle of the colorscale to the value 0.
            This is specially usefull for diverging colormaps such as : BrBG, bwr, coolwarm, seismiv, spectral, etc.
        level : float
            level should be ≥0 and <50. Adjust the lower and upper colorscale to level% and (100-level)% of the data range.
            e.g. if level=1, the colorscale will display 1-99% of the data range
        vmin : float
            Minimum value used for the colorscale
        vmax : flaot
            Maximum value used for the colorscale


        Returns
        -------
        matplotlib.image.AxesImage
            matplolib axis instance returned by imshow

        Examples
        --------
        >>> topo = pySPM.SPM_image(...)
        >>> fig, (ax, ax2) = plt.subplots(2, 3, figsize=(15, 10))
        >>> topo.show(ax=ax[0], cmap='gray', title="color map=\"gray\"")
        >>> topo.show(ax=ax[1], sig=2, title="standard deviation=2")
        >>> topo.show(ax=ax[2], adaptive=True, title="Adaptive colormap")
        >>> topo.show(ax=ax2[0], dmin=4e-8, cmap='gray', title="raise the lowest value for the colormap of +40nm")
        >>> topo.show(ax=ax2[1], dmin=3e-8, dmax=-3e-8, cmap='gray',title="raise lower of +30nm and highest of -30nm")
        >>> topo.show(ax=ax2[2], pixels=True, title="Set axis value in pixels");
        """
        mpl.rc('axes', grid=False)
        
        if ax is None:
            ax = plt.gca()
        if title == None:
            title = u"{0} - {1}".format(self.type, self.channel)
        if wrap is not None:
            title = "\n".join([title[i*wrap:(i+1)*wrap]
                               for i in range(int(len(title)/wrap)+1)])
        unit = self.size['real']['unit']
        sunit = 'afpnum kMGTPE'
        if len(unit) == 1 or unit in ['pixels']:
            isunit = 6
        elif unit[0] in sunit:
            isunit = sunit.find(unit[0])
            unit = unit[1:]
        else:
            isunit = 6
        W = self.size['real']['x']
        H = self.size['real']['y']
        fact = int(np.floor(np.log(W)/np.log(10)/3))
        isunit += fact
        W, H = W/10**(fact*3), H/10**(fact*3)
        if cmap == None:
            cmap = 'gray'
            if unit == 'm' and self.channel == "Topography":
                cmap = 'hot'
        mi, ma = np.nanmin(self.pixels), np.nanmax(self.pixels)            
        if adaptive:
            img = np.asarray(256**2*(self.pixels-mi)/(ma-mi), dtype=np.uint16)
            mi, ma = 0, 1
            img = skimage.exposure.equalize_adapthist(img, clip_limit=0.03)
        else:
            img = mul*self.pixels
            mi *= mul
            ma *= mul
        
        if sig == None:
            vmin = mi+dmin
            vmax = ma+dmax
        else:
            std = np.nanstd(img)
            avg = np.nanmean(img)
            vmin = avg - sig * std
            vmax = avg + sig * std
        if 'level' in kargs:
            if kargs['level'] < 0 or kargs['level']>=50:
                raise ValueError("The level shoud have a value in [0,50)")
            vmax = np.percentile(img, 100-kargs['level'])
            vmin = np.percentile(img, kargs['level'])
            del kargs['level']
        if 'vmin' in kargs:
            vmin = kargs['vmin']
            del kargs['vmin']
        if 'vmax' in kargs:
            vmax = kargs['vmax']
            del kargs['vmax']
        if symmetric:
            vmax = abs(max(vmin,vmax))
            vmin = -vmax
        if not flip:
            if pixels:
                r = ax.imshow(np.flipud(img), cmap=cmap, vmin=vmin, vmax=vmax, **kargs)
            else:
                r = ax.imshow(np.flipud(img), extent=[0, W, 0, H], cmap=cmap, vmin=vmin, vmax=vmax, **kargs)
        else:
            if pixels:
                r = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, **kargs)
            else:
                r = ax.imshow(img, cmap=cmap, extent=[0, W, 0, H], vmin=vmin, vmax=vmax, **kargs)
        if pixels:
            ax.set_xlim((0, self.pixels.shape[1]))
            ax.set_ylim((self.pixels.shape[0], 0))

        if not pixels:
            if isunit != 6:
                u = sunit[isunit]
                if u == 'u':
                    u = '$\\mu$'
                ax.set_xlabel(u'x [{0}{1}]'.format(u, unit))
                ax.set_ylabel(u'y [{0}{1}]'.format(u, unit))
            else:
                ax.set_xlabel(u'x [{0}]'.format(unit))
                ax.set_ylabel(u'y [{0}]'.format(unit))
        if title != None:
            ax.set_title(title)
        return r
        
    def real2px(self, x, y):
        """
        Transform a real (x,y) value in pixels
        Units should be the same as the one plotted by pySPM.SPM_image.show
        """
        return self.real2pixels(x,y)
        
    def real2pixels(self, x, y, float=False):
        """
        Transform a real (x,y) value in pixels
        Units should be the same as the one plotted by pySPM.SPM_image.show
        """
        W = self.size['real']['x']
        fact = int(np.floor(np.log(W)/np.log(10)/3))*3
        if not float:
            px = np.digitize(x, np.linspace(0,self.size['real']['x']/(10**fact),self.pixels.shape[1]), right=True)
            py = np.digitize(y, np.linspace(0,self.size['real']['y']/(10**fact),self.pixels.shape[0]), right=False)
        else:
            px = x*(self.pixels.shape[1]-1)/(self.size['real']['x']/(10**fact))
            py = y*(self.pixels.shape[0]-1)/(self.size['real']['y']/(10**fact))
        return px, py
        
    def px2real(self, x, y):
        """
        Transform  a (x,y) value from pixels to real
        Units are the same as the one plotted by pySPM.SPM_image.show
        """
        W = self.size['real']['x']
        fact = int(np.floor(np.log(W)/np.log(10)/3))*3
        rx = x*self.size['real']['x']/(10**fact)/self.pixels.shape[1]
        ry = (self.pixels.shape[0]-y)*self.size['real']['y']/(10**fact)/self.pixels.shape[0]
        return rx, ry
    
    def circular_profile(self, x0, y0, Ra=1, Rn=0, width=1, N=20, A=0, B=360,\
        cmap='jet', axImg=None, axPolar=None, axProfile=None, plotProfileEvery=1,\
        xtransf=lambda x: x*1e9, ytransf=lambda x:x*1e9,\
        ToFcorr=False, fit=lambda x, *p: p[3]+p[2]*CDF(x, *p[:2]), p0=None,errors=False,bounds=(-np.inf,np.inf), **kargs):
        """
        Create radial profiles from point x0,y0 with length Ra (outer radius) and Rn (negative Radius).
        Start from angle A° to angle B° with N profiles.
        If you want to apply the ToF-correction, please set ToFcorr to the number of scans used to record the ToF-SIMS image.
        Return the fitting uncertainty on sigma if errors is set to True
        The fitting function can be adjusted by fit and the default parameters by p0 which is an array of function where the first parameter passed will be the x-values and the second the y-values.
        """
        from matplotlib import colors, cm
        
        # Create a colormap for each profile
        CM =  plt.get_cmap(cmap) 
        cNorm  = colors.Normalize(vmin=0, vmax=N)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=CM)
        res = []
        cov = []
        angles = []
        assert A<B        
        for i, angle in enumerate(np.linspace(A, B, N)):
            a = np.radians(angle)
            angles.append(a)
            l, p = self.get_profile(
                x0-Rn*np.cos(a),
                y0+Rn*np.sin(a),
                x0+Ra*np.cos(a),
                y0-Ra*np.sin(a),
                ax=axImg, width=width, color=scalarMap.to_rgba(i), **kargs)
            profile = np.mean(p, axis=1)
            if ToFcorr:
                profile = -np.log(1.001-profile/ToFcorr)
            if p0 is None:
                p0 = [l[len(l)//2], 100e-9, np.max(profile)-np.min(profile),np.min(profile) ]
            else:
                for j,p in enumerate(p0):
                    if callable(p):
                        p0[j] = p(l,profile)
            if kargs.get('debug',False):
                print("calculate fit parameters are", p0)
            popt, pcov = scipy.optimize.curve_fit(fit, l , profile, p0)
            res.append(popt)
            cov.append([np.sqrt(abs(pcov[i,i])) for i in range(len(popt))])
            if axProfile and i%plotProfileEvery == 0:
                axProfile.plot(xtransf(l-popt[0]), profile, color=scalarMap.to_rgba(i), linestyle=':')
                axProfile.plot(xtransf(l-popt[0]), fit(l,*popt), color=scalarMap.to_rgba(i))
        # close loop
        if A%360 == B%360:
            angles.append(angles[0])
            res.append(res[0])
            cov.append(cov[0])
        
        # Plot polar
        angles = np.array(angles)
        res = np.array(res)
        cov = np.array(cov)
        fact = 2*np.sqrt(2*np.log(2))
        if axPolar:
            axPolar.plot(angles, ytransf(res[:,1]), color=kargs.get('sig_color','C0'))
            axPolar.plot(angles, ytransf(fact*res[:,1]), color=kargs.get('fwhm_color','C1'))
            if errors:
                axPolar.fill_between(angles, ytransf(res[:,1]-cov[:,1]),ytransf(res[:,1]+cov[:,1]), color=kargs.get('sig_color','C0'), alpha=kargs.get('fillalpha',.5))
                axPolar.fill_between(angles, fact*ytransf(res[:,1]-cov[:,1]),ytransf(res[:,1]+cov[:,1]), color=kargs.get('fwhm_color','C1'), alpha=kargs.get('fillalpha',.5))
            
        return angles, res, cov
        
    def get_profile(self, x1, y1, x2, y2, width=0, ax=None, pixels=True, color='w', axPixels=None, **kargs):
        """
        retrieve the profile of the image between pixel x1,y1 and x2,y2

        Parameters
        ----------
        x1, y1, x2, y2 : ints
            coordinates for the profile
        ax : matplotlib axis
            defines the matplotlib axis on which the position of the profile should be drawn (in not None)
        width : int
            the width of the profile (for averaging/statistics) in pixels
        color : string
            color used to plot the profiles lines
        axPixels : bool
            If True the image plotted in the ax axis is displayed in pixels
        
        Returns
        -------
        x data : 1D numpy array
        profile : 1D numpy array
        """
        if kargs.get('debug',False):
            print("get_profile input coordinates:", x1, y1, x2, y2)
        if axPixels is None:
            axPixels = pixels
        W = self.size['real']['x']
        fact = int(np.floor(np.log(W)/np.log(10)/3))*3
        if not pixels:
            if kargs.get('debug', False):
                print("Image range (real scale):", self.size['real']['x']/(10**fact), self.size['real']['y']/(10**fact))
            x1, y1 = self.real2pixels(x1, y1, float=True)
            x2, y2 = self.real2pixels(x2, y2, float=True)
            y1 = self.pixels.shape[0]-y1
            y2 = self.pixels.shape[0]-y2
            if kargs.get('debug', False):
                print("Pixel coordinates:", x1, y1, x2, y2)
            if not axPixels:
                xvalues, p = get_profile(np.flipud(self.pixels), x1, y1, x2, y2, ax=ax, width=width, color=color,\
                    transx = lambda x: x*(self.size['real']['x']/(10**fact))/self.pixels.shape[1],\
                    transy = lambda x: (self.pixels.shape[0]-x)*(self.size['real']['y']/(10**fact))/self.pixels.shape[0],\
                    **kargs)
            else:
                values, p = get_profile(np.flipud(self.pixels), x1, y1, x2, y2, ax=ax, width=width, color=color, **kargs)
        else:
            if axPixels:
                values, p = get_profile(np.flipud(self.pixels), x1, y1, x2, y2, ax=ax, width=width, color=color, **kargs)
            else:
                values, p = get_profile(np.flipud(self.pixels), x1, y1, x2, y2, ax=ax, width=width, color=color,\
                    transx = lambda x: x*(self.size['real']['x']/(10**fact))/self.pixels.shape[1],\
                    transy = lambda x: (self.pixels.shape[0]-x)*(self.size['real']['y']/(10**fact))/self.pixels.shape[0],\
                    **kargs)

        dx = (x2-x1)*self.size['real']['x']/self.size['pixels']['x']
        dy = (y2-y1)*self.size['real']['y']/self.size['pixels']['y']
        rd = np.sqrt(dx**2+dy**2)
        xvalues = np.linspace(0, rd, len(p))
        return xvalues, p

    def plot_profile(self, x1, y1, x2, y2, width=0, ax=None, pixels=True, img=None, imgColor='w', ztransf=lambda x: x, zunit=None, **kargs):
        """
        Retrieve and plot a profile from an image

        Parameters
        ----------
        x1, y1, x2, y2 : int
            coordinate of the profile in real size or in pixels (if pixels is True)
        width : float
            the width of the profiles in pixels for better statistics
        ax : matplotlib axis
            The axis in which the profile will be plotted
        pixels : bool
            If True the coordinates are given in pixels and not in real units
        img : matplotlib axis
            The axis in which the profile position will be drawn
        imgColor : string
            The color used to display the profile positions
        ztransf : function
            function to transform the profile data. This can be used to scale the data.
            Most profiles are retrieved in 'm' and a 'nm' value can be used by using ztransf=lambda x: x*1e9
        zunit : string
            the zunit name used if ztransft is used
        color : string
            The color of the profile
        col : string
            can be used instead of color
        stdplot : bool
            If True display the ±nσ plots where n is given by the sig parameter
        sig : int
            The number of sigmas used in stdplot
        label : string
            The label used for plotting the profile (useful if you perform a ax.legend() afterwards)

        Returns
        -------
        dictionary : {'plot': matplotlib_plot_instance, 'l': profile_xaxis, 'z': profile_yaxis}

        Examples
        --------
        >>> topo = pySPM.SPM_image(...)
        >>> fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        >>> topo.plot_profile(70, 100, 170, 200, ax=ax[1], img=ax[0], ztransf=lambda x:x*1e9, zunit='nm');
        >>> topo.show(ax=ax[0], pixels=True);
        """
        col = kargs.get('color',kargs.get('col','C0'))
        W = self.size['real']['x']
        fact = int(np.floor(np.log(W)/np.log(10)/3))*3
        if ax == None:
            ax = plt.gca()
        xvalues, p = self.get_profile(x1, y1, x2, y2, width=width, color=imgColor, ax=img, pixels=pixels, **kargs)
        d = np.sqrt((x2-x1)**2+(y2-y1)**2)
        dx = (x2-x1)
        dy = (y2-y1)
        if pixels:
            rd = d
            u = ''
            unit = 'px'
        else:
            unit = self.size['real']['unit']
            sunit = 'afpnum kMGTPE'
            if len(unit) == 1:
                isunit = 6
            elif unit[0] in sunit:
                isunit = sunit.find(unit[0])
                unit = unit[1:]
            else:
                isunit = 6
            isunit += fact//3
            if isunit != 6:
                u = sunit[isunit]
            else:
                u=''
            if u == 'u':
                u = '$\\mu$'
            rd = np.sqrt(dx**2+dy**2)
        xvalues = np.linspace(0, rd, len(p))
        lab = kargs.get("label", "")
        if width < 2:
            profile = ztransf(p)
        else:
            profile = ztransf(np.mean(p, axis=1))
            s = np.std(p)
            if kargs.get('stdplot', True):
                for ns in range(1, kargs.get('sig', 2)+1):
                    ax.fill_between(xvalues, profile-ns*s, profile+ns*s, color=col, alpha=.2, label=[lab+' ($\\sigma,\ldots {}\\sigma$)'.format(kargs.get('sig',2)),None][ns>1])
        
        Plot = ax.plot(xvalues, profile, color=col, linewidth=kargs.get('linewidth',1),linestyle=kargs.get('linestyle','-'), label=lab+[' (mean)',''][width<2])
        if kargs.get('min',False):
            minStyle = kargs.get('minStyle', kargs.get('minmaxStyle', '--'))
            minColor = kargs.get('minColor', kargs.get('minmaxColor', col))
            minMarker = kargs.get('minMarker', kargs.get('minmaxMarker', ''))
            ax.plot(xvalues, np.min(p, axis=1), color=minColor, linewidth=kargs.get('linewidth',1),linestyle=minStyle, marker=minMarker, label=lab+' (min)')
        if kargs.get('max', False):
            maxStyle = kargs.get('maxStyle',kargs.get('minmaxStyle','--'))
            maxColor = kargs.get('maxColor',kargs.get('minmaxColor',col))
            maxMarker = kargs.get('maxMarker',kargs.get('minmaxMarker',''))
            ax.plot(xvalues, np.max(p, axis=1), color=maxColor, linestyle=maxStyle, linewidth=kargs.get('linewidth',1), marker=maxMarker, label=lab+' (max)')
            
        ax.set_xlabel("Distance [{1}{0}]".format(unit, u))
        if zunit is not None:
            ax.set_ylabel("{1} [{0}]".format(zunit, self.channel))
        else:
            ax.set_ylabel("{1} [{0}]".format(self.zscale, self.channel))
       
        return {'plot': Plot, 'l': xvalues, 'z': profile}

    def get_bin_threshold(self, percent, high=True, adaptive=False, binary=True, img=False):
        """
        Threshold the image into binary values
        
        Parameters
        ----------
        percent : float
            The percentage where the thresholding is made
        high : bool
            If high a value of 1 is returned for values > percent
        adaptive : bool
            If True, performs an adaptive thresholding (see skimage.filters.threshold_adaptive)
        binary : bool
            If True return bool data (True/False) otherwise  numeric (0/1)
        img : bool
            If True return a SPM_image otherwise a numpy array
        """
        if adaptive:
            if binary:
                return self.pixels > threshold_local(self.pixels, percent)
            return threshold_local(self.pixels, percent)
        mi = np.min(self.pixels)
        norm = (self.pixels-mi)/(np.max(self.pixels)-mi)
        if high:
            r = norm > percent
        else:
            r = norm < percent
        if not img:
            if binary:
                return r
            return np.ones(self.pixels.shape)*r
        else:
            I = copy.deepcopy(self)
            I.channel = "Threshold from "+I.channel
            if binary:
                I.pixels = r
            else:
                I.pixels = np.ones(self.pixels.shape)*r
            return I

    def spline_offset(self, X, Y, Z=None, inline=True, ax=None, output='img', **kargs):
        """
        subtract a spline interpolated by points corrdinates.
        if Z is None, the image values will be used (default)
        """
        if ax is not None:
            if 'num' in kargs and kargs['num']:
                text_color = 'k'
                if 'text_color' in kargs:
                    text_color = kargs['text_color']
                    del kargs['text_color']
                for i in range(len(X)):
                    l = self.pixels.shape[1]-X[i] < 20
                    ax.annotate(str(i), (X[i], Y[i]), ([
                                5, -5][l], 0), textcoords='offset pixels', va="center", ha=["left", "right"][l], color=text_color)
                del kargs['num']
            ax.plot(X, Y, 'o', **kargs)
        import scipy.interpolate
        T = np.flipud(self.pixels) - np.min(self.pixels)
        if Z is None:
            Z = [T[Y[i], X[i]] for i in range(len(X))]
        x = np.arange(self.pixels.shape[1])
        y = np.arange(self.pixels.shape[0])
        xx, yy = np.meshgrid(x, y)
        I = scipy.interpolate.SmoothBivariateSpline(X, Y, Z)
        z = I.ev(xx, yy)
        if inline:
            self.pixels -= z
            return z
        else:
            if output == 'img':
                New = copy.deepcopy(self)
                New.pixels -= z
                return New
            elif output == 'spline':
                return z
            else:
                raise ValueError(
                    "The output parameter should be either 'img' or 'spline'")

    def get_shadow_mask(self, angle, BIN=None, pb=False):
        """
        If an image is recorded with a beam incident with a certain angle, the topography will shadow the data.
        This function generates the shadow mask for a given topography and a given incident angle.

        Parameters
        ----------
        angle : float
            The incidence angle in degrees
        BIN : numpy array
            Data. If given will move the recorded pixels at the correct x,y positions
        pb : bool
            display a progressbar ?

        Note
        ----
        This function is old, might not be optimized or working properly
        """
        if BIN is not None:
            BIN = BIN*1.0
        slope = np.tan(np.radians(angle))
        neg = False
        if slope < 0:
            neg = True
            slope = -slope
            topo = np.fliplr(self.pixels)
            if BIN is not None:
                BIN = np.fliplr(BIN)
        else:
            topo = self.pixels
        x = np.linspace(0, self.size['real']['x'], self.pixels.shape[1])
        if self.size['real']['unit'] == 'um':
            x *= 1e-6
        elif self.size['real']['unit'] == 'nm':
            x *= 1e-9
        mask = np.zeros(self.pixels.shape)
        AFM_bin_shadow = np.zeros(self.pixels.shape)
        Y = range(self.pixels.shape[0])
        if pb:
            Y = tqdm(Y)
        for yi in Y:
            for xi in range(self.pixels.shape[1]):
                cut = self.pixels.shape[1]-2
                y_ray = slope*(x-x[xi]) + topo[yi, xi]
                while cut > xi and y_ray[cut] > topo[yi, cut]:
                    cut -= 1
                if xi == cut:
                    if BIN is not None:
                        AFM_bin_shadow[yi, xi] = BIN[yi, xi]
                    continue
                # Cut has been found
                if BIN is not None:
                    x1 = x[cut]
                    x2 = x[cut+1]
                    y1 = topo[yi, cut]
                    y2 = topo[yi, cut+1]
                    x0 = x[xi]
                    y0 = topo[yi, xi]
                    if y2 == y1:
                        x_cut = (y1+slope*x0-y0)/slope
                        y_cut = y1
                    else:
                        numerator = x1/(x2-x1)+(y0-slope*x0-y1)/(y2-y1)
                        denominator = 1/(x2-x1)-slope/(y2-y1)
                        x_cut = numerator / denominator
                        y_cut = slope*(x_cut-x0)+y0
                    if x_cut >= x1 and x_cut <= x2:
                        y1 = BIN[yi, cut]
                        y2 = BIN[yi, cut+1]
                        yint = (((y2-y1)/(x2-x1))*(x_cut-x1))+y1
                    else:
                        yint = BIN[yi, xi]
                    AFM_bin_shadow[yi, xi] = yint
                mask[yi, xi] = 1
        if neg:
            mask = np.fliplr(mask)
            AFM_bin_shadow = np.fliplr(AFM_bin_shadow)
        if BIN is not None:
            return (mask, AFM_bin_shadow)
        return mask

    def adjust_position(self, fixed):
        """
        Shift the current pixels to match a fixed image.
        The shift is determined by position where the cross-correlation is maximized.
        """
        adj = copy.deepcopy(self)
        cor = np.fft.fft2(fixed.pixels)
        cor = np.abs(np.fft.ifft2(np.conj(cor) * np.fft.fft2(self.pixels)))
        cor = cor / fixed.pixels.size
        ypeak, xpeak = np.unravel_index(cor.argmax(), cor.shape)
        shift = [-(ypeak-1), -(xpeak-1)]
        adj.pixels = np.roll(self.pixels, shift[0], axis=0)
        adj.pixels = np.roll(adj.pixels, shift[1], axis=1)
        return adj

    def align(self, tform, cut=True):
        """
        Apply an Affine transform on the data

        Parameters
        ----------
        tform : skimage.transform
            the affine transform to perform
        cut : bool
            If True cut the data
        """
        New = copy.deepcopy(self)
        New.pixels = tf.warp(self.pixels, tform, preserve_range=True)
        if not cut:
            return New
        cut = [0, 0] + list(self.pixels.shape)
        if tform.translation[0] >= 0:
            cut[2] -= tform.translation[0]
        elif tform.translation[0] < 0:
            cut[0] -= tform.translation[0]
        if tform.translation[1] >= 0:
            cut[1] += tform.translation[1]
        elif tform.translation[1] < 0:
            cut[3] += tform.translation[1]
        cut = [int(x) for x in cut]
        New.cut(cut, inplace=True)
        return New, cut

    def get_fft(self):
        """
        return the FFT2 transform opf the image
        """
        return np.fft.fftshift(np.fft.fft2(self.pixels))

    def corr_fit2d(self, nx=2, ny=1, poly=False, inline=True, mask=None):
        """
        Subtract a fitted 2D-polynom of nx and ny order from the data

        Parameters
        ----------
        nx : int
            the polynom order for the x-axis
        ny : int
            the polynom order for the y-axis
        poly : bool
            if True the polynom is returned as output
        inline : bool
            create a new object?
        mask : 2D numpy array
            mask where the fitting should be performed
        """
        r, z = fit2d(self.pixels, nx, ny, mask=mask)
        if inline:
            self.pixels -= z
        else:
            N = copy.deepcopy(self)
            N.pixels -= z
            if poly:
                return N, z
            return N
        if poly:
            return z
        return self

    def zero_min(self, inline=True):
        """
        Shift the values so that the minimum becomes zero.
        """
        if inline:
            self.pixels -= np.min(self.pixels)
            return self
        else:
            N = copy.deepcopy(self)
            N.pixels -= np.min(N.pixels)
            return N

    def filter_lowpass(self, p, inline=True):
        """
        Execute a lowpass filter on the data
        """
        F = self.get_fft()
        mask = self.getRmask() < p
        if inline:
            self.pixels = np.real(np.fft.ifft2(np.fft.fftshift(F*mask)))
        else:
            C = copy.deepcopy(self)
            C.pixels = np.real(np.fft.ifft2(np.fft.fftshift(F*mask)))
            return C

    def _resize_infos(self):
        """
        Internal to recalculate the real size when the image is cropped or cut
        """
        self.size['real']['x'] *= self.pixels.shape[1]/self.size['pixels']['x']
        self.size['real']['y'] *= self.pixels.shape[0]/self.size['pixels']['y']
        self.size['pixels']['x'] = int(self.pixels.shape[1])
        self.size['pixels']['y'] = int(self.pixels.shape[0])
        if 'recorded' in self.size:
            self.size['recorded']['real']['x'] \
                *= (self.pixels.shape[1]/self.size['pixels']['x'])
            self.size['recorded']['real']['y'] \
                *= (self.pixels.shape[0]/self.size['pixels']['y'])
            self.size['recorded']['pixels']['x'] = int(self.pixels.shape[1])
            self.size['recorded']['pixels']['y'] = int(self.pixels.shape[0])
        
    def filter_scars_removal(self, thresh=.5, inline=True):
        """
        Filter function to remove scars from images.
        """
        if not inline:
            C = copy.deepcopy(self)
        else:
            C = self
        for y in range(1, self.pixels.shape[0]-1):
            b = self.pixels[y-1, :]
            c = self.pixels[y, :]
            a = self.pixels[y+1, :]
            mask = np.abs(b-a) < thresh*(np.abs(c-a))
            C.pixels[y, mask] = b[mask]
        if not inline:
            return C
        return self

    def cut(self, c, inline=False, pixels=True, **kargs):
        """
        Clip/Crop the image

        Parameters
        ----------
        c : list [llx,lly,urx,ury]
            list of the lowe-left (ll) and upper-right (ur) coordinates
        inline: bool
            perform the transformation inline or produce a new SPM_image?
        pixels : bool
            Are the coordinates given in pixels?

        Returns
        -------
        self if inplace, clipped SPM_image otherwises2 = pySPM.Nanoscan("%s/CyI5b_PCB_ns.xml"%(Path))
        """
        if 'inplace' in kargs:
            inline=kargs['inplace']
        if kargs.get('debug',False):
            print("cut) Input coordinates:", c)
        if not pixels:
            c = [z for s in zip(*self.real2pixels(c[0::2], c[1::2])) for z in s]
            if kargs.get('debug',False):
                print("cut) pixel coordinates:", c)
        if not inline:
            new = copy.deepcopy(self)
            new.pixels = cut(self.pixels, c, **kargs)
            new._resize_infos()
            return new
        else:
            self.pixels = cut(self.pixels, c, **kargs)
            self._resize_infos()
            return self
    
    def zoom(self, zoom_factor, inplace=False, order=3):
        """
        Resize the image to a new pixel size (but keep the real size) by pixel interpolation.

        Parameters
        ----------
        zoom_factor : float
            > 1: up sampling
            < 1: down sampling
        order : int
            The spline interpolation order to use. (default: 3). Use 0 for binary or very sharp images.
        inplace : bool
            create a new image?
        """
        from scipy.ndimage.interpolation import zoom
        if not inplace:        
            new = copy.deepcopy(self)
            new.pixels = zoom(new.pixels, zoom_factor, order=order)
            new.size['pixels']['x'] = new.pixels.shape[1]
            new.size['pixels']['y'] = new.pixels.shape[0]
            return new
        else:
            self.pixels = zoom(self.pixels, zoom_factor, order=order)
            self.size['pixels']['x'] = self.pixels.shape[1]
            self.size['pixels']['y'] = self.pixels.shape[0]
            return self

# Note: The following functions are not part of the SPM_image class.
# All following functions are performed on numpy arrays

def cut(img, c, **kargs):
    """
    Clip / Crop a numpy array

    Parameters
    ----------
    img : 2D numpy array
        The input image array
    c : list [llx, lly, urx, ury]
        the lower-left (ll) and upper-right (ur) coordinates used for the clippings2 = pySPM.Nanoscan("%s/CyI5b_PCB_ns.xml"%(Path))
    """
    if kargs.get('debug',False):
        print("cut in x", c[0], "->", c[2], " - in y", c[1], "->", c[3])
    if c[3] < c[1]:
        c = [c[0],c[3],c[2],c[1]]
    if c[2] < c[0]:
        c = [c[2],c[1],c[0],c[3]]
    if c[2]-c[0] == img.shape[1] and c[3]-c[1] == img.shape[0]:
        raise Exception("Reshaping the same array again?")
    return img[c[1]:c[3], c[0]:c[2]]


def normalize(data, sig=None, vmin=None, vmax=None):
    """
    Normalize the input data. Minimum_value -> 0 and maximum_value -> 1

    Parameters
    ----------
    data : numpy array
        input data
    sig : float or None
        if not None:
            mean(data)-sig*standard_deviation(data) -> 0
            mean(data)+sig*standard_deviation(data) -> 1
    vmin : float or None
        if not None, define the lower bound i.e.  vmin -> 0
    vmax : float or None
        if not None, defines the upper bound i.e. vmax -> 0

    Note
    ----
    All values below the lower bound will be = 0
    and all values above the upper bound will be = 1


    """
    if sig is None:
        mi = np.min(data)
        ma = np.max(data)
    else:
        s = sig*np.std(data)
        mi = np.mean(data)-s
        ma = np.mean(data)+s
    if vmin is not None:
        mi = vmin
    if vmax is not None:
        ma = vmax
    N = (data-mi)/(ma-mi)
    N[N < 0] = 0
    N[N > 1] = 1
    return N


def imshow_sig(img, sig=1, ax=None, **kargs):
    """
    Shortcut to plot a numpy array around it's mean with bounds ±sig sigmas

    Parameters
    ----------
    img : 2D numpy array
        input image to display
    sig : float
        The number of standard-deviation to plot
    ax : matplotlib axis
        matplotlib axis to use. If None, the current axis (plt.gca() will be used).
    **kargs : additional parameters
        will be passed to the imshow function of matplotls2 = pySPM.Nanoscan("%s/CyI5b_PCB_ns.xml"%(Path))ib
    """
    if ax == None:
        fig, ax = plt.subplots(1, 1)
    std = np.std(img)
    avg = np.mean(img)
    vmin = avg - sig * std
    vmax = avg + sig * std
    ax.imshow(img, vmin=vmin, vmax=vmax, **kargs)


def adjust_position(fixed, to_adjust, shift=False):
    """ Shift the current pixels to match a fixed image by rolling the data"""
    adj = copy.deepcopy(to_adjust)
    cor = np.fft.fft2(fixed)
    cor = np.abs(np.fft.ifft2(np.conj(cor) * np.fft.fft2(to_adjust)))
    cor = cor / to_adjust.size
    ypeak, xpeak = np.unravel_index(cor.argmax(), cor.shape)
    shift = [-(ypeak-1), -(xpeak-1)]
    adj = np.roll(to_adjust, shift[0], axis=0)
    adj = np.roll(adj, shift[1], axis=1)
    if shift:
        return adj, shift
    return adj


def tukeyfy(A, alpha, type='default'):
    """
    Apply a Tukey window on the current image

    Parameters
    ----------
    A : 2D numpy array
        input array
    alpha : float
        Size of the Tukey windows in percent of the image (≥0 and ≤1)
    type : string
        if not "default" perform a mean centering (data will blend down to its mean instead of 0)
    """
    tuky = tukeywin(A.shape[0], alpha)
    tukx = tukeywin(A.shape[1], alpha)
    tuk = np.multiply(tukx[:, None].T, tuky[:, None])
    if type is 'default':
        return A * tuk
    avg = np.mean(A)
    return avg+(A-avg) * tuk


def tukeywin(window_length, alpha=0.5):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
    that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
    at \alpha = 0 it becomes a Hann window.

    We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
    output

    Reference
    ---------
    http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html

    '''
    # Special cases
    if alpha <= 0:
        return np.ones(window_length)  # rectangular window
    elif alpha >= 1:
        return np.hanning(window_length)

    # Normal case
    x = np.linspace(0, 1, window_length)
    w = np.ones(x.shape)

    # first condition 0 <= x < alpha/2
    first_condition = x < alpha/2
    w[first_condition] = 0.5 * \
        (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2)))

    # second condition already taken care of

    # third condition 1 - alpha / 2 <= x <= 1
    third_condition = x >= (1 - alpha/2)
    w[third_condition] = 0.5 * \
        (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2)))
    return w


def overlay(ax, mask, color, **kargs):
    """
    Plot an overlay on an existing axis

    Parameters
    ----------
    ax : matplotlib axis
        input axis
    mask : 2D numpy array
        Binary array where a mask should be plotted
    color : string
        The color of the mask to  plot
    **kargs: additional parameters
        passed to the imshow function of matploltib
    """
    m = ma.masked_array(mask, ~mask)
    col = np.array(colors.colorConverter.to_rgba(color))
    I = col[:, None, None].T*m[:, :, None]
    ax.imshow(I, **kargs)


def normP(x, p, trunk=True):
    """
    Normalize the input data accroding to its percentile value.

    Parameters
    ----------
    x : 2D numpy array
        input data
    p : float
        percentile to normalize the data.
        lower bound = p percentile
        upper bound = (100-p) percentile
    trunk : bool
        If True the data are truncated between 0 and 1
    """
    thresh_high = np.percentile(x, 100-p)
    thresh_low = np.percentile(x, p)
    if thresh_low == thresh_high:
        thresh_high = np.max(x)
        thresh_low = np.min(x)
    if thresh_low == thresh_high:
        thresh_high = thresh_low+1
    r = (x-thresh_low)/(thresh_high-thresh_low)
    if trunk:
        r[r < 0] = 0
        r[r > 1] = 1
    return r

def beam_profile(target, source, mu=1e-6, tukey=0, meanCorr=False, source_tukey=None, real=np.abs, **kargs):
    """
    Calculate the PSF by deconvolution of the target
    with the source using a Tikhonov regularization of factor mu.
    """
    if source_tukey is None:
        source_tukey = tukey
    if kargs.get('source_centering', False):
        source = 2*source-1
    if meanCorr:
        target = target-np.mean(target)
    if tukey>0:
        target = tukeyfy(target, tukey)
    if source_tukey>0:
        source = tukeyfy(source, tukey)
    tf = np.fft.fft2(source)
    tf /= np.size(tf)
    recon_tf = np.conj(tf) / (np.abs(tf)**2 + mu)
    return np.fft.fftshift(real(np.fft.ifft2(np.fft.fft2(target) * recon_tf)))/np.size(target)


def beam_profile1d(target, source, mu=1e-6, real=np.abs):
    source = source
    tf = np.fft.fft(source)
    tf /= np.size(tf)
    recon_tf = np.conj(tf) / (np.abs(tf)**2 + mu)
    F = np.fft.fft(target) * recon_tf
    return np.fft.fftshift(real(np.fft.ifft(F))), F


def zoom_center(img, sx, sy=None):
    """
    Zoom by taking the sx × sy central pixels

    Parameters
    ----------
    img : 2D numpy array
        The input data
    sx : int
        The number of pixels along the x-axis to take

    sy : int or None
        The number of pixels alongs the y-axis to take.
        If None take the same value as for sx
    """
    if sy is None:
        sy = sx
    assert type(sx) is int
    assert type(sy) is int
    return img[
        img.shape[0]//2-sy//2: img.shape[0]//2 + sy//2,
        img.shape[1]//2-sx//2: img.shape[1]//2 + sx//2]

def px2real(x, y, size, ext):
    rx = ext[0]+(x/size[1])*(ext[1]-ext[0])
    ry = ext[2]+(y/size[0])*(ext[3]-ext[2])
    return rx, ry


def real2px(x, y, size, ext):
    px = size[1]*(x-ext[0])/(ext[1]-ext[0])
    py = size[0]*(y-ext[2])/(ext[3]-ext[2])
    return px, py


def gaussian(x, mu, sig, A=1):
    """
    Deprecated, please use pySPM.utils.Gauss
    Will be removed in the next revision
    """
    from warnings import warn
    warn("Function pySPM.SPM.gaussian is deprecated. Please use pySPM.utils.Gauss instead")
    return A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def stat(x):
    """
    Quick function to display the min/max, mean and standard-deviation of a given data.

    Note
    ----
    Please use the more adavanced function pySPM.utils.stat_info which also provide information about the quartile and also plot the distribution
    """
    print("Min: {mi:.3f}, Max: {ma:.3f}, Mean: {mean:.3f}, Std: {std:.3f}".format(mi=np.min(x), ma=np.max(x), mean=np.mean(x), std=np.std(x)))


def fit2d(Z0, dx=2, dy=1, mask=None):
    """
    Fit the input data with a 2D polynom of order dx × dy

    Parameters
    ----------
    Z0 : 2D numpy array
        input data
    dx : int
        order of the polynom for the x-axis
    dy : int
        order of the polynom for the y-xis
    mask : 2D numpy array
        Give a mask where True values only will be used to perform the fitting

    Returns
    -------
    numpy array
        fitting parameters
    2D numpy array
        result of the polynom
    """
    x = np.arange(Z0.shape[1], dtype=np.float)
    y = np.arange(Z0.shape[0], dtype=np.float)
    X0, Y0 = np.meshgrid(x, y)
    if mask is not None:
        X = X0[mask]
        Y = Y0[mask]
        Z = Z0[mask]
    else:
        X = X0
        Y = Y0
        Z = Z0
    x2 = X.ravel()
    y2 = Y.ravel()
    A = np.vstack([x2**i for i in range(dx+1)])
    A = np.vstack([A]+[y2**i for i in range(1, dy+1)])
    res = scipy.optimize.lsq_linear(A.T, Z.ravel())
    r = res['x']
    Z2 = r[0]*np.ones(Z0.shape)
    for i in range(1, dx+1):
        Z2 += r[i]*(X0**i)
    for i in range(1, dy+1):
        Z2 += r[dx+i]*(Y0**i)
    return r, Z2


def warp_and_cut(img, tform, cut=True):
    """
    Perform an Affine transform on the input data and cut them if cut=True

    Parameters
    ----------
    img : 2D numpy array
        input data
    tform : skimage.transform
        An Affine fransform to perform on the data
    cut : bool
        Should the data be cutted?
    """
    New = tf.warp(img, tform, preserve_range=True)
    Cut = [0, 0] + list(img.shape)
    if tform.translation[0] >= 0:
        Cut[2] -= tform.translation[0]
    elif tform.translation[0] < 0:
        Cut[0] -= tform.translation[0]
    if tform.translation[1] >= 0:
        Cut[1] += tform.translation[1]
    elif tform.translation[1] < 0:
        Cut[3] += tform.translation[1]
    Cut = [int(x) for x in Cut]
    if cut:
        New = cut(New, Cut)
    return New, Cut


def get_profile(I, x1, y1, x2, y2, width=0, ax=None, color='w', alpha=0, N=None,\
        transx=lambda x: x, transy=lambda x: x, interp_order=1, **kargs):
    """
    Get a profile from an input matrix.
    Low-level function. Doc will come laters2 = pySPM.Nanoscan("%s/CyI5b_PCB_ns.xml"%(Path))
    """
    d = np.sqrt((x2-x1)**2+(y2-y1)**2)
    if N is None:
        N = int(d)+1
    P = []
    dx = -width/2*(y2-y1)/d
    dy = width/2*(x2-x1)/d
    for w in np.linspace(-width/2, width/2, max(1,width)):
        dx = -w*(y2-y1)/d
        dy = w*(x2-x1)/d
        x = np.linspace(x1+dx, x2+dx, N)
        y = np.linspace(y1+dy, y2+dy, N)
        M = scipy.ndimage.interpolation.map_coordinates(I, np.vstack((y, x)), order=interp_order)
        P.append(M)
    if kargs.get('debug',False):
        print("get_profile input coordinates:",x1,y1,x2,y2)
    if not ax is None:
        x1 = transx(x1)
        x2 = transx(x2)
        y1 = transy(y1)
        y2 = transy(y2)
        if kargs.get('debug',False):
            print("Drawing coordinates:",x1,y1,x2,y2)
        dx = -width/2*(y2-y1)/d
        dy = width/2*(x2-x1)/d
        if type(color) in [tuple, list]:
            ax.plot([x1, x2], [y1, y2], color=color, alpha=kargs.get('linealpha',1))
            ax.plot([x1-dx, x1+dx], [y1-dy, y1+dy], color=color, alpha=kargs.get('linealpha',1))
            ax.plot([x2-dx, x2+dx], [y2-dy, y2+dy], color=color, alpha=kargs.get('linealpha',1))
        else:
            ax.plot([x1, x2], [y1, y2], color, alpha=kargs.get('linealpha',1))
            ax.plot([x1-dx, x1+dx], [y1-dy, y1+dy], color, alpha=kargs.get('linealpha',1))
            ax.plot([x2-dx, x2+dx], [y2-dy, y2+dy], color, alpha=kargs.get('linealpha',1))
        if alpha>0:
            import matplotlib.patches
            ax.add_patch(matplotlib.patches.Rectangle(
                (x1+dx,y1+dy), 
                2*np.sqrt(dx**2+dy**2),
                np.sqrt((x2-x1)**2+(y2-y1)**2),
                -np.degrees(np.arctan2(x2-x1,y2-y1)), color=color, alpha=alpha))
    if len(P)==1:
        return np.linspace(0, d, N), P[0]
    return np.linspace(0, d, N), np.vstack(P).T


def dist_v2(img, dx=1, dy=1):
    """
    Return a 2D array with the distance in pixel with the clothest corner of the array.
    """
    x2 = np.arange(img.shape[1])
    x2 = (np.minimum(x2, img.shape[1]-x2) * dx)**2
    y2 = np.arange(img.shape[0])
    y2 = (np.minimum(y2, img.shape[0] - y2) * dy)**2
    X, Y = np.meshgrid(x2, y2)
    return np.sqrt(X+Y)

def generate_k_matrices(A, dx, dy):
    """
    GENERATE_K_MATRICES k-Matrix generation (helper function).
    generates k-matrices for the 2D-channel CHANNEL.
   
    K is a matrix of the same size as the pixel matrix A, containing the real-life frequency distance of each 
    pixel position to the nearest corner of an matrix that is one pixel
    wider/higher.
    KX is of the same size as K and contains the real-life  difference in x-direction of each pixel position to the nearest corner 
    of a matrix that is one pixel wider/higher.
    Similarly, KY is of the  same size as K, containing the real-life difference in y-direction of  each pixel position to the nearest corner of an matrix that is one 
    pixel wider/higher.
    """
    ny, nx = A.shape
    dkx = 2*np.pi/(nx*dx)
    dky = 2*np.pi/(ny*dy)
    
    ky = np.arange(0, ny);
    ky = (np.mod(ky+ny/2, ny) - ny/2) * dky
    
    kx = np.arange(0, nx);
    kx = (np.mod(kx+nx/2, nx) - nx/2) * dkx
    
    kx, ky = np.meshgrid(kx, ky)
    k = dist_v2(A, dkx, dky)
    k[0, 0] = 1.0 # Prevent division by zero error and illegal operand errors. This may be improved...
    return k, kx, ky
        
def mfm_tf(nx, dx, ny, dy, tf_in, derivative=0, transform=0, z=0, A=0, theta=None, phi=None, d=None, delta_w=None):
    """
    Draft for the MFM tf function
    """
    k, kx, ky = generate_k_matrices(tf_in, dx, dy)
    # Distance loss
    tf_out = np.exp(-z*k)
    if d is not None:
        tf_out = tf_out / 2.0
        if not np.isinf(d):
            if d == 0:
                tf_out *= k
            else:
                tf_out *= 1 - np.exp(-d*k)
    if A == 0:
        if transform != 0:
            assert theta is not None
            assert phi is not None
            tf_out *= ((np.cos(theta)+1j*(np.cos(phi)*np.sin(-theta)*kx+np.sin(phi)*np.sin(-theta)*ky)) / k)**transform
        if derivative == 1:
            tf_out *= k
    else:
        pass # TODO
    return tf_out * tf_in

def mfm_inv_calc_flat(img, z, tf_in, thickness=None, delta_w=None, amplitude=0, derivative=0, transform=0, mu=1e-8):
    """
    MFM inv calc function
    """
    theta = np.radians(12)
    phi = np.radians(-90)
    ny, nx = img.shape
    tf = mfm_tf(nx, 1, ny, 1, tf_in, derivative, transform, z, amplitude, theta, phi, thickness, delta_w)
    tf[0,0] = np.real(np.mean(tf))
    recon_tf = np.conj(tf) / (mu+np.abs(tf)**2)
    work_img = np.abs(np.fft.ifft2(np.fft.fft2(img) * recon_tf ))
    return work_img

def getTikTf(Img, mu, tukey=0, source_tukey=0, debug=False, d=200, real=np.real):
    import scipy
    def fit(x, a ,A, bg, x0):
        return bg+(A-bg)*np.exp(-abs(x-x0)/a)
    
    x = np.arange(Img.shape[1])
    y = np.arange(Img.shape[0])
    X, Y = np.meshgrid(x, y)
    x0 = Img.shape[1]/2
    y0 = Img.shape[0]/2
    R = np.sqrt((X-x0)**2+(Y-y0)**2)
    
    Z = beam_profile(Img, Img, mu=mu, tukey=tukey, source_tukey=source_tukey, real=real)
    zoom = zoom_center(Z, d)
    P = zoom[zoom.shape[0]//2, :]
    p0 = (1,np.max(zoom), 0, len(P)/2)
    popt, pcov = scipy.optimize.curve_fit(fit, np.arange(len(P)), P, p0, bounds=((0,0,-np.inf,0),np.inf))
    bg = popt[2]
    a = popt[0]
    if debug:
        return bg+np.exp(-np.abs(R)/a), Z, p0, popt
    return bg+np.exp(-np.abs(R)/a)
    
DEPRECATED_METHODS = {
    'getRowProfile': 'get_row_profile',
    'plotProfile':   'plot_profile',
    'getProfile':    'get_profile',
    'getShadowMask': 'get_shadow_mask',
    'addScale':      'add_scale',
    'getExtent':     'get_extent',
    'getBinThreshold':'get_bin_threshold',
    'corrFit2d':      'corrr_fit2d',
    'getFFT':         'get_fft',
    'filterLowPass':  'filter_lowpass',
    'ResizeInfos':    '_resize_infos'
    }

def method_alias(old_name, new_name):
    def _alias_meth(self, *args, **kargs):
        from warnings import warn
        warn("Function {name} is deprecated. Please use \"{new_name}\" instead".format(name=old_name, new_name=new_name))
        return getattr(SPM_image, new_name)(self, *args, **kargs)
    return _alias_meth
        
def fun_alias(old_name, new_name):
    def _alias_fun(*args, **kargs):
        from warnings import warn
        warn("Function {name} is deprecated. Please use \"{new_name}\" instead".format(name=old_name, new_name=new_name))
        return getattr(SPM, new_name)(*args, **kargs)
    return _alias_fun
        
DEPRECATED_FUNCTIONS = {
        'NormP':'normP',
        'getProfile':'get_profile',
        'BeamProfile':'beam_profile',
        'BeamProfile1D':'beam_profile1d',
        'Normalize':'normalize',
        'Align':'align',
        'ZoomCenter':'zoom_center',
    }
    
for x in DEPRECATED_METHODS:
    setattr(SPM_image, x, method_alias(x, DEPRECATED_METHODS[x]))

thisModule = sys.modules[__name__]
for x in DEPRECATED_FUNCTIONS:
    setattr(thisModule, x, fun_alias(x, DEPRECATED_FUNCTIONS[x]))    
