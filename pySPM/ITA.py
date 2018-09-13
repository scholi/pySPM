# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
This module gives the ability to ready and parse the ToF-SIMS ITA files from iontof.
You can mainly retrieve images and spectra for each channel and scan.
"""

import numpy as np
import struct
import os.path
import zlib
import re
import copy

from .ITM import ITM
from .collection import Collection
from .SPM import SPM_image
from .Block import MissingBlock
from .utils import in_ipynb

class ITA(ITM):
    def __init__(self, filename):
        """
        Open an ITA file.
        
        Parameters
        ----------
        filename : string
            the path of the ita file
        
        Returns
        -------
        Class<ITA>
            ITA Object

        Examples
        --------
        >>> import pySPM
        >>> filename = "myfile.ita"
        >>> A = pySPM.ITA(filename)
        """
        ITM.__init__(self, filename)
        try:
            self.sx = self.root.goto(
                'filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.XSize').getLong()
            self.sy = self.root.goto(
                'filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.YSize').getLong()
        except MissingBlock:
            self.sx = self.size['pixels']['x']
            self.sy = self.size['pixels']['y']
        try:
            # self.Nscan = int(self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data'\
            #    '/ImageStackScans/Image.NumberOfScans').getLong())
            self.Nimg = int(self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data'
                                           '/ImageStackScans/Image.NumberOfImages').getLong())
        except:
            try:
                self.Nimg = self.root.goto('propend/Measurement.ScanNumber').getKeyValue()['int']
            except:
                raise TypeError(
                    "Invalid file format. Maybe the file is corrupted?")

        self.Width = self.root.goto('Meta/SI Image[0]/res_x').getLong()
        self.Height = self.root.goto('Meta/SI Image[0]/res_y').getLong()

        try:
            RAW = zlib.decompress(self.root.goto(
                'Meta/SI Image[0]/intensdata').value)
            data = struct.unpack("<{0}I".format(self.Width*self.Height), RAW)
            self.img = np.array(data).reshape((self.Height, self.Width))
        except:
            import warnings
            warnings.warn("No SI image found. Skipping it.")
            self.img = None

        self.fov = self.root.goto('Meta/SI Image[0]/fieldofview').getDouble()

    def getChannelsByName(self, name, strict=False):
        """
        Retrieve the channels for a given assignment name in the form of a list of dictionaries.
        The output can be formatted in a human readable way with the pySPM.ITA.showChannels function (see examples).

        Parameters
        ----------
        name : string or list of strings
            A regular expression (regex) used for the search
        strict : bool
            If strict is True, the search name won't be treated as a regexp, but rather the whole name should match.

        Returns
        -------
        list
            A list of dictionaries where each dictionary is a description of the selected channel. Which contains:
                - clsid : class ID. A useless information for the end-user
                - desc : a description string encoded in utf16.
                - color : a 32bits color encoding of the peak
                - symbolID : Not used
                - id : The ID of the channel. (The total counts is 0, the sum
                  of the rest 1 and the first peak is 2, ... )
                - SN : an utf16 serial number which is useless for the end-used
                - assign : a utf16 string with the element assignment of the
                  peak (e.g.: CH-, Na+, ...)
                - lmass : a long value indicating the lower mass of the peak (in u)
                - umass : a long value indicating the upper mass of the peak (in u)
                - cmass : a long value indicating the center mass of the peak
   
        Examples
        --------
        >>> A = pySPM.ITA("myfile.ita")
        >>> ch = A.getChannelsByName("C")
        >>> A.showChannels(ch)
        	CH- (), mass: 12.99 - 13.03
            C_2- (), mass: 23.97 - 24.03
            C_2H- (), mass: 24.98 - 25.04
            CN- (), mass: 25.97 - 26.04
            Cl- (), mass: 34.93 - 35.01
            C_2O- (), mass: 39.96 - 40.04
            CHNO- (), mass: 42.97 - 43.05
            CHO_2- (), mass: 44.95 - 45.04
            Cs- (), mass: 132.81 - 133.01
        >>> ch = A.getChannelsByName("C[^a-z]") # Only carbon atoms (meaning that the char after C cannot be lowercase)
        >>> A.showChannels(ch)
        	CH- (), mass: 12.99 - 13.03
            C_2- (), mass: 23.97 - 24.03
            C_2H- (), mass: 24.98 - 25.04
            CN- (), mass: 25.97 - 26.04
            C_2O- (), mass: 39.96 - 40.04
            CHNO- (), mass: 42.97 - 43.05
            CHO_2- (), mass: 44.95 - 45.04
        >>> ch = A.getChannelsByName("CH", True) # Only CH channel and not CHNO and CHO_2
        >>> A.showChannels(ch)
        	CH- (), mass: 12.99 - 13.03
        """
        res = []
        if strict:
            if type(name) in [list, tuple]:
                name = ['^'+n+'[+-]?$' for n in name]
            else:
                name = '^'+name+'[+-]?$'
        if type(name) is not list:
            name = [name]
        for n in name:
            for P in self.peaks:
                p = self.peaks[P]
                ma = re.compile(n, re.U)
                if ma.match(p['assign']['utf16']) or ma.match(p['desc']['utf16']):
                    res.append(p)
        return res

    def showChannels(self, ch):
        """
        Format a list of channels where each channel is represented by a dictionary (like the ones produced by pySPM.ITA.getChannelsByName) to a human readable output.

        Parameters
        ----------
        ch : list
            A list of dictionaries representing the channels

        Returns
        -------
        None
            It will print a list of channels with the assignment, the description in parenthesis followed by the lower - upper mass range.
        """
        for z in ch:
            print("\t{name} ({desc}), mass: {lower:.2f} - {upper:.2f}"
                  .format(desc=z['desc']['utf16'], name=z['assign']['utf16'],
                          lower=z['lmass']['float'], upper=z['umass']['float']))

    def getChannelByMass(self, mass, full=False):
        """
        Retrieves the first channel ID which has a mass range containing a given mass.

        Parameters
        ---------
        mass : int, float
            The mass. If zero, the channel 0 will be returned and correspond to the Total count channel.
        full : bool
            If True, not only the ID is retrieved but the whole dictionary similarly as with pySPM.ITA.getChannelsByName

        Returns
        -------
        int
            The first channel ID containing the mass given in argument. If a mass 0 is given, the output will be 0 which corresponds to the total count channel.
        """
        if mass == 0:
            return 0
        for P in self.peaks:
            p = self.peaks[P]
            
            if p['id']['long'] > 1 and p['lmass']['float'] <= mass and mass <= p['umass']['float']:
                if full:
                    return p
                return p['id']['long']
        raise ValueError('Mass {:.2f} Not Found'.format(mass))

    def getShiftCorrectedImageByName(self, names, **kargs):
        """
        Retrieve the drift corrected (or shift corrected) image for the sum of all channels matching a given name. The shift correction applied is the one saved in the ITA file.

        Parameters
        ---------
        names : string or list of strings
            A channel name of a list of channel names to be summed up

        Returns
        -------
        pySPM.SPM.SPM_image
            The image of the sum of all the selected channels
        list of dictionaries
            The list of all the channels selected. This list can be displayed in a human readable form by the pySPM.ITA.showChannels function

        """
        return self.getSumImageByName(names, Shifts=[(-x,-y) for x,y in self.getSavedShift()],**kargs)
        
    def __getSumImage(self, scans, channels, **kargs):
        """
        An internal function to retrieve the sum of several scans for several channel ID.
        """
        Z = np.zeros((self.sy, self.sx))
        if 'Shifts' in kargs:
            Shifts = kargs['Shifts']
        else:
            Shifts = [(-x,-y) for x,y in self.getSavedShift()]            
        for ch in channels:
            ID = ch['id']['long']
            Z += self.fastGetImage(ID, scans, Shifts)
        return Z
       
    def getSumImageByName(self, names, scans=None, strict=False, prog=False, raw=False, **kargs):
        """
        Retrieve the image for the sum of several scans and channels selected by their channel name.

        Parameters
        ----------
        names : string or list of strings
            Similar as for pySPM.ITA.getChannelsByName
        scans : int, list of int or None
            The list of the scan number to be summed up. For the case of None (default) all the available scans are taken.
        strict : bool
            Is the name selection strict? (see pySPM.ITA.getChannelsByName)
        prog : bool
            If True a progressbar will be displayed to show the summing progress as this might be quite slow.
        raw : bool
            If True a numpy array will be returned instead of a pySPM.SPM.SPM_image
        """
        if scans is None:
            scans = range(self.Nscan)
        if type(scans) == int:
            scans = [scans]
        
        channels = self.getChannelsByName(names, strict)
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm
            scans = tqdm(scans)
        Z = self.__getSumImage(scans, channels)
        if raw:
            return Z, channels
        channel_title = ",".join([z['assign']['utf16'] for z in channels])
        return self.image(np.flipud(Z), channel=channel_title), channels

    def show(self, ax=None):
        """
        Shows the total SI image with the indication of the field of view.

        Parameters
        ----------
        ax : matplotlib axis or None
            The axis in which the image will be shown. If None the current axis will be used (ax = plt.gca())

        Returns
        -------
        None
        """
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.imshow(self.img, extent=(0, self.fov*1e6, 0, self.fov*1e6))
        ax.set_title("Total SI")
        ax.set_xlabel("x [$\mu$m]")
        ax.set_ylabel("y [$\mu$m]")

    def getShiftsByMass(self, masses, centered=True, prog=False, Filter=None):
        """
        Deprecated. A relic function that the developer is not even sure what it was supposed to do ;)
        """
        Shifts = [(0, 0)]
        if Filter is None:
            Filter = lambda z: z
        S0 = Filter(self.getSumImageByMass(masses, 0))
        Y = range(1, self.Nscan)
        if prog:
            from tqdm import tqdm
            Y = tqdm(Y)
        for i in Y:
            S = Filter(self.getSumImageByMass(masses, i))
            Shift = np.real(np.fft.fftshift(np.fft.ifft2(
                np.conj(np.fft.fft2(S0)) * np.fft.fft2(S))))
            cord = np.unravel_index(np.argmax(Shift), S0.shape)
            trans = (cord[1]-S0.shape[1]/2, cord[0]-S0.shape[0]/2)
            Shifts.append(trans)
        if centered:
            avSx = np.round(np.mean([z[0] for z in Shifts]))
            avSy = np.round(np.mean([z[1] for z in Shifts]))
            Shifts = [(z[0]-avSx, z[1]-avSy) for z in Shifts]
        return Shifts

    def getXsectionByMass(self, x1, y1, x2, y2, masses, N=None, prog=False, ax=None, col='w-', **kargs):
        """
        Retrieves a Cross-Section for a given mass along the profile determined by coordinates (x1,y1) and (x2,y2).
        The output is a 2D image where the x-axis correspond to the position along the profile and the y-axis the scan number.

        Parameters
        ----------
        x1 : int
        y1 : int
        x2 : int
        y2 : int
            profile coordinates in pixel: (x1,y1) -> (x2,y2)
        masses : int, float, list of floats
            The masse or list of masses to sum
        N : int or None
            The number of value used along the profile (which will be interpolated).
            None (default) will take the roundest number of values closest to the pixel length of the profile
        prog : bool
            If True display a progressbar
        ax : None or matplotlib axis
            if not None, the axis representing the 2D image can be given in order to display the profile's position
        col : string (matplotlib color format)
            The color of the profile used in case ax is given
        **kargs : arguments
            All supplementary arguments are passed to the pySPM.ITA.getSumImageByMass

        Returns
        -------
        np.ndarray
            2D numpy array containing the sum of all channels. The values are the count number
        """
        if True:
            y1 = self.Height-1-y1
            y2 = self.Height-1-y2
            
        if N is None:
            N = int(np.sqrt((x2-x1)**2+(y2-y1)**2))+1
        x = np.linspace(x1, x2, N)
        y = np.linspace(y1, y2, N)
        out = np.zeros((self.Nscan, N))
        Y = range(self.Nscan)
        if ax is not None:
            ax.plot([x1, x2], [y1, y2], col)
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm as tqdm
            Y = tqdm(Y)
        from scipy.ndimage import map_coordinates
        for s in Y:
            Z = self.getSumImageByMass(masses, s, **kargs)
            P = map_coordinates(Z.pixels, np.vstack((y, x)))
            out[s, :] = P
        return out

    def getAddedImageByName(self, names, strict=False, raw=False, **kargs):
        """
        Retrieve the image for the sum of all scan (precomputed by iontof, but not shift-corrected) for given names

        Parameters
        ----------
        names : string or list of strings
            name of the channel (see pySPM.ITA.getChannelsByName)
        strict : bool
            If True the names are the exact names (see pySPM.ITA.getChannelsByName)
        raw : bool
            If True a 2D numpy array will be returned
        **kargs: supplementary arguments
            passed to pySPM.ITA.getAddedImage

        Returns
        -------
        pySPM.SPM.SPM_image
            Image of the result
        list of dictionaries
            List of all selected peaks used to compute the image.
            Note: Pass this list to pySPM.ITA.showChannels in order to print a human readable representation of it.
        """
        Z = np.zeros((self.sy, self.sx))
        channels = self.getChannelsByName(names, strict)
        for ch in channels:
            ID = ch['id']['long']
            Z += self.getAddedImage(ID, **kargs)
        ch = self.get_masses(channels)
        if raw:
            return Z, ch
        return self.image(np.flipud(Z), channel=",".join([z['assign'] for z in ch])), ch

    def getSavedShift(self):
        """
        getSavedShift returns the shifts saved with the file. Usually this is the shift correction you perform with the IonToF software.

        Returns
        -------
        List of tuples
            each tuple is a (Δx,Δy) in pixels (one for each scan).
        """
        try:
            X = zlib.decompress(self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data'
                                           '/ImageStackScans/ShiftCoordinates/ImageStack.ShiftCoordinates').value)
        except:
            return [(0,0) for x in range(self.Nscan)]
        D = struct.unpack('<'+str(len(X)//4)+'i', X)
        dx = D[::2]
        dy = D[1::2]
        return list(zip(dx, dy))
        
    def getShiftCorrectedImageByMass(self, masses, **kargs):
        """
        Shortcut function for pySPM.ITA.getSumImageByMass using the saved shift corrections.
        """
        return self.getSumImageByMass(masses, Shifts=[(-x,-y) for x,y in self.getSavedShift()], **kargs)
        
    def getSumImageByMass(self, masses, scans=None, prog=False, raw=False, **kargs):
        """
        Similar to pySPM.ITA.getSumImageByName but instead of the names, the mass or list of mass is provided
        see pySPM.ITA.getSumImageByName for more details
        """
        if scans is None:
            scans = range(self.Nscan)
        if type(scans) is int:
            scans = [scans]
        if type(masses) is int or type(masses) is float:
            masses = [masses]
        if prog:
            if in_ipynb():
                try:
                    from tqdm import tqdm_notebook as tqdm
                except:
                    from tqdm import tqdm
            else:
                from tqdm import tqdm
            scans = tqdm(scans, leave=False)
        channels = [self.getChannelByMass(m, full=True) for m in masses]
        Z = self.__getSumImage(scans, channels, **kargs)
        if raw:
            return Z, channels
        channels_name = [["{:.2f}u".format(m['cmass']['float']),m['assign']['utf16']][m['assign']['utf16']!=''] for m in channels]
        return self.image(np.flipud(Z), channel="Masses: "+",".join(channels_name))

    def image(self, I, channel="Unknown", zscale="Counts"):
        """
        Create a pySPM.SPM.SPM_image for a given numpy array with the same real size information as the tof-sims data.

        Parameters
        ----------
        I : numpy 2D array
            A given array
        channel : string
            A channel name describing the image. It will be printed as title when the SPM_image will be displayed.
        zscale : string
            The unit of the zscale. By default it's "Counts"

        Returns
        -------
        pySPM.SPM.SPM_image
            A SPM_image created with the data of a given array.

        Example
        -------
        >>> A = pySPM.ITA("myfile.ita")
        >>> Au,_ = A.getAddedImageByName("Au") # retrieve the gold channel (total counts)
        >>> Au_tofcorr = A.image(-np.log(1-np.fmin(.999, Au.pixels(A.Nscan))), "Au", zscale="yield") # Create a new image with the tof-corrected data
        """
        return SPM_image(I, real=self.size['real'], _type="TOF", zscale=zscale, channel=channel)

    def getAddedImageByMass(self, masses, raw=False, **kargs):
        """
        Retrieve the image for the sum of all scan (precomputed by iontof, but not shift-corrected) for (a) given masse(s)

        Parameters
        ----------
        masses : float or list of float
            mass of the channels to be used
        raw : bool
            If True a 2D numpy array will be returned
        **kargs: supplementary arguments
            passed to pySPM.ITA.getAddedImage

        Returns
        -------
        pySPM.SPM.SPM_image
            Image of the result
        list of dictionaries
            Only returned if raw is True
            List of all selected peaks used to compute the image.
            Note: Pass this list to pySPM.ITA.showChannels in order to print a human readable representation of it.
        """
        if type(masses) is int or type(masses) is float:
            masses = [masses]
        Z = np.zeros((self.sy, self.sx))
        channels = []
        for m in masses:
            ch = self.getChannelByMass(m)
            m = self.get_masses()[ch]
            if m['assign'] != '':
                channels.append(m['assign'])
            else:
                channels.append("{cmass:.2f}u".format(**m))
            Z += self.getAddedImage(ch, **kargs)
        if raw:
            return Z, channels
        return self.image(np.flipud(Z), channel=",".join(channels))

    def getAddedImage(self, channel, **kargs):
        """
        Retrieve the numpy 2D array of a given channel ID for the sum of all scan (precomputed by iontof, but not shift-corrected)
        Note: It is preferable to use the pySPM.ITA.getAddedImageByMass or pySPM.ITA.getAddedImageByName
        """
        assert type(channel) is int
        assert channel >= 0 and channel < self.Nimg
        c = self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded'
                           '/Image['+str(channel)+']/ImageArray.Long')
        D = zlib.decompress(c.value)
        V = np.array(struct.unpack('<'+str(self.sx*self.sy)+'I', D),
                     dtype=np.float).reshape((self.sy, self.sx))
        return V
    
    def fastGetImage(self, channel, scans, Shifts=False, prog=False):
        """
        Retieve a 2D numpy array corresponding to a given channel ID for given scan(s) and return their sum.

        Parameters
        ----------
        channel : int
            The channel ID
        scans: int, list of ints or 1D numpy array
            List of scans
        shifts : False or list of tuples
            List of the shift correction in pixels for ALL the scans ( not only the selected ones).
            If Flase not shift correction is performed
        prog : bool
            Display a progressbar ?

        Returns
        -------
        2D numpy array
            array data of the image
        """
        Z = np.zeros((self.sy, self.sx))
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
                scans = tqdm(scans)
            except:
                warning.warn("tqdm_notebook not available")
                try:
                    from tqdm import tqdm
                    scans = tqdm(scans)
                except:
                    warning.warn("cannot load tqdm library")
            
        im_root =  self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image['+str(channel)+']')
        for scan in scans:
            c = im_root.goto('ImageArray.Long['+str(scan)+']')
            D = zlib.decompress(c.value)
            V = np.array(struct.unpack('<'+str(self.sx*self.sy)+'I', D),
                            dtype=np.float).reshape((self.sy, self.sx))
            if Shifts:
                r = [int(z) for z in Shifts[scan]]
                V = np.roll(np.roll(V, -r[0], axis=1), -r[1], axis=0)
                rx = [max(0,-r[0]), self.sx-max(1,r[0])]
                ry = [max(0,-r[1]), self.sy-max(1,r[1])]
                Z[ry[0]:ry[1],rx[0]:rx[1]] += V[ry[0]:ry[1],rx[0]:rx[1]]
            else:
                Z += V
        return Z
        
    def getImage(self, channel, scan, Shifts=None, ShiftMode='roll', const=0):
        """
        getImage retrieve the image of a specific channel (ID) and a specific scan.

        Parameters
        ----------
        channel : int
            channel ID
        scan : int
            scan number (start with 0)
        Shifts : None or list of tuples
                None: No shift
                list of tuple in the form of where each tuple is in the form (dx,dy) for a given scan
        ShiftMode :  string
            roll : roll the data over. easy but non-physical
            const : replace missing values by a constant (given by argument const)
            NaN : the same as const but with const=numpy.NaN
        const : float
            if ShiftMode is 'const' then this parameter defines the constant used (default 0)
        """
        assert type(channel) is int
        assert type(scan) is int
        assert channel >= 0 and channel < self.Nimg
        assert scan >= 0 and scan < self.Nscan
        c = self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans'
                           '/Image['+str(channel)+']/ImageArray.Long['+str(scan)+']')
        D = zlib.decompress(c.value)
        V = np.array(struct.unpack('<'+str(self.sx*self.sy)+'I', D),
                     dtype=np.float).reshape((self.sy, self.sx))
        if not Shifts is None:
            r = [int(z) for z in Shifts[scan]]
            V = np.roll(np.roll(V, -r[0], axis=1), -r[1], axis=0)
            if ShiftMode == 'const' or ShiftMode == 'NaN':
                if ShiftMode == 'NaN':
                    const = np.nan
                if r[1] < 0:
                    V[:-r[1], :] = const
                elif r[1] > 0:
                    V[-r[1]:, :] = const
                if r[0] < 0:
                    V[:, :-r[0]] = const
                elif r[0] > 0:
                    V[:, -r[0]:] = const
        return V
        
    def showSpectrumAround(self, m0, delta=0.15, sf=None, k0=None, **kargs):
        """
        Display the Spectrum around a given mass.

        Parameters
        ----------
        m0 : float
            The central mass around which the spectrum will be plotted (in u)
        delta : float
            The spectrum will be plotted between m0-delta and m0+delta
        sf : float or None
        k0 : float or None
            sf and k0 are the mass calibration parameters. If None values saved with the file will be used.
        **kargs : supplementary arguments
            Passed to pySPM.utils.showPeak
        """
        polarity = '+'
        if self.getValue('Instrument.Analyzer_Polarity_Switch')['string'] == 'Negative':
            polarity = '-'
        from . import utils
        m, D = self.getSpectrum(sf=sf,k0=k0)
        return utils.showPeak(m,D,m0,delta,polarity=polarity, **kargs)

class ITA_collection(Collection):
    """
    ITA_collection is a super class containing a collection of tof-sims images.
    for details on Collection see pySPM.collection.Collection
    """
    def __init__(self, filename, channels1=None, channels2=None, name=None, mass=False, strict=False):
        """
        Opening a ToF-SIMS ITA file as an image collection

        Parameters
        ----------
        filename : string
            The filename
        channels1 : None or a list of names
        channels2 : None or a list of names
            channels1 and channels2 can be list of names or masses if mass=True
        name : string or None
            Name of the collection. If None, the basename of the filename is used (e.g. path/myfile.ita => name=myfile)
        mass : bool
            if True the channel lists are in mass and not names
        strict : bool
            Is the channel name strict? (see pySPM.ITA.getChannelsByName)

        Returns
        -------
        pySPM.ITA_collection class
        """
        self.ita = ITA(filename)
        self.filename = filename
        self.PCA = None
        if name is None:
            name = os.path.basename(filename)
        self.name = name
        Collection.__init__(self, sx=self.ita.fov, sy=self.ita.fov*self.ita.sy/self.ita.sx,
                            unit='m', name=name)
        self.msg = ""
        if channels1 is None:
            mass = True
            masses = self.ita.get_masses()
            channels1 = [x['cmass'] for x in masses if x['id'] > 1]
        CHS = [channels1]
        if channels2 is not None:
            CHS.append(channels2)
        for channels in CHS:
            if channels is channels2:
                strict = False
            if type(channels) is list:
                for x in channels:
                    if mass:
                        try:
                            I = self.ita.getAddedImageByMass(x)
                            m = masses[2+channels1.index(x)]
                            if m['assign'] != '':
                                self.add(I, m['assign'])
                            else:
                                self.add(I, "{cmass:.2f}u".format(cmass=x))
                        except:
                            pass
                    else:
                        Z, ch = self.ita.getAddedImageByName(x, strict)
                        self.msg += "{0}\n".format(x)
                        for z in ch:
                            self.msg += "\t{name} ({desc}), mass: {lower:.2f} - {upper:.2f}\n"\
                                .format(desc=z['desc'], name=z['assign'],
                                        lower=z['lmass'], upper=z['umass'])
                        self.add(Z, x)
            elif type(channels) is dict:
                for x in channels:
                    if mass:
                        self.add(self.ita.getAddedImageByMass(channels[x]), x)
                    else:
                        Z, ch = self.ita.getAddedImageByName(
                            channels[x], strict)
                        self.msg += "{0}\n".format(x)
                        for z in ch:
                            self.msg += "\t{name} ({desc}), mass: {lower:.2f} - {upper:.2f}\n"\
                                .format(desc=z['desc'], name=z['assign'],
                                        lower=z['lmass'], upper=z['umass'])
                        self.add(Z, x)
            else:
                raise TypeError(
                    "Channels should be a list or a dictionnary. Got {}".format(type(channels)))

    def __getitem__(self, key):
        """
        Retrieve the image of a given channel

        Example
        -------
        >>> A = pySPM.ITA_collection("myfile.ita")
        >>> A['Au-']
        <pySPM.SPM.SPM_image at 0x????????>
        """
        if key not in self.channels:
            return None
        return self.channels[key]

    def runPCA(self, channels=None):
        """
        Perform a Principle Component Analysis (PCA) on the channels

        Parameters
        ----------
        channels : None or list of strings
            List of channels to use for the PCA. If None all channels will be used.
        """
        from pySPM import PCA
        if channels is None:
            channels = self.channels.keys()
        self.PCA = PCA.ITA_PCA(self, channels)

    def showPCA(self, num=None, **kargs):
        """
        Run PCA if not already done and display the PCA images.

        Parameters
        ----------
        num : int or None
            The number of PC component to display. If None display all PAs
        **kargs : additional parameters
            passed to pySPM.PCA.showPCA
        
        Returns
        -------
        None
            Plot num PCA into a 1×num subplots

        """
        if self.PCA is None:
            self.runPCA()
        self.PCA.showPCA(num=num, **kargs)

    def loadings(self, num=None):
        """
        Return a pandas DataFrame with the num first loadings

        Parameters
        ----------
        num : int or None
            The number of PC to use. If None use all PCs

        Note
        ----
        The results can be used in combination with pySPM.PCA.hinton to create nice hinton plots
        >>> col = pySPM.ITA_collection("myfile.ita")
        >>> L = col.loadings(3)
        >>> col.PCA.hincton(matrix=L)
        Display a hinton plot with num lines representing the strength of each loading. Blue means negative loadings and Red means positive ones.
        The size of each square is proportional to the absolute value of each loading.
        """
        if self.PCA is None:
            self.runPCA()
        if num is None:
            return self.PCA.loadings()
        return self.PCA.loadings()[:num]

    def StitchCorrection(self, channel, stitches, gauss=0, debug=False):
        """
        When an image is created by stitching of several images (while moving the stage during the measurement) the resulting image can have several artifacts due to charging.
        The goal of this function is the try to suppress this stitching artifacts by givings a channel name which is known to be homogeneous everywhere

        Parameters
        ----------
        channel : string
            name of a channel with a known homogeneous yield (i.e. where the visible variation of the yield is only due to charging and not to a material density variation
        stitches : list or tuple of two ints
            stitches=(N,M) where N×M is the numer of images stitched
        gauss : float
            if >0 a gauss filter will be applied on the reference image
        debug : bool
            if True returns additionally to the new collection also the reference image

        Returns
        -------
        pySPM.ITA_collection
            A new collection with corrected data

        """
        import copy
        from scipy.ndimage.filters import gaussian_filter
        N = ITA_collection(self.filename, [], name=self.name)
        size = list(self.channels.values())[0].pixels.shape
        S = np.zeros((int(size[0]/stitches[0]), int(size[1]/stitches[1])))
        sy, sx = S.shape
        for i in range(stitches[0]):
            for j in range(stitches[1]):
                S += self.channels[channel].pixels[sy*i:sy*(i+1), sx*j:sx*(j+1)]
        S[S == 0] = 1
        if gauss>0:
            S = gaussian_filter(S, gauss)
        for x in self.channels:
            F = np.zeros(size)
            for i in range(stitches[0]):
                for j in range(stitches[1]):
                    F[sy*i:sy*(i+1), sx*j:sx*(j+1)] = \
                        self.channels[x].pixels[sy*i:sy*(i+1), sx*j:sx*(j+1)]/S
            new_channel = copy.deepcopy(self[x])
            new_channel.pixels = F
            N.add(new_channel, x)
        if debug:
            return N, S
        return N
