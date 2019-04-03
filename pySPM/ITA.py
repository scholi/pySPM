# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
This module gives the ability to ready and parse the ToF-SIMS ITA files from iontof.
You can mainly retrieve images and spectra for each channel and scan.
"""
from __future__ import absolute_import

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
from .utils.misc import deprecated, aliased, alias, PB
import warnings

@aliased
class ITA(ITM):
    def __init__(self, filename, *args, **kargs):
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
        ITM.__init__(self, filename, *args, **kargs)
        try:
            self.sx = self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.XSize').get_ulong()
            self.sy = self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.YSize').get_ulong()
        except MissingBlock:
            self.sx = self.size['pixels']['x']
            self.sy = self.size['pixels']['y']
        try:
            # self.Nscan = int(self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data'\
            #    '/ImageStackScans/Image.NumberOfScans').getLong())
            self.Nimg = int(self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data'
                                           '/ImageStackScans/Image.NumberOfImages').get_ulong())
        except:
            self.Nimg = 0
        self.img = self.get_intensity()
        
        try:
            self.fov = self.root.goto('Meta/SI Image[0]/fieldofview').get_double()
        except MissingBlock:
            self.fov = self.get_value("Registration.Raster.FieldOfView")['float']
            
    @alias("getIntensity")
    def get_intensity(self):
        """
        Retrieve the total Ion image
        """
        try:
            X, Y = self.size['pixels']['x'], self.size['pixels']['y']
            img = self.image(np.flipud(np.array(self.root.goto('Meta/SI Image/intensdata').get_data("f"), dtype=np.float32).reshape((Y, X))), channel="SI count")
        except Exception as e:
            try:
                img = self.get_added_image(0).pixels
            except:
                try:
                    img = self.get_added_image_by_SN(self.get_channel_SN("total"))
                except:
                    import warnings
                    warn.warnings("SI image cannot be retrieved")
                    return None
        return img
        
    def get_channel_SN(self, channel):
        """
        New ITA fileformat assign a serial number (SN) in the form of a UUID for each channel.
        The SN corresponding to a given channel name can be retrieved by this function.

        Parameters
        ----------
        channel : string
            The channel name assigned to a given peak
        """
        for x in  self.root.goto("MassIntervalList"):
            if x.name == 'mi':
                l = x.dict_list()
                if l['assign']['utf16'] == channel or l['desc']['utf16'] == channel:
                    return l['SN']['utf16']

        raise Exception("Channel name \"{channel}\" not found".format(channel=channel))
    
    @alias("getChannelsByName")
    def get_channels_by_name(self, name, strict=False):
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

    @deprecated("showChannels")
    def show_channels(self, ch):
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

    @alias("getChannelByMass")
    def get_channel_by_mass(self, mass, full=False):
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

    @alias("getShiftCorrectedImageByName")
    def get_shift_corrected_image_by_name(self, names, **kargs):
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
        return self.get_sum_image_by_name(names, shifts=[(-x, -y) for x, y in self.get_saved_shift()], **kargs)
        
    def __get_sum_image(self, scans, channels, **kargs):
        """
        An internal function to retrieve the sum of several scans for several channel ID.
        """
        Z = np.zeros((self.sy, self.sx))
        if 'shifts' in kargs:
            shifts = kargs['shifts']
        elif 'Shifts' in kargs:
            shifts = kargs['Shifts']
        else:
            shifts = [(-x,-y) for x,y in self.get_saved_shift()]            
        for ch in channels:
            ID = ch['id']['long']
            Z += self.fast_get_image(ID, scans, shifts)
        return Z

    @alias("getSumImageBySN")
    def get_sum_image_by_sn(self, SN, scans=None, prog=False, raw=False, **kargs):
        """
        Retrieve the image for the sum of several scans for a given channel SN.
        """
        if scans is None:
            scans = range(self.Nscan)
        if type(scans) == int:
            scans = [scans]
        if prog:
            scans= PB(scans)

        Z = np.zeros((self.sy, self.sx))
        for s in scans:
            node = self.root.goto("filterdata/TofCorrection/ImageStack/Reduced Data/Images/{SN}/ScanData/EDROff/{scan}".format(SN=SN, scan=s))
            dat = node.decompress()
            data = struct.unpack("<{}I".format(len(dat)//4), dat)
            Z += np.array(data, dtype=np.float).reshape((self.sy, self.sx))
        if raw:
            return Z
        channel = self.get_channel_by_sn(SN)
        return self.image(np.flipud(Z), channel=channel)

    @alias("getSumImageByName")
    def get_sum_image_by_name(self, names, scans=None, strict=False, prog=False, raw=False, **kargs):
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
        
        channels = self.get_added_image_by_name(names, strict)
        if prog:
            scans = PB(scans)
        Z = self.__get_sum_image(scans, channels)
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

    @alias("getShiftsByMass")
    def get_shifts_by_mass(self, masses, centered=True, prog=False, Filter=None):
        """
        Deprecated. A relic function that the developer is not even sure what it was supposed to do ;)
        """
        Shifts = [(0, 0)]
        if Filter is None:
            Filter = lambda z: z
        S0 = Filter(self.get_added_image_by_mass(masses, 0))
        Y = range(1, self.Nscan)
        if prog:
            Y = PB(Y)
        for i in Y:
            S = Filter(self.get_sum_image_by_mass(masses, i))
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

    @alias("getXsectionByMass")
    def get_xsection_by_mass(self, x1, y1, x2, y2, masses, N=None, prog=False, ax=None, flip=False, col='w-', **kargs):
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
        flip : bool
            Flip the y-coordinates?
        **kargs : arguments
            All supplementary arguments are passed to the pySPM.ITA.getSumImageByMass

        Returns
        -------
        np.ndarray
            2D numpy array containing the sum of all channels. The values are the count number
        """
        y1 = self.sy-1-y1
        y2 = self.sy-1-y2           
        if N is None:
            N = int(np.sqrt((x2-x1)**2+(y2-y1)**2))+1
        x = np.linspace(x1, x2, N)
        y = np.linspace(y1, y2, N)
        out = np.zeros((self.Nscan, N))
        Y = range(self.Nscan)
        if ax is not None:
            if not flip:
                ax.plot([x1, x2], [self.sy-1-y1, self.sy-1-y2], col)
        if prog:
            Y = PB(Y)
        from scipy.ndimage import map_coordinates
        for s in Y:
            Z = self.get_sum_image_by_mass(masses, s, **kargs)
            P = map_coordinates(Z.pixels, np.vstack((y, x)))
            out[s, :] = P
        return out

    @alias("getAddedImageByName")
    def get_added_image_by_name(self, names, strict=False, raw=False, **kargs):
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
        channels = self.get_channels_by_name(names, strict)
        for ch in channels:
            ID = ch['id']['long']
            Z += self.get_added_image(ID, **kargs)
        ch = self.get_masses(channels)
        if raw:
            return Z, ch
        return self.image(np.flipud(Z), channel=",".join([z['assign'] for z in ch])), ch

    @alias("getSavedShift")
    def get_saved_shift(self):
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
        
    @alias("getShiftCorrectedImageByMass")
    def get_shift_corrected_image_by_mass(self, masses, **kargs):
        """
        Shortcut function for pySPM.ITA.get_sum_image_by_mass using the saved shift corrections.
        """
        return self.get_sum_image_by_mass(masses, shifts=[(-x,-y) for x,y in self.get_saved_shift()], **kargs)
        
    @alias("getSumImageByMass")
    def get_sum_image_by_mass(self, masses, scans=None, prog=False, raw=False, **kargs):
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
            scans = PB(scans, leave=False)
        channels = [self.get_channel_by_mass(m, full=True) for m in masses]
        Z = self.__get_sum_image(scans, channels, **kargs)
        if raw:
            return Z, channels
        channels_name = [["{:.2f}u".format(m['cmass']['float']),m['assign']['utf16']][m['assign']['utf16']!=''] for m in channels]
        return self.image(np.flipud(Z), channel="Masses: "+",".join(channels_name))

    @alias("getAddedImageByMass")
    def get_added_image_by_mass(self, masses, raw=False, **kargs):
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
            ch = self.get_channel_by_mass(m)
            m = self.get_masses()[ch]
            if m['assign'] != '':
                channels.append(m['assign'])
            else:
                channels.append("{cmass:.2f}u".format(**m))
            Z += self.get_added_image(ch, **kargs)
        if raw:
            return Z, channels
        return self.image(np.flipud(Z), channel=",".join(channels))

    @alias("get_channel_by_sn","getChannelBySN")
    def get_channel_by_SN(self, SN):
        for node in self.root.goto("MassIntervalList"):
            if node.name == "mi":
                l = node.dict_list()
                if l['SN']['utf16']==SN:
                    name = l['assign']['utf16']
                    if not name:
                        name = l['desc']['utf16']
                    if not name:
                        name = '{:.2f}u'.format(l['cmass']['float'])
                    return name

    @alias("get_added_image_by_sn","getAddedImageBySN")
    def get_added_image_by_SN(self, SN, raw=False):
        """
        New ITA fileformat save images with their respective serial number (SN).
        This function return the image for a given SN.

        Parameters
        ----------

        SN: Serial Number of the channel
        """
        node = self.root.goto("filterdata/TofCorrection/ImageStack/Reduced Data/Images/{SN}/SumImage/EDROff".format(SN=SN))
        dat = node.decompress()
        data = struct.unpack("<{}I".format(len(dat)//4), dat)
        img = np.array(data, dtype=np.float).reshape((self.sy, self.sx))
        if raw:
            return img
        channel = self.get_channel_by_SN(SN)
        return self.image(np.flipud(img), channel=channel)

    @alias("getAddedImage")
    def get_added_image(self, channel, **kargs):
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
    
    @alias("fastGetImage")
    def fast_get_image(self, channel, scans, shifts=False, prog=False, **kargs):
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
        # Old parameter name compatibility
        if 'Shifts' in kargs:
            shifts = kargs.pop("Shifts")
            
        Z = np.zeros((self.sy, self.sx))
        if prog:
            scans = PB(scans)
            
        im_root =  self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image['+str(channel)+']')
        for scan in scans:
            c = im_root.goto('ImageArray.Long['+str(scan)+']')
            V = np.array(c.get_data('I'), dtype=np.float).reshape((self.sy, self.sx))
            if shifts:
                r = [int(z) for z in shifts[scan]]
                V = np.roll(np.roll(V, -r[0], axis=1), -r[1], axis=0)
                rx = [max(0,-r[0]), self.sx-max(1,r[0])]
                ry = [max(0,-r[1]), self.sy-max(1,r[1])]
                Z[ry[0]:ry[1],rx[0]:rx[1]] += V[ry[0]:ry[1],rx[0]:rx[1]]
            else:
                Z += V
        return Z
        
    @alias("getImage")
    def get_image(self, channel, scan, shifts=None, shift_mode='roll', const=0, **kargs):
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
        # Compatibility with old parameter names
        if 'Shifts' in kargs:
            shifts = kargs.pop("Shifts")
        if 'ShiftMode' in kargs:
            shift_mode = kargs.pop("ShiftMode")
            
        assert type(channel) is int
        assert type(scan) is int
        assert channel >= 0 and channel < self.Nimg
        assert scan >= 0 and scan < self.Nscan
        c = self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans'
                           '/Image['+str(channel)+']/ImageArray.Long['+str(scan)+']')
        V = np.array(c.get_data(), dtype=np.float).reshape((self.sy, self.sx))
        if not shifts is None:
            r = [int(z) for z in shifts[scan]]
            V = np.roll(np.roll(V, -r[0], axis=1), -r[1], axis=0)
            if shift_mode == 'const' or shift_mode == 'NaN':
                if shift_mode == 'NaN':
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
        
    @alias("getOperation")
    def get_opertion(self, OpID):
        """
        Test function to retrieve the operations used in the Worksheet.
        """
        Nop = self.root.goto('Presentation/Imaging Worksheet/Worksheet/OPERATIONS/OpCount').get_ulong()
        for i in range(Nop):
            blk = self.root.goto('Presentation/Imaging Worksheet/Worksheet/OPERATIONS/Operation[{}]'.format(i))
            if blk.goto_item('OpID').get_ulong() == OpID:
                return blk
        return None
        
    @alias("showWorksheet")
    def show_worksheet(self, page=0):
        """
        In Dev. function to display the worksheet
        """
        import matplotlib as mpl
        from .utils import sp
        num_pages = self.root.goto('Presentation/Imaging Worksheet/Worksheet/PAGES/COUNT').get_ulong()
        assert page < num_pages
        Nitems = self.root.goto('Presentation/Imaging Worksheet/Worksheet/PAGES/Page[{}]/ItemCount'.format(page)).get_ulong()
        sett = self.root.goto('Presentation/Imaging Worksheet/Worksheet/PAGES/Page[{}]/SETTINGS'.format(page)).dict_list()
        Nx = sett['Xsize']['ulong']
        Ny = sett['Ysize']['ulong']
        items = self.root.goto('Presentation/Imaging Worksheet/Worksheet/PAGES/Page[{}]/Items'.format(page)).get_data()
        ax = sp(len(items))
        IntV = {}
        for x in self.root.goto("MassIntervalList"):
            if x.name == 'mi':
                d = x.dict_list()
                IntV[d['id']['long']] = d['desc']['utf16']+d['assign']['utf16']
        for i, it in enumerate(items):
            blk = self.get_operation(it)
            OPTYPE = blk.goto_item('OPTYPE').get_ulong()
            while OPTYPE !=3:
                OPTYPE = blk.goto_item('OPTYPE').get_ulong()
                if OPTYPE == 4:
                    blk = self.getOperation(blk.goto_item('ArgOpIDs').get_ulong())
                elif OPTYPE==3:
                    palette = np.array(blk.goto_item('BMP-Palette').get_data('B')).reshape((256, 4))
                    B, G, R = palette[:, 0], palette[:, 1], palette[:, 2]
                    dimx = blk.goto('Cache/IImage-Cache-DimX').get_ulong()
                    dimy = blk.goto('Cache/IImage-Cache-DimY').get_ulong()
                    img = np.array(blk.goto('Cache/IImage-Cache-Intensities').get_data('d')).reshape((dimy, dimx))
                    RGB = np.hstack([R[:, None], G[:, None], B[:, None]])/256
                    cm = mpl.colors.ListedColormap(RGB)
                    ax[i].imshow(img, cmap=cm)

    def add_new_images(self, miblock, scans=None, added=None, prog=False, **kargs):
        # Compatibility with old parameter names
        if 'Scans' in kargs:
            scans = kargs.pop("Scans")
        if 'Added' in kargs:
            added = kargs.pop("Added")
            
        assert scans is not None or added is not None
        lvl = 3 # zlib encoding level
        sy, sx = self.size['pixels']['y'], self.size['pixels']['x']
        SN = miblock.goto("SN").get_string()
        if added is None:
            added_img = np.zeros((sy, sx), dtype=np.uint32)
        chID = miblock.goto("id").get_ulong()
        if scans is not None:
            N = self.root.goto("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.NumberOfImages").get_ulong()
        AN = self.root.goto("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image.NumberOfImages").get_ulong()
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "Image.MassIntervalSN", SN.encode('utf8'))
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "Image.XSize", struct.pack("<I", sx))
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "Image.YSize", struct.pack("<I", sy))
        if scans is not None:
            RS = range(self.Nscan)
            if prog:
                RS = PB(RS)
            for i in RS:
                img = np.flipud(scans[i].astype(np.uint32, casting='unsafe'))
                data = zlib.compress(struct.pack("<{}I".format(sx*sy), *np.ravel(img)), level=lvl)
                self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image[{}]".format(N), "ImageArray.Long", data, id=i, _type=128)
                if added is None:
                    added_img += img

        if added is None:
            added = added_img
        else:
            added = np.flipud(added)
        data = zlib.compress(struct.pack("<{}I".format(sx*sy), *np.ravel(added.astype(np.uint32, casting='unsafe'))), level=lvl)
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "ImageArray.Long", data, _type=128)
        
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "Image.PulsesPerPixel", struct.pack("<I", self.spp*self.Nscan))
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "Image.MaxCountsPerPixel", struct.pack("<I", int(np.max(added))))
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "Image.MinCountsPerPixel", struct.pack("<I", int(np.min(added))))
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "Image.TotalCountsDbl", struct.pack("<d", np.sum(added)))
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded/Image[{}]".format(AN), "Image.TotalCounts", struct.pack("<I", int(np.sum(added))))
        
        if scans is not None:
            self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans", "Image.NumberOfImages", struct.pack("<I", N+1))
        self.root.edit_block("filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded", "Image.NumberOfImages", struct.pack("<I", AN+1))
        self.Nimg += 1

@aliased
class ITA_collection(Collection):
    """
    ITA_collection is a super class containing a collection of tof-sims images.
    for details on Collection see pySPM.collection.Collection
    """
    def __init__(self, filename, channels1=None, channels2=None, name=None, mass=False, strict
=False):
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
                            I = self.ita.get_added_image_by_mass(x)
                            m = masses[2+channels1.index(x)]
                            if m['assign'] != '':
                                self.add(I, m['assign'])
                            else:
                                self.add(I, "{cmass:.2f}u".format(cmass=x))
                        except:
                            pass
                    else:
                        Z, ch = self.ita.get_added_image_by_name(x, strict)
                        self.msg += "{0}\n".format(x)
                        for z in ch:
                            self.msg += "\t{name} ({desc}), mass: {lower:.2f} - {upper:.2f}\n"\
                                .format(desc=z['desc'], name=z['assign'],
                                        lower=z['lmass'], upper=z['umass'])
                        self.add(Z, x)
            elif type(channels) is dict:
                for x in channels:
                    if mass:
                        self.add(self.ita.get_added_image_by_mass(channels[x]), x)
                    else:
                        Z, ch = self.ita.get_added_image_by_name(
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

    @alias("runPCA")
    def run_pca(self, channels=None):
        """
        Perform a Principle Component Analysis (PCA) on the channels

        Parameters
        ----------
        channels : None or list of strings
            List of channels to use for the PCA. If None all channels will be used.
        """
        from .PCA import ITA_PCA
        if channels is None:
            channels = self.channels.keys()
        self.PCA = ITA_PCA(self, channels)
    
    @alias("showPCA")
    def show_pca(self, num=None, loadings=True, **kargs):
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
            self.run_pca()
        self.PCA.show_pca(num=num, loadings=loadings, **kargs)

    def loadings(self, num=None, ax=None):
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
            self.run_pca()
        if num is None:
            L = self.PCA.loadings()
        else:
            L = self.PCA.loadings()[:num]
        if ax is not None:
            if ax is True:
                self.PCA.hinton(matrix=L)
            else:
                self.PCA.hinton(matrix=L, ax=ax)
        return L
    
    @deprecated("StitchCorrection")
    def stitch_correction(self, channel, stitches, gauss=0, debug=False):
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
