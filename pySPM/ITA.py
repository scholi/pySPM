"""
This module gives the ability to ready and parse the ITA files.
You can mainly retrieve the images for each channel and scan.
"""

import numpy as np
import struct
import os.path
import zlib
import re

from .ITM import ITM
from .collection import Collection
from .SPM import SPM_image

class ITA(ITM):
    def __init__(self, filename):
        ITM.__init__(self, filename)
        try:
            self.sx = self.root.goto(
                'filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.XSize').getLong()
            self.sy = self.root.goto(
                'filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.YSize').getLong()
        except ValueError:
            self.sx = self.size['pixels']['x']
            self.sy = self.size['pixels']['y']
        try:
            # self.Nscan = int(self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data'\
            #    '/ImageStackScans/Image.NumberOfScans').getLong())
            self.Nimg = int(self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data'
                                           '/ImageStackScans/Image.NumberOfImages').getLong())
        except:
            try:
                self.Nimg = int(self.root.goto('propend/Measurement.ScanNumber').getKeyValue(16)['Value'])
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
        res = []
        if strict:
            name = '^'+name+'[+-]?$'
        if type(name) is not list:
            name = [name]
        for n in name:
            for P in self.peaks:
                p = self.peaks[P]
                ma = re.compile(n, re.U)
                if ma.match(p[b'assign']['utf16']) or ma.match(p[b'desc']['utf16']):
                    res.append(p)
        return res

    def showChannels(self, ch):
        for z in ch:
            print("\t{name} ({desc}), mass: {lower:.2f} - {upper:.2f}"
                  .format(desc=z[b'desc']['utf16'], name=z[b'assign']['utf16'],
                          lower=z[b'lmass']['float'], upper=z[b'umass']['float']))

    def getChannelByMass(self, mass, full=False):
        if mass == 0:
            return 0
        for P in self.peaks:
            p = self.peaks[P]
            if p[b'id']['long'] > 1 and p[b'lmass']['float'] <= mass and mass <= p[b'umass']['float']:
                if full:
                    return p
                return p[b'id']['long']
        raise ValueError('Mass {:.2f} Not Found'.format(mass))

    def getShiftCorrectedImageByName(self, names, **kargs):
        return self.getSumImageByName(names, Shifts=[(-x,-y) for x,y in self.getSavedShift()],**kargs)
        
    def __getSumImage(self, scans, channels, **kargs):
        Z = np.zeros((self.sy, self.sx))
        if 'Shifts' in kargs:
            Shifts = kargs['Shifts']
        else:
            Shifts = [(-x,-y) for x,y in self.getSavedShift()]            
        for ch in channels:
            ID = ch[b'id']['long']
            Z += self.fastGetImage(ID, scans, Shifts)
        return Z
       
    def getSumImageByName(self, names, scans=None, strict=False, prog=False, raw=False, **kargs):
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
        channel_title = ",".join([z[b'assign']['utf16'] for z in channels])
        return self.image(np.flipud(Z), channel=channel_title), channels

    def show(self, ax=None):
        """
        Shows the total SI image with the indication of the field of view.
        ax (=None): if you provide an ax argument, the image can be plottet in the axis of your choice
        """
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.imshow(self.img, extent=(0, self.fov*1e6, 0, self.fov*1e6))
        ax.set_title("Total SI")
        ax.set_xlabel("x [$\mu$m]")
        ax.set_ylabel("y [$\mu$m]")

    def getShiftsByMass(self, masses, centered=True, prog=False, Filter=None):
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
            from tqdm import tqdm
            Y = tqdm(Y)
        from scipy.ndimage import map_coordinates
        for s in Y:
            Z = self.getSumImageByMass(masses, s, **kargs)
            P = map_coordinates(Z.pixels, np.vstack((y, x)))
            out[s, :] = P
        return out

    def getAddedImageByName(self, names, strict=False, raw=False, **kargs):
        Z = np.zeros((self.sy, self.sx))
        channels = self.getChannelsByName(names, strict)
        for ch in channels:
            ID = ch[b'id']['long']
            Z += self.getAddedImage(ID, **kargs)
        ch = self.get_masses(channels)
        if raw:
            return Z, ch
        return self.image(np.flipud(Z), channel=",".join([z['assign'] for z in ch])), ch

    def getSavedShift(self):
        """
        getSavedShift returns the shifts saved with the file. Usually this is the shift correction you perform with the IonToF software.
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
        return self.getSumImageByMass(masses, Shifts=[(-x,-y) for x,y in self.getSavedShift()], **kargs)
        
    def getSumImageByMass(self, masses, scans=None, prog=False, raw=False, **kargs):
        if scans is None:
            scans = range(self.Nscan)
        if type(scans) is int:
            scans = [scans]
        if type(masses) is int or type(masses) is float:
            masses = [masses]
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm
            scans = tqdm(scans)
        channels = [self.getChannelByMass(m, full=True) for m in masses]
        Z = self.__getSumImage(scans, channels)
        if raw:
            return Z, channels
        channels_name = [["{:.2f}u".format(m[b'cmass']['float']),m[b'assign']['utf16']][m[b'assign']['utf16']!=''] for m in channels]
        return self.image(np.flipud(Z), channel="Masses: "+",".join(channels_name))

    def image(self, I, channel="Unknown"):
        return SPM_image(I, real=self.size['real'], _type="TOF", zscale='Counts', channel=channel)

    def getAddedImageByMass(self, masses, raw=False, **kargs):
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
        assert type(channel) is int
        assert channel >= 0 and channel < self.Nimg
        c = self.root.goto('filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScansAdded'
                           '/Image['+str(channel)+']/ImageArray.Long')
        D = zlib.decompress(c.value)
        V = np.array(struct.unpack('<'+str(self.sx*self.sy)+'I', D),
                     dtype=np.float).reshape((self.sy, self.sx))
        return V
    
    def fastGetImage(self, channel, scans, Shifts=False, prog=False):
        Z = np.zeros((self.sy, self.sx))
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                warning.warn("tqdm_notebook not available")
                try:
                    from tqdm import tqdm
                except:
                    warning.warn("cannot load tqdm library")
            scans = tqdm(scans)
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
        
    def getImage(self, channel, scan, Shifts=None, ShiftMode='roll', **kargs):
        """
        getImage retrieve the image of a specific channel (ID) and a specific scan.
        channel: channel ID
        scan: scan nummber (start with 0)
        Shifts: None=No shift, otherwise provide an array of tuple ((x,y) shift for each scan)
        ShiftMode:  * roll (roll the data over. easy but unphysical)
                    * const (replace missing values by a constant. given by argument const)
                    * NaN (the same as const with const=NaN)
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
                    kargs['const'] = np.nan
                if 'const' not in kargs:
                    raise KeyError('Missing argument const')
                if r[1] < 0:
                    V[:-r[1], :] = kargs['const']
                elif r[1] > 0:
                    V[-r[1]:, :] = kargs['const']
                if r[0] < 0:
                    V[:, :-r[0]] = kargs['const']
                elif r[0] > 0:
                    V[:, -r[0]:] = kargs['const']
        return V


class ITA_collection(Collection):
    def __init__(self, filename, channels1=None, channels2=None, name=None, mass=False, strict=False):
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
        if key not in self.channels:
            return None
        return self.channels[key]

    def runPCA(self, channels=None):
        from pySPM import PCA
        if channels is None:
            channels = self.channels.keys()
        self.PCA = PCA.ITA_PCA(self, channels)

    def showPCA(self, num=None, **kargs):
        if self.PCA is None:
            self.runPCA()
        self.PCA.showPCA(num=num, **kargs)

    def loadings(self, num=None):
        if self.PCA is None:
            self.runPCA()
        if num is None:
            return self.PCA.loadings()
        return self.PCA.loadings()[:num]

    def StitchCorrection(self, channel, stitches, gauss=0, debug=False):
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
