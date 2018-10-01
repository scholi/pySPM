import zlib
import struct
import numpy as np
from .Block import Block
from . import utils

class ITAX:
    def __init__(self, filename):
        self.f = open(filename, 'rb')
        self.f.read(8)
        self.root = Block(self.f)
        mo = self.root.goto('CommonDataObjects/MeasurementOptions').getList()[0]
        self.meas_options = self.root.goto('CommonDataObjects/MeasurementOptions/{name}[{id}]/pickle'.format(**mo)).unpickle()
        self.size = {'pixels':dict(x=self.meas_options['raster_resolution'], y=self.meas_options['raster_resolution']), 'real':dict(x=self.meas_options['raster_fov'], y=self.meas_options['raster_fov'], unit='m')}

    def getSnapshots(self):
        """
        Retrieve teh RGB images of the camera which are captured
        Return a list of array for all the snapshots
        """
        snapshots = []
        for x in self.root.goto('CommonDataObjects/Attachments'):
            for y in x.getList():
                if y['name'] == 'Video Snapshot':
                    self.f.seek(y['bidx'])
                    blk = Block(self.f)
                    sx = blk.goto('res_x').getLong()
                    sy = blk.goto('res_y').getLong()
                    raw = blk.goto("imagedata").value
                    data = zlib.decompress(raw)
                    I = np.flipud(np.array(struct.unpack("<"+str(3*sx*sy)+"B", data)).reshape((sy, sx, 3)))
                    snapshots.append(I)
                    del blk
        return snapshots
    
    def showSnapshots(self):
        """
        Retrieve and plot all recorded snapshots (camera images)
        """
        from .utils import sp
        s = self.getSnapshots()
        ax = sp(len(s))
        for i, S in enumerate(s):
            ax[i].imshow(S)
            
    def get_mass_cal(self):
        k0 = self.root.goto("CommonDataObjects/DataViewCollection/*/properties/Context.MassScale.K0", lazy=True).getKeyValue()['float']
        sf = self.root.goto("CommonDataObjects/DataViewCollection/*/properties/Context.MassScale.SF", lazy=True).getKeyValue()['float']
        return sf, k0

    def getSpectrum(self, sf=None, k0=None, time=False, **kargs):
        """
        Retrieve the spectrum in a similar way as for ITA file
        """
        slen = self.root.goto("CommonDataObjects/DataViewCollection/*/sizeSpectrum").getLong()
        raw = self.root.goto("CommonDataObjects/DataViewCollection/*/dataSource/simsDataCache/spectrum/correctedData").value
        spectrum = np.array(struct.unpack("<"+str(slen)+"d", raw))
        CH = 2*np.arange(slen)        
        if time:
            return CH, spectrum
        if sf is None:
            sf = self.root.goto("CommonDataObjects/DataViewCollection/*/properties/Context.MassScale.SF", lazy=True).getKeyValue()['float']
        if k0 is None:
            k0 = self.root.goto("CommonDataObjects/DataViewCollection/*/properties/Context.MassScale.K0", lazy=True).getKeyValue()['float']
        m = utils.time2mass(CH, sf, k0)
        return m, spectrum
    
    def getProfile(self, name):
        """
        retrieve the depth profile for a given channel name
        """
        SN = None
        for x in self.root.goto("CommonDataObjects/MeasurementOptions/*/massintervals"):
            if x.name == 'mi':
                v = x.dictList()
                lab = v['assign']['utf16'] or v['desc']['utf16']
                if lab == name:
                    SN = v['SN']['utf16']
                    break
        if SN is None:
            raise Exception("Profile \"{}\" not found".format(name))
        path = "CommonDataObjects/DataViewCollection/*/dataSource/simsDataCache/{SN}/profile".format(SN=SN)
        raw = self.root.goto(path, lazy=True).decompress()
        return struct.unpack("<"+str(len(raw)//8)+"d", raw)