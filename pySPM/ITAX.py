import zlib
import struct
import numpy as np
from .Block import Block

class ITAX:
    def __init__(self, filename):
        self.f = open(filename, 'rb')
        self.f.read(8)
        self.root = Block(self.f)
        self.meas_options = self.root.goto('CommonDataObjects/MeasurementOptions/*/pickle').unpickle()
        self.size = {'pixels':dict(x=self.meas_options['raster_resolution'], y=self.meas_options['raster_resolution']), 'real':dict(x=self.meas_options['raster_fov'], y=self.meas_options['raster_fov'], unit='m')}

    def getSnapshots(self):
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
        from .utils import sp
        s = self.getSnapshots()
        ax = sp(len(s))
        for i,S in enumerate(s):
            ax[i].imshow(S)
