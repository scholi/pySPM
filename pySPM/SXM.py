# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
This module handle sxm file format used by Nanonis instruments
"""

import os
# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

import struct
import numpy as np
from .SPM import SPM_image

class SXM:
    def __init__(self, filename):
        self.filename = filename
        assert os.path.exists(self.filename)
        self.f = open(self.filename, 'rb')
        l = ''
        key = ''
        self.header = {}
        while l!=b':SCANIT_END:':
            l = self.f.readline().rstrip()
            if l[:1]==b':':
                key = l.split(b':')[1].decode('ascii')
                self.header[key] = []
            else:
                if l: # remove empty lines
                    self.header[key].append(l.decode('ascii').split())
        while self.f.read(1)!=b'\x1a':
            pass
        assert self.f.read(1)==b'\x04'
        assert self.header['SCANIT_TYPE'][0][0] in ['FLOAT','INT','UINT','DOUBLE']
        self.data_offset = self.f.tell()
        self.size = dict(pixels={
                'x': int(self.header['SCAN_PIXELS'][0][0]),
                'y': int(self.header['SCAN_PIXELS'][0][1])
            },real={
                'x': float(self.header['SCAN_RANGE'][0][0]),
                'y': float(self.header['SCAN_RANGE'][0][1]),
                'unit': 'm'
            })

    def list_channels(self):
        print("Channels")
        print("========")
        h = self.header['DATA_INFO'][0]
        i = h.index('Name')
        for z in self.header['DATA_INFO'][1:]:
            print("  - "+z[i])

    def get_channel(self, name, direction='forward',corr=None):
        chID = 0
        zscale = ''
        for x in self.header['DATA_INFO'][1:]:
            if x[1]==name:
                if x[3]=='both' and direction=='backward':
                    chID += 1
                if x[3]== 'both' or direction == x[3]:
                    break
                return None
                zscale = x[2]
            elif x[3] == 'both':
                chID += 2
            else:
                chID += 1
                    
        size = self.size['pixels']['x']*self.size['pixels']['y']
        self.f.seek(self.data_offset+chID*size*4)
        data = np.array(struct.unpack('<>'['MSBFIRST'==self.header['SCANIT_TYPE'][0][1]]+str(size)+{'FLOAT':'f','INT':'i','UINT':'I','DOUBLE':'d'}[self.header['SCANIT_TYPE'][0][0]],self.f.read(4*size))).reshape((self.size['pixels']['y'],self.size['pixels']['x']))
        return SPM_image(
            channel=name,
            BIN=data,
            real=self.size['real'],
            _type = 'Nanonis SXM',
            zscale = zscale,
            corr=corr)
