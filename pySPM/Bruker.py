# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Module to handle SPM images recorded by a Bruker AFM
"""

import struct
import re
import numpy as np
import pySPM
import warnings


class Bruker:
    """
    Class to handle SPM images recorded by a Bruker AFM
    """

    def __init__(self, path):
        self.path = path
        self.file = open(self.path, 'rb')
        self.layers = []
        self.scanners = []
        mode = ''
        while True:
            line = self.file.readline().rstrip().replace(b'\\', b"")
            if line == b'*Ciao image list':
                self.layers.append({})
                mode = 'Image'
            elif line == b'*Scanner list':
                self.scanners.append({})
                mode = 'Scanner'
            else:
                args = line.split(b": ")
                if len(args) > 1:
                    if mode == 'Image':
                        self.layers[-1][args[0]] = args[1:]
                    elif mode == 'Scanner':
                        self.scanners[-1][args[0]] = args[1:]
                if line == b"*File list end":
                    break

    def _get_raw_layer(self, i):
        """
        Internal function to retrieve raw data of a layer
        """
        off = int(self.layers[i][b'Data offset'][0])
        cols = int(self.layers[i][b'Number of lines'][0])
        rows = int(self.layers[i][b'Samps/line'][0])
        byte_length = int(self.layers[i][b'Data length'][0])
        length = rows*cols
        bpp = byte_length // length
        byte_length = length * bpp
        
        self.file.seek(off)
        return np.array(
            struct.unpack("<"+str(length)+{2:'h',4:'i',8:'q'}[bpp], self.file.read(byte_length)),
            dtype='float64').reshape((cols, rows))

    def list_channels(self, encoding='latin1'):
        print("Channels")
        print("========")
        for x in [z[b'@2:Image Data'][0] for z in self.layers]:
            print("\t"+x.decode(encoding))

    def get_channel(self, channel="Height Sensor", backward=False, corr=None, debug=False, encoding='latin1', lazy=True):
        """
        Load the SPM image contained in a channel
        """
        for i in range(len(self.layers)):
            layer_name = self.layers[i][b'@2:Image Data'][0].decode(encoding)
            result = re.match(
                r'([^ ]+) \[([^\]]*)\] "([^"]*)"', layer_name).groups()
            if result[2] == channel:
                if debug:
                    print("channel "+channel+" Found!")
                bck = False
                if self.layers[i][b'Line Direction'][0] == b'Retrace':
                    bck = True
                if bck == backward:
                    if debug:
                        print("Direction found")
                    var = self.layers[i][b'@2:Z scale'][0].decode(encoding)
                    if debug:
                        print("@2:Z scale",var)
                    if '[' in var:
                        result = re.match(r'[A-Z]+\s+\[([^\]]+)\]\s+\(-?[0-9\.]+ .*?\)\s+(-?[0-9\.]+)\s+(.*?)$', var).groups()
                        if debug:
                            print(result)
                        scale = float(result[1])/65536.0
                        
                        result2 = self.scanners[0][b'@'+result[0].encode(encoding)][0].split()
                        if debug:
                            print("result2", result2)
                        scale2 = float(result2[1])
                        if len(result2)>2:
                            zscale = result2[2]
                        else:
                            zscale = result2[0]
                        if b'/V' in zscale:
                            zscale = zscale.replace(b'/V',b'')
                        if debug:
                            print("scale: {:.3e}".format(scale))
                            print("scale2: {:.3e}".format(scale2))
                            print("zscale: "+str(zscale))
                        var = self.layers[i][b'@2:Z offset'][0].decode(encoding)
                        result = re.match(r'[A-Z]+\s+\[[^\]]+\]\s+\(-?[0-9\.]+ .*?\)\s+(-?[0-9\.]+)\s+.*?$', var).groups()
                        offset = float(result[0])
                    else:
                        result = re.match(r'[A-Z]+ \(-?[0-9\.]+ [^\)]+\)\s+(-?[0-9\.]+) [\w]+', var).groups()
                        scale = float(result[0])/65536.0
                        scale2 = 1
                        zscale = b'V'
                        result = re.match(r'[A-Z]+ \(-?[0-9\.]+ .*?\)\s+(-?[0-9\.]+) .*?', self.layers[i][b'@2:Z offset'][0].decode(encoding)).groups()
                        offset = float(result[0])
                    data = self._get_raw_layer(i)*scale*scale2

                    scan_size = self.layers[i][b'Scan Size'][0].split()
                    if scan_size[2][0] == 126:
                        scan_size[2] = b'u'+scan_size[2][1:]
                    size = {
                        'x': float(scan_size[0]),
                        'y': float(scan_size[1]),
                        'unit': scan_size[2].decode(encoding)}
                    image = pySPM.SPM_image(
                        channel=[channel,'Topography'][channel=='Height Sensor'],
                        BIN=data,
                        real=size,
                        _type='Bruker AFM',
                        zscale=zscale.decode(encoding),
                        corr=corr)
                    return image
        if lazy:
            return self.get_channel(channel=channel,backward=not backward, corr=corr, debug=debug, encoding=encoding, lazy=False)
        raise Exception("Channel {} not found".format(channel))
