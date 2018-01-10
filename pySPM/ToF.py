# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
This is an old module which is deprecated now. It is kept in the project as it can be used to read BIF6 and BIF3D images.
"""

import numpy as np
import struct
import os
import re
import matplotlib.pyplot as plt
import pySPM

Elements = {
    'H': {1: (1.00782503223, .999885), 2: (2.01410177812, .000115)},
    'Li': {6: (6.0151228874, .0759), 7: (7.0160034366, .9241)},
    'Be': {9: (9.012183065, 1)},
    'B': {10: (10.01293695, .199), 11: (11.00930536, .801)},
    'C': {12: (12, .9893), 13: (13.00335483507, .01107)},
    'N': {14: (14.00307400443, .99636), 15: (15.00010889888, .00364)},
    'O': {16: (15.99491461957, .99757), 17: (16.99913175650, .00038), 18: (17.99915961286, .00205)},
    'F': {19: (18.9984036273, 1)},
    'Ag': {107: (106.9050916, .51839), 109: (108.9047553, .48161)},
    'Au': {197: (196.96656879, 1)},
    'Si': {28: (27.97692653465, .92223), 29: (28.97649466490, .04685), 30: (29.973770136, 0.3092)},
    'Ti': {46: (45.95262772, .0825), 47: (46.95175879, 0.0744), 48: (47.94794198, .7372), 49: (48.94786568, .0541), 50: (49.94478689, .0518)},
}


def getSpecElt(Elts):
    if len(Elts) == 1:
        return Elements[Elts[0]]
    tm = {}
    r1 = Elements[Elts[0]]
    r2 = getSpecElt(Elts[1:])
    for x in r1:
        for y in r2:
            if x+y not in tm:
                tm[x+y] = [r1[x][0]+r2[y][0], r1[x][1]*r2[y][1]]
            else:
                tm[x+y][1] += r1[x][1]*r2[y][1]
    return tm


def SplitElts(elt):
    elts = []
    for x, i in re.findall('([A-Z][a-z]?)([0-9]*)', elt):
        if i == '':
            i = 1
        else:
            i = int(i)
        elts += [x for k in range(i)]
    return elts


def showElts(Elts):
    r = getSpecElt(Elts)
    plt.bar(r.keys(), [z[1] for z in r.values()])
    plt.show()


class BIF6:
    def __init__(self, filename):
        self.f = open(filename, 'rb')
        self.header = struct.unpack('xx4s5H', self.f.read(16))
        self.size = (self.header[2], self.header[3])
        self.N = self.size[0]*self.size[1]
        self.cat = []
        for i in range(self.header[1]):
            self.f.seek(16+i*(4*self.N+16))
            self.cat.append(struct.unpack('4f', self.f.read(16)))

    def getImgID(self, ID):
        assert ID >= 0 and ID <= self.header[1]
        self.f.seek(32+ID*(4*self.N+16))
        return np.array(struct.unpack(str(self.N)+'I', self.f.read(4*self.N))).reshape(self.size)

    def getImgMass(self, masses, raw=False):
        if type(masses)is float or type(masses) is int:
            masses = [masses]
        SUM = None
        for i, x in enumerate(self.cat):
            for m in masses:
                if m >= x[0]-.5 and m <= x[1]+.5:
                    if SUM is None:
                        SUM = self.getImgID(i)
                    else:
                        SUM += self.getImgID(i)
        if raw:
            return SUM
        if SUM is None:
            return None
        return pySPM.SPM_image(np.flipud(SUM))

    def getImgElt(self, elt):
        r = {}
        A = getSpecElt(SplitElts(elt))
        for k in A:
            I = self.getImgMass(k)
            if I is not None:
                r[k] = {
                    'data': I,
                    'mass': A[k][0],
                    'abund': A[k][1]}
        return r

    def showImgElt(self, elt, size=10, abundCorr=False, tot=True):
        A = self.getImgElt(elt)
        ks = [z for z in A if A[z]['data'] != None]
        fig, ax = plt.subplots((len(ks)+1)//4+1, 4,
                               figsize=(10*((len(ks)+1)//4+1), 10))
        mi = np.min(A[ks[0]]['data'].pixels)
        ma = np.max(A[ks[0]]['data'].pixels)
        A[ks[0]]['CT'] = ma
        for i, k in enumerate(ks[1:]):
            M = np.min(A[k]['data'].pixels)
            mi = min(mi, M)
            A[k]['minCT'] = M
            M = np.max(A[k]['data'].pixels)
            ma = max(ma, M)
            A[k]['CT'] = M
        ks.sort()
        di = 0
        if tot:
            SUM = sum([A[k]['data'].pixels for k in ks])
            ax[0][0].imshow(SUM, vmin=0)
            ax[0][0].set_title("Total")
            di = 1
        for i, k in enumerate(ks):
            if abundCorr:
                ax[(i+di)//4][(i+di) % 4].imshow(A[k]
                                                 ['data'].pixels/A[k]['abund'], vmin=0, vmax=ma)
            else:
                A[k]['data'].show(ax=ax[(i+di)//4][(i+di) %
                                                   4], vmin=0, vmax=ma)
            ax[(i+di)//4][(i+di) % 4].set_title(
                'mass: {mass:.3} - abundancy: {abund:.3f} - CT: {CT}'.format(**A[k]))

    def __exit__(self, exc_type, exc_value, traceback):
        self.f.close()


class BIF3D:

    def __init__(self, path):
        self.Peaks = {}
        self.RPeaks = {}
        self.Path = os.path.dirname(path)
        self.Basename = os.path.basename(path)
        self.size = None
        for x in os.listdir(self.Path):
            if x[:len(self.Basename)] == self.Basename:
                if x[-6:] == ".BIF3D":
                    s = self.Basename+r' \(([0-9]+)\)(?: - (.*?))?\.BIF3D'
                    r = re.search(s, x)
                    ID = int(r.group(1))
                    self.Peaks[ID] = r.group(2)
                    if r.group(2) != None:
                        self.RPeaks[r.group(2)] = ID
        self.data = {}

    def getIDs(self, channels):
        assert type(channels) in [str, int, list, tuple]
        if type(channels) == int:
            return channels
        if type(channels) == str:
            assert channels in self.RPeaks.keys()
            return self.RPeaks[channels]
        if type(channels) == list or type(channels) == tuple:
            ID = []
            for x in channels:
                ID.append(self.getIDs(x))
            return ID

    def getChannel(self, channel):
        k = self.getIDs(channel)
        if k in self.data:
            return self.data[k]
        v = self.Peaks[k]
        if v != None:
            path = '{path}{sep}{basename} ({ID}) - {Peak}.BIF3D'
        else:
            path = '{path}{sep}{basename} ({ID}).BIF3D'
        f = open(path.format(path=self.Path, sep=os.sep,
                             basename=self.Basename, ID=k, Peak=v), 'rb')
        A = f.read()
        f.close()
        size = struct.unpack("II", A[32:40])
        if self.size == None:
            self.size = size
        else:
            assert self.size == size
        return np.array(struct.unpack("{0}d".format(self.size[0]*self.size[1]), A[640:])).reshape(self.size)

    def loadChannel(self, channel):
        k = self.getIDs(channel)
        if not k in self.data:
            self.data[k] = self.getChannel(channel)

    def loadChannels(self, channels):
        if not type(channels) in [list, tuple]:
            self.loadChannel(channels)
        else:
            for x in channels:
                self.loadChannel(x)

    def getChannels(self, *channels):
        k = self.getIDs(channels)
        self.loadChannels(k)
        if type(k) == int:
            return self.data[k]
        D = np.zeros(self.size)
        for x in k:
            D += self.data[x]
        return D

    def listChannels(self):
        return self.Peaks

    def show(self, ax, *channels):
        D = self.getChannels(*channels)
        ax.imshow(D)
        if len(channels) > 1:
            ax.set_title(','.join([str(z) for z in channels]))
        else:
            ax.set_title(str(channels[0]))
