from pyTOF import Block
import numpy as np
import struct
import os.path
import zlib
import re

class ITM:
    def __init__(self, filename):
        """
        Create the ITM object out of the filename.
        Note that this works for all .ITA,.ITM, .ITS files as they have the same structure
        """
        self.filename = filename
        assert os.path.exists(filename)
        self.f = open(self.filename, 'rb')
        self.Type = self.f.read(8)
        assert self.Type == b'ITStrF01'
        self.root = Block.Block(self.f)

    def getSize(self):
        """
        Return a dict of the size of the image
        """
        d=self.root.goto('LateralShiftCorrection').dictList()
        return {
            'pixels':{
                'X':d[b'ImageStack.Raster.Resolution.X']['long'],
                'Y':d[b'ImageStack.Raster.Resolution.Y']['long']},
            'real':{
                'X':d[b'ImageStack.FieldOfView.X']['float'],
                'Y':d[b'ImageStack.FieldOfView.Y']['float']},
            'Scans':d[b'ImageStack.NumberOfShiftCoordinates']['long']}

    def getIntensity(self):
            """
            Retieve the total Ion image
            """
            S = self.getSize()
            X,Y = S['pixels']['X'],S['pixels']['Y']
            return np.array(struct.unpack('<'+str(X*Y)+'I',zlib.decompress(self.root.goto('Meta/SI Image/intensdata').value))).reshape((Y,X))

    def getValues(self, pb=False):
        """
        Beta function: Retieve a list of the values
        """
        Vals={}
        List=self.root.goto('propend').getList()
        if pb:
            import tqdm
            List=tqdm.tqdm(List)
        for l in List:
            Node = self.root.goto('propend').gotoItem(l['name'],l['id'])
            r = Node.getKeyValue(16)
            del Node
            S=Vals
            K=r['Key'].split('.')
            for k in K:
                if k not in S: S[k]={}
                S=S[k]
            S['value']=r['SVal']
        return Vals
        Data = self.root.goto('rawdata')
        List=Data.getList()
        if pb:
            List=tqdm.tqdm(List)
        for l in List:
            if l['name']==b'  20':
                Node = Data.gotoItem(l['name'],l['id'])
                r = Node.getKeyValue()
                del Node
                S=Vals
                K=r['Key'].split('.')
                for k in K:
                    if k not in S: S[k]={}
                    S=S[k]
                S['value']=r['SVal']
        return Vals

    def showValues(self, pb=False):
        from pyTOF import GUI
        GUI.ShowValues(self.getValues(pb))
        
    def getSpectrum(self):
        """
        Retieve a mass,spectrum array
        """
        RAW = zlib.decompress(self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/IITFSpecArray/CorrectedData').value)
        D = np.array(struct.unpack("<{0}f".format(len(RAW)//4), RAW))
        ch = np.arange(len(D))
        V = self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/IMassScaleSFK0')
        sf = V.goto('sf').getDouble()
        k0 = V.goto('k0').getDouble()
        chW = V.goto('channelwidth').getDouble()*1e-6
        masses = ((ch-k0/2)/(sf/2))**2
        return masses,D

    def showSpectrum(self, low=0, high=None,ax=None, log = False):
        """
        Plot the (summed) spectrum
        low and high: mass boundary of the plotted data
        ax: matplotlib axis
        log: plot the log of the intensity if True
        """
        m, s = self.getSpectrum()
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        if high is None:
            high = m[-1]
        mask = np.logical_and(m >= low, m <= high)
        M = m[mask]
        S = s[mask]
        if log:
            S = np.log(S)
        ax.plot(M, S)
        self.getMassInt()
        
        for P in self.peaks:
            p = self.peaks[P]
            c = p[b'cmass']['float']
            mask = (m>=p[b'lmass']['float'])*(m<=p[b'umass']['float'])
            if c >= low and c <= high:
                i = np.argmin(abs(m-c))
                ax.axvline(m[i], color='red')
                ax.fill_between(m[mask], *ax.get_ylim(), color='red', alpha=.2)
                ax.annotate(p[b'assign']['utf16'], (m[i],ax.get_ylim()[1]), (2, -10), va='top', textcoords='offset points')
            
    def getMassInt(self):
        """
        retrieve the peak list as a dictionnary
        """
        R = [z for z in self.root.goto('MassIntervalList').getList() if z['name'].decode() == 'mi']
        N = len(R)
        self.peaks = {}
        for x in R:
            try:
                X = self.root.goto('MassIntervalList/mi['+str(x['id'])+']')
                d = X.dictList()
                self.peaks[d[b'id']['long']] = d
            except ValueError:
                pass
                
    def showMassInt(self):
        """
        Display the peak list (assignment name with lower, center and upper mass)
        """
        self.getMassInt()
        for P in self.peaks:
            p = self.peaks[P]
            print(p[b'id']['long'], p[b'desc']['utf16'], p[b'assign']['utf16'], p[b'lmass']['float'], p[b'cmass']['float'], p[b'umass']['float'])
            
    def showStage(self, ax = None, markers=False):
        """
        Display an image of the stage used
        ax: maplotlib axis to be ploted in. If None the current one (or new) will be used
        markers: If True will display on the map the Position List items.
        """
        W = self.root.goto('SampleHolderInfo/bitmap/res_x').getLong()
        H = self.root.goto('SampleHolderInfo/bitmap/res_y').getLong()
        
        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 1, figsize=(W*10/H,10))

        Dat = zlib.decompress(self.root.goto('SampleHolderInfo/bitmap/imagedata').value)
        I = 255*np.array(struct.unpack("<"+str(W*H*3)+"B", Dat)).reshape((H, W, 3))
        ax.imshow(I);
        if markers:
            X = self.root.goto('Meta/SI Image[0]/stageposition_x').getDouble()
            Y = self.root.goto('Meta/SI Image[0]/stageposition_y').getDouble()

            def toXY(xy, W, H):
                sx = 23
                sy = 23
                return (913+sx*xy[0], 1145+sy*xy[1])
    
            for x in self.root.goto('SampleHolderInfo/positionlist'):
                if x.name == b'shpos':
                    y = pickle.loads(x.goto('pickle').value)
                    pos = toXY((y['stage_x'],y['stage_y']), W, H)
                    if pos[0] >= 0 and pos[0] < W and pos[1] >= 0 and pos[1] < H:
                            ax.annotate(y['name'], xy=pos, xytext=(-15, -25), textcoords='offset points', arrowprops=dict(arrowstyle='->', facecolor='black'));
                pos = toXY((X, Y), W, H)
                ax.plot(pos[0], pos[1], 'xr');
        ax.set_xlim((0, W))
        ax.set_ylim((0, H));
                
    def showPeaks(self):
        self.getMassInt()
        for p in self.peaks:
            P = {k.decode('utf8'):self.peaks[p][k] for k in self.peaks[p]}
            print("{0}) {peaklabel}".format(p,**P))