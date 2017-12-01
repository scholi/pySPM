from pySPM import Block, utils
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
        
        The ITM has the speciality to contain the "rawdata" block. Which contains a lot of sub-blocks all having a name which concists of spaces followed by numbers.
            Blocks: '  20': Contains parameters of the ToF-SIMS recorder during the measurement such as the Emission Current and the Suppressor Voltages.
                    All the parameters saves have the sames names that the ones found in propend and propstart.
                    '   2': Is the start block. Contains no info
                    '   3': It's the end block. Contains a byte (unknown meaning)
                    '   6': Uint32 indicating the scan number
                    '   7': uint64 telling the pixel id of the next '  14' block. The id is the scanNumber*(width*height)+height*row+col
                    '  14': binary block compressed with zlib. The raw data are encoded by a suite of uint32.
                            The sequence is as follow. The first 4-btyes should end with \xC0, the the second with \xD0, the third to \x40, the a suite of uint32 which don't end with \xc0
                            and the sequence starts again. Each sequence tell:
                                The first uint32 terminating by \xC0 tells the x coordinates of the pixel location (just replace the \xC0 byte by 0 and read it as a uint32)
                                The second uint32 terminating by\cD0 telld the y coord. of the pixel
                                The third terminating by \x40 tells the pixel number (increasing monotonically. I know iontof like to write 10 times the same information everywhere)
                                The rest are the detected peaks measured in channel unit for that specific pixel. In order to get the mass see the channel2mass function
        """
        self.filename = filename
        assert os.path.exists(filename)
        self.f = open(self.filename, 'rb')
        self.Type = self.f.read(8)
        assert self.Type == b'ITStrF01'
        self.root = Block.Block(self.f)
        d = self.root.goto('Meta/SI Image').dictList()
        self.size = {
            'pixels': {
                'x': d[b'res_x']['long'],
                'y': d[b'res_y']['long']},
            'real': {
                'x': d[b'fieldofview']['float'],
                'y': d[b'fieldofview']['float']*d[b'res_y']['long']/d[b'res_x']['long'],
                'unit': 'm'}}
        try:
            self.size['Scans'] = \
                self.root.goto(
                    'filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.NumberOfScans').getULong()
        except:
            pass
        self.peaks = {}
        self.MeasData = {}
        try:
            self.Nscan = self.root.goto('propend/Measurement.ScanNumber').getKeyValue()['int']
        except:
            self.Nscan = None
            
        try:
            R = [z for z in self.root.goto('MassIntervalList').getList() if z[
                'name'].decode() == 'mi']
            N = len(R)
            for x in R:
                try:
                    X = self.root.goto('MassIntervalList/mi['+str(x['id'])+']')
                    d = X.dictList()
                    self.peaks[d[b'id']['long']] = d
                except ValueError:
                    pass
        except:
            pass
    
    def getPeakList(self, name):
        """
        Retrieve extra MassIntervalList (MIL) by it's name. Those are saved and visible in iontof software in the spectrum window.
        """
        PeakList = []
        for x in self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/ExtraMILs'):
            if x.goto('Name').getString() == name:
                for y in [z for z in x.getList() if z['name'].decode() == 'mi']:
                    d = x.goto('mi[{}]'.format(y['id'])).dictList()
                    PeakList.append({key.decode('utf8'):d[key] for key in d})
                return PeakList
        return None
        
    def showPeakList(self, name):
        for m in self.getPeakList(name):
            print("{id[long]}: ({desc[utf16]}) [{assign[utf16]}] {lmass[float]:.2f}u - {umass[float]:.2f}u (center: {cmass[float]:.2f}u)".format(**m))
        
    def get_summary(self):
        """
        Retrieve a summary of the important data concerning the measurement
        """
        values = self.getValues(start=False,names=["Instrument.SputterGun.Energy","Instrument.SputterGun.Species",
            "Instrument.LMIG.Extractor","Instrument.LMIG.Lens_Source","Instrument.Analyzer.ExtractionDelay",
            "Analysis.AcquisitionTime","Analysis.SputterTime","Analysis.TotalScans","Analysis.TotalTime"])
        def Get(v,k):
            return v.get(k,['Unknown'])[0]
            
        return {
            'pixels': self.size['pixels'],
            'fov': self.root.goto('Meta/SI Image[0]/fieldofview').getDouble(),
            'LMIG': {
                'Extractor': Get(values,"Instrument.LMIG.Extractor"),
                'Lens_Source': Get(values,"Instrument.LMIG.Lens_Source")},
            'ExtractionDelay': Get(values,"Instrument.Analyzer.ExtractionDelay"),
            'SputterSpecies': values.get('Instrument.SputterGun.Species',['Off'])[0],
            'SputterEnergy': values.get('Instrument.SputterGun.Energy',['Off'])[0],
            'AnalysisTime': Get(values,"Analysis.AcquisitionTime"),
            'SputterTime': Get(values,"Analysis.SputterTime"),
            'Scans': values.get("Analysis.TotalScans",[self.Nscan])[0],
            'TotalTime':  Get(values, "Analysis.TotalTime"),
            'peaks': self.get_masses(),
            }
        
    def getIntensity(self):
        """
        Retieve the total Ion image
        """
        X, Y = self.size['pixels']['x'], self.size['pixels']['y']
        return np.array(struct.unpack('<'+str(X*Y)+'I',
                                      zlib.decompress(self.root.goto('Meta/SI Image/intensdata').value))).reshape((Y, X))

    def get_LMIG_info(self):
        rs = self.getValues(start=True)
        re = self.getValues()
        Val = [["Parameter name", "Value at start", "Value at end"]]
        for i, x in enumerate(rs):
            Val.append([x, rs[x], re[x]])

    def getValues(self, pb=False, start=True,end=True,names=[], startsWith="", nest=False, hidePrefix=True):
        """
        Beta function: Retieve a list of the values
        """
        Vals = {}
        startEnd = []
        if start:
            startEnd.append(True)
        if end:
            startEnd.append(False)
        for start in startEnd:
            List = self.root.goto(['propend', 'propstart'][start]).getList()
            if pb:
                import tqdm
                List = tqdm.tqdm(List)
            for l in List:
                Node = self.root.goto(['propend', 'propstart'][
                                      start]).gotoItem(l['name'], l['id'])
                r = Node.getKeyValue()
                del Node
                S = Vals
                if r['key'] in names or ( names==[] and r['key'].startswith(startsWith) ):
                    if hidePrefix:
                        key_name = r['key'][len(startsWith):]
                    else:
                        key_name = r['key']
                    if nest:
                        K = key_name.split('.')
                        for k in K:
                            if k not in S:
                                S[k] = {}
                            S = S[k]
                        S['value @'+['end', 'start'][start]] = r['string']
                    else:
                        if key_name in Vals:
                            Vals[key_name].append(r['string'])
                        else:
                            Vals[key_name] = [r['string']]
                            
        return Vals

    def autoMassCal(self, pos=True, debug=False, Range=5000, FittingPeaks = ['C','CH','CH2','CH3','Na'], sf=None, k0=None):
        """
        perform an auto callibration for positive spectrum. (in test, might be unreliable)
        """
        from scipy.optimize import curve_fit
        TimeWidth = 1e10*self.root.goto('propend/Instrument.LMIG.Chopper.Width').getKeyValue()['float']
        t, D = self.getSpectrum(time=True)
        N = np.prod(list(self.size['pixels'].values())) * self.Nscan
        mask = D > N*0.01/TimeWidth
        times_start = t[1:][np.nonzero(mask[1:]*(~mask[:-1]))]
        times_end = t[np.nonzero(mask[:-1]*(~mask[1:]))]
        times = (times_start + times_end)/2
        tH = times[0]
        mask = (t>=tH-TimeWidth)*(t<=tH+TimeWidth)
        tH = t[mask][np.argmax(D[mask])]
        if sf is None or k0 is None:
            sf=1e5
            for i in range(3):
                k0 = tH-sf*np.sqrt(utils.getMass('H+'))
                m = utils.time2mass(t, sf=sf, k0=k0)
                mP = utils.time2mass(times, sf=sf, k0=k0)
                t0 = times[np.argmin(np.abs(mP-12))]
                t1 = times[np.argmin(np.abs(mP-12))+1]
                sf = np.sqrt((t1-k0)**2 - (t0-k0)**2)
                k0 = tH-sf*np.sqrt(utils.getMass('H+'))
        ts = []
        for x in [utils.mass2time(utils.getMass(x+'+'), sf=sf, k0=k0) for x in FittingPeaks]:
            mask = (t>=(x-Range))*(t<=(x+Range))
            t_peak = t[mask][np.argmax(D[mask])]
            ts.append(t_peak)
        ms = [utils.getMass(x+'+') for x in FittingPeaks]
        if debug:
            return utils.fitSpectrum(ts, ms), ts, ms
        return utils.fitSpectrum(ts, ms)
    
    def showValues(self, pb=False, gui=False, **kargs):
        html = True
        if 'html' in kargs:
            html = kargs['html']
            del kargs['html']
        if gui:
            from pySPM import GUI
            Vals = self.getValues(pb, nest=True, **kargs)
            GUI.ShowValues(Vals)
        else:
            Vals = self.getValues(pb, **kargs)
            Table = [["Parameter Name", "Value @start", "Value @end"]]
            for x in Vals:
                Table.append(tuple([x]+Vals[x]))
            if not html:
                print(utils.aa_table(Table, header=True))
            else:
                from IPython.core.display import display, HTML
                res = utils.html_table(Table, header=True)
                display(HTML(res))

    def channel2mass(self, channels, sf=None, k0=None, binning=1):
        """
        Calculate the mass from the time of flight
        The parameters sf and k0 will be read from the file (calibration used by iontof) in case they are None
        The binning just multiply the channels in order to get the correct answer (used for the ITA files for the moment)
        default values for uncalibrated spectrum are sf = 100'000 &  k0 = 0
        The time can be retrieved by knowing the channel width (typically 100ps). So time = channel * channelWidth
        The mass = (2*q*U/L**2)*(channel*channelWidth)**2
        L is the path length, thus twice the columns length and is about 2m I think
        For an extractor of 2keV and L=2m, 1u = (time[s]*310620.843175[u**.5/s])**2 = [channel/sf]**2 => sf is about 32000 for k0 = 0
        Caution. If data are binned, the channel number should be multiplied by the binning factor in order to use the same sf and k0 factor!
        """
        if sf is None or k0 is None:
            V = self.root.goto(
                'filterdata/TofCorrection/Spectrum/Reduced Data/IMassScaleSFK0')
        if sf is None:
            for id in [x['id'] for x in V.getList() if x['name']==b'sf']:
                try:
                    sf = V.gotoItem('sf',id).getDouble()
                    break
                except:
                    pass
        if k0 is None:
            for id in [x['id'] for x in V.getList() if x['name']==b'k0']:
                try:
                    k0 = V.gotoItem('k0',id).getDouble()
                    break
                except:
                    pass
        return ((binning*channels-k0)/(sf))**2

    def getSpectrum(self, sf=None, k0=None, time=False):
        """
        Retieve a mass,spectrum array
        """
        RAW = zlib.decompress(self.root.goto(
            'filterdata/TofCorrection/Spectrum/Reduced Data/IITFSpecArray/CorrectedData').value)
        D = np.array(struct.unpack("<{0}f".format(len(RAW)//4), RAW))
        ch = 2*np.arange(len(D))
        if time:
            return ch, D
        return self.channel2mass(ch, sf=sf, k0=k0), D
    
    def getMeasData(self, name='Instrument.LMIG.Emission_Current', prog=False, debug=False):
        """
        Allows to recover the data saved during the measurements.
        This function is like getValues, but instead of giving the values at the beginning and at the end, 
        it track the changes of them during the measurement.
        """
        if name in self.MeasData:
            return self.MeasData[name]
        rawdata = self.root.goto('rawdata')
        L = rawdata.getList()
        i = 1
        while L[-i]['name'] !=  b'  20':
            i += 1
        max_index = L[-i]['id']
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm as tqdm
            T = tqdm(L)
        else:
            T=L
        for i,elt in enumerate(T):
            if elt['name'] != b'  20':
                continue
            idx = elt['bidx']
            self.f.seek(idx)
            child = Block.Block(self.f)
            r = child.getKeyValue(0)
            if not r['key'] in self.MeasData:
                self.MeasData[r['key']] = []
            self.MeasData[r['key']].append((idx,r['float']))
        if name in self.MeasData:
            return self.MeasData[name]
        else:
            raise KeyError(name)
    
    def showMeasData(self, name='Instrument.LMIG.Emission_Current', prog=False, ax=None, mul=1, scans=2, **kargs):
        t = self.getMeasData('Measurement.AcquisitionTime')
        S = self.getMeasData("Measurement.ScanNumber")
        idx = [x[0] for x in t]
        time = [x[1] for x in t]
        ScanData = [x[1] for x in S]
        ScanIdx = [x[0] for x in S]
        s = np.interp(ScanIdx,idx,time)
        
        Meas = self.getMeasData(name,prog=prog)
        
        MeasData = [x[1] for x in Meas]
        MeasIdx = [x[0] for x in Meas]
        t = np.interp(MeasIdx,idx, time)
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        p = ax.plot(t, np.array(MeasData)*mul, **kargs)
        ax.set_xlabel("Time [s]");
        ax.set_ylabel(name)
        if scans:
            assert type(scans) is int
            lim = ax.get_ylim()
            axs = ax.twiny()
            axs.set_xticks(.5*s[:-1:scans]+.5*s[1::scans])
            axs.set_xticklabels([str(i+1) for i in range(0,self.Nscan,scans)])
            colors = [i%2 for i in range(0,self.Nscan,scans)]
            for i,tick in enumerate(axs.xaxis.get_ticklabels()):
                tick.set_color(["black","green"][colors[i]])
            axs.set_xlim(ax.get_xlim())
            for i in range(1,self.Nscan-1,2):
                ax.fill_between([s[i],s[i+1]],*lim,color='green',alpha=.1)
            axs.set_xlabel("Scan number")
            axs.set_xlim(ax.get_xlim())
            axs.set_ylim(*lim)
        return p

    def showSpectrum(self, low=0, high=None, sf=None, k0=None, ax=None, log=False, showPeaks=True):
        """
        Plot the (summed) spectrum
        low and high: mass boundary of the plotted data
        ax: matplotlib axis
        log: plot the log of the intensity if True
        """
        m, s = self.getSpectrum(sf=sf,k0=k0)
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
        self.get_masses()
        if showPeaks:
            for P in [x for x in self.peaks if self.peaks[x][b'desc']['utf16'] not in ['total','sum of rest']]:
                p = self.peaks[P]
                c = p[b'cmass']['float']
                mask = (m >= p[b'lmass']['float'])*(m <= p[b'umass']['float'])
                if c >= low and c <= high:
                    i = np.argmin(abs(m-c))
                    ax.axvline(m[i], color='red')
                    ax.fill_between(m[mask], *ax.get_ylim(), color='red', alpha=.2)
                    ax.annotate(p[b'assign']['utf16'], (m[i], ax.get_ylim()[
                                1]), (2, -10), va='top', textcoords='offset points')

    def get_masses(self, mass_list=None):
        """
        retrieve the peak list as a dictionnary
        """
        if mass_list is None:
            masses = self.peaks
        else:
            masses = mass_list
        result = []
        for P in masses:
            if mass_list is None:
                p = self.peaks[P]
            else:
                p = P
            result.append({
                'id': p[b'id']['long'],
                'desc': p[b'desc']['utf16'],
                'assign': p[b'assign']['utf16'],
                'lmass': p[b'lmass']['float'],
                'cmass': p[b'cmass']['float'],
                'umass': p[b'umass']['float']})
        return result

    def getRawSpectrum(self, scans=None, ROI=None, FOVcorr=True, deadTimeCorr=True, **kargs):
        """
        Reconstruct the spectrum from RAW data.
        scans: List of scans to use. if None all scans are used (default)
        ROI: Region Of Interest. It's an image of the same size as the measurement having a value of True for pixels to be taken into accoun.
        FOVcorr: Correction for the primary time of flight variation
        
        Δt = (√2/2)∙x∙√(mp/(2∙E)) where x is the x-corrdinate (in m), mp, the primary ion mass (in kg) and E the primary energy
        as E=½∙mp∙v² ⇒ v = √(2E/mp)
        
        Note: if E is given in eV, it should be multiplied by the electron charge.
        
        deadTimeCorr: Dead time correction (Poisson statistics)
        
        Icorr(c) = -N*log(1-I(c)/N'(c))
        where N = Nscan*size['x']*size['y']
        N'(c) = N - sum_{c'=c-ct}^{c-1}I(c')
        
        see Ref. T. Stephan, J. Zehnpfenning and A. Benninghoven, J. vac. Sci. A 12 (2), 1994
 
        """

        gun = self.root.goto('propend/Instrument.PrimaryGun.Species').getKeyValue()['string'] # Primary Gun Species (Bi1,Bi3,Bi3++)
        nrj = self.root.goto('propend/Instrument.PrimaryGun.Energy').getKeyValue()['float'] # Primary ion energy (in eV)
        dx = self.size['real']['x']/self.size['pixels']['x'] # distance per pixel
        
        # Calculate the mass of the primary ion
        if '+' in gun:
            mp = utils.getMass(gun)
        else:
            mp = utils.getMass(gun+'+')
            
        # Perform the time of flight correction?
        if FOVcorr:
            DT = dx*(1/5e-11)*.5*np.sqrt(2)*np.sqrt((1e-3*mp/utils.NA)/(2*nrj*utils.qe)) # delta time in channel per pixel. The 1e10 is the channelwidth (100ps)
            # sqrt(2)/2 is from the sin(45°), nrj=E=.5*mp*v^2
        else:
            DT = 0
            
        if scans is None:
            scans = range(self.Nscan)
            
        assert hasattr(scans, '__iter__')
        
        # Allocate vector for the spectrum
        number_channels = int(round(self.root.goto('propend/Measurement.CycleTime').getKeyValue()['float']\
            / self.root.goto('propend/Registration.TimeResolution').getKeyValue()['float']))
        multi_roi = False
        if ROI is None or not(type(ROI)==list or type(ROI)==tuple):
            Spectrum = np.zeros(number_channels, dtype=np.float32)
        else:
            multi_roi = True
            Spectrum = np.zeros((number_channels,len(ROI)), dtype=np.float32)
        # Display a progress bar?
        if kargs.get('prog',False):
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm as tqdm
            T = tqdm(scans)
        else:
            T = scans
                
        for s in T:
            Data = self.getRawData(s)[2]
            for xy in Data:
                dt = DT*(self.size['pixels']['x']/2-xy[0]) # time correction for the given x coordinate (in channel number)
                ip = int(dt)
                fp = dt%1
                if multi_roi:
                    li = []
                    for k,R in enumerate(ROI):
                        if R[xy[1],xy[0]]:
                            li.append(k)
                    for x in Data[xy]:
                        for k in li:
                            Spectrum[x-ip,k] += (1-fp)
                            Spectrum[x-ip-1,k] += fp
                else:
                    if ROI is None or ROI[xy[1],xy[0]]:
                        for x in Data[xy]:
                            Spectrum[x-ip] += (1-fp)
                            Spectrum[x-ip-1] += fp
        sf = kargs.get('sf',self.root.goto('MassScale/sf').getDouble())
        k0 = kargs.get('k0',self.root.goto('MassScale/k0').getDouble())
        masses = self.channel2mass(np.arange(number_channels),sf=sf,k0=k0)
        if deadTimeCorr:
            dt = 2000 # 100ns*(1ch/50ps) = 2000 channels
            N = self.Nscan*self.size['pixels']['x']*self.size['pixels']['y'] # total of count events
            Np = np.zeros(Spectrum.shape)
            if multi_roi:
                for i in range(len(ROI)):
                    Np[:,i] = N-np.convolve(Spectrum[:,i], np.ones(dt-1,dtype=int), 'full')[:-dt+2]
            else:
                Np = N-np.convolve(Spectrum, np.ones(dt-1,dtype=int), 'full')[:-dt+2]
            Np[Np==0] = 1
            Spectrum = -N*np.log(1-Spectrum/Np)
        return masses, Spectrum
        
    def getRawData(self, scan=0):
        """
        Function which allows you to read and parse the raw data.
        With this you are able to reconstruct the data.
        Somehow the number of channel is double in the raw data compared
        to the compressed version saved in the ITA files.
        scan: The scan number. Start at 0
        """
        assert scan < self.Nscan
        found = False
        RAW = b''
        list_ = self.root.goto('rawdata').getList()
        startFound = False
        for x in list_:
            if x['name'] == b'   6':
                if not startFound:
                    startFound = x['id'] == scan
                else:
                    break
            elif startFound and x['name'] == b'  14':
                self.root.f.seek(x['bidx'])
                child = Block.Block(self.root.f)
                RAW += zlib.decompress(child.value)
        Blocks = {}
        _Block = []
        PixelList = []
        PixelOrder = np.zeros(
            (self.size['pixels']['y'], self.size['pixels']['x']))
        i = 0
        while i < len(RAW):
            b = RAW[i:i+4]
            if b[3:4] == b'\xc0':
                if len(_Block):
                    Blocks[(x, y)] = _Block
                b = b[:3] + b'\x00'
                x = struct.unpack('<I', b)[0]
                i += 4
                b = RAW[i:i+4]
                if b[3:4] != b'\xd0':
                    raise TypeError("Expecting a D0 block at {}".format(i+3))
                b = b[:3] + b'\x00'
                y = struct.unpack('<I', b)[0]
                i += 4
                b = RAW[i:i+4]
                if b[3] < 64:
                    raise TypeError("Expecting a 40 or higher block at {}, got {:02x}".format(i+3,b[3]))
                b = b[:3] + bytes([b[3]&16])
                _Block = []
                id = struct.unpack('<I', b)[0]
                PixelOrder[y, x] = id
                PixelList.append((x,y))
            else:
                _Block.append(struct.unpack('<I', b)[0])
            i += 4
        return PixelList, PixelOrder, Blocks

    def show_masses(self, mass_list=None):
        """
        Display the peak list (assignment name with lower, center and upper mass)
        """
        for m in self.get_masses():
            if mass_list is None or m['id'] in [z['id'] for z in mass_list]:
                print(
                    "{id}: ({desc}) [{assign}] {lmass:.2f}u - {umass:.2f}u (center: {cmass:.2f}u)".format(**m))

    def showStage(self, ax=None, markers=False):
        """
        Display an image of the stage used
        ax: maplotlib axis to be ploted in. If None the current one (or new) will be used
        markers: If True will display on the map the Position List items.
        """
        W = self.root.goto('SampleHolderInfo/bitmap/res_x').getLong()
        H = self.root.goto('SampleHolderInfo/bitmap/res_y').getLong()

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 1, figsize=(W*10/H, 10))

        Dat = zlib.decompress(self.root.goto(
            'SampleHolderInfo/bitmap/imagedata').value)
        I = 255*np.array(struct.unpack("<"+str(W*H*3)+"B", Dat)
                         ).reshape((H, W, 3))
        ax.imshow(I)
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
                    pos = toXY((y['stage_x'], y['stage_y']), W, H)
                    if pos[0] >= 0 and pos[0] < W and pos[1] >= 0 and pos[1] < H:
                        ax.annotate(y['name'], xy=pos, xytext=(-15, -25), textcoords='offset points',
                                    arrowprops=dict(arrowstyle='->', facecolor='black'))
                pos = toXY((X, Y), W, H)
                ax.plot(pos[0], pos[1], 'xr')
        ax.set_xlim((0, W))
        ax.set_ylim((0, H))

    def showPeaks(self):
        self.getMassInt()
        for p in self.peaks:
            P = {k.decode('utf8'): self.peaks[p][k] for k in self.peaks[p]}
            print("{0}) {peaklabel}".format(p, **P))
            
    def modify_block_and_export(self, path, new_data, output, **kargs):
        self.root.modify_block_and_export(path, new_data, output, **kargs)
