# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

from __future__ import absolute_import
from . import Block
import numpy as np
import struct
import os.path
import zlib
import re
import os
from .utils.misc import deprecated, aliased, alias

class InvalidRAWdataformat(Exception):
    def __init__(self, block, msg):
        self.block = block
        self.msg = msg
        
    def __str__(self):
        return "Invalid RAW dataformat seen in block "+self.block.path+self.block.name+' : '+self.msg

@aliased
class ITM:
    def __init__(self, filename, debug=False, readonly=False):
        """
        Create the ITM object out of the filename.  Note that this works for
        all .ITA,.ITM, .ITS files as they have the same structure
        
        The ITM has the specialty to contain the "rawdata" block. Which
        contains a lot of sub-blocks all having a name which consists of spaces
        followed by numbers.
        Blocks: '  20': Contains parameters of the ToF-SIMS recorder during the
        measurement such as the Emission Current and the Suppressor Voltages.
        All the parameters saves have the sames names that the ones found in
        propend and propstart. 
        '   2': Is the start block. Contains no info
        '   3': It's the end block. Contains a byte (unknown meaning)
        '   6': Uint32 indicating the scan number
        '   7': uint64 telling the pixel id of the next
        '  14' block. The id is the scanNumber*(width*height)+height*row+col
        '  14': binary block compressed with zlib. The raw data are encoded by
        a suite of uint32. The sequence is as follow. The first 4-btyes should
        end with \xC0, the second with \xD0, the third to \x40, the a suite
        of uint32 which don't end with \xc0 and the sequence starts again. Each
        sequence tell: The first uint32 terminating by \xC0 tells the x
        coordinates of the pixel location (just replace the \xC0 byte by 0 and
        read it as a uint32) The second uint32 terminating by\cD0 tells the y
        coord.  of the pixel. The third terminating by \x40 tells the pixel
        number (increasing monotonically. I know iontof like to write 10 times
        the same information everywhere). The rest are the detected peaks
        measured in channel unit for that specific pixel. In order to get the
        mass see the channel2mass function
        """
        self.filename = filename
        if not os.path.exists(filename):
            print("ERROR: File \"{}\" not found".format(filename))
            raise FileNotFoundError
        if readonly:
            self.f = open(self.filename, 'rb')
        else:
            self.f = open(self.filename, 'r+b')
        self.Type = self.f.read(8)
        assert self.Type == b'ITStrF01'
        self.root = Block.Block(self.f)
        try:
            d = self.root.goto('Meta/SI Image').dictList()
            self.size = {
                'pixels': {
                    'x': d['res_x']['long'],
                    'y': d['res_y']['long']},
                'real': {
                    'x': d['fieldofview']['float'],
                    'y': d['fieldofview']['float']*d['res_y']['long']/d['res_x']['long'],
                    'unit': 'm'}}
        except:
            s = self.getValue('Registration.Raster.Resolution')['int']
            fov = self.getValue("Registration.Raster.FieldOfView")['float']
            self.size = {'pixels':dict(x=s,y=s),'real':dict(x=fov,y=fov,unit='m')}
        try:
            self.size['Scans'] = \
                self.root.goto(
                    'filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.NumberOfScans').getULong()
        except:
            pass
        self.polarity = self.getValue("Instrument.Analyzer_Polarity_Switch")['string']
        self.peaks = {}
        self.MeasData = {}
        self.rawlist = None
        try:
            self.Nscan = self.root.goto("filterdata/TofCorrection/ImageStack/Reduced Data/NumberOfScans").getLong()
        except:
            try:
                self.Nscan = self.root.goto('propend/Measurement.ScanNumber').getKeyValue()['int']
            except:
                self.Nscan = None
        self.spp = self.root.goto("propend/Registration.Raster.ShotsPerPixel").getKeyValue()['int']
        try:
            R = [z for z in self.root.goto('MassIntervalList').getList() if z[
                'name'] == 'mi']
            N = len(R)
            for x in R:
                try:
                    X = self.root.goto('MassIntervalList/mi['+str(x['id'])+']')
                    d = X.dictList()
                    self.peaks[d['id']['long']] = d
                except ValueError:
                    pass
        except Exception as e:
            if debug:
                raise e
    
    @alias("get_peak_list")
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
    
    @alias("show_peak_list")
    def showPeakList(self, name):
        for m in self.getPeakList(name):
            print("{id[long]}: ({desc[utf16]}) [{assign[utf16]}] {lmass[float]:.2f}u - {umass[float]:.2f}u (center: {cmass[float]:.2f}u)".format(**m))
    
    @alias("get_property_trend")
    def getPropertyTrend(self, name):
        """
        Sometimes values might be saved in ITA files under PropertyTrend.
        You can recover them here
        """
        for x in self.root.goto("PropertyTrends"):
            if (x.name=='PropertyTrend' and x.goto("Trend.Name").getString() == name) or x.name == name:
                N = x.goto('Trend.Data.NumberEntries').getLong()
                data = x.goto("Trend.Data").value
                dat = np.array([struct.unpack('<4d',data[32*i:32*(i+1)])[2:4] for i in range(N)])
                return dat
        return None
        
    def get_summary(self, numeric=False):
        """
        Retrieve a summary of the important data concerning the measurement
        """
        from .utils import time2hms, funit
        def Get(k, default='Unknown', numeric=numeric):
            try:
                v = self.getValue(k)
                if len(v['string'].split())==1:
                    return v['string']
                unit = v['string'].split()[-1]
                if unit[0] in 'afpnumkMGE':
                    unit = unit[1:]
                if not numeric:
                    if unit == 's' and v['float']>60:
                        return time2hms(v['float'])
                    r = funit(v['float'], unit)
                    if int(r['value'])==r['value']:
                        return "{value:.0f} {unit}".format(**r)
                    return "{value:.2f} {unit}".format(**r)
                
                return v['string']
            except:
                return default
            
        return {
            'pixels': self.size['pixels'],
            'fov': self.root.goto('Meta/SI Image[0]/fieldofview').getDouble(),
            'LMIG': {
                'Extractor': Get("Instrument.LMIG.Extractor"),
                'Lens_Source': Get("Instrument.LMIG.Lens_Source")},
            'ExtractionDelay': Get("Instrument.Analyzer.ExtractionDelay"),
            'SputterSpecies': Get('Instrument.SputterGun.Species','Off'),
            'SputterEnergy': Get('Instrument.SputterGun.Energy','Off'),
            'AnalysisTime': Get("Analysis.AcquisitionTime"),
            'SputterTime': Get("Analysis.SputterTime"),
            'Scans':  Get("Analysis.TotalScans",self.Nscan),
            'TotalTime':  Get( "Analysis.TotalTime"),
            'peaks': self.get_masses(),
            'polarity': Get("Instrument.Analyzer_Polarity_Switch"),
            'CycleTime': Get( "Measurement.CycleTime"),
            'UpperMass': Get( "Measurement.UpperMass"),
            'LMIGDropouts': Get( "Measurement.LMIGDropouts"),
            'ShotsPerPixel' : Get( "Registration.Raster.ShotsPerPixel"),
            'RasterMode': Get("Registration.Raster.Mode")
            }

    def show_summary(self, fig=None, plot=True, **kargs):
        from . import funit
        s = self.get_summary(numeric=False)
        print("""
        Analysis time: {AnalysisTime}
        Delayed extraction: {ExtractionDelay}
        LMIG's Lens source: {LMIG[Lens_Source]}
        Sputter species: {SputterSpecies} @ {SputterEnergy}
        Field of View: {Fov[value]:.2f} {Fov[unit]}
        Pixel size: {pixels[x]} × {pixels[y]}
        Polarity: {polarity}
        Pixel Size: {pxs[value]:.2f} {pxs[unit]}
        Shots per pixel: {ShotsPerPixel}
        Raster mode: {RasterMode}
        Number of scans: {Nscan}
        Cycle time: {CycleTime} (Upper mass: {UpperMass})""".format(Nscan=self.Nscan, Fov=funit(s['fov'],'m'), pxs=funit(s['fov']/s['pixels']['x'], 'm'), **s))
        print("Peaks:")
        for x in s['peaks']:
            if x['assign']:
                print("\t▸ {assign} ({lmass:.3f} - {umass:.3f}u)".format(**x))
            elif x['desc']:
                print("\t▸ {desc} ({lmass:.3f} - {umass:.3f}u)".format(**x))
            else:
                print("\t▸ {cmass:.2f}u ({lmass:.3f} - {umass:.3f}u)".format(**x))
        if plot:
            import matplotlib.pyplot as plt
            import matplotlib as mpl
            if fig is None:
                fig = plt.figure(figsize=kargs.get("figsize",(21,10)))
            
            SI = self.getIntensity()
            Snapshot = self.getSnapshot()
            EC = self.getPropertyTrend("Instrument.LMIG.Emission_Current")
            Supp = self.getPropertyTrend("Instrument.LMIG.Suppressor")
            Press = self.getPropertyTrend("Instrument.VCU.Pressure.Main")
            N = (SI is not None)+(Snapshot is not None)+1+(EC is not None or Supp is not None or Press is not None)
            gs = mpl.gridspec.GridSpec(2, N)
            
            index = 0
            if SI is not None:
                ax = plt.subplot(gs[0, index])
                desc = self.root.goto('Meta/SI Image/description').getString()
                SI.show(ax=ax, **{k:kargs[k] for k in kargs if k not in ['high', 'low']})
                index += 1
            if Snapshot is not None:
                ax = plt.subplot(gs[0, index])
                desc = self.root.goto('Meta/Video Snapshot/description').getString()
                ax.imshow(Snapshot)
                ax.set_title(desc)
                index += 1
            axStage = plt.subplot(gs[0, index+1])
            axSpectra = plt.subplot(gs[1,:])
            if EC is not None or Supp is not None or Press is not None:
                ax = plt.subplot(gs[0, index])
                from mpl_toolkits.axes_grid1 import make_axes_locatable
                divider = make_axes_locatable(ax)
                index2 = 0
                from .utils import  s2hms
                if Press is not None:
                   t, tunit = s2hms(Press[:,0])
                   ax.plot(t, Press[:,1]*1e6, 'C2')
                   ax.set_xlabel("Time [{}]".format(tunit))
                   ax.set_ylabel("Pressure ($\cdot 10^{-8}$) [mbar]")
                   index2 += 1
                   if index2%2 == 0:
                    ax.yaxis.set_label_position("right")
                if EC is not None:
                   axc = divider.append_axes("bottom", size=1.2, sharex=ax)
                   t, tunit = s2hms(EC[:,0])
                   axc.plot(t, EC[:,1]*1e6, 'C0')
                   axc.set_xlabel("Time [{}]".format(tunit))
                   axc.set_ylabel("Emission Current [$\mu$A]")
                   index2 += 1
                   if index2%2 == 0:
                    axc.yaxis.set_label_position("right")

                if Supp is not None:
                   axb = divider.append_axes("bottom", size=1.2, sharex=ax)
                   t, tunit = s2hms(Supp[:,0])
                   axb.plot(t, Supp[:,1], 'C1')
                   axb.set_xlabel("Time [{}]".format(tunit))
                   axb.set_ylabel("LMIG Suppressor [V]")
                   index2 += 1
                   if index2%2 == 0:
                    axb.yaxis.set_label_position("right")
                index += 1
            
            self.showStage(ax=axStage, markers=True)
            self.showSpectrum(low=kargs.get('low',0), high=kargs.get('high', None), ax=axSpectra)
            
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
        from .SPM import SPM_image
        return SPM_image(I, real=self.size['real'], _type="TOF", zscale=zscale, channel=channel)     
       
    @alias("get_intensity")
    def getIntensity(self):
        """
        Retieve the total Ion image
        """
        try:
            X, Y = self.size['pixels']['x'], self.size['pixels']['y']
            img = self.image(np.flipud(np.array(self.root.goto('Meta/SI Image/intensdata').getData("f"), dtype=np.float32).reshape((Y, X))), channel="SI count")
        except Exception as e:
            try:
                img = self.getAddedImage(0).pixels
            except:
                try:
                    img = self.getAddedImageBySN(self.get_channel_SN("total"))
                except:
                    import warnings
                    warn.warnings("SI image cannot be retrieved")
                    return None
        return img

    def get_LMIG_info(self):
        rs = self.getValues(start=True)
        re = self.getValues()
        Val = [["Parameter name", "Value at start", "Value at end"]]
        for i, x in enumerate(rs):
            Val.append([x, rs[x], re[x]])

    @alias("get_value")
    def getValue(self, name, end=True):
        return self.root.goto('prop{}/{name}'.format(['start','end'][end],name=name)).getKeyValue()
    
    @alias("get_values")
    def getValues(self, pb=False, start=False,end=True,names=[], startsWith="", nest=False, hidePrefix=True, numeric=False):
        """
        Beta function: Retrieve a list of the values
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
                if numeric:
                    value = r['float']
                else:
                    value = r['string']
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
                        S['value @'+['end', 'start'][start]] = value
                    else:
                        if key_name in Vals:
                            Vals[key_name].append(value)
                        else:
                            Vals[key_name] = [value]
        return Vals

    @alias("auto_mass_cal")
    def autoMassCal(self, t=None, S=None, pos=True, debug=False, error=False, Range=5000, fitting_peaks = ['C','CH','CH2','CH3','Na'], sf=None, k0=None, apply=False, **kargs):
        """
        perform an auto calibration for spectrum. (in test, might be unreliable)
        """
        if 'FittingPeaks' in kargs:
            fitting_peaks = kargs['FittingPeaks']
        from .utils import get_mass, time2mass, fitSpectrum, mass2time
        from scipy.optimize import curve_fit
        TimeWidth = 1e10*self.root.goto('propend/Instrument.LMIG.Chopper.Width').getKeyValue()['float']
        if t is None or S is None:
            t, S = self.getSpectrum(time=True)
        N = np.prod(list(self.size['pixels'].values())) * self.Nscan
        mask = S > N*0.01/TimeWidth
        times_start = t[1:][np.nonzero(mask[1:]*(~mask[:-1]))]
        times_end = t[np.nonzero(mask[:-1]*(~mask[1:]))]
        times = (times_start + times_end)/2
        tH = times[0]
        mask = (t>=tH-TimeWidth)*(t<=tH+TimeWidth)
        tH = t[mask][np.argmax(S[mask])]
        if sf is None or k0 is None:
            sf=1e5
            for i in range(3):
                k0 = tH-sf*np.sqrt(get_mass('H'))
                m = time2mass(t, sf=sf, k0=k0)
                mP = time2mass(times, sf=sf, k0=k0)
                t0 = times[np.argmin(np.abs(mP-12))]
                t1 = times[np.argmin(np.abs(mP-12))+1]
                sf = np.sqrt((t1-k0)**2 - (t0-k0)**2)
                k0 = tH-sf*np.sqrt(get_mass('H'))
        ts = []
        for x in [mass2time(get_mass(x), sf=sf, k0=k0) for x in fitting_peaks]:
            mask = (t>=(x-Range))*(t<=(x+Range))
            t_peak = t[mask][np.argmax(S[mask])]
            ts.append(t_peak)
        ms = [get_mass(x) for x in fitting_peaks]
        sf, k0, dsf, dk0 = fitSpectrum(ts, ms, error=True)
        if apply:
            self.setSF(sf)
            self.setK0(k0)
        if debug:
            return sf, k0, dsf, dk0, ts, ms
        if error:
            return sf, k0, dsf, dk0
        return sf, k0
    
    @alias("show_values")
    def showValues(self, pb=False, gui=False, **kargs):
        from .utils import html_table, aa_table
        html = True
        if 'html' in kargs:
            html = kargs['html']
            del kargs['html']
        if gui:
            from pySPM.tools import values_display
            Vals = self.getValues(pb, nest=True, **kargs)
            values_display.ShowValues(Vals)
        else:
            Vals = self.getValues(pb, **kargs)
            Table = [["Parameter Name", "Value @start", "Value @end"]]
            for x in Vals:
                Table.append(tuple([x]+Vals[x]))
            if not html:
                print(aa_table(Table, header=True))
            else:
                from IPython.core.display import display, HTML
                res = html_table(Table, header=True)
                display(HTML(res))
                
    def get_mass_cal(self, alt=False):
        try:
            if alt:
                V = self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/IMassScaleSFK0')
                sf = V.goto('sf',lazy=True).getDouble()
                k0 = V.goto('k0',lazy=True).getDouble()
            else:
                sf = self.root.goto('MassScale/sf').getDouble()
                k0 = self.root.goto('MassScale/k0').getDouble()
        except:
            import warnings
            warnings.warn("Failed to get sf,k0, find alternative")
            if alt:
                sf = self.root.goto('MassScale/sf').getDouble()
                k0 = self.root.goto('MassScale/k0').getDouble()
            else:
                V = self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/IMassScaleSFK0')
                sf = V.goto('sf',lazy=True).getDouble()
                k0 = V.goto('k0',lazy=True).getDouble()
        return sf, k0

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
        
        For information the channel width is 50ps and is strangely not saved in the file.
        """
        if sf is None or k0 is None:
            sf0, k00 = self.get_mass_cal()
        if sf is None:
            sf = sf0
        if k0 is None:
            k0 = k00
        return ((binning*channels-k0)/(sf))**2

    @alias("get_spectrum")
    def getSpectrum(self, sf=None, k0=None, time=False, error=False, **kargs):
        """
        Retieve a mass,spectrum array
        This only works for .ita and .its files.
        For this reason it is implemented in the itm class.
        """
        RAW = zlib.decompress(self.root.goto(
            'filterdata/TofCorrection/Spectrum/Reduced Data/IITFSpecArray/'+['CorrectedData','Data'][kargs.get('uncorrected',False)]).value)
        D = np.array(struct.unpack("<{0}f".format(len(RAW)//4), RAW))
        ch = 2*np.arange(len(D))
        if time:
            return ch, D
        m = self.channel2mass(ch, sf=sf, k0=k0)
        if error:
            Dm = 2*np.sqrt(m)*np.sqrt(Dk0**2+m*Dsf**2)/sf
            return m, D, Dm
        return m, D
    
    @alias("get_meas_data")
    def getMeasData(self, name='Instrument.LMIG.Emission_Current', prog=False, debug=False):
        """
        Allows to recover the data saved during the measurements.
        This function is like getValues, but instead of giving the values at the beginning and at the end, 
        it track the changes of them during the measurement.
        """
        if name in self.MeasData:
            return self.MeasData[name]
        self.rawdata = self.root.goto('rawdata')
        L = self.rawdata.getList()
        i = 1
        while L[-i]['name'] !=  '  20':
            i += 1
        max_index = L[-i]['id']
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm as tqdm
            if prog is True:
                T = tqdm(L)
            else:
                T = prog(L)
        else:
            T=L
        for i,elt in enumerate(T):
            if elt['name'] != '  20':
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
    
    def show_stability(self, ax=None, prog=False):
        from .utils.plot import DualPlot
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        self.showMeasData(ax=ax, scans=3, mul=1e6, prog=prog);
        axb = DualPlot(ax)
        self.showMeasData("Instrument.LMIG.Suppressor", ax=axb, color='orange', scans=False, prog=prog);
        
    @alias("show_meas_data")
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
                ax.fill_between([s[i], s[i+1]],*lim,color='green',alpha=.1)
            axs.set_xlabel("Scan number")
            axs.set_xlim(ax.get_xlim())
            axs.set_ylim(*lim)
        return p
        
    @deprecated("SpectraPerPixel")
    def spectra_per_pixel(self, pixel_aggregation=None, peak_lim=0, scans=None, prog=False, safe=True, FOVcorr=True, smooth=False):
        """
        This function return a 2D array representing the spectra per pixel. The first axis correspond to each aggregated pixel and the second axis the spectral time.
        In order to keep the 2D array small enough the spectra are filtered in order to keep only strictly positive values (or larger than peak_lim).
        
        Parameters
        ----------
        pixel_aggregation: int
            the number of pixels aggregation. pixel_aggregation 
        
        """
        if pixel_aggregation is None:
            pixel_aggregation = max(1,int(self.size['pixels']['x']//64))
                 
        from .utils import get_mass, constants as const
        gun = self.root.goto('propend/Instrument.PrimaryGun.Species').getKeyValue()['string'] # Primary Gun Species (Bi1,Bi3,Bi3++)
        
        # if the + is missing in the name, add it
        if gun[-1]!='+':
            gun += '+'
            
        Q = gun.count('+') # number of charge
        
        nrj = self.root.goto('propend/Instrument.PrimaryGun.Energy').getKeyValue()['float'] # Primary ion energy (in eV)
        dx = self.size['real']['x']/self.size['pixels']['x'] # distance per pixel
        
        # Calculate the mass of the primary ion
        mp = get_mass(gun)
        
        if FOVcorr:
            DT = dx*(1/5e-11)*.5*np.sqrt(2)*np.sqrt((1e-3*mp/const.NA)/(Q*2*nrj*const.qe)) # delta time in channel per pixel. The 5e-11 is the channelwidth (50ps)
            # sqrt(2)/2 is from the sin(45°), nrj=E=.5*mp*v^2
        else:
            DT = 0
            
        if peak_lim is None:
            peak_lim = 0
            
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm
            PB = tqdm(total=3, postfix={'task':"Calculating total spectrum"})
            IT = lambda x: tqdm(x, leave=False)
        else:
            IT = lambda x: x
            
        pixel_size = int(self.size['pixels']['x']*self.size['pixels']['y']//pixel_aggregation**2)
        channels = round(self.getValue("Measurement.CycleTime")['float']/self.getValue("Registration.TimeResolution")['float'])
        
        # calculate total spectra
        m = np.zeros(channels)
        
        if scans is None:
            scans = range(self.Nscan)
        for scan in IT(scans):
            r = self.getRawData(scan)
            for p in r[2]:
                x, y = p
                for t in r[2][p]:
                    dt = DT*(self.size['pixels']['x']/2-x) # time correction for the given x coordinate (in channel number)
                    ip = int(dt)
                    fp = dt%1
                    m[t-ip] += (1-fp)
                    m[t-ip-1] += fp
                    
        # calculate the extreme cases
        max_time = np.nonzero(m)[0][-1]
        max_count = np.max(m)
        
        if prog:
            PB.update(1)
            PB.set_postfix({'task':"calculating aggregated spectrum"})

        # Select the peaks which are higher that peak_lim counts
        t = np.arange(channels)
        tx = t[m>peak_lim] # reduced time vector
        mx = m[m>peak_lim] # reduced mass vector
        rev = {x:-1 for x in range(max_time)}
        for k in range(len(tx)):
            rev[tx[k]] = k
        
        if safe:
            import psutil
            free_ram = psutil.virtual_memory().free
            if pixel_size*tx.size*8>= free_ram:
                raise Exception("""You don't have sufficient free RAM to perform this operation.
                Free RAM: {ram:.1f}Mb
                Number of pixels: {Npix}
                Spectrum size [value>{peak_lim}] : {ss}
                It is advised that you clean up memory or use a higher pixel_aggregation value.
                You can force the execution of this command by using the argument safe=False.
                """.format(ram=free_ram/1024**2, peak_lim=peak_lim, Npix=pixel_size, ss=pixel_size*tx.size*4/1024**2))
                
        size = (pixel_size, tx.size)
        spec = np.zeros(size, dtype='float32')
        for scan in IT(scans):
            r = self.getRawData(scan)
            for p in r[2]:
                x, y = p
                i = (self.size['pixels']['x']//pixel_aggregation)*(y//pixel_aggregation)+x//pixel_aggregation

                for t in r[2][p]:
                    dt = DT*(self.size['pixels']['x']/2-x) # time correction for the given x coordinate (in channel number)
                    ip = int(dt)
                    fp = dt%1
                    j1 = rev.get(t-ip, 0)
                    j2 = rev.get(t-ip-1, 0)
                    if j1>0:
                        spec[i][j1] += 1-fp
                    if j2>0:
                        spec[i][j2] += fp
                        
        if prog:
            PB.update(1)
            PB.set_postfix({'task':'smooth spectra'})

        if smooth:
            if smooth is True:
                smooth = (51, 3)
            from scipy.signal import savgol_filter
            for i in IT(range(pixel_size)):
                s = np.zeros(channels)
                s[m>peak_lim] = spec[i]
                s = savgol_filter(s, *smooth)
                spec[i] = s[m>peak_lim]
            
        if prog:
            PB.update(1)
            PB.set_postfix({'task':'done'})
            
        return m>peak_lim, spec
            
    @alias("show_spectrum")
    def showSpectrum(self, low=0, high=None, sf=None, k0=None, ax=None, log=False, showPeaks=True, **kargs):
        """
        Plot the (summed) spectrum
        low and high: mass boundary of the plotted data
        ax: matplotlib axis
        log: plot the log of the intensity if True
        """
        m, s = self.getSpectrum(sf=sf,k0=k0,**kargs)
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        if high is None:
            high = m[-1]
        mask = np.logical_and(m >= low, m <= high)
        M = m[mask]
        S = s[mask]
        if log:
            ax.log = True
            S[S>=1] = np.log10(S[S>=1])
            S[S<1] = 0
        ax.plot(M, S, **kargs)
        self.get_masses()
        if showPeaks:
            ymax = ax.get_ylim()[1]
            labels = []
            pos = []
            index = 0
            for P in [x for x in self.peaks if self.peaks[x]['desc']['utf16'] not in ['total','sum of rest']]:
                p = self.peaks[P]
                c = p['cmass']['float']
                mask = (m >= p['lmass']['float'])*(m <= p['umass']['float'])
                if c >= low and c <= high:
                    i = np.argmin(abs(m-c))
                    pos.append(m[i])
                    #ax.fill_between(m[mask], 0, ymax, color=['red','green','blue'][index%3], alpha=.2)
                    index += 1
                    labels.append(p['assign']['utf16'])
            from .utils import put_Xlabels
            put_Xlabels(ax, pos, labels);
            ax.set_xlabel("Mass [u]")
            ax.set_ylabel("Total counts [-]");
            
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
                'id': p['id']['long'],
                'desc': p['desc']['utf16'],
                'assign': p['assign']['utf16'],
                'lmass': p['lmass']['float'],
                'cmass': p['cmass']['float'],
                'umass': p['umass']['float']})
        return result

    @alias("get_raw_spectrum")
    def getRawSpectrum(self, scans=None, ROI=None, FOVcorr=True, deadTimeCorr=True, **kargs):
        """
        Reconstruct the spectrum from RAW data.
        scans: List of scans to use. if None all scans are used (default)
        ROI: Region Of Interest. It can be:
            1) None: All pixels are used
            2) A list of images. The pixels will be added to the spectrum i if the value of the i-th image is True for that pixel
            3) An image with integer value. The current pixels will be added to the i-th spectrum if the value of ROI at that pixel is i.
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
        from .utils import get_mass, constants as const
        if ROI is None and 'roi' in kargs:
            ROI = kargs['roi']
            
        gun = self.root.goto('propend/Instrument.PrimaryGun.Species').getKeyValue()['string'] # Primary Gun Species (Bi1,Bi3,Bi3++)
        
        # if the + is missing in the name, add it
        if gun[-1]!='+':
            gun += '+'
            
        Q = gun.count('+') # number of charge
        
        nrj = self.root.goto('propend/Instrument.PrimaryGun.Energy').getKeyValue()['float'] # Primary ion energy (in eV)
        dx = self.size['real']['x']/self.size['pixels']['x'] # distance per pixel
        
        # Calculate the mass of the primary ion
        mp = get_mass(gun)
            
        # Perform the time of flight correction?
        if FOVcorr:
            DT = dx*(1/5e-11)*.5*np.sqrt(2)*np.sqrt((1e-3*mp/const.NA)/(Q*2*nrj*const.qe)) # delta time in channel per pixel. The 5e-11 is the channelwidth (50ps)
            # sqrt(2)/2 is from the sin(45°), nrj=E=.5*mp*v^2
        else:
            DT = 0
            
        if scans is None:
            scans = range(self.Nscan)
            
        assert hasattr(scans, '__iter__')
        
        # Allocate vector for the spectrum
        number_channels = int(round(self.root.goto('propend/Measurement.CycleTime').getKeyValue()['float']\
            / self.root.goto('propend/Registration.TimeResolution').getKeyValue()['float']))
        
        if kargs.get('prog',False):
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm as tqdm
            T = tqdm(scans)
        else:
            T = scans
            
        if ROI is None:
            Spectrum = np.zeros(number_channels, dtype=np.float32)
            for s in T:
                Data = self.getRawData(s)[2]
                if kargs.get('prog', False):
                    LData = tqdm(Data, leave=False)
                else:
                    LData = Data
                for xy in LData:
                    dt = DT*(self.size['pixels']['x']/2-xy[0]) # time correction for the given x coordinate (in channel number)
                    ip = int(dt)
                    fp = dt%1
                    for x in Data[xy]:
                        Spectrum[x-ip] += (1-fp)
                        Spectrum[x-ip-1] += fp
        elif type(ROI) is np.ndarray:
            assert np.min(ROI)>=0
            Spectrum = np.zeros((number_channels,np.max(ROI)+1), dtype=np.float32)
            for s in T:
                Data = self.getRawData(s)[2]
                if kargs.get('prog', False):
                    LData = tqdm(Data, leave=False)
                else:
                    LData = Data
                for xy in LData:
                    dt = DT*(self.size['pixels']['x']/2-xy[0]) # time correction for the given x coordinate (in channel number)
                    ip = int(dt)
                    fp = dt%1
                    id = ROI[xy[1], xy[0]]
                    for x in Data[xy]:
                        Spectrum[x-ip, id] += (1-fp)
                        Spectrum[x-ip-1, id] += fp
        elif type(ROI) in [list, tuple]:
            multi_roi = True
            Spectrum = np.zeros((number_channels, len(ROI)), dtype=np.float32)
            for s in T:
                Data = self.getRawData(s)[2]
                for xy in Data:
                    dt = DT*(self.size['pixels']['x']/2-xy[0]) # time correction for the given x coordinate (in channel number)
                    ip = int(dt)
                    fp = dt%1
                    li = []
                    for k,R in enumerate(ROI):
                        if R[xy[1], xy[0]]:
                            li.append(k)
                    for x in Data[xy]:
                        for k in li:
                            Spectrum[x-ip, k] += (1-fp)
                            Spectrum[x-ip-1, k] += fp
        sf = kargs.get('sf', self.root.goto('MassScale/sf').getDouble())
        k0 = kargs.get('k0', self.root.goto('MassScale/k0').getDouble())
        masses = self.channel2mass(np.arange(number_channels), sf=sf, k0=k0)
        if deadTimeCorr:
            dt = 1300 # 65ns*(1ch/50ps) = 1300 channels
            N = self.Nscan*self.size['pixels']['x']*self.size['pixels']['y'] # total of count events
            Np = np.zeros(Spectrum.shape)
            if Spectrum.ndim>1:
                for i in range(Spectrum.shape[1]):
                    Np[:,i] = N-np.convolve(Spectrum[:, i], np.ones(dt-1, dtype=int), 'full')[:-dt+2]
            else:
                Np = N-np.convolve(Spectrum, np.ones(dt-1,dtype=int), 'full')[:-dt+2]
            Np[Np==0] = 1
            Spectrum = -N*np.log(1-Spectrum/Np)
        if kargs.get('time', False):
            return np.arange(number_channels), Spectrum
        return masses, Spectrum

    @alias("get_raw_raw_data")
    def getRawRawData(self, scan=0):
        assert scan < self.Nscan
        found = False
        RAW = b''
        if self.rawlist is None:
            self.rawlist = self.root.goto('rawdata').getList()
        startFound = False
        for x in self.rawlist:
            if x['name'] == '   6':
                if not startFound:
                    startFound = x['id'] == scan
                else:
                    break
            elif startFound and x['name'] == '  14':
                self.root.f.seek(x['bidx'])
                child = Block.Block(self.root.f)
                RAW += zlib.decompress(child.value)
        if type(RAW) is str:
            return bytearray(RAW)
        return RAW

    @alias("get_raw_data")
    def getRawData(self, scan=0):
        """
        Function which allows you to read and parse the raw data.
        With this you are able to reconstruct the data.
        Somehow the number of channel is double in the raw data compared
        to the compressed version saved in the ITA files.
        scan: The scan number. Start at 0
        """
        RAW = self.getRawRawData(scan)
        Blocks = {}
        _Block = []
        PixelList = []
        PixelOrder = np.zeros((self.size['pixels']['y'], self.size['pixels']['x']))
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
                b = bytearray(RAW[i:i+4])
                if b[3] < 64:
                    raise TypeError("Expecting a 40 or higher block at {}, got {:02x}".format(i+3,b[3]))
                b = b[:3] + bytearray([b[3]&16])
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

    @alias("show_stage")
    def showStage(self, ax=None, markers=False, plist=False):
        """
        Display an image of the stage used
        ax: maplotlib axis to be ploted in. If None the current one (or new) will be used
        markers: If True will display on the map the Position List items.
        """
        import pickle
        import matplotlib as mpl
        W = self.root.goto('SampleHolderInfo/bitmap/res_x').getLong()
        H = self.root.goto('SampleHolderInfo/bitmap/res_y').getLong()

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 1, figsize=(W*10/H, 10))

        Dat = zlib.decompress(self.root.goto(
            'SampleHolderInfo/bitmap/imagedata').value)
        I = np.array(struct.unpack("<"+str(W*H*3)+"B", Dat), dtype=np.uint8).reshape((H, W, 3))
        ax.imshow(I)
        if markers:
            X = self.root.goto('Meta/SI Image[0]/stageposition_x').getDouble()
            Y = self.root.goto('Meta/SI Image[0]/stageposition_y').getDouble()

            def toXY(xy, W, H):
                sx = 23
                sy = 23
                return (913+sx*xy[0], 1145+sy*xy[1])
            if plist:
                for x in self.root.goto('SampleHolderInfo/positionlist'):
                    if x.name == 'shpos':
                        y = pickle.loads(x.goto('pickle').value)
                        pos = toXY((y['stage_x'], y['stage_y']), W, H)
                        if pos[0] >= 0 and pos[0] < W and pos[1] >= 0 and pos[1] < H:
                            ax.annotate(y['name'], xy=pos, xytext=(-15, -25), textcoords='offset points',
                                        arrowprops=dict(arrowstyle='->', facecolor='black'))
            pos = toXY((X, Y), W, H)
            ax.plot(pos[0], pos[1], 'xr')
            ll = toXY((X-self.fov*5e2, Y-self.fov*5e3), W, H)
            ur = toXY((X+self.fov*5e2, Y+self.fov*5e3), W, H)
            ax.add_patch(mpl.patches.Rectangle(ll, ur[0]-ll[0], ur[1]-ll[1], ec='lime', fill=False))
        ax.set_xlim((0, W))
        ax.set_ylim((0, H))
        ax.set_xticks([])
        ax.set_yticks([])

    @alias("get_snapshot")
    def getSnapshot(self):
        """
        Return the video snapshot
        """
        try:
            dl = self.root.goto('Meta/Video Snapshot').dictList()
            sx = dl['res_x']['long']
            sy = dl['res_y']['long']
            img = np.array(self.root.goto('Meta/Video Snapshot/imagedata').getData('B')).reshape((sy, sx, 3))
            return  img
        except Exception as e:
            return None
    @alias("show_peaks")
    def showPeaks(self):
        for p in self.peaks:
            P = self.peaks[p]
            label = P['assign']['utf16']
            if label:
                print("{0}) {peaklabel}".format(p, peaklabel=label))
            else:
                print("{0}) {cmass:.2f}u [{lmass:.2f}u-{umass:.2f}u]".format(p, cmass=P['cmass']['float'],lmass=P['lmass']['float'],umass=P['umass']['float']))
            
    def modify_block_and_export(self, path, new_data, output, reload=True, **kargs):
        self.root.modify_block_and_export(path, new_data, output, **kargs)
        if reload:
            self.f.close()
            self.f.open(path, "rb+")
            self.f.read(8)
            self.root = Block.Block(self.f)

    def reconstruct(self, channels, scans=None, sf=None, k0=None, prog=False, time=False):
        """
        Reconstruct an Image from a raw spectra by defining the lower and upper mass
        channels: list of (lower_mass, upper_mass)
        upper_mass: upper mass of the peak
        scans: The list of the scans to take into account (if None, all scans are taken)
        sf/k0: mass calibration. If none take the saved values
        time: If true the upper/lower_mass will be understood as time value
        prog: If True display a progressbar with tqdm
        """
        from .utils import mass2time
        from . import SPM_image
        assert hasattr(channels, '__iter__')
        if not hasattr(channels[0], '__iter__'):
            channels = [channels]
        for c in channels:
            assert len(c)==2
        if scans is None:
            scans = range(self.Nscan)
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm as tqdm
            scans = tqdm(scans)
        left = np.array([x[0] for x in channels])
        right = np.array([x[1] for x in channels])
        if not time:
            if sf is None or k0 is None:
                sf, k0 = self.get_mass_cal()
            left = mass2time(left, sf=sf, k0=k0)
            right = mass2time(right, sf=sf, k0=k0)
        Counts = [np.zeros((self.size['pixels']['x'],self.size['pixels']['y'])) for x in channels]
        for s in scans:
            PixelList, PixelOrder, Data = self.getRawData(s)
            for xy in Data:
                for i,ch in enumerate(channels):
                    Counts[i][xy[1],xy[0]] += np.sum((Data[xy]>=left[i])*(Data[xy]<=right[i]))
        res = [SPM_image(C, real=self.size['real'], _type='TOF', channel="{0[0]:.2f}{unit}-{0[1]:.2f}{unit}".format(channels[i],unit=["u", "s"][time]), zscale="Counts") for i,C in enumerate(Counts)]
        if len(res) == 1:
            return res[0]
        return res
        
    def create_new_miblock(self, assign, lmass=None, umass=None, cmass=None, desc="", _uuid=None, **kargs):
        import uuid
        if _uuid is None:
            _uuid = str(uuid.uuid4())
        if _uuid[0] != '{':
            _uuid = '{'+_uuid+'}'
        _uuid = _uuid.upper()
        if lmass is None:
            lmass = self.channel2mass(0)
        if umass is None:
            up = int(np.round(self.getValue("Measurement.CycleTime")['float']/self.getValue("Registration.TimeResolution")['float'], 0))
            umass = self.channel2mass(up)
        new_index = max([x.head['ID'] for x in self.root.goto("MassIntervalList") if x.name == 'mi'])+1
        new_index2 = self.root.goto("Measurement Options/massintervals/NextId").getLong()
        new_id = max([x.goto("id").getLong() for x in self.root.goto("MassIntervalList") if x.name == 'mi'])+1
        new_id2 = max([x.goto("id").getLong() for x in self.root.goto("Measurement Options/massintervals") if x.name == 'mi'])+1
        NID = self.root.goto("MassIntervalList/NextId", lazy=True)
        NID.rewrite(struct.pack("<I", NID.getLong()+1))
        defaults = {
            'clsid': ('raw', b'\x85\x1b\xa9\xe4hZ\xc8M\x93\xe0\x7f\xc0\xf8\xe2\xbb\xd9'),
            'color': ('raw', b'\xff\x00\x00\x00'),
            'symbolID': ('I',8),
            'RSF': ('d', 0),
            'label': ('B', 1),
            'peaklabel': ('I', 0),
            'edrfl': ('I', 0),
            'attr.ConsiderBurstMode': ('B', 1),
            'attr.Normalized':('B', 0),
            'attr.PoissonCorrectable':('B', 1),
            'attr.Sticky':('B', 0),
            'attr.Visible':('B', 1)
            }
        p = dict(new_index=new_index,new_index2=new_index2)
        self.root.edit_block("Measurement Options/massintervals/mi[{new_index2}]".format(**p), "id", struct.pack("<I", new_id))
        self.root.edit_block("MassIntervalList/mi[{new_index}]".format(**p), "id", struct.pack("<I", new_id))
        for path in ["Measurement Options/massintervals/mi[{new_index2}]","MassIntervalList/mi[{new_index}]"]:
            self.root.edit_block(path.format(**p), "desc", desc.encode('utf16')[2:])
            self.root.edit_block(path.format(**p), "SN", _uuid.encode('utf16')[2:])
            self.root.edit_block(path.format(**p), "assign", assign.encode('utf16')[2:])
            self.root.edit_block(path.format(**p), "lmass", struct.pack("<d", lmass))
            self.root.edit_block(path.format(**p), "umass", struct.pack("<d", umass))
            if cmass is None:
                cmass = (lmass+umass)/2
            self.root.edit_block(path.format(**p), "cmass", struct.pack("<d", cmass))
            self.root.edit_block(path.format(**p), "desc", desc.encode('utf16')[2:])
            for key in defaults:
                data = kargs.get(key, defaults[key][1])
                fmt = defaults[key][0]
                if fmt == 'utf16':
                    data = data.encode('utf16')[2:]
                elif fmt=='raw':
                    pass
                else:
                    data = struct.pack("<"+fmt, data)
                self.root.edit_block(path.format(**p), key, data)
        blk = self.root.goto("MassIntervalList/mi[{new_index}]".format(**p))
        d = blk.dictList()
        self.peaks[d['id']['long']] = d
        return blk
      
    def __del__(self):
        if self.f.writable():
            self.f.flush()
            os.fsync(self.f)
        self.f.close()

    @alias("set_k0")
    def setK0(self, k0):
        import struct
        b = self.root.goto("MassScale/k0")
        buffer = struct.pack("<d", k0)
        b.rewrite(buffer);

    @alias("set_sf")
    def setSF(self, sf):
        import struct
        b = self.root.goto("MassScale/sf")
        buffer = struct.pack("<d", sf)
        b.rewrite(buffer);
