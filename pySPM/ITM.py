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
from .utils.misc import deprecated, aliased, alias, PB
from warnings import warn

class InvalidRAWdataformat(Exception):
    def __init__(self, block, msg):
        self.block = block
        self.msg = msg
        
    def __str__(self):
        return "Invalid RAW dataformat seen in block "+self.block.path+self.block.name+' : '+self.msg

@aliased
class ITM:
    def __init__(self, filename, debug=False, readonly=True, precond=False, label=None):
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
        
        Parameters
        ----------
        filename : string
            Path of the ITA/ITM/ITS file
        debug : bool
            if True display some debug message.
        readonly : bool
            The pySPM library can now EDIT ITM/ITA files. In case you want to avoid that, please set readonly=True.
            This might be also useful in case the file is open by another program which locks the file.
        precond : bool
            If True will run the preconditioner (adjust k0 so that H peak is correct and adjust scaling factor on the H peak)
        """
        self.filename = filename
        if label is None:
            self.label = os.path.basename(filename)
        else:
            self.label = label
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
            d = self.root.goto('Meta/SI Image').dict_list()
            self.size = {
                'pixels': {
                    'x': d['res_x']['long'],
                    'y': d['res_y']['long']},
                'real': {
                    'x': d['fieldofview']['float'],
                    'y': d['fieldofview']['float']*d['res_y']['long']/d['res_x']['long'],
                    'unit': 'm'}}
        except:
            s = self.get_value('Registration.Raster.Resolution')['int']
            fov = self.get_value("Registration.Raster.FieldOfView")['float']
            self.size = {'pixels':dict(x=s,y=s),'real':dict(x=fov,y=fov,unit='m')}
        try:
            self.size['Scans'] = \
                self.root.goto(
                    'filterdata/TofCorrection/ImageStack/Reduced Data/ImageStackScans/Image.NumberOfScans').get_ulong()
        except:
            pass
        self.polarity = self.get_value("Instrument.Analyzer_Polarity_Switch")['string']
        self.peaks = {}
        self.meas_data = {}
        self.rawlist = None
        try:
            self.Nscan = self.root.goto("filterdata/TofCorrection/ImageStack/Reduced Data/NumberOfScans").getLong()
        except:
            try:
                self.Nscan = self.root.goto('propend/Measurement.ScanNumber').get_key_value()['int']
            except:
                self.Nscan = None
        self.spp = self.root.goto("propend/Registration.Raster.ShotsPerPixel").get_key_value()['int']
        try:
            R = [z for z in self.root.goto('MassIntervalList').get_list() if z['name'] == 'mi']
            N = len(R)
            for x in R:
                try:
                    X = self.root.goto('MassIntervalList/mi['+str(x['id'])+']')
                    d = X.dict_list()
                    self.peaks[d['id']['long']] = d
                except ValueError:
                    pass
        except Exception as e:
            if debug:
                raise e
        self.sf, self.k0 = self.get_mass_cal()
        self.scale = 1
        if precond:
            self.precond()
    
    @alias("getPeakList")
    def get_peak_list(self, name):
        """
        Retrieve extra MassIntervalList (MIL) by it's name. Those are saved and visible in iontof software in the spectrum window.
        """
        PeakList = []
        for x in self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/ExtraMILs'):
            if x.goto('Name').get_string() == name:
                for y in [z for z in x.get_list() if z['name'].decode() == 'mi']:
                    d = x.goto('mi[{}]'.format(y['id'])).dict_list()
                    PeakList.append({key.decode('utf8'):d[key] for key in d})
                return PeakList
        return None
    
    @alias("showPeakList")
    def show_peak_list(self, name):
        for m in self.get_peak_list(name):
            print("{id[long]}: ({desc[utf16]}) [{assign[utf16]}] {lmass[float]:.2f}u - {umass[float]:.2f}u (center: {cmass[float]:.2f}u)".format(**m))
    
    @alias("getPropertyTrend")
    def get_property_trend(self, name):
        """
        Sometimes values might be saved in ITA files under PropertyTrend.
        You can recover them here
        """
        for x in self.root.goto("PropertyTrends"):
            if (x.name=='PropertyTrend' and x.goto("Trend.Name").get_string() == name) or x.name == name:
                N = x.goto('Trend.Data.NumberEntries').get_long()
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
            'fov': self.root.goto('Meta/SI Image[0]/fieldofview').get_double(),
            'Floodgun': Get("Instrument.Timing.Floodgun"),
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
        Floodgun: {Floodgun}
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
            
            SI = self.get_intensity()
            Snapshot = self.get_snapshot()
            EC = self.get_property_trend("Instrument.LMIG.Emission_Current")
            Supp = self.get_property_trend("Instrument.LMIG.Suppressor")
            Press = self.get_property_trend("Instrument.VCU.Pressure.Main")
            N = (SI is not None)+(Snapshot is not None)+1+(EC is not None or Supp is not None or Press is not None)
            gs = mpl.gridspec.GridSpec(2, N)
            
            index = 0
            if SI is not None:
                ax = plt.subplot(gs[0, index])
                desc = self.root.goto('Meta/SI Image/description').get_string()
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
            
            self.show_stage(ax=axStage, markers=True)
            self.show_spectrum(low=kargs.get('low',0), high=kargs.get('high', None), ax=axSpectra)
            
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
       
    @alias("getIntensity")
    def get_intensity(self):
        """
        Retrieve the total Ion image
        """
        X, Y = self.size['pixels']['x'], self.size['pixels']['y']
        img = self.image(np.flipud(np.array(self.root.goto('Meta/SI Image/intensdata').get_data("f"), dtype=np.float32).reshape((Y, X))), channel="SI count")
        return img

    def get_LMIG_info(self):
        rs = self.get_values(start=True)
        re = self.get_values()
        Val = [["Parameter name", "Value at start", "Value at end"]]
        for i, x in enumerate(rs):
            Val.append([x, rs[x], re[x]])

    @alias("getValue")
    def get_value(self, name, end=True):
        return self.root.goto('prop{}/{name}'.format(['start','end'][end],name=name)).get_key_value()
    
    @alias("getValues")
    def get_values(self, prog=False, start=False, end=True, names=[], startsWith="", nest=False, hidePrefix=True, numeric=False):
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
            List = self.root.goto(['propend', 'propstart'][start]).get_list()
            if prog:
                List = PB(List)
            for l in List:
                Node = self.root.goto(['propend', 'propstart'][
                                      start]).goto_item(l['name'], l['id'])
                r = Node.get_key_value()
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

    @alias("autoMassCal")
    def auto_mass_cal(self, t=None, S=None, pos=True, debug=False, error=False, Range=5000, fitting_peaks = ['C','CH','CH2','CH3','Na'], apply=False, **kargs):
        """
        perform an auto calibration for spectrum. (in test, might be unreliable)
        """
        # old parameter name compatibility
        if 'FittingPeaks' in kargs:
            fitting_peaks = kargs['FittingPeaks']
        if type(fitting_peaks) is str:
            fitting_peaks = fitting_peaks.split(",")
        negative = self.polarity=='Negative'
        fitting_peaks = [x+[['+','-'][negative],'']['+' in x or '-' in x] for x in fitting_peaks]
        from .utils import get_mass, time2mass, fit_spectrum, mass2time
        from scipy.optimize import curve_fit
        time_width = 1e10*self.root.goto('propend/Instrument.LMIG.Chopper.Width').get_key_value()['float']
        if t is None or S is None:
            t, S = self.get_spectrum(time=True)
        N = np.prod(list(self.size['pixels'].values())) * self.Nscan
        mask = S > N*0.01/time_width
        times_start = t[1:][np.nonzero(mask[1:]*(~mask[:-1]))]
        times_end = t[np.nonzero(mask[:-1]*(~mask[1:]))]
        times = (times_start + times_end)/2
        tH = times[0]
        mask = (t>=tH-time_width)*(t<=tH+time_width)
        tH = t[mask][np.argmax(S[mask])]
        if not 'sf' in kargs:
            sf = self.sf
        else:
            sf = kargs.pop('sf')
        if not 'k0' in kargs:
            k0 = self.k0
        else:
            k0 = kargs.pop('k0')
            
        if sf is None or k0 is None:
            sf = 72000
            for i in range(3):
                k0 = tH-sf*np.sqrt(get_mass('H'+['+','-'][self.polarity=='Negative']))
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
        sf, k0, dsf, dk0 = fit_spectrum(ts, ms, error=True)
        if apply:
            self.set_sf(sf)
            self.set_k0(k0)
        self.sf = sf
        self.k0 = k0
        if debug:
            return sf, k0, dsf, dk0, ts, ms
        if error:
            return sf, k0, dsf, dk0
        return sf, k0
    
    @alias("showValues")
    def show_values(self, pb=False, gui=False, **kargs):
        from .utils import html_table, aa_table
        html = True
        if 'html' in kargs:
            html = kargs['html']
            del kargs['html']
        if gui:
            from pySPM.tools import values_display
            Vals = self.get_values(pb, nest=True, **kargs)
            values_display.show_values(Vals)
        else:
            Vals = self.get_values(pb, **kargs)
            Table = [["Parameter Name", "Value @start", "Value @end"]]
            for x in Vals:
                Table.append(tuple([x]+Vals[x]))
            if not html:
                print(aa_table(Table, header=True))
            else:
                from IPython.core.display import display, HTML
                res = html_table(Table, header=True)
                display(HTML(res))
                
    def reset_mass_cal(self, alt=False):
        self.sf, self.k0 = self.get_mass_cal(alt=alt)
        
    def get_mass_cal(self, alt=False):
        try:
            if not alt:
                V = self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/IMassScaleSFK0')
                sf = V.goto('sf',lazy=True).get_double()
                k0 = V.goto('k0',lazy=True).get_double()
            else:
                sf = self.root.goto('MassScale/sf').get_double()
                k0 = self.root.goto('MassScale/k0').get_double()
        except:
            import warnings
            warnings.warn("Failed to get sf,k0, find alternative")
            if not alt:
                sf = self.root.goto('MassScale/sf').get_double()
                k0 = self.root.goto('MassScale/k0').get_double()
            else:
                V = self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/IMassScaleSFK0')
                sf = V.goto('sf',lazy=True).get_double()
                k0 = V.goto('k0',lazy=True).get_double()
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
        
        For information the channel width is 50ps and can be retrieved by pySPM.ITM.get_value("Registration.TimeResolution")
        """
        if sf is None or k0 is None:
            sf0, k00 = self.sf, self.k0
        if sf is None:
            sf = sf0
        if k0 is None:
            k0 = k00
        return ((binning*channels-k0)/(sf))**2

    def shift_sf(self, dsf):
        """
        This function adjust the value of sf by adding the value given as parameter and recalculate k0 such that the H peak does not move.
        """
        from .utils import get_mass
        self.k0 = self.k0-dsf*np.sqrt(get_mass("H+"))
        self.sf += dsf
        
    def rescale(self, elt, amp=10000, delta=.02):
        from .utils.fit import peak_fit
        m, s = self.get_spectrum(scale=1)
        p0 = peak_fit(m, s, elt, delta=delta)
        self.scale = amp / p0[2]
        return p0
        
    def adjust_sf(self, elt, delta=0.02):
        from .utils import LG, get_mass
        from .utils.fit import peak_fit
        m, s = self.get_spectrum()
        mH = get_mass("H"+"+-"[self.polarity=='Negative'])
        p0 = peak_fit(m, s, mH, delta=.02)
        tH = 2*np.argmin(abs(m-p0[0]))
        
        mX = get_mass(elt)
        p0 = peak_fit(m, s, mX, delta=delta)
        tX = 2*np.argmin(abs(m-p0[0]))
        self.sf = (tX-tH)/(np.sqrt(mX)-np.sqrt(mH))
        self.k0 = tH-self.sf*np.sqrt(mH)
    
    def precond(self, amp=10000, do_mass_cal=False):
        """
        This is a pre-conditioner function which adjust k0 so that the H peak is exactly centered on the correct mass and it adjusts the scale factor so that the H peak is exactly equivalent to amp (10'000 by default).
        
        This function is useful in case people want to compare several measurements together.
        
        Parameters
        ----------
        amp : float, int
            The amplitude of the H peak after scaling
        apply : bool
            Appy the changes of k0 and sf to the file
        """
        from .utils import get_mass, LG
        from .utils.fit import peak_fit
        
        m, s  = self.getSpectrum()
        mask = (m>.5)*(m<1.5)
        amp0 = np.max(s[mask])
        mH = get_mass("H"+"+-"[self.polarity=='Negative'])
        if do_mass_cal or amp0 < max(100, .05*self.Nscan*self.size['pixels']['x']*self.size['pixels']['y']):
            warn("the initial mass calibration seems to be wrong, let's try to perform it from scratch")
            self.sf, self.k0 = self.auto_mass_cal(sf=72000, k0=0)
            m, s  = self.getSpectrum()
            mask = (m>.98)*(m<1.02)
            amp0 = np.max(s[mask])
        try:
            p0 = peak_fit(m, s, mH)
        except:
            if not do_mass_cal:
                return self.precond(amp=amp, apply=apply, do_mass_cal=True)
        tH = 2*np.argmin(abs(m-p0[0]))
        self.k0 = tH - self.sf*np.sqrt(mH)
        if amp is not None:
            self.scale = amp / p0[2]
    
    @alias("getSpectrum")
    def get_spectrum(self, sf=None, k0=None, scale=None, time=False, error=False, **kargs):
        """
        Retieve a mass,spectrum array
        This only works for .ita and .its files.
        For this reason it is implemented in the itm class.
        """
        RAW = zlib.decompress(self.root.goto(
            'filterdata/TofCorrection/Spectrum/Reduced Data/IITFSpecArray/'+['CorrectedData','Data'][kargs.get('uncorrected',False)]).value)
        if scale is None:
            scale = self.scale
        D = scale*np.array(struct.unpack("<{0}f".format(len(RAW)//4), RAW))
        ch = 2*np.arange(len(D)) # We multiply by two because the channels are binned.
        if time:
            return ch, D
        m = self.channel2mass(ch, sf=sf, k0=k0)
        if error:
            Dm = 2*np.sqrt(m)*np.sqrt(Dk0**2+m*Dsf**2)/sf
            return m, D, Dm
        return m, D
    
    @alias("getMeasData")
    def get_meas_data(self, name='Instrument.LMIG.Emission_Current', prog=False, debug=False):
        """
        Allows to recover the data saved during the measurements.
        This function is like getValues, but instead of giving the values at the beginning and at the end, 
        it track the changes of them during the measurement.
        """
        if name in self.meas_data:
            return self.meas_data[name]
        self.rawdata = self.root.goto('rawdata')
        L = self.rawdata.get_list()
        i = 1
        while L[-i]['name'] !=  '  20':
            i += 1
        max_index = L[-i]['id']
        if prog:
            T = PB(L)
        else:
            T = L
        for i, elt in enumerate(T):
            if elt['name'] != '  20':
                continue
            idx = elt['bidx']
            self.f.seek(idx)
            child = Block.Block(self.f)
            r = child.get_key_value(0)
            if not r['key'] in self.meas_data:
                self.meas_data[r['key']] = []
            self.meas_data[r['key']].append((idx, r['float']))
        if name in self.meas_data:
            return self.meas_data[name]
        else:
            raise KeyError(name)
    
    def show_stability(self, ax=None, prog=False):
        from .utils.plot import dual_plot
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        self.show_meas_data(ax=ax, scans=3, mul=1e6, prog=prog);
        axb = dual_plot(ax)
        self.show_meas_data("Instrument.LMIG.Suppressor", ax=axb, color='orange', scans=False, prog=prog);
        
    def show_meas_data(self, name='Instrument.LMIG.Emission_Current', prog=False, ax=None, mul=1, scans=2, **kargs):
        t = self.get_meas_data('Measurement.AcquisitionTime')
        S = self.get_meas_data("Measurement.ScanNumber")
        idx = [x[0] for x in t]
        time = [x[1] for x in t]
        ScanData = [x[1] for x in S]
        ScanIdx = [x[0] for x in S]
        s = np.interp(ScanIdx,idx,time)
        
        Meas = self.get_meas_data(name, prog=prog)
        
        MeasData = [x[1] for x in Meas]
        MeasIdx = [x[0] for x in Meas]
        t = np.interp(MeasIdx, idx, time)
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
            axs.set_xticklabels([str(i+1) for i in range(0, self.Nscan,scans)])
            colors = [i%2 for i in range(0, self.Nscan, scans)]
            for i, tick in enumerate(axs.xaxis.get_ticklabels()):
                tick.set_color(["black", "green"][colors[i]])
            axs.set_xlim(ax.get_xlim())
            for i in range(1, self.Nscan-1,2):
                ax.fill_between([s[i], s[i+1]], *lim, color='green', alpha=.1)
            axs.set_xlabel("Scan number")
            axs.set_xlim(ax.get_xlim())
            axs.set_ylim(*lim)
        return p
        
    @alias("showSpectrumAround")
    def show_spectrum_around(self, m0, delta=None, amp_scale=1, sf=None, k0=None, **kargs):
        """
        Display the Spectrum around a given mass.

        Parameters
        ----------
        m0 : float
            The central mass around which the spectrum will be plotted (in u)
        delta : float
            The spectrum will be plotted between m0-delta and m0+delta
        sf : float or None
        k0 : float or None
            sf and k0 are the mass calibration parameters. If None values saved with the file will be used.
        **kargs : supplementary arguments
            Passed to pySPM.utils.showPeak
        """
        polarity = '+'
        if self.get_value('Instrument.Analyzer_Polarity_Switch')['string'] == 'Negative':
            polarity = '-'
        from . import utils
        m, D = self.get_spectrum(sf=sf, k0=k0)
        if 'label' not in kargs:
            kargs['label'] = self.label
        return utils.show_peak(m, D*amp_scale, m0, delta, polarity=polarity, sf=sf, k0=k0, **kargs)
        
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
            pixel_aggregation = max(1, int(self.size['pixels']['x']//64))
                 
        from .utils import get_mass, constants as const, closest_arg
        gun = self.root.goto('propend/Instrument.PrimaryGun.Species').get_key_value()['string'] # Primary Gun Species (Bi1,Bi3,Bi3++)
        
        # if the + is missing in the name, add it
        if gun[-1]!='+':
            gun += '+'
            
        Q = gun.count('+') # number of charge
        
        nrj = self.root.goto('propend/Instrument.PrimaryGun.Energy').get_key_value()['float'] # Primary ion energy (in eV)
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
            pb = PB(total=3, postfix={'task':"Calculating total spectrum"})
            IT = lambda x: PB(x, leave=False)
        else:
            IT = lambda x: x
            
        pixel_size = ((self.size['pixels']['x']+pixel_aggregation-1)//pixel_aggregation)*((self.size['pixels']['y']+pixel_aggregation-1)//pixel_aggregation)
        channels = round(self.get_value("Measurement.CycleTime")['float']/self.get_value("Registration.TimeResolution")['float'])
        
        # calculate total spectra
        m = np.zeros(channels)
        
        if scans is None:
            scans = range(self.Nscan)
        dts = DT*(self.size['pixels']['x']/2-np.arange(self.size['pixels']['x'])) # time correction for the given x coordinate (in channel number)
        for scan in IT(scans):
            raw = self.get_raw_raw_data(scan)
            rawv = struct.unpack('<{}I'.format(len(raw)//4), raw)
            i = 0
            while i<len(rawv):
                b = rawv[i]
                if b& 0xc0000000:
                    x = b & 0x0fffffff
                    dt = dts[x]
                    #fp = dt%1
                    ip = int(dt)
                    i += 3
                else:
                    m[b-ip] += 1#(1-fp)
                    #m[b-ip-1] += fp
                    i += 1
                    
        # calculate the extreme cases
        max_time = np.nonzero(m)[0][-1]
        max_count = np.max(m)
        
        if prog:
            pb.update(1)
            pb.set_postfix({'task':"calculating aggregated spectrum"})

        # Select the peaks which are higher that peak_lim counts
        t = np.arange(channels)
        tx = t[m>peak_lim] # reduced time vector
        mx = m[m>peak_lim] # reduced mass vector
        rev = {x:-1 for x in range(max_time+1)}
        for k in range(len(tx)):
            rev[tx[k]] = k
        
        if safe:
            import psutil
            free_ram = psutil.virtual_memory().free
            if pixel_size*tx.size*4>= free_ram:
                raise Exception("""You don't have sufficient free RAM to perform this operation.
                Free RAM: {ram:.1f}Mb
                Number of pixels: {Npix}
                Spectrum size [value>{peak_lim}] : {tx} elements
                Array size: {N} elements = {ss}Mb
                It is advised that you clean up memory or use a higher pixel_aggregation value.
                You can force the execution of this command by using the argument safe=False.
                """.format(ram=free_ram/1024**2, peak_lim=peak_lim, Npix=pixel_size, tx=tx.size, N=pixel_size*tx.size, ss=pixel_size*tx.size*4/1024**2))
                
        size = (pixel_size, tx.size)
        spec = np.zeros(size, dtype='float32')
        for scan in IT(scans):
            raw = self.get_raw_raw_data(scan)
            rawv = struct.unpack('<{}I'.format(len(raw)//4), raw)
            k = 0
            while k<len(rawv):
                b = rawv[k]
                if b& 0xc0000000:
                    x = b & 0x0fffffff
                    y = rawv[k+1] & 0x0fffffff
                    dt = dts[x]
                    #fp = dt%1
                    ip = int(dt)
                    k += 3
                    i = (self.size['pixels']['x']//pixel_aggregation)*(y//pixel_aggregation)+x//pixel_aggregation
                else:
                    j1 = rev[b-ip]
                    if j1<0:
                        j1 = closest_arg(tx, b-ip)
                    #j2 = closest_arg(tx, b-ip-1)
                    spec[i][j1] += 1
                    #if j2>0:
                    #    spec[i][j2] += fp
                    k += 1   
        if prog:
            pb.update(1)
            pb.set_postfix({'task':'smooth spectra'})

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
            pb.update(1)
            pb.set_postfix({'task':'done'})
            
        return m>peak_lim, spec
            
    @alias("showSpectrum")
    def show_spectrum(self, low=0, high=None, sf=None, k0=None, ax=None, log=False, show_peaks=False, **kargs):
        """
        Plot the (summed) spectrum
        low and high: mass boundary of the plotted data
        ax: matplotlib axis
        log: plot the log of the intensity if True
        """
        # old notation compatibility
        if 'showPeaks' in kargs:
            warn("The parameter showPeaks is deprecated. Please use show_peaks")
            show_peaks = kargs.pop("showPeaks")
            
        m, s = self.get_spectrum(sf=sf, k0=k0, **kargs)
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
        if 'label' not in kargs:
            kargs['label'] = self.label
        ax.plot(M, S, **kargs)
        self.get_masses()
        if show_peaks:
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
                'umass': p['umass']['float'],
                'SN': p['SN']['utf16']})
        return result

    @alias("getRawSpectrum")
    def get_raw_spectrum(self, scans=None, ROI=None, FOVcorr=True, deadTimeCorr=True, **kargs):
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
            
        gun = self.root.goto('propend/Instrument.PrimaryGun.Species').get_key_value()['string'] # Primary Gun Species (Bi1,Bi3,Bi3++)
        
        # if the + is missing in the name, add it
        if gun[-1]!='+':
            gun += '+'
            
        Q = gun.count('+') # number of charge
        
        nrj = self.root.goto('propend/Instrument.PrimaryGun.Energy').get_key_value()['float'] # Primary ion energy (in eV)
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
        number_channels = int(round(self.root.goto('propend/Measurement.CycleTime').get_key_value()['float']\
            / self.root.goto('propend/Registration.TimeResolution').get_key_value()['float']))
        
        if kargs.get('prog',False):
            T = PB(scans)
        else:
            T = scans
        if kargs.get('debug', False):
            import time
            t0 = time.time()
        dts = DT*(self.size['pixels']['x']/2-np.arange(self.size['pixels']['x'])) # time correction for the given x coordinate (in channel number)
        if ROI is None:
            Spectrum = np.zeros(number_channels, dtype=np.float32)
            for s in T:
                raw = self.get_raw_raw_data(s)
                rawv = struct.unpack('<{}I'.format(len(raw)//4), raw)
                i = 0
                while i < len(rawv):
                    b = rawv[i]
                    if b & 0xc0000000:
                        x = b & 0x0fffffff
                        dt = dts[x]
                        ip = int(dt)
                        fp = dt%1
                        i += 3
                    else:
                        Spectrum[b-ip] += (1-fp)
                        Spectrum[b-ip-1] += fp
                        i += 1
        elif type(ROI) is np.ndarray:
            assert np.min(ROI)>=0
            Spectrum = np.zeros((number_channels, np.max(ROI)+1), dtype=np.float32)
            for s in T:
                raw = self.get_raw_raw_data(s)
                rawv = struct.unpack('<{}I'.format(len(raw)//4), raw)
                i = 0
                while i < len(rawv):
                    b = rawv[i]
                    if b & 0xc0000000:
                        x = b & 0x0fffffff
                        y = rawv[i+1] & 0x0fffffff
                        dt = dts[x]
                        fp = dt%1
                        ip = int(dt)
                        id = ROI[y, x]
                        i += 3
                    else:
                        Spectrum[b-ip, id] += (1-fp)
                        Spectrum[b-ip-1, id] += fp
                        i += 1
        elif type(ROI) in [list, tuple]:
            multi_roi = True
            Spectrum = np.zeros((number_channels, len(ROI)), dtype=np.float32)
            for s in T:
                raw = self.get_raw_raw_data(s)
                rawv = struct.unpack('<{}I'.format(len(raw)//4), raw)
                i = 0
                while i < len(rawv):
                    b = rawv[i]
                    if b & 0xc000000000:
                        x = b & 0x0fffffff
                        y = rawv[i+1] & 0x0fffffff
                        dt = dts[x]
                        fp = dt%1
                        ip = int(dt)
                        li = []
                        for k, R in enumerate(ROI):
                            if R[y,x]:
                                li.append(k)
                        i += 3
                    else:
                        for k in li:
                                Spectrum[b-ip, k] += (1-fp)
                                Spectrum[b-ip-1, k] += fp
                        i += 1
        if kargs.get('debug', False):
            t1 = time.time()
            print("Sepctra calc. time: ", t1-t0)
            t0 = t1
        sf = kargs.get('sf', self.root.goto('MassScale/sf').get_double())
        k0 = kargs.get('k0', self.root.goto('MassScale/k0').get_double())
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
        if kargs.get('debug', False):
            t1 = time.time()
            print("Dead time correction time: ", t1-t0)
            t0 = t1
        if kargs.get('time', False):
            return np.arange(number_channels), Spectrum
        return masses, Spectrum

    @alias("getRawRawData")
    def get_raw_raw_data(self, scan=0):
        assert scan < self.Nscan
        found = False
        RAW = b''
        if self.rawlist is None:
            self.rawlist = self.root.goto('rawdata').get_list()
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

    def get_pixel_order(self, scan=0):
        raw = self.get_raw_raw_data(scan)
        rawv = struct.unpack('<{}I'.format(len(raw)//4), raw)
        pixel_order = np.zeros((self.size['pixels']['y'], self.size['pixels']['x']))
        i = 0
        while i < len(rawv):
            b = rawv[i]
            if b & 0xc0000000:
                x = b & 0x0fffffff
                y = rawv[i+1] & 0x0fffffff
                id = rawv[i+2] & 0x3fffffff
                pixel_order[y,x] = id
                i += 3
            else:
                i += 1
        return pixel_order
        
    @alias("getRawData")
    def get_raw_data(self, scan=0):
        """
        Function which allows you to read and parse the raw data.
        It return a dictionary of list where each key is the pixel position (x,y) and the value is the list of all times (channel number).
        """
        raw = self.get_raw_raw_data(scan)
        rawv = struct.unpack('<{}I'.format(len(raw)//4), raw)
        blocks = {}
        i = 0
        _block = []
        x,y = -1,-1
        while i < len(rawv):
            b = rawv[i]
            if b & 0xc0000000:
                blocks[(x,y)] = _block # Somehow it's faster to append the data to a list first and then attribute it to the dict
                x = b & 0x0fffffff
                y = rawv[i+1] & 0x0fffffff
                _block = []
                i += 3
            else:
                _block.append(b)
                i += 1
        del blocks[(-1,-1)]
        return blocks

    def show_masses(self, mass_list=None):
        """
        Display the peak list (assignment name with lower, center and upper mass)
        """
        for m in self.get_masses():
            if mass_list is None or m['id'] in [z['id'] for z in mass_list]:
                print(
                    "{id}: ({desc}) [{assign}] {lmass:.2f}u - {umass:.2f}u (center: {cmass:.2f}u)".format(**m))

    @alias("showStage")
    def show_stage(self, ax=None, markers=False, plist=False):
        """
        Display an image of the stage used
        ax: maplotlib axis to be ploted in. If None the current one (or new) will be used
        markers: If True will display on the map the Position List items.
        """
        import pickle
        import matplotlib as mpl
        W = self.root.goto('SampleHolderInfo/bitmap/res_x').get_ulong()
        H = self.root.goto('SampleHolderInfo/bitmap/res_y').get_ulong()

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 1, figsize=(W*10/H, 10))

        Dat = zlib.decompress(self.root.goto(
            'SampleHolderInfo/bitmap/imagedata').value)
        I = np.array(struct.unpack("<"+str(W*H*3)+"B", Dat), dtype=np.uint8).reshape((H, W, 3))
        ax.imshow(I)
        if markers:
            X = self.root.goto('Meta/SI Image[0]/stageposition_x').get_double()
            Y = self.root.goto('Meta/SI Image[0]/stageposition_y').get_double()

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

    @alias("getSnapshot")
    def get_snapshot(self):
        """
        Return the video snapshot
        """
        try:
            dl = self.root.goto('Meta/Video Snapshot').dict_list()
            sx = dl['res_x']['ulong']
            sy = dl['res_y']['ulong']
            img = np.array(self.root.goto('Meta/Video Snapshot/imagedata').get_data('B')).reshape((sy, sx, 3))
            return  img
        except Exception as e:
            return None
            
    @alias("showPeaks")
    def show_peaks(self):
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
            scans = PB(scans)
        left = np.array([x[0] for x in channels])
        right = np.array([x[1] for x in channels])
        if not time:
            if sf is None or k0 is None:
                sf, k0 = self.get_mass_cal()
            left = mass2time(left, sf=sf, k0=k0)
            right = mass2time(right, sf=sf, k0=k0)
        Counts = [np.zeros((self.size['pixels']['x'], self.size['pixels']['y'])) for x in channels]
        for s in scans:
            Data = self.get_raw_data(s)
            for xy in Data:
                for i,ch in enumerate(channels):
                    Counts[i][xy[1],xy[0]] += np.sum((Data[xy]>=left[i])*(Data[xy]<=right[i]))
        res = [SPM_image(C, real=self.size['real'], _type='TOF', channel="{0[0]:.2f}{unit}-{0[1]:.2f}{unit}".format(channels[i],unit=["u", "s"][time]), zscale="Counts") for i,C in enumerate(Counts)]
        if len(res) == 1:
            return res[0]
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
            up = int(np.round(self.get_value("Measurement.CycleTime")['float']/self.get_value("Registration.TimeResolution")['float'], 0))
            umass = self.channel2mass(up)
        new_index = max([x.head['ID'] for x in self.root.goto("MassIntervalList") if x.name == 'mi'])+1
        new_index2 = self.root.goto("Measurement Options/massintervals/NextId").get_ulong()
        new_id = max([x.goto("id").get_ulong() for x in self.root.goto("MassIntervalList") if x.name == 'mi'])+1
        new_id2 = max([x.goto("id").get_ulong() for x in self.root.goto("Measurement Options/massintervals") if x.name == 'mi'])+1
        NID = self.root.goto("MassIntervalList/NextId", lazy=True)
        NID.rewrite(struct.pack("<I", NID.get_ulong()+1))
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
        p = dict(new_index=new_index, new_index2=new_index2)
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
        d = blk.dict_list()
        self.peaks[d['id']['long']] = d
        return blk
      
    def __del__(self):
        if self.f.writable():
            self.f.flush()
            os.fsync(self.f)
        self.f.close()

    @alias("setK0")
    def set_k0(self, k0):
        import struct
        try:
            b = self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/IMassScaleSFK0/k0')
        except:
            b = self.root.goto("MassScale/k0")
        buffer = struct.pack("<d", k0)
        b.rewrite(buffer);

    @alias("setSF")
    def set_sf(self, sf):
        import struct
        try:
            b = self.root.goto('filterdata/TofCorrection/Spectrum/Reduced Data/IMassScaleSFK0/sf')
        except:
            b = self.root.goto("MassScale/sf")
        buffer = struct.pack("<d", sf)
        b.rewrite(buffer);
