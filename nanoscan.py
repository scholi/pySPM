import os
import base64
import xml.etree.ElementTree as ET
import struct
import numpy as np
import pySPM.SPM
from pySPM.SPM import SPM_image, funit

def getCurve(filename, channel='Normal Deflection', backward=False):
    """
    function to retrieve data which are not in the form of images.
    This is typically used for 1D channel where the normal deflection is recorded while z is swept.
    """
    tree = ET.parse(filename)
    root = tree.getroot()
    namespace = {'spm':'http://www.nanoscan.ch/SPM'}
    RAW = root.findall("spm:vector/spm:contents/spm:direction/spm:vector/" \
        "spm:contents/spm:name[spm:v='{direction}']/../spm:channel/spm:vector/" \
        "spm:contents/spm:name[spm:v='{channel}']/../spm:data/spm:v" \
        .format(direction=['forward', 'backward'][backward], \
        channel=channel), namespace)[0].text
    start = float(root.findall("spm:vector/spm:contents/spm:axis/spm:vector/" \
        "spm:contents/spm:start/spm:vector/spm:v", namespace)[0].text)
    stop = float(root.findall("spm:vector/spm:contents/spm:axis/spm:vector/" \
        "spm:contents/spm:stop/spm:vector/spm:v", namespace)[0].text)
    unit = root.findall("spm:vector/spm:contents/spm:axis/spm:vector/spm:contents" \
        "/spm:unit/spm:v", namespace)[0].text
    BIN = base64.b64decode(RAW)
    N = len(BIN)
    vals = np.array(struct.unpack("<"+str(N//4)+"f", BIN))
    x = np.linspace(start, stop, len(vals))
    return x, vals

class Nanoscan(SPM_image):    
    def __init__(self, filename=None, channel='Topography', backward=False, **kargs):
        if not os.path.exists(filename):
            raise IOError('File "{0}" Not Found'.format(filename))
        if filename[-4:] != '.xml':
            raise TypeError("Nanoscan files should be xml files")
        self.filename = filename
        tree = ET.parse(filename)
        self.root = tree.getroot()
        
        if self.root.tag == "{http://www.nanoscan.ch/SPM}scan":
            namespaces = {'spm':"http://www.nanoscan.ch/SPM"}
            _type = "Nanoscan"
            try:
                RAW = self.root.findall("spm:vector//spm:direction/spm:vector/spm:contents" \
                    "/spm:name[spm:v='%s']/../spm:channel//spm:contents/spm:name[spm:v='%s']" \
                    "/../spm:data/spm:v"%(["forward", "backward"][backward], channel), \
                    namespaces)[0].text
            except:
                raise 'Channel {0} in {1} scan not found'.format(channel, direction)
                return
            pixel_size = [int(z.text) for z in self.root.findall("spm:vector/spm:contents/spm:size/spm:contents//spm:v", namespaces)]
            self.fbPath = "spm:vector/spm:contents/spm:instrumental_parameters/spm:contents/spm:z_control/spm:contents"
            
            uval = float(self.root.findall(".//spm:area//spm:contents/spm:size/spm:contents/spm:fast_axis/spm:v",namespaces)[0].text)
            udispu = self.root.findall(".//spm:area//spm:contents/spm:display_unit/spm:v",namespaces)[0].text
            udisps = float(self.root.findall(".//spm:area/spm:contents/spm:display_scale/spm:v",namespaces)[0].text)
            uname = self.root.findall(".//spm:area/spm:contents/spm:unit/spm:v",namespaces)[0].text
            x = funit(uval*udisps, udispu)
            uval = float(self.root.findall(".//spm:area//spm:contents/spm:size/spm:contents/spm:slow_axis/spm:v",namespaces)[0].text)
            y = funit(uval*udisps, udispu)
          
            BIN = base64.b64decode(RAW)
            recorded_size = len(BIN)/4
            size = {'pixels':{ \
                'x':pixel_size[0], \
                'y':pixel_size[1] \
            },'real':{ \
                'unit':x['unit'], \
                'x':x['value'], \
                'y':y['value'], \
            }}
            size['recorded'] = {\
                'pixels':{\
                    'y':int(recorded_size/pixel_size[0]),\
                    'x':pixel_size[0]}}
            size['recorded']['real'] = { \
                'x':size['real']['x'], \
                'y':size['real']['y']*size['recorded']['pixels']['y']/float(size['pixels']['y'])}

            image_array = np.array(struct.unpack("<%if"%(recorded_size),BIN)).reshape( \
                (size['recorded']['pixels']['y'],size['recorded']['pixels']['x']))
            
        elif self.root.tag == "channel_list": # ToF-SIMS data (old and no more used. Kept for old script compatibility)
            _type = "ToF-SIMS (xml)"
            channel = "Counts"
            x   = int(self.root.findall("./channel/axis[name='x']/count")[0].text)
            y   = int(self.root.findall("./channel/axis[name='y']/count")[0].text)
            RAW = self.root.findall("./channel/pixels")[0].text
            BIN = base64.b64decode(RAW)
            image_array = np.array(struct.unpack("<%if"%(x*y),BIN)).reshape(x, y)
            size = {
                'pixels':{ \
                    'x':x, \
                    'y':y \
                }, 'real':{ \
                    'unit':'m', \
                    'x':float(self.root.findall("./channel/axis[name='x']/variable/extent")[0].text), \
                    'y':float(self.root.findall("./channel/axis[name='y']/variable/extent")[0].text)}, \
                'recorded':{'real':{ \
                    'unit':'m', \
                    'x':float(self.root.findall("./channel/axis[name='x']/variable/extent")[0].text), \
                    'y':float(self.root.findall("./channel/axis[name='y']/variable/extent")[0].text)}}}
        SPM_image.__init__(self, image_array, channel=channel, **kargs)
        self.type = _type
        self.size = size
        self.direction = ['forward', 'backward'][backward]
        self.namespaces = namespaces
        
    def get_scanspeed(self):
        print(self.direction,self.namespaces)
        return {'value':float(self.root.findall("spm:vector//spm:direction/spm:vector/" \
                    "spm:contents/spm:name[spm:v='%s']/../spm:point_interval/spm:v"%(self.direction), self.namespaces)[0].text) * self.size['pixels']['x'], \
                'unit': self.root.findall("spm:vector//spm:direction/spm:vector/" \
                    "spm:contents/spm:name[spm:v='%s']/../spm:point_interval_unit/spm:v"%(self.direction), self.namespaces)[0].text}

    def get_feedback(self):
            self.feedback = {'channel':\
                self.root.findall('{0}/spm:z_feedback_channel/spm:v'.format(self.fbPath), namespaces)[0].text}
            self.feedback['P'] = {\
                'value':float(self.root.findall('{0}/spm:proportional_z_gain/spm:v'\
                    .format(self.fbPath), self.namespaces)[0].text),
                'unit':self.root.findall('{0}/spm:proportional_z_gain_unit/spm:v'\
                    .format(self.fbPath), self.namespaces)[0].text}
            self.feedback['I'] = {\
                'value':float(self.root.findall('{0}/spm:integral_z_time/spm:v'\
                    .format(fbPath), self.namespaces)[0].text),
                'unit':self.root.findall('{0}/spm:integral_z_time_unit/spm:v'\
                    .format(fbPath), self.namespaces)[0].text}
            if self.feedback['channel'] == 'df':
                self.feedback['channel'] = u'?f'