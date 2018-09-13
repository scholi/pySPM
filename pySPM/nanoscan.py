# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

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
    namespace = {'spm': 'http://www.nanoscan.ch/SPM'}
    RAW = root.findall("spm:vector/spm:contents/spm:direction/spm:vector/"
                       "spm:contents/spm:name[spm:v='{direction}']/../spm:channel/spm:vector/"
                       "spm:contents/spm:name[spm:v='{channel}']/../spm:data/spm:v"
                       .format(direction=['forward', 'backward'][backward],
                               channel=channel), namespace)[0].text
    start = float(root.findall("spm:vector/spm:contents/spm:axis/spm:vector/"
                               "spm:contents/spm:start/spm:vector/spm:v", namespace)[0].text)
    stop = float(root.findall("spm:vector/spm:contents/spm:axis/spm:vector/"
                              "spm:contents/spm:stop/spm:vector/spm:v", namespace)[0].text)
    unit = root.findall("spm:vector/spm:contents/spm:axis/spm:vector/spm:contents"
                        "/spm:unit/spm:v", namespace)[0].text
    BIN = base64.b64decode(RAW)
    N = len(BIN)
    vals = np.array(struct.unpack("<"+str(N//4)+"f", BIN))
    x = np.linspace(start, stop, len(vals))
    return x, vals


class Nanoscan():

    def __init__(self, filename=None):
        if not os.path.exists(filename):
            raise IOError('File "{0}" Not Found'.format(filename))
        if filename[-4:] != '.xml':
            raise TypeError("Nanoscan files should be xml files")
        self.filename = filename
        tree = ET.parse(filename)
        self.root = tree.getroot()

        if self.root.tag == "{http://www.nanoscan.ch/SPM}scan":
            self.namespaces = {'spm': "http://www.nanoscan.ch/SPM"}
            self.type = "Nanoscan"
            self.fbPath = "spm:vector/spm:contents/spm:instrumental_parameters/spm:contents/spm:z_control/spm:contents"
            self.pixel_size = [int(z) for z in self.__grab(
                "spm:vector/spm:contents/spm:size/spm:contents//spm:v")]
            uval = float(self.__grab(
                ".//spm:area//spm:contents/spm:size/spm:contents/spm:fast_axis/spm:v"))
            udispu = self.__grab(
                ".//spm:area//spm:contents/spm:display_unit/spm:v")
            udisps = float(self.__grab(
                ".//spm:area/spm:contents/spm:display_scale/spm:v"))
            uname = self.__grab(".//spm:area/spm:contents/spm:unit/spm:v")
            x = funit(uval*udisps, udispu)
            uval = float(self.__grab(
                ".//spm:area//spm:contents/spm:size/spm:contents/spm:slow_axis/spm:v"))
            y = funit(uval*udisps, udispu)
            self.size = {
                'unit': x['unit'],
                'x': x['value'],
                'y': y['value']}
            try:
                self.feedback = {'channel': self.__grab(
                    '{0}/spm:z_feedback_channel/spm:v'.format(self.fbPath))}
                self.feedback['P'] = {
                    'value': float(self.__grab('{0}/spm:proportional_z_gain/spm:v'.format(self.fbPath))),
                    'unit': self.__grab('{0}/spm:proportional_z_gain_unit/spm:v'.format(self.fbPath))}
                self.feedback['I'] = {
                    'value': float(self.__grab('{0}/spm:integral_z_time/spm:v'.format(self.fbPath))),
                    'unit': self.__grab('{0}/spm:integral_z_time_unit/spm:v'.format(self.fbPath))}
                if self.feedback['channel'] == 'df':
                    self.feedback['channel'] = u'Δf'
            except:
                self.feedback = {}
            self.scan_speed = {
                z: {
                    'value': float(self.__grab("spm:vector//spm:direction/spm:vector/spm:contents/spm:name[spm:v='{dir}']/../spm:point_interval/spm:v".format(dir=z))) * self.pixel_size[0],
                    'unit': self.__grab("spm:vector//spm:direction/spm:vector/spm:contents/spm:name[spm:v='{dir}']/../spm:point_interval_unit/spm:v".format(dir=z))} for z in ['forward', 'backward']}
        else:
            raise TypeError(
                "Unknown or wrong data type. Expecting a valid Nanoscan xml")
    
    def list_channels(self):
        """
        Printout the list of stored channels
        """
        for d in ['forward','backward']:
            print(d)
            print("="*len(d))
            for z in self.root.findall("spm:vector//spm:direction/spm:vector/spm:contents/spm:name[spm:v='{}']/../spm:channel//spm:contents/spm:name/spm:v".format(d), self.namespaces):
                print("  - "+z.text)
            print()

    def get_channel(self, channel='Topography', backward=False, corr=None):
        try:
            RAW = self.__grab("spm:vector//spm:direction/spm:vector/spm:contents"
                              "/spm:name[spm:v='{direction}']/../spm:channel//spm:contents/spm:name[spm:v='{channel}']"
                              "/../spm:data/spm:v".format(direction=["forward", "backward"][backward], channel=channel))
        except:
            raise 'Channel {0} in {1} scan not found'.format(
                channel, direction)
            return None

        BIN = base64.b64decode(RAW)
        recorded_length = len(BIN)/4

        py = int(recorded_length/self.pixel_size[0])
        recorded_size = {
            'x': self.size['x'],
            'y': self.size['y']*py/float(self.pixel_size[1]),
            'unit': self.size['unit']}

        image_array = np.array(struct.unpack("<%if" % (recorded_length), BIN)).reshape(
            (py, self.pixel_size[0]))
        return SPM_image(image_array, channel=channel, _type=self.type, real=recorded_size, corr=corr)

    def __grab(self, path):
        result = [z.text for z in self.root.findall(path, self.namespaces)]
        if len(result) == 1:
            result = result[0]
        return result

    def arraySummary(self):
        from pySPM.utils import htmlTable
        res = [y.format(**self.__dict__) for y in
               ["{filename}", "{pixel_size[0]}×{pixel_size[1]}", "{size[x][value]}×{size[y][value]} {size[x][unit]}",
                "{scan_speed[forward][value]} {scan_speed[forward][unit]}",
                "{feedback[channel]}", "{P[value]:.2f} {P[unit]}", "{I[value]:.2f} {I[unit]}"]]

    def getSummary(self):
        x = funit(self.size['x'], self.size['unit'])
        y = funit(self.size['y'], self.size['unit'])
        P = funit(self.feedback['P'])
        I = funit(self.feedback['I'])
        return u"""Feedback: {feedback[channel]} : P:{P[value]}{P[unit]} : I:{I[value]}{I[unit]}
Size: {pixel_size[0]}×{pixel_size[1]} pixels = {x[value]:.3} {x[unit]}×{y[value]:.3} {y[unit]}
Scan Speed: {scanSpeed[value]}{scanSpeed[unit]}/line""".format(
            x=x, y=y, P=P, I=I,
            feedback=self.feedback, pixel_size=self.pixel_size, size=self.size,
            scanSpeed=self.scan_speed['forward'])

    @staticmethod
    def show_dir_summary(path):
        from pySPM.utils import htmlTable
        res = [["Filename", "pixel size", "real size",
                "scan_speed", "feedback", "P", "I"]]
        for x in os.listdir(path):
            try:
                A = pySPM.Nanoscan(path+x)
                res.append([y.format(f=os.path.basename(A.filename), **A.__dict__) for y in
                        ["{f}", "{pixel_size[0]}×{pixel_size[1]}", "{size[x]}×{size[y]} {size[unit]}", "{scan_speed[forward][value]} {scan_speed[forward][unit]}",
                         "{feedback[channel]}", "{feedback[P][value]:.2f} {feedback[P][unit]}", "{feedback[I][value]:.2f} {feedback[I][unit]}"]])
            except:
                print("Cannot read image \""+x+"\" skipping it")
        htmlTable(res, header=True)
