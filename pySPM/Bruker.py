# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Module to handle SPM images recorded by a Bruker AFM
"""

import contextlib
import re
import struct

import numpy as np

import pySPM


class Bruker:
    """
    Class to handle SPM images recorded by a Bruker AFM
    """

    def __init__(self, path):
        self.path = path
        with open(self.path, "rb") as file:
            self.layers = []
            self.scanners = []
            mode = ""
            while True:
                line = file.readline().rstrip().replace(b"\\", b"")
                if line == b"*Ciao image list":
                    self.layers.append({})
                    mode = "Image"
                elif line == b"*Scanner list":
                    self.scanners.append({})
                    mode = "Scanner"
                elif line.startswith(b"*EC"):
                    mode = "EC"
                else:
                    args = line.split(b": ")
                    if len(args) > 1:
                        if mode == "Image":
                            self.layers[-1][args[0]] = args[1:]
                        elif mode == "Scanner":
                            self.scanners[-1][args[0]] = args[1:]
                    if line == b"*File list end":
                        break

    def _get_bpp(self, i):
        int(self.layers[i][b"Data offset"][0])
        cols = int(self.layers[i][b"Number of lines"][0])
        rows = int(self.layers[i][b"Samps/line"][0])
        byte_length = int(self.layers[i][b"Data length"][0])
        length = rows * cols
        bpp = byte_length // length
        return bpp

    def _get_raw_layer(self, i, debug=False, mock_data=False):
        """
        Internal function to retrieve raw data of a layer
        """
        off = int(self._get_layer_val(i, "Data offset"))
        if debug:
            print("RAW offset: ", off)
        cols, rows = self._get_res(i)
        byte_length = int(self._get_layer_val(i, "Data length"))
        length = rows * cols
        bpp = byte_length // length
        byte_length = length * bpp
        if mock_data:
            return np.zeros((rows, cols))

        with open(self.path, "rb") as file:
            file.seek(off)
            return np.array(
                struct.unpack(
                    "<" + str(length) + {2: "h", 4: "i", 8: "q"}[bpp],
                    file.read(byte_length),
                ),
                dtype="float64",
            ).reshape((rows, cols))

    def list_channels(self, encoding="latin1"):
        print("Channels")
        print("========")

        for layer in self.layers:
            with contextlib.suppress(KeyError):
                print(layer[b"@2:Image Data"][0].decode(encoding))
        for layer in self.layers:
            with contextlib.suppress(KeyError):
                print(layer[b"@3:Image Data"][0].decode(encoding) + " (MFM)")

    def _get_layer_val(self, index: int, name: str, first=True):
        lname = name.lower()
        lname2 = name[0] + lname[1:]
        if name.encode() in self.layers[index]:
            val = self.layers[index][name.encode()]
        elif lname.encode() in self.layers[index]:
            val = self.layers[index][lname.encode()]
        elif lname2.encode() in self.layers[index]:
            val = self.layers[index][lname2.encode()]
        if first:
            return val[0]
        return val

    def _get_res(self, layer_index):
        row_key = (
            "Valid data len X"
            if b"Valid data len X" in self.layers[layer_index]
            else "Number of lines"
        )
        col_key = (
            "Valid data len Y"
            if b"Valid data len Y" in self.layers[layer_index]
            else "Samps/line"
        )
        xres = int(self._get_layer_val(layer_index, row_key))
        yres = int(self._get_layer_val(layer_index, col_key))
        return xres, yres

    def _get_layer_size(self, layer_index, encoding, debug=False):
        scan_size = self._get_layer_val(layer_index, "Scan Size").split()
        xres, yres = self._get_res(layer_index)
        if scan_size[2][0] == 126:
            scan_size[2] = b"u" + scan_size[2][1:]
        size = {}
        size = {
            "x": float(scan_size[0]),
            "y": float(scan_size[1]) * yres / xres,
            "unit": scan_size[2].decode(encoding),
        }
        if debug:
            print("scan size", scan_size)

        return size

    def get_channel(
        self,
        channel="Height Sensor",
        backward=False,
        corr=None,
        debug=False,
        encoding="latin1",
        lazy=True,
        mfm=False,
        mock_data=False,
    ):
        """
        Load the SPM image contained in a channel
        """
        for i in range(len(self.layers)):
            if mfm:
                # Handle case of MFM, where image data is stored at layer 3 on the backward channel
                backward = True
                _type = "Bruker MFM"
                try:
                    layer_name = self._get_layer_val(i, "@3:Image Data").decode(
                        encoding
                    )
                except KeyError:
                    continue
            else:
                _type = "Bruker AFM"
                try:
                    layer_name = self._get_layer_val(i, "@2:Image Data").decode(
                        encoding
                    )
                except KeyError:
                    continue
            result = re.match(r'([^ ]+) \[([^]]*)] "([^"]*)"', layer_name).groups()
            if result[2] == channel:
                if debug:
                    print("channel " + channel + " Found!")
                bck = False
                if self._get_layer_val(i, "Line Direction") == b"Retrace":
                    bck = True
                if bck == backward:
                    if debug:
                        print("Direction found")
                    var = self._get_layer_val(i, "@2:Z scale").decode(encoding)
                    if debug:
                        print("@2:Z scale", var)
                    if "[" in var:
                        result = re.match(
                            r"[A-Z]+\s+\[([^]]+)]\s+\(-?[0-9.]+ .*?\)\s+(-?[0-9.]+)\s+(.*?)$",
                            var,
                        ).groups()
                        if debug:
                            print(result)
                        bpp = int(self._get_layer_val(i, "Bytes/pixel"))
                        if debug:
                            print("BPP", bpp)
                        # scale = float(result[1])
                        scale = float(result[1]) / 256**bpp

                        result2 = self.scanners[0][b"@" + result[0].encode(encoding)][
                            0
                        ].split()
                        if debug:
                            print("result2", result2)
                        scale2 = float(result2[1])
                        zscale = result2[2] if len(result2) > 2 else result2[0]
                        if b"/V" in zscale:
                            zscale = zscale.replace(b"/V", b"")
                        if debug:
                            print(f"scale: {scale:.3e}")
                            print(f"scale2: {scale2:.3e}")
                            print("zscale: " + str(zscale))
                        var = self._get_layer_val(i, "@2:Z offset").decode(encoding)
                        result = re.match(
                            r"[A-Z]+\s+\[[^]]+]\s+\(-?[0-9.]+ .*?\)\s+(-?[0-9.]+)\s+.*?$",
                            var,
                        ).groups()
                        offset = float(result[0])
                    else:
                        if debug:
                            print("mode 2")
                        result = re.match(
                            r"[A-Z]+ \(-?[0-9.]+ [^)]+\)\s+(-?[0-9.]+) [\w]+", var
                        ).groups()
                        scale = float(result[0]) / 65536.0
                        scale2 = 1
                        zscale = b"V"
                        result = re.match(
                            r"[A-Z]+ \(-?[0-9.]+ .*?\)\s+(-?[0-9.]+) .*?",
                            self._get_layer_val(i, "@2:Z offset").decode(encoding),
                        ).groups()
                        offset = float(result[0])
                    if debug:
                        print("Offset:", offset)
                    data = (
                        self._get_raw_layer(i, debug=debug, mock_data=mock_data)
                        * scale
                        * scale2
                    )
                    xres, yres = self._get_res(i)
                    if debug:
                        print("xres/yres", xres, yres)
                    scan_size = self._get_layer_val(i, "Scan Size").split()
                    aspect_ratio = [
                        float(x)
                        for x in self._get_layer_val(i, "Aspect Ratio").split(b":")
                    ]
                    if debug:
                        print("aspect ratio", aspect_ratio)
                        print("scan size", scan_size)
                    if scan_size[2][0] == 126:
                        scan_size[2] = b"u" + scan_size[2][1:]
                    size = {
                        "x": float(scan_size[0]),
                        "y": float(scan_size[1]) * yres / xres,
                        "unit": scan_size[2].decode(encoding),
                    }
                    image = pySPM.SPM_image(
                        channel=[channel, "Topography"][channel == "Height Sensor"],
                        BIN=data,
                        real=size,
                        _type=_type,
                        zscale=zscale.decode(encoding),
                        corr=corr,
                    )
                    return image
        if lazy:
            return self.get_channel(
                channel=channel,
                backward=not backward,
                corr=corr,
                debug=debug,
                encoding=encoding,
                lazy=False,
                mock_data=mock_data,
            )
        raise Exception(f"Channel {channel} not found")
