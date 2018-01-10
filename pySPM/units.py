# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Deprecated module. Will be removed in final version
"""

import math

mag = "afpnum1kMGTPE"
SI_units = ['m', 'kg', 'K', 'A', 's', 'cd', 'mol']
units_conv = {'Hz': 's-1', 'N': 'kg*m*s-2',
              'Pa': 'kg*m-1*s-2', 'J': 'kg*m2*s-2', 'V': 'kg*m2*s-3*A-1'}


class unit:

    def __init__(self, name, value=1):
        if name[0] in mag and len(name) > 1:
            self.imag = mag.index(name[0])
            self.name = name[1:]
        else:
            self.imag = mag.index('1')
            self.name = name
        self.value = value

    def __add__(self, b):
        if self.name != b.name:
            raise TypeError
        if self.imag > b.imag:
            value = self.value + b.value * 10**(b.imag-a.imag)
        else:
            value = self.value * 10**(a.imag-b.imag) + b.value
        return unit(self.name, value)

    def __str__(self):
        return "{value} {unit}".format(**self)

    def __repr__(self):
        shift = math.floor(math.log10(self.value)/3)
        prefix = ""
        return "{prefix}{value} {unit}".format(unit=self.name, value=self.value)
