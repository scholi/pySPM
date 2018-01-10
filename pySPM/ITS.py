# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

from pySPM import Block, utils
import numpy as np
import struct
import os.path
import zlib
import re

class InvalidRAWdataformat(Exception):
    def __init__(self, block, msg):
        self.block = block
        self.msg = msg
        
    def __str__(self):
        return "Invalid RAW dataformat seen in block "+self.block.parent+'/'+self.block.name+' : '+self.msg
    
class ITS:
    def __init__(self, filename, debug=False):
        """
        ITS
        """
        self.filename = filename
        assert os.path.exists(filename)
        self.f = open(self.filename, 'rb')
        self.Type = self.f.read(8)
        assert self.Type == b'ITStrF01'
        self.root = Block.Block(self.f)