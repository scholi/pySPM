"""
Module to handle block type used in iontof file formats ITA,ITM,ITS, etc...
"""

import sys
import binascii
import struct
import numpy as np


class Block:
    """
    Class to handle a iontof-Block
    One iontof file ITA,ITM,ITS contains a lot of Blocks forming a hierarchical structure.
    Each Block can have children (sub-Blocks) and values (data)
    """
    def __init__(self, fp):
        """
        Init the class
        fp: file pointer (the one created by open(...) of an ITA,ITM,ITS, etc... file pointing at the beginning of a block
        
        Each block start with one byte of type followed by 4 bytes that should always be \x19\x00\x00\x00 (all those 5 bytes are saved in self.Type)
        Then follows 5 uint32: length, z, u ,x ,y
            length: The length of the block's name
            z: Not used. My guess is that it's the block ID. So it starts at 0 and is increased for all following Blocks of the same name.
                We usually find the ID from the children's list (see below) and this information is never used as it's redundant.
            u: The number of children / sub-blocks. Might be = 0 even if the black has children. Check the value L (defined below) if so
            x: The length of the block's value
            y: Redundant. Seems to be always = x
        Then follow length-bytes representing the name of the block
        Then follow x-bytes forming the value of the block
        
        Blocks of types \x01\x19\x00\x00\x00 and \x03\x19\x00\x00\x00 are blocks that contains sub-blocks. There is no big difference between the two. I guess that types \x01 is the first one and type \x03 are the continuation blocks
            Those block have a value which starts with 41-bytes.
                3 uint32 -> (length, nums, ID).
                    length: We actually don't need it. It's a redundant information. That is the length of the sub-headers. (It stop just before the sub-blocks names)
                    nums: redundant info. Not used
                    ID: Block ID
                9 unknown padding bytes
                1 uint32 -> L
                    The variable u (see above) contains the number of children. If u == 0, then L will tell the correct number of children
                8 unknown padding bytes
                1 uint64 -> NextBlock
                    Big blocks can be chunked in several ones. NextBlock tells the position in the file of the next chunk. If = 0, then it's the last chunk
            Then 33 bytes for each sub-block follows:
                3 uint32 -> index, slen, id
                    index: The position of the sub-block name in the header
                    slen: The length of the sub-block name (which is store later). So basically the sub-block name is: Block.value[index:index+slen]
                    id: start at 0 and increase monotonically for sub-blocks having the same name
                4 unknown padding bytes
                2 uint64 -> blen, bidx
                    blen: Block length
                    bidx: Position of the Block in the file
            All the names of all sub-blocks follows (concatenated). You need to use their name length (see slen above) in order to chunk them properly
            You should then go to the position NextBlock in the file and read the block their. It's the continuation of the current one.
            
        Blocks of types \x00\x19\x00\x00\s00 have no children and only a value, usually formed by a string (UTF-16), an int or a float.
        The type of the value was not discovered to be written in the Block. The user should deduce it depending on the block's name and location.
        
        Blocks of type \x80\x19\x00\x00\x00 have it's values compressed with zlib (notice the \x78\x5e at the beginning which is typical for zlib encoded data)
            Once decompressed it usually store an array in binary.
        """
        self.f = fp
        self.offset = self.f.tell()
        self.Type = self.f.read(5)
        if self.Type[1:] != b'\x19\x00\x00\x00':
            raise ValueError('Wrong block type ({Type}) found @{pos}'\
                .format(pos=self.offset, Type=binascii.hexlify(self.Type[1:])))
        if len(self.Type) < 5:
            raise ValueError('EOF reached. Block cannot be read')
        self.head = dict(zip(['length', 'z', 'u', 'x', 'y'], \
            struct.unpack('<5I', self.f.read(20))))
        self.name = self.f.read(self.head['length'])
        self.value = self.f.read(self.head['x'])
        self.List = None
        self.iterP = 0

    def getName(self):
        """
        Return the name of the Block
        """
        return self.name

    def getList(self):
        """
        Return the list of the sub-blocks (chrildren) of the current Block.
        """
        if not self.Type[0:1] in [b'\x01', b'\x03']:
            return []
        if self.List is None:
            return self.createList()
        return self.List

    def createList(self, limit=None):
        """
        Generate a list (self.List) containing all the children (sub-blocks) of the current Block
        """
        offset = self.offset
        self.List = []
        while True:
            self.f.seek(offset)
            head = dict(zip(['length', 'z', 'u', 'x', 'y'], \
                struct.unpack('<5x5I', self.f.read(25))))
            name = self.f.read(head['length'])
            data = self.f.tell()
            length, nums, ID, L, NextBlock = struct.unpack('<III9xI8xQ', self.f.read(41))
            N = head['u']
            if N == 0:
                N = nums
            for i in range(N):                
                self.f.seek(data+42+33*i)
                S = dict(\
                    zip(['index', 'slen', 'id', 'blen', 'bidx'],\
                    struct.unpack('<III4xQQ', self.f.read(32))))
                self.f.seek(data+S['index'])
                S['name'] = self.f.read(S['slen'])
                self.List.append(S)
            if NextBlock == 0:
                break
            offset = NextBlock
        return self.List
            
    def getString(self):
        """
        Decode the value of the Block to UTF-16 (standard for all strings in this fileformat)
        """
        return self.value.decode('utf16')

    def dictList(self):
        """
        Return a dictionnary of the value decoded with various formats (raw, long, float, utf16)
        As the type of the data is not known, this function is very helpful for debugging purpose
        """
        d = {}
        for i, l in enumerate(self.getList()):
            self.f.seek(l['bidx'])
            child = Block(self.f)
            if child.Type[0:1] == b'\x00':
                value = binascii.hexlify(child.value)
                d[child.name] = {'raw':value}
                if len(child.value) == 4:
                    d[child.name]['long'] = child.getLong()
                elif len(child.value) == 8:
                    d[child.name]['float'] = child.getDouble()
                    d[child.name]['long'] = child.getLongLong()
                if len(child.value)%2 == 0:
                    d[child.name]['utf16'] = child.value.decode('utf16', "ignore")
            del child
        return d

    def showList(self):
        """
        Show a list of all the chrildren (sub-blocks) of the current Block.
        It will also display the value/data of all the children (if any)
        """
        print('List of', len(self.getList()), 'elements. Type:', self.subType)
        for i, l in enumerate(self.List):
            self.f.seek(l['bidx'])
            other = ''
            try:
                child = Block(self.f)
                if child.Type[0:1] == b'\x00':
                    if len(child.value) == 4:
                        vL = child.getLong()
                        Dtype = 'long'
                    elif len(child.value) == 8:
                        vL = child.getDouble()
                        Dtype = 'double'
                        other += ' = '+str(child.getLongLong())+" (long64)"
                    elif len(child.value)%2 == 0:
                        vL = child.value.decode('utf16', "ignore")
                        if len(vL) > 20:
                            vL = vL[:20]+'...'
                        Dtype = 'UTF-16'
                    elif len(child.value) == 2:
                        vL = child.getShort()
                        Dtype = 'short'
                    elif len(child.value) == 1:
                        vL = child.getByte()
                        Dtype = 'byte'
                    else:
                        vL = '???'
                        Dtype = '???'
                    value = binascii.hexlify(child.value)
                    if len(value) > 16:
                        value = value[:16]+b'...'
                    print(u"{name} ({id}) <{blen}> @{bidx}, value = {value} (hex) = {vL} ({Dtype}){other}"\
                        .format(value=value, vL=vL, Dtype=Dtype, other=other, **l))
                else:
                    print("{name} ({id}) [{T}] <{blen}> @{bidx}".format(T=child.Type[0], **l))
                del child
            except ValueError:
                pass

    def __iter__(self):
        """
        Return an iterator over all the children of the current block)
        """
        self.pointer = 0
        return self

    def __next__(self):
        L = self.getList()
        if self.pointer >= len(L):
            raise StopIteration
        it = L[self.pointer]
        self.pointer += 1
        return self.gotoItem(it['name'], it['id'])

    def gotoItem(self, name, idx=0):
        """
        Return a new Block instance of a child of the current Block
        name: name of the children's block
        """
        Idx = self.getIndex(name, idx)
        self.f.seek(Idx)
        return Block(self.f)

    def getIndex(self, name, idx=0):
        """
        Get the index of the children having a name=name.
        This function is more intended for internal usage.
        You are encouraged to use the function _goto_ instead
        
        If more than one children have the same name,
        the second one can by retrieved by idx=1, the third with idx=2, etc.
        """
        if type(name) is str:
            name = name.encode()
        for l in self.getList():
            if l['name'] == name and l['id'] == idx:
                return l['bidx']
        raise ValueError('Item "{name}" (index={index}) not found!'.format(name=name, index=idx))

    def goto(self, path):
        """
        Return a sub Block having a specific path
        path: path is similar to filepath.
        The block X conatained in B which is itself contained in A,
        can be retrieved with the path: A/B/X
        if the block B has several children having the same name,
        A/B/X[n] will return the n-th child (note that 0 is the first child)
        """
        s = self
        for p in path.split('/'):
            idx = 0
            if '[' in p and p[-1] == ']':
                i = p.index('[')
                idx = int(p[i+1:-1])
                p = p[:i]
            s = s.gotoItem(p, idx)
        return s

    def getLongLong(self):
        """
        Decode the value as an 64-Integer
        """
        return struct.unpack('<q', self.value)[0]

    def getDouble(self):
        """
        Decode the value as a 64-float (Double)
        """
        return struct.unpack('<d', self.value)[0]

    def getShort(self):
        """
        Decode the value as an 16-Integer (Short)
        """
        return struct.unpack('<h', self.value)[0]

    def getByte(self):
        """
        Decode the value as a 1-Byte
        """
        return struct.unpack('<B', self.value)[0]

    def getULong(self):
        """
        Decode the value as an unsigned 32-Integer (Long)
        """
        return struct.unpack('<I', self.value)[0]

    def getLong(self):
        """
        Decode the value as an 32-Integer (Long)
        """
        return struct.unpack('<i', self.value)[0]

    def getKeyValue(self, offset=0):
        """
        Return a dictionnary of key/values pairs of the data
        Note that the function has no idea if the data are stored as so.
        """
        L = struct.unpack("<I", self.value[offset:offset+4])[0]
        Key = self.value[offset+4:offset+4+L].decode('utf16', 'ignore')
        Value = struct.unpack("<10xd", self.value[offset+4+L:offset+22+L])[0]
        L2 = struct.unpack("<I", self.value[offset+22+L:offset+26+L])[0]
        SVal = self.value[offset+26+L:offset+26+L+L2].decode('utf16', 'ignore')
        return {'Key':Key, 'Value':Value, 'SVal':SVal}

    def show(self, maxlevel=3, level=0, All=False, out=sys.stdout, digraph=False, parent=None, ex=None):
        """
        Display the children of the current Block (recursively if maxlevel > 1)
        Very usefull for debugging purpose and looking for the path of valuable data.
        out: file instance to write the results (default terminal)
        digraph: if True return a digraph (http://www.graphviz.org/) representation of the children
        level: internal variable used to call the function recursively.
        ex: execute function
        """
        if not ex is None:
            ex(self)
        if parent == None:
            parent = self.name.decode('utf8')
        if digraph and level == 0:
            out.write('digraph {{\n graph [nodesep=.1 rankdir=LR size="10,120"]\n'.format(root=parent))
        for l in self.getList():
            if l['id'] == 0 or All:
                if digraph:
                    out.write('"{parent}-{name}" [label="{name}"]\n"{parent}" -> "{parent}-{name}"\n'\
                        .format(parent=parent, name=l['name'].decode('utf8')))
                else:
                    if ex is None:
                        out.write("{tab}{name} ({id}) @{bidx}\n".format(tab="\t"*level, **l))
                if level < maxlevel:
                    try:
                        self.gotoItem(l['name'], l['id'])\
                            .show(maxlevel, level+1, All=All, out=out, digraph=digraph\
                            , parent=parent+'-'+l['name'].decode('utf8'), ex=ex)
                    except:
                        pass
        if digraph and level == 0:
            out.write('}')

    def getIndexes(self, key):
        r = []
        for x in self.getList():
            if x['name'] == key.decode('utf8'):
                r.append(x['id'])
        return r
