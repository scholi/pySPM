# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Module to handle block type used in iontof file formats ITA,ITM,ITS, etc...
"""

from __future__ import absolute_import
import sys
import binascii
import struct
import os
from .utils import dec_debug, do_debug
from .utils.misc import deprecated, aliased, alias

class MissingBlock(Exception):
    def __init__(self, parent, name, index):
        self.block_name = parent.path+parent.name+'/'+name
        self.index = index
        
    def __str__(self):
        return "Missing block \"{name}\" with index {index}".format(name=self.block_name,index=self.index)
        
@aliased
class Block:
    """
    Class to handle a iontof-Block
    One iontof file ITA,ITM,ITS contains a lot of Blocks forming a hierarchical structure.
    Each Block can have children (sub-Blocks) and values (data).
    
    Note: This class was created by reverse engineering on the fileformat of iontof and is most probably not 100% accurate.
    Nevertheless is works in very good agreement with the developer's data.
    """
    def __init__(self, fp, parent=None):
        """
        Init the class
        fp: file pointer (the one created by open(...) of an ITA,ITM,ITS, etc... file pointing at the beginning of a block
        
        Each block start with one byte of type followed by 4 bytes that should always be \x19\x00\x00\x00 (all those 5 bytes are saved in self.Type)
        Note: the value \x19\x00\x00\x00 is the unit32 for 25 which is the pre-header length of the block.
        
        Then follows 5 uint32: slen, ID, N ,length1, length2
            slen: The length of the block's name
            ID: Block ID. Start at 0 and is increased monotonically for each blocks of the same name with the same parent. We usually find the ID from the children's list (see below) and this information is never used as it's redundant.
            N: The number of children / sub-blocks
            length1: The length of the block's value
            length2: Redundant. Seems to be always = x. We suspect that length2 is the length till the next following block in the file. Length1 and length2 might be different in case a long block is modified with a smaller content. Length1 should thus math the new data length and length2 remain as large as the old data.
        Then follow slen-bytes representing the name of the block
        Then follow length1-bytes forming the value of the block
        
        Blocks of types \x01\x19\x00\x00\x00 and \x03\x19\x00\x00\x00 are blocks that contains sub-blocks. There is no big difference between the two. Types \x01 are the first one and type \x03 are the continuation blocks (see NextBlock information below)
            Those block have a value which starts with 41-bytes.
                2 uint32 -> (length, nums).
                    length: We actually don't need it. It's a redundant information. That is the length of the sub-headers. (It stop just before the sub-blocks names)
                    nums: nums is an indication on the number of children the block can contains. It was found that the block size is = 53*nums
                5 bytes: type of the NextBlock (usually 00 00 00 00 00 or 03 19 00 00 00)
                5 uint32 -> a,b,c,d,e
                    a,b,c,d,e are the NextBlock header (a=slen, b=ID, ...)
                1 uint64 -> NextBlock
                    Big blocks can be chunked in several ones. NextBlock tells the position in the file of the next chunk. If = 0, then it's the last chunk
            Then 33 bytes for each sub-block follows:
                1 byte: Type of the child
                3 uint32 -> index, slen, id
                    index: The position of the sub-block name in the header
                    slen: The length of the sub-block name (which is store later). So basically the sub-block name is: Block.value[index:index+slen]
                    id: start at 0 and increase monotonically for sub-blocks having the same name
                1 uint32: 1 if the child contains data (we believe)?
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
        self.parent = parent
        self.path = ''
        if parent is not None:
            if parent.parent is None:
                self.path = '/'
            else:
                self.path = parent.path+parent.name+'/'
        self.offset = self.f.tell()
        self.Type = self.f.read(5)
        if self.Type[1:] != b'\x19\x00\x00\x00':
            raise ValueError('Wrong block type ({Type}) found @{pos}'\
                .format(pos=self.offset, Type=binascii.hexlify(self.Type[1:])))
        if len(self.Type) < 5:
            raise ValueError('EOF reached. Block cannot be read')
        self.head = dict(zip(['name_length', 'ID', 'N', 'length1', 'length2'], \
            struct.unpack('<5I', self.f.read(20))))
        self.name = self.f.read(self.head['name_length']).decode('ascii')
        self.value = self.f.read(self.head['length1'])
        self.List = None
        self.iterP = 0
    
    def inc(self):
        """
        Increment the value (interpreted as unsigned in32) of a block by one and write the changes to the file
        """
        self.rewrite(struct.pack("<I", self.getLong()+1))
        
    def add_child(self, blk):
        """
        Add a new child to a given block. /!\ will overwrite the ITA/ITM file.
        """
        import os
        assert self.Type[0] in [1,3]
        
        length, nums, NextBlock = struct.unpack('<II25xQ', self.value[:41])
        if NextBlock != 0:
            return self.goto_next_block().add_child(blk)
        
        # Check the available size of the block
        header_length = 41+33*self.head['N']
        
        children_names_length = 0
        last_id = None
        lowest_index = struct.unpack("<I", self.value[:4])[0]
        children_names_length = self.head['length1']-lowest_index
        free_space = lowest_index-header_length
        new_name = blk.name not in self
        
        new_subblock = False
        if new_name:
            if free_space < 33+len(blk.name):
                new_subblock = True
            index = lowest_index-len(blk.name)
        else:
            if free_space < 33:
                new_subblock = True
            else:
                index = 0
                for i in range(self.head['N']):
                    self.f.seek(self.offset+25+len(self.name)+41+33*i)
                    s = struct.unpack("<B4I2Q", self.f.read(33))
                    if s[2] == len(blk.name):
                        self.f.seek(self.offset+25+len(self.name)+s[1])
                        if blk.name == self.f.read(s[2]).decode('utf8'):
                            index = s[1]
                            break
                assert index > 0
        if new_subblock:
            raise Exception("Block is too small to fit the data.")
        self.f.seek(self.offset+25+len(self.name)+header_length)
        self.f.write(struct.pack("<B4I2Q", blk.Type[0], index, len(blk.name), blk.head['ID'], [0,1][blk.Type[0] in [0,128]], blk.head['length1'], blk.offset))
        self.f.seek(self.offset+25+len(self.name)+index)
        self.f.write(blk.name.encode('utf8'))
        self.f.seek(self.offset+13)
        self.head['N'] += 1
        self.f.write(struct.pack("<I", self.head['N']))
        self.f.seek(self.offset+25+len(self.name))
        if new_name:
            self.f.write(struct.pack("<I", struct.unpack("<I", self.value[:4])[0]-len(blk.name)))
        
        self.refresh()
        return blk
        
    def refresh(self):
        """
        Reload self.value from the file.
        This function is useful in case a block was overwritten
        """
        self.List = None
        self.f.seek(self.offset+self.head['name_length']+25)
        self.value = self.f.read(self.head['length1'])
        
    def edit_child(self, old_block, new_block, debug=False):
        """
        Edit the children list of a given block.
        """
        if not self.Type[0] in [1,3]:
            raise Exception("The current block should be of type 01 or 03 (folder)")
        
        header_length = 41+33*self.head['N']
        children_names_length = 0
        last_id = None
        lowest_index = struct.unpack("<I", self.value[:4])[0]
        children_names_length = self.head['length1']-lowest_index
        found = False
        for i in range(self.head['N']):
            self.f.seek(self.offset+25+len(self.name)+41+33*i)
            entry = list(struct.unpack("<B4I2Q", self.f.read(33)))
            if entry[6]==old_block.offset:
                found = True
                self.f.seek(self.offset+25+len(self.name)+41+33*i)
                entry[6] = new_block.offset
                entry[5] = new_block.head['length1']
                if do_debug(debug):
                    print("Update inode")
                self.f.write(struct.pack("<B4I2Q", *entry))
                break
        if not found:
            raise Exception('Child {} not found in {}'.format(old_block.name, self.name))
        self.refresh()
        
    def create_dir(self, name, children=[], nums=100, assign=True, id=0):
        """
        Create a directory in the ITStr format.
        
        Parameters
        ----------
        name: string
            name of the directory
        children: list of blocks
            list of the childrens' block which are added
        nums: int
            this number will give the approximate number of children one can create. The size of the block will be 53*nums
        assign: bool
            If True, will add the created dir as a child of self
        id: int
            Blocks having the same name should have a different id.
        """
        size = 53*nums
        value = struct.pack("<2IB6IQ{}x".format(size-41), size, nums,  0,  0, 0, 0, 0, 0, 0, 0)
        blk = self.create_block(name, value, _type=1, id=id)
        for x in children:
            blk.add_child(x)
        if assign:
            self.add_child(blk)
        return blk
    
    def edit_block(self, path, name, value, id=0, _type=0, force=False, debug=False):
        """
        This function will go to a given path and create all necessary "folder" (block of type 1).
        It will then either create a new child with a given name and value if it does not exists and if so it will edit it.
        Please note that for safety the new value should be of the exact same length than the existing one.
        If this is not the case, you should consider using the modify_block_and_export() method.
        You can also use the force=True parameter to force to edit a block when its content is not the same. Be careful, because this will actually keep the old data in the file (but won't be accessible or seen anymore). This means that the size of your ita can grow quickly if you perform a lot of edits...
        """
        self.f.seek(self.offset)
        parent = self
        if path != '':
            for p in path.split('/'):
                idx = 0
                if '[' in p and p[-1] == ']':
                    i = p.index('[')
                    idx = int(p[i+1:-1])
                    p = p[:i]
                if p is '*':
                    e = parent.get_list()[idx]
                    p = e['name']
                    idx = e['id']
                try:
                    parent = parent.goto_item(p, idx)
                except MissingBlock:
                    parent = parent.create_dir(p, children=[], id=idx)
        try:
            if do_debug(debug):
                print("Accessing block \"{}\"[{}]".format(name, id))
            child = parent.goto_item(name, id)
            if child.head['length1'] == len(value):
                if do_debug(debug):
                    print("rewrite block")
                child.rewrite(value, debug=dec_debug(debug))
                return child
            elif force:
                blk = parent.create_block(name, value, id=id, _type=_type, debug=dec_debug(debug))
                parent.edit_child(child, blk, debug=dec_debug(debug))
            else:
                raise Exception("Use the force=True parameter if you wish to replace an existing block with another data size")
        except MissingBlock:
            if do_debug(debug):
                print("create new block")
            return parent.add_child(parent.create_block(name, value, id=id, _type=_type))
    
    def create_block(self, name, value, id=0, _type=0, debug=False):
        """
        Create a new block and write it at the end of the file
        """
        if do_debug(debug):
            print("Creating new block \"{}\" of size {}".format(name, len(value)))
        if type(name) is str:
            name = name.encode('utf8')
        self.f.seek(0, 2) # goto end of file
        offset = self.f.tell()
        size = len(value)
        slen = len(name)
        self.f.write(struct.pack("<B6I", _type, 25, slen, id, 0, size, size))
        self.f.write(name)
        self.f.write(value)
        self.f.seek(offset)
        return Block(self.f)
        
    @deprecated("DepthFirstSearch")
    def depth_first_search(self, callback=None, filter=lambda x: True, func=lambda x: x):
        """
        Perform a depth first search on the blocks.
        pass each block where filter(block) is true to a callback
        Also apply a function func to the result and create a list for the result
        """
        res = []
        if filter(self):
            res += [func(self)]
            if callback is not None:
                callback(self)
        if self.Type[0] in [1,3]:
            for x in self:
                res += x.depth_first_search(callback=callback, filter=filter, func=func)
        return res
    
    @deprecated("getName")
    def get_name(self):
        """
        Return the name of the Block
        """
        return self.name

    @deprecated("getList")
    def get_list(self):
        """
        Return the list of the sub-blocks (children) of the current Block.
        """
        if not self.Type[0:1] in [b'\x01', b'\x03']:
            return []
        if self.List is None:
            return self.create_list()
        return self.List
    
    @deprecated("gotoFollowingBlock")
    def goto_following_block(self):
        offset = self.offset+25+self.head['name_length']+self.head['length1']
        if offset < os.fstat(self.f.fileno()).st_size:
            self.f.seek(offset)
            return Block(self.f, parent=None)
        return None

    #deprecated("gotoNextBlock")
    def goto_next_block(self):
        offset = self.offset
        self.f.seek(offset)
        head = dict(zip(['name_length', 'ID', 'N', 'length1', 'length2'], struct.unpack('<5x5I', self.f.read(25))))
        name = self.f.read(head['name_length'])
        length, nums, NextBlock = struct.unpack('<II25xQ', self.f.read(41))
        if NextBlock==0:
            return None
        self.f.seek(NextBlock)
        return Block(self.f, parent=self.parent)
    
    def getNthChild(self, n=0):
        L = self.get_list()
        assert len(L)>n
        return self.goto_item(L[n]['name'],L[n]['id'])

    @deprecated("createList")
    def create_list(self, limit=None, debug=False):
        """
        Generate a list (self.List) containing all the children (sub-blocks) of the current Block
        """
        length, nums, next_block = struct.unpack('<II25xQ', self.value[:41])
        self.nums = nums
        offset = self.offset
        self.List = []
        while True:
            self.f.seek(offset)
            data = self.f.read(25)
            if len(data)<25:
                raise Exception('Children of {} @{} cannot be read. Data might be corrupted!'.format(self.name, offset))
            head = dict(zip(['name_length', 'ID', 'N', 'length1', 'length2'], \
                struct.unpack('<5x5I', data)))
            name = self.f.read(head['name_length'])
            data = self.f.tell()
            length, nums, next_block = \
                struct.unpack('<II25xQ', self.f.read(41))
            N = head['N']
            ## The following is commented as believed to be erroneous
            #if N == 0:
            #    N = nums
            for i in range(N):                
                self.f.seek(data+42+33*i)
                S = dict(\
                    zip(['index', 'slen', 'id', 'blen', 'bidx'],\
                    struct.unpack('<III4xQQ', self.f.read(32))))
                self.f.seek(data+S['index'])
                S['name'] = self.f.read(S['slen']).decode('ascii')
                self.List.append(S)
            if next_block == 0:
                break
            offset = next_block
        return self.List
               
    @deprecated("getString")
    def get_string(self):
        """
        Decode the value of the Block to UTF-16 (standard for all strings in this fileformat)
        """
        return self.value.decode('utf16')

    @deprecated("dictList")
    def dict_list(self):
        """
        Return a dictionary of the value decoded with various formats (raw, long, float, utf16)
        As the type of the data is not known, this function is very helpful for debugging purpose
        """
        d = {}
        for i, l in enumerate(self.get_list()):
            self.f.seek(l['bidx'])
            child = Block(self.f, parent=self)
            if child.Type[0:1] == b'\x00':
                value = binascii.hexlify(child.value)
                d[child.name] = {'raw':value}
                if len(child.value) == 4:
                    d[child.name]['long'] = child.get_long()
                    d[child.name]['ulong'] = child.get_ulong()
                elif len(child.value) == 8:
                    d[child.name]['float'] = child.get_double()
                    d[child.name]['long'] = child.get_longlong()
                if len(child.value)%2 == 0:
                    d[child.name]['utf16'] = child.value.decode('utf16', "ignore")
            del child
        return d
    
    @deprecated("showList")
    def show_list(self):
        """
        Show a list of all the children (sub-blocks) of the current Block.
        It will also display the value/data of all the children (if any)
        """
        print('List of', len(self.get_list()))
        for i, l in enumerate(self.List):
            self.f.seek(l['bidx'])
            other = ''
            try:
                child = Block(self.f, parent=self)
                if child.Type[0:1] == b'\x00':
                    if len(child.value) == 4:
                        vL = child.get_long()
                        Dtype = 'long'
                    elif len(child.value) == 8:
                        vL = child.get_double()
                        Dtype = 'double'
                        other += ' = '+str(child.get_longlong())+" (long64)"
                    elif len(child.value) == 2:
                        vL = child.get_short()
                        Dtype = 'short'
                    elif len(child.value) == 1:
                        vL = child.get_byte()
                        Dtype = 'byte'
                    else:
                        vL = '???'
                        Dtype = '???'
                    
                    value = binascii.hexlify(child.value)
                    if len(value) > 16:
                        value = value[:16]+b'...'
                    if len(child.value)%2 == 0:
                        vS = child.value.decode('utf16', "ignore")
                        if len(vS) > 20:
                            vS = vS[:20]+'...'
                        print(u"{name} ({id}) <{blen}> @{bidx}, value = {value} (hex) = \"{vS}\" (UTF-16)= {vL} ({Dtype}){other}"\
                            .format(value=value, vL=vL, Dtype=Dtype, other=other, vS=vS, **l))
                    else:
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
        L = self.get_list()
        if self.pointer >= len(L):
            raise StopIteration
        it = L[self.pointer]
        self.pointer += 1
        return self.goto_item(it['name'], it['id'])

    @deprecated("gotoItem")
    def goto_item(self, name, idx=0, lazy=False):
        """
        Return a new Block instance of a child of the current Block
        name: name of the children's block
        """
        Idx = self.get_index(name, idx, lazy=lazy)
        self.f.seek(Idx)
        return Block(self.f, parent=self)

    @deprecated("getIndex")
    def get_index(self, name, idx=0, lazy=False):
        """
        Get the index (here to understand the offset in bytes from the beginning of the file till the block of interest) of the children having a name=name.
        This function is more intended for internal usage.
        You are encouraged to use the function _goto_ instead
        
        If more than one children have the same name,
        the second one can by retrieved by idx=1, the third with idx=2, etc.
        
        Sometimes the id does not start with 0, but with random high values.
        Instead of looking at the correct id, you can use lazy=True with idx=0 in order to fetch the first one saved.
        """
        if type(name) is bytes:
            name = name.decode('ascii')
        i = 0
        if name is '*':
            return self.get_list()[idx]['bidx']
        for l in self.get_list():
            if l['name'] == name:
                if (lazy and i==idx) or (not lazy and l['id'] == idx):
                    return l['bidx']
                i+=1
        raise MissingBlock(self, name, idx)
       
    def goto(self, path, lazy=False):
        """
        Return a sub Block having a specific path
        path: path is similar to filepath.
        The block X contained in B which is itself contained in A,
        can be retrieved with the path: A/B/X
        if the block B has several children having the same name,
        A/B/X[n] will return the n-th child (note that 0 is the first child)$
        
        As the id sometimes start at weird values and we just want the first ones
        saved whatever its id is, we can use the lazy=True argument.
        """
        if path.startswith('/'):
            path = path[1:]
        if path == '':
            return self
        self.f.seek(self.offset)
        s = Block(self.f)
        for p in path.split('/'):
            idx = 0
            if '[' in p and p[-1] == ']':
                i = p.index('[')
                idx = int(p[i+1:-1])
                p = p[:i]
            if p is '*':
                e = s.get_list()[idx]
                p = e['name']
                idx = e['id']
            s = s.goto_item(p, idx, lazy=lazy)
        return s

    @deprecated("getLongLong")
    def get_longlong(self):
        """
        Decode the value as an 64-Integer
        """
        return struct.unpack('<q', self.value)[0]

    @deprecated("getDouble")
    def get_double(self):
        """
        Decode the value as a 64-float (Double)
        """
        return struct.unpack('<d', self.value)[0]

    @deprecated("getShort")
    def get_short(self):
        """
        Decode the value as an 16-Integer (Short)
        """
        return struct.unpack('<h', self.value)[0]

    @deprecated("getByte")
    def get_byte(self):
        """
        Decode the value as a 1-Byte
        """
        return struct.unpack('<B', self.value)[0]
        
    def get_bytes(self):
        """
        Decode the values as unsigned bytes
        """
        return struct.unpack('<{}B'.format(len(self.value)), self.value)

    @deprecated("getULong")
    def get_ulong(self):
        """
        Decode the value as an unsigned 32-Integer (Long)
        """
        return struct.unpack('<I', self.value)[0]

    @deprecated("getLong")
    def get_long(self):
        """
        Decode the value as an 32-Integer (Long)
        """
        return struct.unpack('<i', self.value)[0]
    
    def unpickle(self):
        import pickle
        return pickle.loads(self.value)
        
    @deprecated("getKeyValue")
    def get_key_value(self, offset=16):
        """
        Return a dictionnary of key/values pairs of the data
        Note that the function has no idea if the data are stored as so.
        """
        L = struct.unpack("<I", self.value[offset:offset+4])[0]
        Key = self.value[offset+4:offset+4+L].decode('utf16', 'ignore')
        int_value, float_value = struct.unpack("<2xqd", self.value[offset+4+L:offset+22+L])
        L2 = struct.unpack("<I", self.value[offset+22+L:offset+26+L])[0]
        SVal = self.value[offset+26+L:offset+26+L+L2].decode('utf16', 'ignore')
        return {'key':Key, 'float':float_value, 'int':int_value,'string':SVal}

    def __contains__(self, name):
        return name in [x['name'] for x in self.get_list()]
        
    def show(self, maxlevel=3, level=0, all=False, out=sys.stdout, digraph=False, parent=None, ex=None, **kargs):
        """
        Display the children of the current Block (recursively if maxlevel > 1)
        Very useful for debugging purpose and looking for the path of valuable data.
        out: file instance to write the results (default terminal)
        digraph: if True return a digraph (http://www.graphviz.org/) representation of the children
        level: internal variable used to call the function recursively.
        ex: execute function
        """
        if 'All' in kargs:
            all = kargs.pop("All")
            
        if not ex is None:
            ex(self)
        if parent == None:
            parent = self.name
        if digraph and level == 0:
            out.write('digraph {{\n graph [nodesep=.1 rankdir=LR size="10,120"]\n'.format(root=parent))
        for l in self.get_list():
            if l['id'] == 0 or all:
                if digraph:
                    out.write('"{parent}-{name}" [label="{name}"]\n"{parent}" -> "{parent}-{name}"\n'\
                        .format(parent=parent, name=l['name'].decode('utf8')))
                else:
                    if ex is None:
                        out.write("{tab}{name} ({id}) @{bidx}\n".format(tab="\t"*level, **l))
                if level < maxlevel:
                    try:
                        self.goto_item(l['name'], l['id'])\
                            .show(maxlevel, level+1, all=all, out=out, digraph=digraph\
                            , parent=parent+'-'+l['name'], ex=ex)
                    except:
                        pass
        if digraph and level == 0:
            out.write('}')

    @deprecated("getIndexes")
    def get_indexes(self, key, debug=False):
        if type(key) is str:
            key = key.encode('utf8')
        r = []
        for x in self.get_list():
            if debug:
                print(x['name'],key)
            if x['name'] == key:
                r.append(x['id'])
        return r

    def decompress(self):
        import zlib
        return zlib.decompress(self.value)

    @deprecated("getData")
    def get_data(self, fmt="I", decompress=True):
        if decompress:
            raw = self.decompress()
        else:
            raw = self.value
        L = len(raw)//int(struct.calcsize(fmt))
        return struct.unpack("<"+str(L)+fmt, raw)
    
    def rewrite(self, content, debug=False):
        assert(len(content)) == self.head['length1']
        if do_debug(debug):
            print("Rewriting block \"{}\"".format(self.name))
        # set pointer at beginning of data
        self.f.seek(self.offset+25+self.head['name_length'])
        self.f.write(content)
        self.refresh()
            
    def modify_block_and_export(self, path, new_data, output, debug=False, prog=False, lazy=False):
        assert not os.path.exists(output) # Avoid to erase an existing file. Erase it outside the library if needed.
        out = open(output, 'wb')
        out.write(b'ITStrF01')
        block = self.goto(path, lazy=lazy)
        block_offset = block.offset
        length_diff = len(new_data)-len(block.value)
        self.f.seek(8)
        FILE_SIZE =  os.fstat(self.f.fileno()).st_size
        if prog:
            try:
                from tqdm import tqdm_notebook as tqdm
            except:
                from tqdm import tqdm as tqdm
            T = tqdm(total=FILE_SIZE)
        debug_msg = ["FileSize: "+str(FILE_SIZE)]
        if debug:
            print("File Size",FILE_SIZE)
        curr = 8
        if prog:
            T.update(8)
        while self.f.tell() < FILE_SIZE:
            debug_msg = debug_msg[-30:]
            if prog:
                ncurr = self.f.tell()
                T.update(ncurr-curr)
                curr = ncurr
            debug_msg.append("Current position: @"+str(self.f.tell()))
            try:
                current = Block(self.f) # Here we don't care about the parent argument. It is used only for debug purpose anyway.
            except Exception as ex:
                print("Error found! Debug info")
                for x in debug_msg:
                    print("\t"+x)
                raise ex
            self.f.seek(current.offset)
            curr_block_length = current.head['length1'] + current.head['name_length'] + 25
            debug_msg.append('Block Name: "{}" / length: {}'.format(current.name, curr_block_length))
            if current.offset == block_offset: # Found the block to change
                debug_msg.append("Block to change FOUND!")
                out.write(self.f.read(5)) # Write block type
                out.write(struct.pack("<5I",block.head['name_length'],block.head['ID'],block.head['N'], length_diff+block.head['length1'], length_diff+block.head['length2']))
                self.f.read(20) # Skip header
                out.write(self.f.read(block.head['name_length'])) # copy block name
                self.f.read(block.head['length1']) # skip data in source
                out.write(new_data) # write new_data
            elif current.Type[0] in [1,3]: # found a container, check for references to block after the modified block
                debug_msg.append("Block container found. Checking children...")
                out.write(self.f.read(25)) # copy header
                out.write(self.f.read(current.head['name_length'])) # copy block name
                SubHeader = list(struct.unpack('<2I5s5IQ', self.f.read(41) )) # read sub-block header
                if SubHeader[8] > block_offset: # Is the nextbloxk after the modified block? Yes => Adjust the offset position
                    SubHeader[8] += length_diff
                out.write(struct.pack('<2I5s5IQ', *SubHeader )) # write sub-block header
                N = current.head['N']
                #if N == 0:
                #    N = SubHeader[1]
                for i in range(N):
                    X, index, slen, id,Y, blen, bidx = struct.unpack('<B4I2Q', self.f.read(33))
                    if bidx == block_offset: # If the children block is the modified block, adjust length
                        blen = len(new_data)
                    elif bidx > block_offset: # If the children is after the modifien block, adjust its offset
                        bidx += length_diff
                    out.write(struct.pack('<B4I2Q',X, index, slen, id, Y,blen, bidx)) # write child info
                # Write the extra bytes used by iontof which seems to be useless as well as the childrens' name
                delta = curr_block_length - (self.f.tell() - current.offset) # number of bytes remaining till the end of the block
                out.write(self.f.read(delta))
            else:
                debug_msg.append("Data Block found. Copy data without check...")
                out.write(self.f.read(curr_block_length))
        if prog:
            T.update(FILE_SIZE-curr)
            T.close()
        out.close()
