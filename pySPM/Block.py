"""
Module to handle block type used in iontof file formats ITA,ITM,ITS, etc...
"""

import sys
import binascii
import struct
import numpy as np
import os

class MissingBlock(Exception):
    def __init__(self, parent, name, index):
        self.block_name = parent.parent+'/'+name
        self.index = index
        
    def __str__(self):
        return "Missing block \"{name}\" with index {index}".format(name=self.block_name,index=self.index)
        
class Block:
    """
    Class to handle a iontof-Block
    One iontof file ITA,ITM,ITS contains a lot of Blocks forming a hierarchical structure.
    Each Block can have children (sub-Blocks) and values (data).
    
    Note: This class was created by reverse engineering on the fileformat of iontof and is most probably not 100% accurate. Nevertheless is works perfercly with our data upto now.
    """
    def __init__(self, fp, parent=''):
        """
        Init the class
        fp: file pointer (the one created by open(...) of an ITA,ITM,ITS, etc... file pointing at the beginning of a block
        
        Each block start with one byte of type followed by 4 bytes that should always be \x19\x00\x00\x00 (all those 5 bytes are saved in self.Type)
        Note: the value \x19\x00\x00\x00 is the unit32 for 25 which is the pre-header length of the block.
        
        Then follows 5 uint32: length, z, u ,x ,y
            length: The length of the block's name
            z: Block ID. Start at 0 and is increased monotonically for each blocks of the same name with the same parent. We usually find the ID from the children's list (see below) and this information is never used as it's redundant.
            u: The number of children / sub-blocks. Might be = 0 even if the block has children. Check the value L (defined below) if so
            x: The length of the block's value
            y: Redundant. Seems to be always = x
        Then follow length-bytes representing the name of the block
        Then follow x-bytes forming the value of the block
        
        Blocks of types \x01\x19\x00\x00\x00 and \x03\x19\x00\x00\x00 are blocks that contains sub-blocks. There is no big difference between the two. I guess that types \x01 is the first one and type \x03 are the continuation blocks
            Those block have a value which starts with 41-bytes.
                2 uint32 -> (length, nums).
                    length: We actually don't need it. It's a redundant information. That is the length of the sub-headers. (It stop just before the sub-blocks names)
                    nums: The variable u (see above) contains the number of children. If u ==0, then nums will tell the correct number of children
                5 bytes: type (usually 00 00 00 00 00 or 03 19 00 00 00)
                5 uint32 -> a,b,L,d,e
                    a,b,d,e are unknown
                    L seems to give information on the number of children
                1 uint64 -> NextBlock
                    Big blocks can be chunked in several ones. NextBlock tells the position in the file of the next chunk. If = 0, then it's the last chunk
            Then 33 bytes for each sub-block follows:
                1 byte: spacing (usually = 0 or 1)
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
        self.parent = parent
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
        
    def gotoNextBlock(self):
        offset = self.offset
        self.f.seek(offset)
        head = dict(zip(['name_length', 'ID', 'N', 'length1', 'length2'], struct.unpack('<5x5I', self.f.read(25))))
        name = self.f.read(head['name_length'])
        length, nums, NextBlock = struct.unpack('<II25xQ', self.f.read(41))
        self.f.seek(NextBlock)
        return Block(self.f, parent=self.parent)
            
    def createList(self, limit=None, debug=False):
        """
        Generate a list (self.List) containing all the children (sub-blocks) of the current Block
        """
        length, nums, NextBlock = struct.unpack('<II25xQ', self.value[:41])
        self.nums = nums
        offset = self.offset
        self.List = []
        while True:
            self.f.seek(offset)
            head = dict(zip(['name_length', 'ID', 'N', 'length1', 'length2'], \
                struct.unpack('<5x5I', self.f.read(25))))
            name = self.f.read(head['name_length'])
            data = self.f.tell()
            length, nums, NextBlock = \
                struct.unpack('<II25xQ', self.f.read(41))
            N = head['N']
            ## The following is commented as believed to be eronous
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
            child = Block(self.f, parent=[self.parent+'/'+self.name,'/'][self.parent==''])
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
        print('List of', len(self.getList()))
        for i, l in enumerate(self.List):
            self.f.seek(l['bidx'])
            other = ''
            try:
                child = Block(self.f, parent=[self.parent+'/'+self.name,'/'][self.parent==''])
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

    def gotoItem(self, name, idx=0, lazy=False):
        """
        Return a new Block instance of a child of the current Block
        name: name of the children's block
        """
        Idx = self.getIndex(name, idx, lazy=lazy)
        self.f.seek(Idx)
        return Block(self.f, parent=[self.parent+'/'+self.name,'/'][self.parent==''])

    def getIndex(self, name, idx=0, lazy=False):
        """
        Get the index of the children having a name=name.
        This function is more intended for internal usage.
        You are encouraged to use the function _goto_ instead
        
        If more than one children have the same name,
        the second one can by retrieved by idx=1, the third with idx=2, etc.
        
        Sometimes the id does not start with 0, but with random high values. Instead of looking at the correct id, you can sue lazy=True with idx=0 in order to fetch the first one saved.
        """
        if type(name) is bytes:
            name = name.decode('ascii')
        i=0
        for l in self.getList():
            if l['name'] == name:
                if (lazy and i==idx) or (not lazy and l['id'] == idx):
                    return l['bidx']
                i+=1
        raise MissingBlock(self,name,idx)

    def goto(self, path, lazy=False):
        """
        Return a sub Block having a specific path
        path: path is similar to filepath.
        The block X conatained in B which is itself contained in A,
        can be retrieved with the path: A/B/X
        if the block B has several children having the same name,
        A/B/X[n] will return the n-th child (note that 0 is the first child)$
        
        As the id sometimes start at weird values and we just want the first ones saved whatever its id is, we can use the lazy=True argument.
        """
        if path == '':
            return self
        s = self
        for p in path.split('/'):
            idx = 0
            if '[' in p and p[-1] == ']':
                i = p.index('[')
                idx = int(p[i+1:-1])
                p = p[:i]
            s = s.gotoItem(p, idx, lazy=lazy)
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
    
    def unpickle(self):
        import pickle
        return pickle.loads(self.value)
        
    def getKeyValue(self, offset=16):
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
            parent = self.name
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

    def getIndexes(self, key, debug=False):
        if type(key) is str:
            key = key.encode('utf8')
        r = []
        for x in self.getList():
            if debug:
                print(x['name'],key)
            if x['name'] == key:
                r.append(x['id'])
        return r
        
    def modify_block_and_export(self, path, new_data, output, debug=False, prog=False, lazy=False):
        assert not os.path.exists(output) # Avoid to erase an existing file. Erase it outside the library if needed.
        out = open(output,'wb')
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
            debug_msg.append('Block Name: "{}" / length: {}'.format(current.name.decode('utf8'), curr_block_length))
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
