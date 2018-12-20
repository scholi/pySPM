# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

import pickle
import zipfile
import os
import shutil
import tempfile

"""
This small utility can save and load python objects to a single file.
It uses ZIP DEFLATE compression method in order to keep the file small.
It works in a similar way as numpy.savez_compressed and numpy.load, but
just work with any python objects and thus do not generate a numpy.array
for every object. Several objects can be added to the same file with either one call of save or several.
The new data will just be appended to your file. It also support the update of the file (in this case the data file is copied in a new one).

the default file extension is pkz (for Pickle Zip)
"""

data_path = '.'

def set_datapath(path):
    global data_path
    data_path = path

def findPKZ(filename):
    if os.path.splitext(filename)[1]=='':
        filename += '.pkz'
    p2 = os.path.join(data_path, filename)
    if os.path.exists(p2):
        return p2
    if os.path.exists(filename):
        return filename
    if not os.path.exists(filename):
        raise IOError("File \"{}\" not found".format(filename))
    return filename
        
def inarxiv(filename, obj)        :
    if os.path.splitext(filename)[1]=='':
        filename += '.pkz'
    filename = os.path.join(data_path, filename)
    out = zipfile.ZipFile(filename, 'a', zipfile.ZIP_DEFLATED)
    file_list = out.namelist()
    return obj in file_list
    
def save(filename, *objs, **obj):
    """
    save python objects to file
    """
    if os.path.splitext(filename)[1]=='':
        filename += '.pkz'
    filename = os.path.join(data_path, filename)
    out = zipfile.ZipFile(filename, 'a', zipfile.ZIP_DEFLATED)
    file_list = out.namelist()
    update = []
    for i,o in enumerate(objs):
        obj[i] = o
        
    # List objects which are already present in the file (so which will be updated)
    for k in obj:
        if k in file_list:
            update.append(k)
            
    if len(update) == 0 : # No files are updated, just append the new data to the file
        for k in obj:
            out.writestr(k, pickle.dumps(obj[k], pickle.HIGHEST_PROTOCOL))
    else:
        # close the current file, copy it to a temp file and reopen the original file in write mode (will thus override all data)
        out.close()
        ft, temp = tempfile.mkstemp()
        shutil.copy(filename, temp)
        out = zipfile.ZipFile(filename , 'w', zipfile.ZIP_DEFLATED)
        
        # Copy all old data which are not updated and clean tempfile
        old = zipfile.ZipFile(temp, 'r')
        for k in [x for x in file_list if x not in update]:
            out.writestr(k, old.read(k))
        old.close()
        os.fdopen(ft).close()
        os.remove(temp)
        
        # Write all new objects
        for k in obj:
            out.writestr(k, pickle.dumps(obj[k], pickle.HIGHEST_PROTOCOL))
    out.close()
    
def load(filename, *keys):
    """
    load python objects from saved pkz files
    """
    filename = findPKZ(filename)
    if len(keys) == 0:
        keys = ['0']
    elif len(keys) == 1 and ',' in keys[0]:
        keys = keys[0].split(',')
    f = zipfile.ZipFile(filename, 'r')
    
    res = []
    for key in keys:
        if type(key) is int:
            key = str(key)
        try: 
            raw = f.read(key)
        except Exception as e:
            raise KeyError("There is no key {} found in the data {}".format(key, filename))
        try:
            res.append(pickle.loads(raw))
        except:
            raise Exception("Cannot pickle recorded data for key {}".format(key))
        
    f.close()
    if len(keys)==1:
        return res[0]
    return res
    
class loader:
    """
    This class act as a dictionary and read default value for compressed pkz files.
    Values can be set to new values without changing the content of the saved data.
    Only retrieved data are kept in memory and not all the data.
    """
    def __init__(self, filename):
        self.filename = findPKZ(filename)
        self.local = {}
    
    def __iter__(self):
        f = zipfile.ZipFile(self.filename, 'r')
        self.keys = set(f.namelist()+list(self.local.keys()))
        f.close()
        return self
        
    def __next__(self):
        if len(self.keys) == 0:
            raise StopIteration()
        else:
            return self.keys.pop()
        
    def __getitem__(self, key):
        if not key in self.local:
            f = zipfile.ZipFile(self.filename, 'r')
            self.local[key] = pickle.loads(f.read(key))
            f.close()
        return self.local[key]
    
    def __setitem__(self, key, value):
        self.local[key] = value
        
    def __delitem__(self, key):
        delf.local.delitem(key)
        
class BidirData:
    """
    This class act as a dictionary and read default value for compressed pkz files.
    Values can be set to new values without changing the content of the saved data.
    Only retrieved data are kept in memory and not all the data.
    """
    def __init__(self, filename):
        try:
            self.filename = findPKZ(filename)
        except IOError:
            if os.path.splitext(filename)[1]=='':
                filename += '.pkz'
            self.filename = os.path.join(data_path, filename)
            out = zipfile.ZipFile(filename, 'a', zipfile.ZIP_DEFLATED)
            out.close()
        self.local = {}
    
    def __iter__(self):
        f = zipfile.ZipFile(self.filename, 'r')
        self.keys = set(f.namelist()+list(self.local.keys()))
        f.close()
        return self
        
    def __next__(self):
        if len(self.keys) == 0:
            raise StopIteration()
        else:
            return self.keys.pop()
        
    def __getitem__(self, key):
        if not key in self.local:
            f = zipfile.ZipFile(self.filename, 'r')
            self.local[key] = pickle.loads(f.read(key))
            f.close()
        return self.local[key]
    
    def __setitem__(self, key, value):
        self.local[key] = value
        save(self.filename, **{key:value})
        
    def __delitem__(self, key):
        delf.local.delitem(key)
        
    def keys(self):
        f = zipfile.ZipFile(self.filename, 'r')
        keys = f.filelist
        f.close()
        return [x.filename for x in keys]
        