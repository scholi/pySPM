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
for every object.

the default file extension is pkz (for Pickle Zip)
"""

def save(filename, *objs, **obj):
    """
    save python objects to file
    """
    if os.path.splitext(filename)[1]=='':
        filename += '.pkz'
    out = zipfile.ZipFile(filename, 'a', zipfile.ZIP_DEFLATED)
    file_list = out.namelist()
    update = []
    for i,o in enumerate(objs):
        obj[i] = o
    for k in obj:
        if k in file_list:
            update.append(k)
    if len(update) == 0 :
        for k in obj:
            out.writestr(k, pickle.dumps(obj[k], pickle.HIGHEST_PROTOCOL))
    else:
        out.close()
        _, temp = tempfile.mkstemp()
        shutil.copy(filename, temp)
        out = zipfile.ZipFile(filename, 'w', zipfile.ZIP_DEFLATED)
        old = zipfile.ZipFile(temp, 'r')
        for k in [x for x in file_list if x not in update]:
            out.writestr(k, old.read(k))
        os.remove(temp)
        for k in update:
            out.writestr(k, pickle.dumps(obj[k], pickle.HIGHEST_PROTOCOL))
    out.close()
    
def load(filename, key='0'):
    """
    load python objects from saved pkz files
    """
    if os.path.splitext(filename)[1]=='' and not os.path.exists(filename):
        filename += '.pkz'
    assert os.path.exists(filename)
    if type(key) is int:
        key = str(key)
    f = zipfile.ZipFile(filename, 'r')
    return pickle.loads(f.read(key))
    
class loader:
    """
    This class act as a dictionary and read default value for compressed pkz files.
    Values can be set to new values without changing the content of the saved data.
    Only retrieved data are kept in memory and not all the data.
    """
    def __init__(self, filename):
        if os.path.splitext(filename)[1]=='' and not os.path.exists(filename):
            filename += '.pkz'
        assert os.path.exists(filename)
        self.f = zipfile.ZipFile(filename, 'r')
        self.local = {}
    
    def __iter__(self):
        self.keys = set(self.f.namelist()+list(self.local.keys()))
        return self
        
    def __next__(self):
        if len(self.keys) == 0:
            raise StopIteration()
        else:
            return self.keys.pop()
        
    def __getitem__(self, key):
        if not key in self.local:
            self.local[key] = pickle.loads(self.f.read(key))
        return self.local[key]
    
    def __setitem__(self, key, value):
        self.local[key] = value
        
    def __delitem__(self, key):
        delf.local.delitem(key)