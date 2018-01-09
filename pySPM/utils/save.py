import pickle
import zipfile

def save(filename, *objs, **obj):
    out = zipfile.ZipFile(filename, 'w', zipfile.ZIP_DEFLATED)
    for i,o in enumerate(objs):
        obj[i] = o
    for k in obj:
        out.writestr(k, pickle.dumps(obj[k] , pickle.HIGHEST_PROTOCOL))
    out.close()
    
def load(filename, key='0'):
    if type(key) is int:
        key = str(key)
    f = zipfile.ZipFile(filename, 'r')
    return pickle.loads(f.read(key))