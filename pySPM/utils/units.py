import sys
import math
import numpy as np
from . import constants as const

_SI_units = ['kg','m','s','A','K','cd','mol']

_units = {
    'V':{'kg':1,'m':2,'s':-3,'A':-1,'K':0,'cd':0,'mol':0},
    'C':{'kg':0,'m':0,'s':1,'A':1,'K':0,'cd':0,'mol':0},
    'N':{'kg':1,'m':1,'s':-2,'A':0,'K':0,'cd':0,'mol':0},
    'J':{'kg':1,'m':2,'s':-2,'A':0,'K':0,'cd':0,'mol':0},
    'W':{'kg':1,'m':2,'s':-3,'A':0,'K':0,'cd':0,'mol':0},
    'Pa':{'kg':1,'m':-1,'s':-2,'A':0,'K':0,'cd':0,'mol':0},
    'Ω':{'kg':1,'m':2,'s':-3,'A':-2,'K':0,'cd':0,'mol':0}
    }

_units_scale_conversion = {
    'eV': (const.qe,'V'),
    'bar':(1e5,'Pa'),
    'atm':(101325,'Pa')
    }

_unit_shift_conversion = {
    '°C':(273.15, 'K'),
    }
   
_unit_scale = {
    'y':-24,'z':-21,'a':-18,'f':-15,'p':-12,'n':-9,'u':-6,'m':-3,'c':-2,'d':-1,
    'h':2,'k':3,'M':6,'G':9,'T':12,'P':15,'E':18,'Z':21,'Y':24,'':0}
    
_unit_revert_scale = {}
import re

_SI_units = ['kg','m','s','A','K','cd','mol']

_units = {
    'V':{'kg':1,'m':2,'s':-3,'A':-1},
    'C':{'s':1,'A':1},
    'N':{'kg':1,'m':1,'s':-2},
    'J':{'kg':1,'m':2,'s':-2},
    'W':{'kg':1,'m':2,'s':-3},
    'Pa':{'kg':1,'m':-1,'s':-2},
    }

_unit_scale = {
    'y':-24,'z':-21,'a':-18,'f':-15,'p':-12,'n':-9,'u':-6,'m':-3,'c':-2,'d':-1,
    'h':2,'k':3,'M':6,'G':9,'T':12,'P':15,'E':18,'Z':21,'Y':24,'':0}

def par_parse(s):
    parsed = []
    count = 0
    opening = None
    closing = 0
    for i,x in enumerate(s):
        if x is '(':
            if opening is None:
                opening = i
            count += 1
        elif x is ')':
            count -= 1
            if count==0 and opening is not None:
                parsed += [s[closing:opening], par_parse(s[opening+1:i])]
                closing = i+1
                opening = None
    if closing < len(s):
        parsed.append(s[closing:])
    return parsed

def op_parse(s):
    r = []
    for x in s:
        if type(x) is list:
            r.append(op_parse(x))
        else:
            r += [x for x in re.split(r'(\*|/|\^)', x) if len(x)>0]
    return r

def parse(unit):
    sp = par_parse(unit)
    sp = op_parse(sp)
    sp = u_parse(sp)
    sp = op_exec(sp)
    return sp

def num_parse(s):
    pass  

def u_parse(s):
    if type(s) is list:
        sub = [u_parse(y) for y in s]
        return sub
    for x in '*/^':
        if x in s:
            return s
    result = None
    if re.match(r'\-?[0-9]+(\.[0-9]+)?', s):
        result = unit({}, float(s))
    elif s in _SI_units:
        result = unit(s)
    elif s in _units:
        result =  unit(_units[s])
    elif s[0] in _unit_scale:
        if s[1:] in _SI_units:
            result = unit(s[1:], 10**(_unit_scale[s[0]]))
        elif s[1:] in _units:
            result = unit(_units[s[1:]], 10**(_unit_scale[s[0]]))
    elif len(s) == 2 and s[1] == 'g' and x[0] in _unit_scale:
        result = unit('kg',10**(_unit_scale[s[0]]-3))
    elif s == 'g':
        result = unit('kg',1e-3)
    elif s in _units_scale_conversion:
        u = _units_scale_conversion[s]
        result = unit(u[1], u[0])
    return result

def op_exec(s):
    s = [op_exec(x) if type(x) is list else x for x in s]
    while '^' in s:
        i = s.index('^')
        a = s[i-1]
        b = s[i+1]
        s = s[:i-1]+[a**b]+s[i+2:]
    while '/' in s:
        i = s.index('/')
        s = s[:i-1]+[s[i-1]/s[i+1]]+s[i+2:]
    while '*' in s:
        i = s.index('*')
        s = s[:i-1]+[s[i-1]*s[i+1]]+s[i+2:]
    return s[0]  

class unit(object):
    def __init__(self, u, value=1):
        self.value = 1
        if type(u) is str:
            if u in _SI_units:
                self.units = {u:1}
            else:
                p = parse(u)
                self.units = p.units
                self.value = p.value
        elif type(u) is dict:
            self.units = {x: u[x] for x in u}
        else:
            raise TypeError(type(u))
        self.value *= value
    
    def __mul__(self, b):
        value = self.value
        units = self.units
        if isinstance(b, unit):
            for x in b.units:
                units[x] = units.get(x,0)+b.units[x]
            value *= b.value
            return unit(units, value)
        return unit(self.units, self.value*b)
    
    __rmul__ = __mul__
    
    def __div__(self, b):
        value = self.value
        units = self.units
        if isinstance(b, unit):
            for x in b.units:
                units[x] = units.get(x,0)-b.units[x]
            value /= b.value
            return unit(units, value)
        return unit(self.units, self.value/b)
       
    def __rdiv__(self, b):
        value = 1/self.value
        units = {x: -self.units[x] for x in self.units}
        if isinstance(b, unit):
            for x in b.units:
                units[x] = units.get(x,0)+b.units[x]
            value *= b.value
            return unit(units, value)
        return unit(units, b/self.value)
       
    def __pow__(self, n):
        if isinstance(n, unit):
            assert n.units == {}
            return unit({x: n.value*self.units[x] for x in self.units}, self.value**n.value)    
        return unit({x: n*x for x in self.units}, self.value**n)
    
    __truediv__ = __div__
    __rtruediv__ = __rdiv__
    
    def __repr__(self):
        if self.units == {}:
            return str(self.value)
        u = '*'.join(['{}^{}'.format(x, self.units[x]) if self.units[x] is not 1 else x for x in self.units if self.units[x]!=0])
        if self.value == 1:
            return u
        return '{:.3e}*'.format(self.value)+u
    
    __str__ = __repr__
    
class SIunit(np.ndarray):
    def __new__(cls, input_array, u={}):
        obj = np.asarray(input_array).view(cls)
        obj.unit = unit(u)
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.unit = getattr(obj, 'unit', {})
    
    def __repr__(self):
        return np.ndarray.__repr__(self)[:-1]+', unit='+str(self.unit)+')'
    
    def __mul__(self, b):
        r = np.ndarray.__mul__(self, b)
        if isinstance(b, SIunit):
            r.unit = r.unit * b.unit
        else:
            r.unit = self.unit
        return r

    def __add__(self, b):
        r = np.ndarray.__add__(self, b)
        if isinstance(b, SIunit):
            for x in self.units:
                assert r.units[x] == b.units[x]
            for x in b.units:
                assert r.units[x] == b.units[x]
            return r