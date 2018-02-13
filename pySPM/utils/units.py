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


        
for x in _unit_scale:
    for y in range(_unit_scale[x],_unit_scale[x]+3):
        _unit_revert_scale[y] = x
        
class UnitMissmatch(Exception):
    pass

class InvalidUnit(Exception):
    pass

def isunitless(units):
    for x in units:
        if units[x]!=0:
            return False
    return True

def SIvalues(x):
    if isinstance(x, SIunit):
        return x.value
    if isinstance(x, np.ndarray):
        return np.array([k.value for k in x])
    if hasattr(x, "__getitem__"):
        return [k.value for k in x]

def parse_unit(unit):
    import re
    # TODO

class SIunit(object):
    def __init__(self, value=1, **kargs):        
        self.units = {u:0 for u in _SI_units}
        self.value = value
        for x in kargs:
            if x in _SI_units:
                self.units[x] += kargs[x]
            elif x in _units:
                self.units = {y: self.units[y]+kargs[x]*_units[x][y] for y in self.units}
            elif x[0] in _unit_scale and x[1:] in _SI_units:
                self.value *= 10**(kargs[x]*_unit_scale[x[0]])
                self.units[x[1:]] += kargs[x]
            elif x[0] in _unit_scale and x[1:] in _units:
                self.value *= 10**(kargs[x]*_unit_scale[x[0]])
                self.units = {y: self.units[y]+kargs[x]*_units[x][y] for y in self.units}
            elif x[-1] == 'g' and x[:-1] in _unit_scale:
                self.units['kg'] += kargs[x]
                self.value *= 10**(kargs[x]*(_unit_scale[x[:-1]]-3))
            else:
                raise InvalidUnit
                
    def convert(self, unit):
        if unit in _units:
            return self.value
        if unit in _unit_shift_conversion:
            return SIunit(self.value-_unit_shift_conversion[unit][0], _unit_shift_conversion[unit][1]).convert(unit)
        if unit in _units_scale_conversion:
            return SIunit(self.value/_units_scale_conversion[unit][0], _unit_shift_conversion[unit][1]).convert(unit)
        
    def __add__(self, B):
        if isinstance(B, np.ndarray):
            if not isinstance(B[0], SIunit):
                raise TypeError
            for b in B:
                for x in self.units:
                    if self.units[x] != b.units[x]:
                        raise UnitMissmatch
            return np.array([SIunit(self.value+b.value, self.units) for b in B])

        for x in self.units:
            if self.units[x] != B.units[x]:
                raise UnitMissmatch
        return SIunit(self.value+B.value, **self.units)
        
        
    def __sub__(self, B):
        if isinstance(B, np.ndarray):
            if not isinstance(B[0], SIunit):
                raise TypeError
            for b in B:
                for x in self.units:
                    if self.units[x] != b.units[x]:
                        raise UnitMissmatch
            return np.array([SIunit(self.value-b.value, self.units) for b in B])

        for x in self.units:
            if self.units[x] != B.units[x]:
                raise UnitMissmatch
        return SIunit(self.value-B.value, **self.units)
       
    def __mul__(self, B):
        if isinstance(B, SIunit):
            new_units = {x:self.units[x]+B.units[x] for x in self.units}
            if isunitless(new_units):
                return self.value*B.value
            return SIunit(self.value*B.value, **new_units)
        elif isinstance(B, np.ndarray):
            if isinstance(B[0], SIunit):
                new_units = {x:self.units[x]+B[0].units[x] for x in self.units}
                if isunitless(new_units):
                    np.array([self.value*x.value for x in B])
                return np.array([SIunit(self.value*x.value, **new_units) for x in B])
            else:
                return np.array([SIunit(self.value*x, **self.units) for x in B])            
        return SIunit(self.value*B, **self.units)
        
    __rmul__ = __mul__
    
    
    def __div__(self, B):
        if isinstance(B, SIunit):
            new_units = {x:self.units[x]-B.units[x] for x in self.units}
            if isunitless(new_units):
                return self.value/B.value
            return SIunit(self.value/B.value, **new_units)
        elif isinstance(B, np.ndarray):
            return (1/B)*SIunit(self.value, **self.units)
        elif isinstance(B, (int, float, complex)):
            return SIunit(self.value/B, **self.units)
            
    __truediv__ = __div__
    
    def __rdiv__(self, B):
        if isinstance(B, np.ndarray):
            return B*SIunit(1/self.value, **{x:-self.units[x] for x in self.units})
        
    def __pow__(self, N):
        return SIunit(self.value**N, **{x:self.units[x]*N for x in self.units})
    
    def get_unit(self):
        if self.units in _units.values():
            return list(_units.keys())[list(_units.values()).index(self.units)]
        return '*'.join([(['{unit}','{unit}^{N}'][self.units[x]!=1]).format(unit=x,N=self.units[x]) for x in self.units if self.units[x]!=0])
        
    def __repr__(self):
        return '{value:.3e} ('.format(value=self.value)+self.get_unit()+')'
    
    def __str__(self):
        return 'SIunit({value:.3e},{unit})'.format(value=self.value, unit=self.get_unit())
        
    def get_scaled(self):
        if self.value == 0:
            return (0, self.get_unit())
            
        f = math.log10(self.value)
        if f<-24:
            f = -24
        elif f>24:
            f = 24
        unit_prefix = _unit_revert_scale[f]
        prefix_exp = _unit_scale[unit_prefix]
        val = self.value * 10**(-prefix_exp)
        return (val, unit_prefix)
            
    def sqrt(self):
        return SIunit(self.value**.5,**{x:self.units[x]/2 for x in self.units})
    
    def log(self):
        for x in self.units:
            if self.units[x] !=0:
                raise InvalidUnit
        return SIunit(math.log(self.value))

        
consts = {
    'qe': SIunit(const.qe, C=1),
    'NA': SIunit(const.NA, mol=-1),
    'me': SIunit(const.me*1e-3/const.NA, kg=1)
    }