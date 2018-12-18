# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Handle elements to calculate mass, abundance, etc.
"""

from __future__ import absolute_import

import sqlite3
import os
import re
from .constants import me
from .misc import deprecated

def formulafy(x):
    """
    Convert the input string x to a Latex representation of the chemical formula.
    e.g. ^13CH4 -> $^{13}CH_{4}$
    This can be used in matplotlib to display nicely a chemical formula.
    """
    import re
    return '$'+re.sub('([a-zA-Z])_?([0-9]+)', r'\1_{\2}', re.sub(r'\^([0-9]+)', r'^{\1}', re.sub('([\\+-]+)$', r'^{\1}', x)))+'$'
    
def get_peaklist(nominal_mass, negative=False):
    """
    Retrieve all elements from the database which have a given nominal_mass
    
    Parameters
    ----------
    nominal_mass : int or tuple/list
        If int provide the nominal mass of the elements
        If tuple/list return all elements bounded between a lower (first element) and upper (second eleemnt) nominal mass.
    negative : boolean
        If True will return the elements in the database for negative polarity
    
    Examples
    --------
    get_peaklist(12) will return all elements with a nominal mass of 12u
    get_peaklist((12,15)) will return all elements with a nominal mass between 12 and 15
    """
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    if type(nominal_mass) in [float, int]:
        c.execute("SELECT Formula from fragments where NominalMass={nm} and (Polarity is NULL or Polarity=={pol})".format(nm=nominal_mass,pol=[1,-1][negative]))
    elif type(nominal_mass) in [tuple, list]:
        c.execute("SELECT Formula from fragments where NominalMass>={nm[0]} and NominalMass<={nm[1]} and (Polarity is NULL or Polarity=={pol})".format(nmL=nominal_mass,pol=[1,-1][negative]))
    return [str(x[0]) for x in c.fetchall()]

def get_properties(elt):
    """
    Retrieve properties of an element such as its density, lattice constant, etc. as a dictionary
    If no information is found in the database then at least Element, Z and mass will be returned
    Examples
    --------
    > get_properties("C")
    {'Element': 'C', 'Z': 6, 'mass': 12.010735896764459}
    
    > get_properties("Si")
    {'Element': 'Si',
    'Q': 0.78,
    'Z': 14,
    'density': 5e+22,
    'ion_nrj': 8.15,
    'lattice_const': 5.431,
    'lattice_type': 'diamond_fcc',
    'mass': 28.085498706418875,
    'sublim_nrj': 4.63,
    'work_func': 4.85}
    """
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT * from properties where Element='{}'".format(elt))
    res = c.fetchone()
    r = {}
    if res:
        for idx, col in enumerate(c.description):
            r[col[0]] = res[idx]
    else:
        res = c.execute("SELECT mass,abund,symbol,Z from elements where symbol='{}'".format(elt))
        if res:
            row = res.fetchall()
            r['Element'] = row[0][2]
            r['Z'] = row[0][3]
            r['mass'] = sum([x[0]*x[1] for x in row])
    return r
    
    
def get_mass(elt, mz=True):
    """
    Calculate the exact molar mass of a givent chemical compound.
    
    Parameters
    ----------
    elt : string
        Element given by its chemical formula
    mz : boolean
        Will output the mass per ion.
        e.g. This means that if the formula is Si++, then the answer will be half the mass of Si
    
    Notes
    -----
    "Si+" will subtract the mass of one electron compared to "Si"
    "Si-" will add the mass of one electron compared to "Si"
    """
    if type(elt) is Molecule:
        elt = str(elt)
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__, "../..")), "data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    m = 0
    for A,x,n in re.findall('(?:\\^([0-9]+))?([A-Z][a-z]?)_?([0-9]*)',elt):
        if n == '':
            n = 1
        else:
            n = int(n)
        if A == '':
            c.execute("SELECT mass from elements where symbol='{sym}' and abund=(SELECT max(abund) from elements where symbol='{sym}')".format(sym=x))
        else:
            c.execute("SELECT mass from elements where symbol='{sym}' and A={A}".format(sym=x, A=A))
        res = c.fetchone()
        if res is None:
            raise Exception("Cannot fetch mass of {}".format(x))
        m += n*res[0]
    m -= me*elt.count('+')
    m += me*elt.count('-')
    Z = max(1, abs(elt.count('+')-elt.count('-')))
    if mz:
        return m/Z
    return m
    
def get_main_isotope(elt):
    """
    Return the A value of the main isotope of a given element.
    e.g. get_main_isotope("Si") will return 28
    """
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT A from elements where symbol='{sym}' and abund=(SELECT max(abund) from elements where symbol='{sym}')".format(sym=elt))
    res = c.fetchone()
    if res is None:
        raise Exception("Cannot fetch element '{sym}'".format(sym=elt))
    return res[0]

def is_main_isotope(elt, A):
    return get_main_isotope(elt)==A

def get_isotopes_of_element(elt):
    """
    Get all the isotopes for a given element
    """
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT symbol, A, abund from elements where symbol='{sym}'".format(sym=elt))
    return c.fetchall()

def _formula2dict(elt, iso=True):
    """
    Convert a formula into a dictionary.
    """
    elts = {}
    for A,x,n in re.findall('(?:\\^([0-9]+))?([A-Z][a-z]?)_?([0-9]*)', elt):
        if A == '':
            A = get_main_isotope(x)
        else:
            A = int(A)
        if n == '':
            n = 1
        else:
            n = int(n)
        if iso:
            if (A,x) not in elts:
                elts[(A,x)] = n
            else:
                elts[(A,x)] += n
        else:
            if x not in elts:
                elts[x] = n
            else:
                elts[x] += n
    return elts

def is_fragment_of(fragment, element):
    f = _formula2dict(fragment, iso=False)
    e = _formula2dict(element, iso=False)
    for x in f:
        if x not in e or e[x]<f[x]:
            return False
    return True

def simplify_formula(elt, debug=False):
    elts = _formula2dict(elt)
    return _dict2formula(elts)
    
def _dict2formula(elts, debug=False):
    if debug:
        print(elts)
    order = "CHSNO"
    res = ""
    for o in order:
        filt = [x[0] for x in elts if x[1]==o]
        miso = get_main_isotope(o)
        filt2 = [x for x in filt if x!=miso]
        for A in filt2:
            n = elts[(A, o)]
            if n==0:
                continue
            fmt = "^{A}{sym}"
            if n>1:
                fmt += "{n}"
            if debug:
                print("-",fmt, A,o,n)
            res += fmt.format(A=A, sym=o, n=n)
        if (miso,o) in elts:
            n = elts[(miso,o)]
            if n==0:
                break
            fmt = "{sym}"
            if n>1:
                fmt += "{n}"
            if debug:
                print("*",fmt, A,o,n)
            res += fmt.format(sym=o, n=n)
    for A, sym in sorted([x for x in elts if x[1] not in order], key=lambda x: x[1]):
        n = elts[(A,sym)]
        if n==0:
            continue
        if is_main_isotope(sym, A):
            fmt = "{sym}"
        else:
            fmt = "^{A}{sym}"
        if n>1:
            fmt += "{n}"
        res += fmt.format(A=A, sym=sym, n=n)
    return res 

def _get_isotopes_elt(n, x, min_abund=0):
    """
    get the isotopes of n atoms x
    """
    res = {'':1}
    iso = get_isotopes_of_element(x)
    for i in range(n):
        r = []
        Nres = {}
        for x in res:
            for y in iso:
                X = simplify_formula(x+'^{1}{0}'.format(*y))
                ab = res[x]*y[2]
                if X not in Nres:
                    if n*ab>=min_abund:
                        Nres[X] = ab
                else:
                    Nres[X] += ab
        res = Nres
    return [(x, res[x]) for x in res]

def get_isotopes(elt, min_abund=0):
    charges = '+'*elt.count('+')+'-'*elt.count('-')
    res = {'':1}
    for A,x,n in re.findall('(?:\\^([0-9]+))?([A-Z][a-z]?)_?([0-9]*)',elt):
        if n == '':
            n = 1
        else:
            n = int(n)
        R = _get_isotopes_elt(n, x, min_abund=min_abund)
        Nres = {}
        for x0 in res:
            for x, ab0 in R:
                ab = res[x0]*ab0
                if ab >= min_abund:
                    Nres[x0+x] = ab
        res = Nres
    R = [(x+charges, get_mass(x+charges), res[x]) for x in res]
    return R

def get_abund(elt):
    from .math import perm, prod
    if type(elt) is Molecule:
        elt = str(elt)
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    abund = 1
    elts = [(int([A,get_main_isotope(x)][A=='']), x , int(n+['', '1'][n==''])) for A, x, n in re.findall('(?:\\^([0-9]+))?([A-Z][a-z]?)_?([0-9]*)', elt)]
    ns = {}
    
    for x in elts:
        if x[1] not in ns:
            ns[x[1]] = []
        ns[x[1]].append(x[2])
        
    permutations = prod([perm(ns[k]) for k in ns])
    
    for A, x, n in elts:
        if A == '':
            c.execute("SELECT abund from elements where symbol='{sym}' and abund=(SELECT max(abund) from elements where symbol='{sym}')".format(sym=x))
        else:
            c.execute("SELECT abund from elements where symbol='{sym}' and A={A}".format(sym=x, A=A))
        res = c.fetchone()
        if res is None:
            raise Exception("Cannot fetch mass of {}".format(x))
        abund *= (res[0])**n
    return permutations*abund

@deprecated("getOrganicAt")
def get_organic_at(m0):
    import numpy as np
    res = []
    for C in range(1, 1+m0//12):
       for N in range(1+(m0-12*C)//15):
        for S in range(1+(m0-12*C-15*N)//32):
            for O in range(1+(m0-12*C-15*N-32*S)//16):
                H = m0-12*C-15*N-32*S-16*O
                elt = ''
                for e in 'CNSOH':
                    if locals()[e]>0:
                        elt += e
                        if locals()[e]>1:
                            elt += str(locals()[e])
                if int(np.round(get_mass(elt), 0))==m0:
                    res.append(elt)
    return res
    
def elts_substract(elt1, elt2):
    r = _formula2dict(elt1)
    d = _formula2dict(elt2)
    for x in d:
        if x not in r:
            raise Exception("Cannot subtract a molecule with more atoms than the main one")
        r[x] -= d[x]
        
def _dict_add(d1, d2):
    r = {k:d1[k] for k in d1}
    for k in d2:
        r[k] = r.get(k, 0) + d2[k]
    return r
    
def _elts_add(eltsa, eltsb):
    r = []
    for x,mx in eltsa:
        for y,my in eltsb:
            r.append((_dict_add(x, y), mx+my))
    return r

@deprecated("eltsNM")
def elts_nm(elts, NM):
    import re
    r = [({}, 0)]
    res = []
    if type(elts) is str:
        elts = re.compile(r'((?:^[0-9]+)?[A-Z][a-z]?(?:_[0-9]+)?)').findall(elts)
    me = [get_mass(x) for x in elts]
    while r:
        rnew = _elts_add(r, list(zip([_formula2dict(x) for x in elts], me)))
        r = [x for x in rnew if x[1]<NM+.5]
        res += [_dict2formula(x[0]) for x in r if x[1]>NM-.5 if _dict2formula(x[0]) not in res]
        res = list(set(res))
    return res
    
class Molecule:    
    def __init__(self, formula):
        self.atoms = _formula2dict(formula)
    
    def formula(self):
        return _dict2formula(self.atoms)
        
    def __add__(self, b):
        r = Molecule(self.formula())
        for x in b.atoms:
            if x in r.atoms:
                r.atoms[x] += b.atoms[x]
            else:
                r.atoms[x] = b.atoms[x]
        return r
        
    __radd__ = __add__
    
    def __mul__(self, n):
        assert type(n) is int
        r = Molecule(self.formula())
        for x in r.atoms:
            r.atoms[x] *= n
        return r
        
    __rmul__ = __mul__
    
    def __sub__(self, b):
        r = Molecule(self.formula())
        for x in b.atoms:
            n = b.atoms[x]
            if x not in r.atoms or r.atoms[x]<n:
                raise Exception("Cannot subtract a larger molecule from a smaller one")
            r.atoms[x] -= n
        return r
    
    def inc(self, elt='H'):
        r = Molecule(self.formula())
        k = list(_formula2dict(elt).keys())[0]
        if k in r.atoms:
            r.atoms[k] += 1
        else:
            r.atoms[k] = 1
        return r
        
    def dec(self, elt='H'):
        r = Molecule(self.formula())
        k = list(_formula2dict(elt).keys())[0]
        if k not in r.atoms or r.atoms[k]<1:
            raise Exception("Cannot decrement molecule below 0 atoms")
        r.atoms[k] -= 1
        return r
            
    def mass(self):
        return get_mass(self.formula())
        
    def abund(self):
        return get_abund(self.formula())
        
    def __repr__(self):
        return self.formula()
        
    def __str__(self):
        return self.formula()

H = Molecule('H')