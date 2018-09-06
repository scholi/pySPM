# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Handle elements to calculate mass, abundance, etc.
"""

import sqlite3
import os
import re

from .constants import me

def get_peaklist(nominal_mass, negative=False):
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT Formula from fragments where NominalMass={nm} and (Polarity is NULL or Polarity=={pol})".format(nm=nominal_mass,pol=[1,-1][negative]))
    return [str(x[0]) for x in c.fetchall()]

def get_mass(elt):
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
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
    return m
    
def get_main_isotope(elt):
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT A from elements where symbol='{sym}' and abund=(SELECT max(abund) from elements where symbol='{sym}')".format(sym=elt))
    res = c.fetchone()
    if res is None:
        raise Exception("Cannot fetch element '{sym}'".format(sym=elt))
    return res[0]

def is_main_isotope(elt, A):
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT A from elements where symbol='{sym}' and abund=(SELECT max(abund) from elements where symbol='{sym}')".format(sym=elt))
    res = c.fetchone()
    if res is None:
        raise Exception("Cannot fetch element '{sym}'".format(sym=elt))
    return res[0]==A


def get_isotopes_of_element(elt):
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT symbol, A, abund from elements where symbol='{sym}'".format(sym=elt))
    return c.fetchall()

def formula2dict(elt, iso=True):
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
    f = formula2dict(fragment, iso=False)
    e = formula2dict(element, iso=False)
    for x in f:
        if x not in e or e[x]<f[x]:
            return False
    return True

def simplify_formula(elt, debug=False):
    elts = formula2dict(elt)
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

def __get_isotopes_elt(n, x, min_abund=0):
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
    res = {'':1}
    for A,x,n in re.findall('(?:\\^([0-9]+))?([A-Z][a-z]?)_?([0-9]*)',elt):
        if n == '':
            n = 1
        else:
            n = int(n)
        R = __get_isotopes_elt(n, x, min_abund=min_abund)
        Nres = {}
        for x0 in res:
            for x, ab0 in R:
                ab = res[x0]*ab0
                if ab >= min_abund:
                    Nres[x0+x] = ab
        res = Nres
    R = [(x,get_mass(x),res[x]) for x in res]
    return R

def get_abund(elt):
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    abund = 1
    for A,x,n in re.findall('(?:\\^([0-9]+))?([A-Z][a-z]?)_?([0-9]*)',elt):
        if n == '':
            n = 1
        else:
            n = int(n)
        if A == '':
            c.execute("SELECT abund from elements where symbol='{sym}' and abund=(SELECT max(abund) from elements where symbol='{sym}')".format(sym=x))
        else:
            c.execute("SELECT abund from elements where symbol='{sym}' and A={A}".format(sym=x, A=A))
        res = c.fetchone()
        if res is None:
            raise Exception("Cannot fetch mass of {}".format(x))
        abund *= (res[0])**n
    return abund

def getOrganicAt(m0):
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
