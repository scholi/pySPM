# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Handle elements to calculate mass, abundance, etc.
"""

import sqlite3
import os
import re

from .constants import me

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

def simplify_formula(elt, debug=False):
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
        if (A,x) not in elts:
            elts[(A,x)] = n
        else:
            elts[(A,x)] += n
    if debug:
        print(elts)
    order = "CHSNO"
    res = ""
    for o in order:
        filt = [x[0] for x in elts if x[1]==o]
        miso = get_main_isotope(o)
        filt2 = [x for x in filt if x!=miso]
        for A in filt2:
            n = elts[(A,o)]
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

def get_isotopes(elt):
    res = {'':1}
    for A,x,n in re.findall('(?:\\^([0-9]+))?([A-Z][a-z]?)_?([0-9]*)',elt):
        if n == '':
            n = 1
        else:
            n = int(n)
        iso = get_isotopes_of_element(x)
        for i in range(n):
            r = []
            for x in res:
                    for y in iso:
                        r.append(simplify_formula(x+'^{1}{0}'.format(*y)))
            res = {}
            for x in r:
                if x not in res:
                    res[x] = 1
                else:
                    res[x] += 1
    return [(x,get_mass(x),res[x]*get_abund(x)) for x in res]

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
