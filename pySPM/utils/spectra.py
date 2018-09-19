# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Helper functions to handle spectras.
"""

from.misc import dec_debug, do_debug

def get_substance_peaks(substance, negative=True):
    import os
    import sqlite3
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__,"../..")),"data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT Peaks.Fragment from Peaks where Peaks.Substance==(SELECT ID from substance where Name LIKE '%{name}%') and Polarity{pol}=0".format(name=substance,pol='><'[negative]))
    return [x[0] for x in c.fetchall()]

        
def formulafy(x):
    import re
    return '$'+re.sub('([a-zA-Z])_?([0-9]+)',r'\1_{\2}',re.sub(r'\^([0-9]+)',r'^{\1}',re.sub('([+-])$',r'^\1',x)))+'$'

def showPeak(m,D,m0, delta=0.15, errors=False, dm0=0, dofit=False, showElts=True, debug=False, Aredux=1, label=None, include=[], include_only=None, exclude=[], polarity="+", colors='rgb', pretty=True, formula=True, **kargs):
    """
    given masses m and Spectrum D, zoom around m0 with Δm=delta.
    Will perform a peak-fitting if dofit is True
    """
    from . import LG, get_mass, get_peaklist
    from scipy.optimize import curve_fit
    import numpy as np
    import copy
    import matplotlib.pyplot as plt
    if do_debug(debug):
        import time
        t0 = time.time()
    if type(include) is str:
        include=[include]
    if type(exclude) is str:
       exclude = [exclude]
    mask = (m>=m0-delta)*(m<=m0+delta)
    negative = False
    if polarity in ['-','Negative','negative','neg','Neg','NEG']:
        negative = True
    if include_only is not None:
        if type(include_only) is str:
            E = [include_only]
        else:
            E = include_only
    else:
        E = get_peaklist(int(round(m0)), negative)
        E = [x for x in E if x not in exclude] + include
    E = list(set(E))
    if formula:
        E_labels = [formulafy(x) for x in E]
    else:
        E_labels = E
    if do_debug(debug):
        print("Elements:",", ".join(E))
    mp = m[mask][np.argmax(D[mask])] # mass of max peak
    dm = dm0
    if dm is None:
        dm=0
    if len(E)>0:
        i = np.argmin(abs(np.array([get_mass(x+'+') for x in E if type(x) is str]+[x for x in E if type(x) is float])-mp)) # which element is the closest to max_peak
        if dm0 is None:
            dm = mp-get_mass(E[i]+'+')
    p0 = [1,dm]+[0,0]*len(E) # delta m is estimnated from the deviation of the highest peak
    m0s = [get_mass(x+'+') for x in E]
    Et = copy.deepcopy(E) # copy element list
    if do_debug(debug):
        print(" ; ".join(E))
    Dt = np.copy(D[mask])
    mt = m[mask]
    
    while len(Et)>0:
        mp = mt[np.argmax(Dt)] # mass of max peak
        ms = [get_mass(x+'+') for x in Et]
        i = np.argmin(abs(ms-mp))
        idx = E.index(Et[i])
        j = np.argmin(abs(mt-ms[i]-dm))
        if do_debug(debug):
            print("Max element:",Et[i],mp,ms[i],Dt[j], idx)
        p0[2+2*idx] = kargs.get('sig0',0.002)
        p0[3+2*idx] = Aredux*Dt[j]
        Dt -= LG(mt,ms[i],p0[2+2*idx] , Amp=p0[3+2*idx], lg=0)
        Dt[Dt<0] = 0
        del Et[i]
    
    def fit(x,*p):
        y = x*0
        for i in range((len(p)-2)//2):
            x0 = m0s[i]+p[1]
            y += LG(x, x0, p[2+2*i], Amp=p[3+2*i], asym=p[0], lg=0)
        return y
        
    ax = kargs.pop('ax',plt.gca())
    
    fit_type = None
    if do_debug(debug):
        print("p0",p0)
        t1 = time.time()
        print("setup time: ",t1-t0)
    if dofit:
        try:
            popt, pcov = curve_fit(fit, m[mask], D[mask], p0=p0,
                    bounds=(
                        [1/kargs.get('asym_max', 10),-0.015]+[0,0]*((len(p0)-1)//2),
                        [kargs.get('asym_max', 10),0.015]+[kargs.get('sig_max', 0.01),np.inf]*((len(p0)-1)//2))
                        )
            fit_type = 0
        except:
            p0[1] = 0
            try:
                popt, pcov = curve_fit(fit, m[mask], D[mask], p0=p0,
                        bounds=(
                        [1/kargs.get('asym_max', 10),-0.015]+[0,0]*((len(p0)-1)//2),
                        [kargs.get('asym_max', 10), 0.015]+[kargs.get('sig_max', 0.01),np.inf]*((len(p0)-1)//2))
                        )
                fit_type = 1
            except Exception as e:
                fit_type = 2
                if do_debug(debug):
                    raise e
                popt = p0
                pcov = np.zeros((len(p0),len(p0)))
                if ax is not None:
                    for x in ['right','left','top','bottom']:
                        ax.spines[x].set_color('red')
    else:
        popt = p0
        popt[1] = 0
        pcov = np.zeros((len(p0),len(p0)))
    if ax is not None:
        if label is None:
            p = ax.plot(m[mask]-popt[1],D[mask])
        else:
            p = ax.plot(m[mask]-popt[1],D[mask], label=label)
    res = {}
    err = np.sqrt(np.diag(pcov))
    for i in range((len(popt)-1)//2):
        Y = LG(m[mask], m0s[i], popt[2+2*i], popt[3+2*i],lg=0, asym=popt[0])
        Area = popt[2*i+3]*popt[2*i+2]*np.sqrt(2*np.pi)*(.5+.5*popt[0])
        Area_err = np.sqrt(
            (err[2*i+3]*popt[2*i+2]*np.sqrt(2*np.pi)*(.5+.5*popt[0]))**2+
            (err[2*i+2]*popt[2*i+3]*np.sqrt(2*np.pi)*(.5+.5*popt[0]))**2+
            (err[0]*popt[2*i+2]*popt[2*i+3]*np.sqrt(2*np.pi)*.5)**2
            )
        if not errors:
            res[E_labels[i]] = {
                'm0': m0s[i],
                'mass': m0s[i]+popt[1],
                'Area' : Area,
                'Amp' : popt[2*i+3],
                'sig' : popt[2*i+2],
                'assym' : popt[0],
                'dm': popt[1]*1e6/m0s[i],
                'fit': fit_type
                }
        else:
            res[E_labels[i]] = {
                'm0': m0s[i],
                'mass': (m0s[i]+popt[1],err[1]),
                'Area' : (Area,Area_err),
                'Amp' : (popt[2*i+3],err[2*i+3]),
                'sig' : (popt[2*i+2],err[2*i+2]),
                'assym' : (popt[0],err[0]),
                'dm': (popt[1]*1e6/m0s[i],err[1]*1e6/m0s[i]),
                'fit': fit_type
                }
    if showElts and ax is not None:
        if do_debug(debug):
            print("put labels")
            t2 = time.time()
            print("fitting time: ",t2-t1)
        from . import put_Xlabels
        if dofit and ax is not None:
            ax.plot(m[mask],D[mask],color=p[0].get_color(), alpha=.1)
            if pretty:
                put_Xlabels(ax, m0s, ["{0}: {res[Area]:.2f}".format(E,res=res[E]) for E in res], color=colors, debug=dec_debug(debug), **kargs)
            resO = [(res[x]['m0'],res[x]) for x in res]
            resO.sort(key=lambda x: x[0])
            for i,(_,r) in enumerate(resO):
                col = colors[i%len(colors)]
                Y = LG(m[mask], r['m0'], r['sig'], r['Amp'],lg=0, asym=r['assym'])
                ax.plot(m[mask], Y, '--', color=col);
        elif pretty:
            put_Xlabels(ax, m0s, E_labels, color=colors, debug=dec_debug(debug), **kargs)
        if not pretty:
            P = list(zip(m0s, E))
            P.sort(key=lambda x: x[0])
            y = ax.get_ylim()[1]
            for i,(mi,Ei) in enumerate(P):
                col = colors[i%len(colors)]
                ax.axvline(mi, color=col, alpha=.5)
                ax.annotate(Ei, (mi,y), (5,0), rotation=90, va='bottom', ha='left', textcoords='offset pixels')

    if do_debug(debug):
        print(popt)
    if (dofit or do_debug(debug)) and ax is not None:
        ax.plot(m[mask]-popt[1], fit(m[mask], *popt), 'r:');
    if do_debug(debug):
        t3 = time.time()
        print("labeling time: ",t3-t2)
    return res

