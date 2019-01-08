# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Helper functions to handle spectras.
"""

from .misc import dec_debug, do_debug, alias

def get_substance_peaks(substance, negative=True):
    import os
    import sqlite3
    DB_PATH = os.path.join(os.path.abspath(os.path.join(__file__, "../..")), "data", "elements.db")
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT Peaks.Fragment from Peaks where Peaks.Substance==(SELECT ID from substance where Name LIKE '%{name}%') and Polarity{pol}=0".format(name=substance, pol='><'[negative]))
    return [x[0] for x in c.fetchall()]
        
def get_dm(m, sf, k0, dsf, dk0):
    import numpy as np
    return 2*np.sqrt(m)*np.sqrt((dk0**2/(sf**2))+m*(dsf**2/(sf**2)))

@alias("showPeak")
def show_peak(m, D, m0, delta=None, errors=False, dm0=0, dofit=False, show_elts=True,
    debug=False, Aredux=1, label=None, include=None, include_only=None, exclude=[],
    polarity="+", colors='rgb', pretty=True, formula=True, auto_scale=True, fakefit=False, zero_axis=True, **kargs):
    """
    Plot a spectrum (given by its mass vector m and intensity vector D) around a mass m0 (Δm=delta).
    
    Parameters
    ----------
    m : numpy.ndarray
        mass vector
    D : numpy.ndarray
        Intensity vector
    m0 : float
        Central mass of the plot's x-axis
    delta: None or float
        the half-span of the plot's x-axis.
        If None will try to determine automatically the peaks region
    errors : boolean
        If true will include the fitting errors (if dofit is True)
    dofit : boolean
        If True will fit the peaks for the elements found
    show_elts : boolean
        If True will display a list of all elements for the plotted region. The basic elements comes from a database, but they can be adjusted with the include, include_only and exclude parameters.
    include : list
        List of elements to be included to the plotting elements
    include_only : list
        Only those elements will be displayed and none from the database
    exclude : list
        Exclude some elements from display
    fakefit : boolean
        If the fit fails, fakefit can be used in order to see and adjust the first fitting guess
    dm0 : float
        The mass shift used for the fitting first guess.
    Aredux : float
        The amplitude reduction used to the first fitting guess.
    label : string
        Plot label. Will be displayed if you use the legend() function from matplotlib.
    polarity : string ("+" or "-")
        will switch the elements chosen from the database
    colors : string
        A series of letters defining matplotlib color used to display the elements.
    pretty : boolean
        If True will display the elements in a pretty way such that the elements text are not overlapping (or at least will avoid it as much as possible)
    formula : boolean
        Will convert the formula in a pretty manner. Exponents and indices will be written correctly
    auto_scale : boolean
        If True will try to rescale the plot in order to have all labels inside.
    zero_axis : boolean or dict
        Will draw an horizontal line for the zero intensities.
        If a dict, the parameters of the lines can be tuned. See the parameters passed to axhline from matplotlib
    """
    from . import LG, get_mass, get_peaklist
    from .elts import formulafy
    from scipy.optimize import curve_fit
    import numpy as np
    import copy
    import matplotlib.pyplot as plt

    if "showElts" in kargs:
        from warnings import warn
        warn("Parameter showElts is deprecated. Please use show_elts")
        show_elts = kargs.pop("showElts")
        
    if include is None:
        include = []
        
    if type(include) is str:
        include = include.split(',')
        
    assert type(include) is list

    if type(m0) is str:
        include.append(m0)
        m0 = get_mass(m0)

    if do_debug(debug):
        import time
        t0 = time.time()
    
    if type(exclude) is str:
       exclude = exclude.split(',')

    if fakefit:
        dofit = True
   
    if delta is None:
        delta0 = kargs.pop('delta0', .5)
        mask = (m>=m0-delta0)*(m<=m0+delta0)
        thresh = max(5, np.max(D[mask])*kargs.get('min_perc', .01))
        li = np.argmax(mask)+np.argmax(D[mask]> thresh)
        ui = np.argmax(mask)+np.sum(mask)-np.argmax(D[mask][::-1]>thresh)
        upper_mass = m[ui]+.01
        lower_mass = m[li]-.01
        
        delta = (upper_mass-lower_mass)/2.0+.01
    else:
        lower_mass = m0 - delta
        upper_mass = m0 + delta
    if do_debug(debug):
        print("Lower mass: {:.2f}\nUpper mass: {:.2f}".format(lower_mass, upper_mass))
    mask = (m>=lower_mass)*(m<=upper_mass)
    
    negative = False
    if polarity in ['-', 'Negative', 'negative', 'neg', 'Neg', 'NEG']:
        negative = True
    if show_elts:
        if include_only is not None:
            if type(include_only) is str:
                E = include_only.split(",")
            else:
                E = include_only
        else:
            E = [x for NM in range(int(round(lower_mass)), int(round(upper_mass))+1) for x in get_peaklist(NM, negative)]
            E = [x for x in E if x not in exclude] + include
        E = list(set([x+[['+', '-'][negative], '']['+' in x or '-' in x] for x in E]))
    else:
        E = []
    if negative:
        E = [x+['-', '']['-' in x] for x in E]
    if do_debug(debug):
        print("Selected Elements: "+", ".join(E))
    m0s = [get_mass(x) for x in E]
    E = [x for x, y in zip(E, m0s) if y>=lower_mass and y<=upper_mass]
    m0s = [get_mass(x) for x in E]
    if formula:
        E_labels = [formulafy(x) for x in E]
    else:
        E_labels = E
    if do_debug(debug):
        print("Displayed elements:", ", ".join(E))
    
    Dt = np.copy(D[mask])
    mt = m[mask]
    if do_debug(debug):
        print(" ; ".join(E))
    ax = kargs.pop('ax', plt.gca())
    if dofit or fakefit:
        mp = m[mask][np.argmax(D[mask])] # mass of max peak
        dm = dm0
        if dm is None:
            dm = 0
        if len(E)>0:
            i = np.argmin(abs(np.array([get_mass(x) for x in E if type(x) is str]+[x for x in E if type(x) is float])-mp)) # which element is the closest to max_peak
            if dm0 is None:
                dm = mp-get_mass(E[i])
        p0 = [kargs.get('asym0', 1), dm]+[0, 0]*len(E) # delta m is estimated from the deviation of the highest peak
        Et = copy.deepcopy(E) # copy element list    
    
        while len(Et)>0:
            mp = mt[np.argmax(Dt)] # mass of max peak
            ms = [get_mass(x) for x in Et]
            i = np.argmin(abs(ms-mp))
            idx = E.index(Et[i])
            j = np.argmin(abs(mt-ms[i]-dm))
            if do_debug(debug):
                print("Max element:", Et[i], mp, ms[i], Dt[j], idx)
            p0[2+2*idx] = kargs.get('sig0', 0.002)
            p0[3+2*idx] = Aredux*Dt[j]
            Dt -= LG(mt, ms[i], p0[2+2*idx], amp=p0[3+2*idx], asym=p0[0], lg=0)
            Dt[Dt<0] = 0
            del Et[i]
        
        def fit(x,*p):
            y = x*0
            for i in range((len(p)-2)//2):
                x0 = m0s[i]+p[1]
                y += LG(x, x0, p[2+2*i], amp=p[3+2*i], asym=p[0], lg=0)
            return y
            
        fit_type = None
        if dofit:
            try:
                assert not fakefit
                popt, pcov = curve_fit(fit, m[mask], D[mask], p0=p0,
                        bounds=(
                            [1/kargs.get('asym_max', 10), -0.015]+[0, 0]*((len(p0)-1)//2),
                            [kargs.get('asym_max', 10), 0.015]+[kargs.get('sig_max', 0.01), np.inf]*((len(p0)-1)//2))
                            )
                fit_type = 0
            except:
                p0[1] = 0
                try:
                    assert not fakefit
                    popt, pcov = curve_fit(fit, m[mask], D[mask], p0=p0,
                        bounds=(
                            [1/kargs.get('asym_max', 4),-0.015]+[0, 0]*((len(p0)-1)//2),
                            [kargs.get('asym_max', 4), 0.015]+[kargs.get('sig_max', 0.01), np.inf]*((len(p0)-1)//2))
                            )
                    fit_type = 1
                except Exception as e:
                    fit_type = 2
                    if do_debug(debug):
                        raise e
                    popt = p0
                    pcov = np.zeros((len(p0), len(p0)))
                    if ax is not None:
                        for x in ['right', 'left', 'top', 'bottom']:
                            ax.spines[x].set_color('red')
        else:
            popt = p0
            popt[1] = 0
            pcov = np.zeros((len(p0), len(p0)))
        if fakefit:
            popt[1] = dm0
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
    else:
        popt = [0, 0]
        res = None
    if ax is not None:
        if label is None:
            p = ax.plot(m[mask]-popt[1], D[mask], color=kargs.get('curve_color', None))
        else:
            p = ax.plot(m[mask]-popt[1], D[mask], label=label, color=kargs.get('curve_color', None))
    
    if show_elts and ax is not None:
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
            
        if 'dsf' in kargs and 'dk0' in kargs and 'sf' in kargs and 'k0' in kargs:
            sf, dsf, k0, dk0 = [kargs[x] for x in 'sf,dsf,k0,dk0'.split(',')]
            P = list(zip(m0s, E))
            P.sort(key=lambda x: x[0])
            for i,(mi,Ei) in enumerate(P):
                col = colors[i%len(colors)]
                dmi = get_dm(mi, sf, k0, dsf, dk0)
                ax.axvline(mi-dmi, color=col, alpha=.5, linestyle=':')
                ax.axvline(mi+dmi, color=col, alpha=.5, linestyle=':')
            
        if not pretty:
            P = list(zip(m0s, E_labels))
            P.sort(key=lambda x: x[0])
            y = ax.get_ylim()[1]
            for i,(mi,Ei) in enumerate(P):
                col = colors[i%len(colors)]
                ax.axvline(mi, color=col, alpha=.5)
                ax.annotate(Ei, (mi,y), (5,0), rotation=90, va='top', ha='left', textcoords='offset pixels', color=col)
        if auto_scale:
            fig = plt.gcf()
            renderer = fig.canvas.get_renderer()
            ymax = ax.get_ylim()[1]
            ymin = ax.get_ylim()[0]
            for child in ax.get_children():
                bbox = child.get_window_extent(renderer)
                bbox_data = bbox.transformed(ax.transData.inverted())
                if ymax<bbox_data.ymax:
                    ymax = bbox_data.ymax
            ax.set_ylim((ymin, ymax))
            ax.autoscale_view()    
            
    # Add the Intensity=0 line in case there is negative intensities (happens when plotting PCA components)
    if zero_axis:
        if type(zero_axis) is dict:
            ax.achline(0, **zero_axis)
        else:
            ax.axhline(0, color='k', alpha=.5, lw=.5)
        
    if dofit and ax is not None:
        ax.plot(m[mask]-popt[1], fit(m[mask], *popt), 'r:');
    return res


def plot_isotopes(elt, amp=None, main=None, ax=None, sig=None, asym=None, lg=0, limit=1, color='C1', show_elts=False, debug=False, dm=0, **kargs):
    """
    plot the isotopes of a given element on a spectral profile plot
    """
    # Old parameter name compatibility
    if 'Amp' in kargs:
        from warnings import warn
        warn("Parameter Amp is deprecated. Please use amp in order to set the amplitude!")
        amp = kargs.pop("Amp")
    if 'showElts' in kargs:
        from warnings import warn
        warn("Parameter showElts is deprecated. Please use show_elts")
        show_elts = kargs.pop("showElts")
    from . import get_isotopes, LG, get_mass, get_abund, Molecule
    import matplotlib.pyplot as plt
    import numpy as np
    import re
    from scipy.optimize import curve_fit
    if type(elt) is Molecule:
        elt = str(elt)
    if ax is None:
        ax = plt.gca()
    if main is None:
        main = ax
    main_isotope = re.sub('^([0-9]+)', '', elt)
    main_iso = (main_isotope, get_mass(main_isotope), get_abund(main_isotope))
    L = main.lines[0]
    m, y = L.get_xdata(), L.get_ydata()
    m -= dm
    if hasattr(ax, 'log') and ax.log:
        y = 10**y
    mask = np.abs(m-main_iso[1])<.1
    if amp is None:
        i = np.argmin(np.abs(m-main_iso[1]))
        amp = np.max(y[i-10:i+10])/main_iso[2]
    if sig is None and asym is None:
        def fit(x, s, a):
            return LG(x, main_iso[1], s, amp=amp*main_iso[2], lg=lg, asym=a)
        (sig, asym),_ = curve_fit(fit, m[mask], y[mask], (0.005, 1))
    elif sig is None:
        def fit(x, s):
            return LG(x, main_iso[1], s, amp=amp*main_iso[2], lg=lg, asym=asym)
        s,_ = curve_fit(fit, m[mask], y[mask], (0.005))
        sig = s[0]
    elif asym is None:
        def fit(x, a):
            return LG(x, main_iso[1], sig, amp=amp*main_iso[2], lg=lg, asym=a)
        a, _ = curve_fit(fit, m[mask], y[mask], (1))
        asym = a[0]
    if debug:
        print(sig, asym)
    isos = get_isotopes(elt, min_abund=limit/amp)
    L = ax.lines[0]
    m, y = L.get_xdata(), L.get_ydata()
    if hasattr(ax, 'log') and ax.log:
        y = 10**y
    r = main.get_xlim()
    s = m*0
    
    for iso in isos:
        s += LG(m, iso[1], sig, amp=amp*iso[2], lg=lg, asym=asym)
        if show_elts:
            ax.annotate(iso[0], (iso[1], amp*iso[2]), (0,5), textcoords='offset pixels', ha='center')
    if hasattr(ax, 'log') and ax.log:
        s[s>=1] = np.log10(s[s>=1])
        s[s<1] = 0
    ax.plot(m+dm, s, color, **kargs);
    return m, s
    
