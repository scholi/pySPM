def showPeak(m,D,m0, delta=0.15, ax=None, dm0=None, dofit=False, showElts=True, debug=False, sig0=0.002, Aredux=1,sf=None,k0=None, label=None):
    """
    gives masses m and Spectrum D, zoom around m0 with Δm=delta.
    Will perform a peak-fitting if dofit is True
    """
    from . import LG, getMass, elts
    from scipy.optimize import curve_fit
    import numpy as np
    import copy
    import matplotlib.pyplot as plt
    
    mask = (m>=m0-delta)*(m<=m0+delta)
    if int(round(m0)) in elts:
        E = elts[int(round(m0))]
    else:
        E=[]
    if type(E) is not list:
        E=[E]
    mp = m[mask][np.argmax(D[mask])] # mass of max peak
    i = np.argmin(abs(np.array([getMass(x+'+') for x in E if type(x) is str]+[x for x in E if type(x) is float])-mp)) # which element is the clothest to max_peak
    if dm0 is None:
        dm = mp-getMass(E[i]+'+')
    else:
        dm = dm0
    p0 = [1,dm,.5,.5]+[0,0]*len(E) # delta m is estimnated from the deviation of the highest peak
    m0s = [getMass(x+'+') for x in E]
    Et = copy.deepcopy(E) # copy element list
    if debug:
        print(" ; ".join(E))
    Dt = np.copy(D[mask])
    mt = m[mask]
    
    while len(Et)>0:
        mp = mt[np.argmax(Dt)] # mass of max peak
        ms = [getMass(x+'+') for x in Et]
        i = np.argmin(abs(ms-mp))
        idx = E.index(Et[i])
        j = np.argmin(abs(mt-ms[i]-dm))
        if debug:
            print("Max element:",Et[i],mp,ms[i],Dt[j], idx)
        p0[4+2*idx] = sig0
        p0[5+2*idx] = Aredux*Dt[j]
        Dt -= LG(mt,ms[i],p0[4+2*idx] ,p0[5+2*idx])
        Dt[Dt<0] = 0
        del Et[i]
    
    def fit(x,*p):
        y = x*0
        for i in range((len(p)-4)//2):
            x0 = m0s[i]+p[1]
            y += LG(x, x0, p[4+2*i], Amp=p[5+2*i],lg=p[2], asym=p[0])
        return y
        
    if ax is None:
        ax = plt.gca()
    if label is None:
        ax.plot(m[mask],D[mask])
    else:
        ax.plot(m[mask],D[mask], label=label)
    for i in range((len(p0)-3)//2):
        ax.axvline(m0s[i], color='r', alpha=.1, linestyle=':');
    if dofit:
        try:
            popt, pcov = curve_fit(fit, m[mask], D[mask], p0=p0,
                    bounds=(
                        [0,-0.015,0,0]+[0,0]*((len(p0)-3)//2),
                        [np.inf,0.015,1,1]+[0.01,np.inf]*((len(p0)-3)//2))
                        )
        except:
            p0[1] = 0
            try:
                popt, pcov = curve_fit(fit, m[mask], D[mask], p0=p0,
                        bounds=(
                        [0,-0.015,0,0]+[0,0]*((len(p0)-3)//2),
                        [np.inf,0.015,1,1]+[0.01,np.inf]*((len(p0)-3)//2))
                        )
            except Exception as e:
                if debug:
                    raise e
                popt=p0
                for x in ['right','left','top','bottom']:
                    ax.spines[x].set_color('red')
    else:
        popt = p0
    res = {}
    
    for i in range((len(popt)-3)//2):
        Y = LG(m[mask], m0s[i]+popt[1], popt[4+2*i], popt[5+2*i],lg=popt[2], asym=popt[0])
        Area = popt[2*i+5]*popt[2*i+4]*np.sqrt(2*np.pi)*(.5+.5*popt[0])
        res[E[i]] = {
            'm0': m0s[i],
            'mass': m0s[i]+popt[1],
            'Area' : Area,
            'Amp' : popt[2*i+5],
            'sig' : popt[2*i+4],
            'assym' : popt[0],
            'dm': popt[1]*1e6,
            'lgL': popt[2],
            'lgR': popt[3]
            }
        if showElts:
            ax.axvline(m0s[i]+popt[1], color='r', alpha=.3);
            ax.annotate(E[i], (m0s[i]+popt[1]/2, ax.get_ylim()[1]), rotation=90, va='top');
        if dofit:
            ax.plot(m[mask], Y,'--');
            ax.annotate("{:.2f}".format(Area), (m0s[i]+popt[1], popt[5+2*i]/2))
    if debug:
        print(popt)
    if dofit or debug:
        ax.plot(m[mask], fit(m[mask], *popt), 'r:');
    return res