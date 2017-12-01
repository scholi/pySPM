def showPeak(m,D,m0, delta=0.15, errors=False,dm0=None, dofit=False, showElts=True, debug=False, Aredux=1,label=None, include=[], exclude=[], **kargs):
    """
    given masses m and Spectrum D, zoom around m0 with Δm=delta.
    Will perform a peak-fitting if dofit is True
    """
    from . import LG, getMass, elts
    from scipy.optimize import curve_fit
    import numpy as np
    import copy
    import matplotlib.pyplot as plt
    
    if type(include) is str:
        include=[include]
    if type(exclude) is str:
       exclude = [exclude]
    mask = (m>=m0-delta)*(m<=m0+delta)
    if int(round(m0)) in elts:
        E = elts[int(round(m0))]
    else:
        E=[]
    if type(E) is not list:
        E=[E]
    E = [x for x in E if x not in exclude] + include
    if debug:
        print("Elements:",", ".join(E))
    mp = m[mask][np.argmax(D[mask])] # mass of max peak
    dm = dm0
    if dm is None:
        dm=0
    if len(E)>0:
        i = np.argmin(abs(np.array([getMass(x+'+') for x in E if type(x) is str]+[x for x in E if type(x) is float])-mp)) # which element is the clothest to max_peak
        if dm0 is None:
            dm = mp-getMass(E[i]+'+')
    p0 = [1,dm]+[0,0]*len(E) # delta m is estimnated from the deviation of the highest peak
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
        
    ax = kargs.get('ax',plt.gca())
    if label is None:
        ax.plot(m[mask],D[mask])
    else:
        ax.plot(m[mask],D[mask], label=label)
    for i in range((len(p0)-1)//2):
        ax.axvline(m0s[i], color='r', alpha=.3, linestyle=':');
    fit_type = None
    if debug:
        print("p0",p0)
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
                if debug:
                    raise e
                popt = p0
                pcov = np.zeros((len(p0),len(p0)))
                for x in ['right','left','top','bottom']:
                    ax.spines[x].set_color('red')
    else:
        popt = p0
        pcov = np.zeros((len(p0),len(p0)))
    res = {}
    err = np.sqrt(np.diag(pcov))
    for i in range((len(popt)-1)//2):
        Y = LG(m[mask], m0s[i]+popt[1], popt[2+2*i], popt[3+2*i],lg=0, asym=popt[0])
        Area = popt[2*i+3]*popt[2*i+2]*np.sqrt(2*np.pi)*(.5+.5*popt[0])
        Area_err = np.sqrt(
            (err[2*i+3]*popt[2*i+2]*np.sqrt(2*np.pi)*(.5+.5*popt[0]))**2+
            (err[2*i+2]*popt[2*i+3]*np.sqrt(2*np.pi)*(.5+.5*popt[0]))**2+
            (err[0]*popt[2*i+2]*popt[2*i+3]*np.sqrt(2*np.pi)*.5)**2
            )
        if not errors:
            res[E[i]] = {
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
            res[E[i]] = {
                'm0': m0s[i],
                'mass': (m0s[i]+popt[1],err[1]),
                'Area' : (Area,Area_err),
                'Amp' : (popt[2*i+3],err[2*i+3]),
                'sig' : (popt[2*i+2],err[2*i+2]),
                'assym' : (popt[0],err[0]),
                'dm': (popt[1]*1e6/m0s[i],err[1]*1e6/m0s[i]),
                'fit': fit_type
                }

        if showElts:
            ax.axvline(m0s[i]+popt[1], color='r', alpha=.5);
            ax.annotate(E[i], (m0s[i]+popt[1], ax.get_ylim()[1]),(3,-1),textcoords='offset pixels', rotation=90, va='top',ha='left');
        if dofit:
            ax.plot(m[mask], Y, '--');
            ax.annotate("{:.2f}".format(Area), (m0s[i]+popt[1], popt[2+2*i]/2))
    if debug:
        print(popt)
    if dofit or debug:
        ax.plot(m[mask], fit(m[mask], *popt), 'r:');
    return res
    