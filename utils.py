Elts={
    'B':[10,11],
    'C':12,
    'N':14,
    'O':16,
    'OH':17,
    'F':19,
    'Na':23,
    'Mg':24,
    'CN':26,
    'Al':27,
    'Si':28,
    'P':31,
    'S':32,
    'SH':33,
    'Cl':[35,37],
    'TiO':[64,65],
    'TiO2':[80,81],
    'Cs':133,
    'Au':197
    }

def show_table(t):
    from IPython.core.display import display, HTML
    display(HTML(html_table(t)))
    
def html_table(t, header=False):
    S = "<table>"
    if header:
        S += "<tr>"+"".join(["<th>{0}</th>".format(x) for x in t[0]])+"</tr>"
        t = t[1:]
    for x in t:
        S += "<tr>"+"".join(["<td>{0}</td>".format(y) for y in x])+"</tr>"
    S += "</table>"
    return S
    
def aa_table(t, header=False):
    """
    print a list of list in a nice ascii-art
    """
    Ncols=len(t[0])
    Lcol=[0]*Ncols
    for x in t:
        for i in range(Ncols):
            Lcol[i]=max(Lcol[i],len(repr(x[i])))
    if header:
        print("  ".join([u"{: <"+str(Lcol[i]+4)+"}" for i in range(Ncols)]).format(*t[0]))
        print("="*sum(Lcol))
        t=t[1:]
    for j,x in enumerate(t):
        print("  ".join([u"{:"+['.','_'][j%2]+"<"+str(Lcol[i]+4)+"}" for i in range(Ncols)]).format(*x))
    

def dict_update(d, u):
    import collections
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = dict_update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d
    
def fact(x):
    import math
    f=[]
    i=2
    while True:
        while x%i == 0:
            f.append(i)
            x /= i
        i += 1
        if x == 1:
            return f