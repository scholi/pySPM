def dec_debug(debug):
    if debug>0:
        return debug-1
    if debug<0:
        return debug+1
    return False

def do_debug(debug):
    return debug==-1 or debug>0

