def hot2val(rgb):
    if type(rgb) in [list, tuple]:
        r,g,b = rgb
    else:
        r = rgb[:,:,0]
        g = rgb[:,:,1]
        b = rgb[:,:,2]
    A = 0.365079
    B = 0.7460321   
    return A*(r-0.0416)/0.9584+(B-A)*g+(1-B)*b
    