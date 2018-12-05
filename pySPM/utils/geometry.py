class Point:
    def __init__(self, xy, y=None):
        if y is None:
            assert type(xy) in [list, tuple]
            assert len(xy)==2
            self.x = xy[0]
            self.y = xy[1]
        else:
            self.x = xy
            self.y = y
            
    def __add__(self, other):
        assert isinstance(other, Point)
        return Point(self.x+other.x, self.y+other.y)
    
    def __sub__(self, other):
        assert isinstance(other, Point)
        return Point(self.x-other.x, self.y-other.y)
    
    def __mul__(self, other):
        assert isinstance(other, Point)
        return Point(self.x*other.x, self.y*other.y)
        
    def __div__(self, other):
        assert isinstance(other, Point)
        return Point(self.x/other.x, self.y/other.y)

class Bbox:
    def __init__(self, *args, **kargs):
        if len(args)==0:
            self.left = kargs['left']
            self.right = kargs['right']
            self.top = kargs['top']
            self.bottom = kargs['bottom']
        elif len(args)==1:
            assert type(args[0]) is dict
            self.left = args[0]['left']
            self.right = args[0]['right']
            self.top = args[0]['top']
            self.bottom = args[0]['bottom']
        elif len(args)==3:
            assert type(args[0]) in [list, tuple]
            self.left = args[0][0]
            self.bottom = args[0][1]
            self.right = self.left + args[1]
            self.top = self.bottom + args[2]
        else:
            assert len(args)==4
            self.left = args[0]
            self.right = args[1]
            self.top = args[2]
            self.bottom = args[3]
    
    def __repr__(self):
        return "Bbox ({left},{bottom}) -> ({right},{top})".format(**self.__dict__)
        
    def is_overlapping(self, other):
        return (other.left<self.right) and (other.bottom<self.top) and (other.right>self.left) and (other.top>self.bottom)
        
    def overlap(self, other):
        if not self.is_overlapping(other):
            return Bbox(0,0,0,0)
        left = max(self.left, other.left)
        right = min(self.right, other.right)
        top = min(self.top, other.top)
        bottom = max(self.bottom, other.bottom)
        return Bbox(left, right, top , bottom)
    
    def show(self, ax=None, **kargs):
        from .plot import pixel2img
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        
        if ax is None:
            ax = plt.gca()
            
        pos = pixel2img((self.left, self.bottom), ax=ax)
        ur = pixel2img((self.right, self.top), ax=ax)
        width = ur[0]-pos[0]
        height = ur[1]-pos[1]
        fill = kargs.pop("fill", None)
        color = kargs.pop("color", 'r')
        ax.add_patch(mpl.patches.Rectangle(pos, width, height, fill=fill, color=color, **kargs))