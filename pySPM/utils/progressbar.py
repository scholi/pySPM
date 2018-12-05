import sys
class Progressbar:
    def __init__(self, iterator=None, total=None, length=80):
        self.it = iterator
        if total is None and iterator is not None:
            try:
                self.total = len(iterator)
            except Exception as e:
                self.total = None
        self.length = length
        
    def __iter__(self):
        self.elapsed = 0
        self.update()
        for x in self.it:
            yield x
            self.elapsed += 1
            self.update()
            
    def __repr__(self):
        if self.total is not None:
            p = int(self.length*self.elapsed/self.total)
            tot = self.total
        else:
            p = 0
            tot = '???'
        return "|"+"="*p+" "*(self.length-p)+"| ({}/{})".format(self.elapsed, tot)
    
    def update(self):
        sys.stderr.write(self.__repr__()+'\r')