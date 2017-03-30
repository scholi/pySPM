from pySPM.SPM import SPM_image
import pySPM
import copy
import matplotlib.pyplot as plt
import re
import numpy as np

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
	
class collection:
	"""Class to handle a collection of SPM images"""
	
	def __init__(self, sx=None,sy=None,unit='px', name='RawData', cls=None):
		"""
			Create a new collection.
			You should provide a size.
				sx for the x-axis
				sy for the y-axis
				and unit for the name of the dimensional units
		"""
		
		if isinstance(cls, collection):
			self.size = cls.size
			self.name = cls.name
			
		if sy is None and sx is not None:
			sy=sx
			
		self.name=name
		self.CH = {}
		self.size = {'x':sx,'y':sy,'unit':unit}
	
	def add(self, Img, name):
		if len(self.CH)==0 and self.size['unit']=='px':
			self.size['x']=len(Img[0])
			self.size['y']=len(Img)
		if name in self.CH:
			raise KeyError('The channel {} is already present in the collection. Please delete it before')
			return
		self.CH[name]=Img
	
	def __getitem__(self, key):
		if key not in self.CH: return None
		return SPM_image(_type=self.name,BIN=self.CH[key],real=self.size,channel=key)
	
	def __setitem__(self, key, value):
		self.add(value, key)
		
	def show(self, ax=None,channels=None, cmap='hot', **kargs):
		if channels is None:
			channels = list(self.CH.keys())
		N=len(channels)
		channels.sort(key=natural_keys)
		if ax is None:
			if N==4:
				fig, ax = plt.subplots(2,2,figsize=(20,self[channels[0]].pixels.shape[0]*20/self[channels[0]].pixels.shape[1]))
			else:
				fig, ax = plt.subplots((N-1)//3+1,min(3,N),figsize=(20,((N-1)//3+1)*20/min(3,N)))
		for i,x in enumerate(channels):
			self[x].show(ax=ax.ravel()[i],cmap=cmap,**kargs)
		plt.tight_layout()
	
	def getMultiVariate(self, channels=None):
		import pandas as pd
		if channels is None:
			channels = self.CH.keys()
		data = pd.DataFrame({k:self.CH[k].ravel() for k in channels})
		return data
	
	def overlay(self, chs, cols=[[0,.5,.3],[1,0,0]],ax=None, sig = None, vmin=None, vmax=None):
		assert len(chs)>=2
		assert len(cols)==len(chs)
		Data = [pySPM.SPM.Normalize(self[ch].pixels,sig=sig,vmin=vmin,vmax=vmax) for ch in chs]
		C    = [np.stack([Data[i]*cols[i][j] for j in range(3)],axis=2) for i in range(len(chs))]
		if ax is None:
			ax = plt.gca()
		ax.imshow(np.sum(C,axis=0))
		return C

		
	def StitchCorrection(self, channel, stitches):
		N = copy.deepcopy(self)
		del N.CH
		N.CH = {}
		size = list(self.CH.values())[0]
		S=np.zeros((size[0]/stitches[0],size[1]/stitches[1]))
		for i in range(stitches[0]):
			for j in range(stitches[1]):
				S+=self.CH[channel][128*i:128*(i+1),128*j:128*(j+1)]
		S[S==0]=1
		for x in self.CH:
			F = np.zeros(size)
			for i in range(stitches[0]):
				for j in range(stitches[1]):
					F[128*i:128*(i+1),128*j:128*(j+1)]=self.CH[x][128*i:128*(i+1),128*j:128*(j+1)]/T
			N.add(F,x)
		return N

def Tsign(p1,p2,p3):
	return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

def PointInTriangle (pt,v1,v2,v3):
    b1 = Tsign(pt, v1, v2) < 0
    b2 = Tsign(pt, v2, v3) < 0
    b3 = Tsign(pt, v3, v1) < 0
    return (b1 == b2) * (b2 == b3)

def overlayTriangle(chs, cols,ax=None,size=512):
	assert len(chs)==3
	assert len(cols)==3
	if ax is None:
		ax = plt.gca()
	p=.9
	r=.8
	RGB=[np.zeros((size,size)) for i in range(3)]
	D=2*r*p*np.sin(np.radians(120))

	x=np.linspace(-.7,1.1,size)
	y=np.linspace(-1,1,size)
	X,Y=np.meshgrid(x,y)
	centers=[(r*p*np.cos(np.radians(120*i)),r*p*np.sin(np.radians(120*i))) for i in range(3)]
	dist=[np.sqrt((X-centers[i][0])**2+(Y-centers[i][1])**2) for i in range(3)]

	for j in range(3):
		RGB[j]=sum([cols[i][j]*np.maximum((D-dist[i])/D,0) for i in range(3)])
		ax.annotate(chs[j],(r*np.cos(np.radians(120*j)),r*np.sin(np.radians(120*j))),color='w',fontsize=20,va="center",ha="center")
		RGB[j][PointInTriangle([X,Y],*centers)==0]=0
	A=np.stack(RGB,axis=2)
	ax.imshow(A,extent=[x[0],x[-1],y[0],y[-1]])
	ax.set_xticks([])
	ax.set_yticks([])
