from pySPM.SPM import SPM_image
import pandas as pd
import copy
import matplotlib.pyplot as plt

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
		if ax is None:
			if N==4:
				fig, ax = plt.subplots(2,2,figsize=(20,self[channels[0]].pixels.shape[0]*20/self[channels[0]].pixels.shape[1]))
			else:
				fig, ax = plt.subplots((N-1)//3+1,min(3,N),figsize=(20,((N-1)//3+1)*20/min(3,N)))
		for i,x in enumerate(channels):
			self[x].show(ax=ax.ravel()[i],cmap=cmap,**kargs)
		plt.tight_layout()
	
	def getMultiVariate(self, channels=None):
		if channels is None:
			channels = self.CH.keys()
		data = pd.DataFrame({k:self.CH[k].ravel() for k in channels})
		return data
	
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