from pySPM.SPM import SPM_image

class collection:
	"""Class to handle a collection of SPM images"""
	
	def __init__(self, sx=None,sy=None,unit='px', name='RawData'):
		"""
			Create a new collection.
			You should provide a size.
				sx for the x-axis
				sy for the y-axis
				and unit for the name of the dimensional units
		"""
		
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
		
	def show(self, ax, channels, cmap='hot'):
		for i,x in enumerate(channels):
			self[x].show(ax=ax[i],cmap=cmap)
			