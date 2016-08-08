import numpy as np
import struct
import os
import re

class ToF:
	def __init__(self, path):
		self.Peaks = {}
		self.RPeaks = {}
		self.Path = os.path.dirname(path)
		self.Basename = os.path.basename(path)
		self.size = None
		for x in os.listdir(self.Path):
			if x[:len(self.Basename)]==self.Basename:
				if x[-6:]==".BIF3D":
					s=self.Basename+r' \(([0-9]+)\)(?: - (.*?))?\.BIF3D'
					r=re.search(s,x)
					ID=int(r.group(1))
					self.Peaks[ID]=r.group(2)
					if r.group(2)!=None:
						self.RPeaks[r.group(2)]=ID
		self.data={}
		
	def getIDs(self, channels):
		assert type(channels) in [str,int,list,tuple]
		if type(channels)==int: return channels
		if type(channels)==str:
			assert channels in self.RPeaks.keys()
			return self.RPeaks[channels]
		if type(channels)==list or type(channels)==tuple:
			ID=[]
			for x in channels:
				ID.append(self.getIDs(x))
			return ID
				
	def getChannel(self, channel):
		k=self.getIDs(channel)
		if k in self.data: return self.data[k]
		v=self.Peaks[k]
		if v!=None:
			path='{path}{sep}{basename} ({ID}) - {Peak}.BIF3D'
		else:
			path='{path}{sep}{basename} ({ID}).BIF3D'
		f=open(path.format(path=self.Path,sep=os.sep,basename=self.Basename,ID=k,Peak=v),'rb')
		A=f.read()
		f.close()
		size=struct.unpack("II",A[32:40]) 
		if self.size == None:
			self.size = size
		else:
			assert self.size==size
		return np.array(struct.unpack("{0}d".format(self.size[0]*self.size[1]),A[640:])).reshape(self.size)
	
	def loadChannel(self,channel):
		k = self.getIDs(channel)
		if not k in self.data:
			self.data[k]=self.getChannel(channel)
			
	def loadChannels(self, channels):
		if not type(channels) in [list,tuple]:
			self.loadChannel(channels)
		else:
			for x in channels:
				self.loadChannel(x)
	
	def getChannels(self, *channels):
		k=self.getIDs(channels)
		self.loadChannels(k)
		if type(k)==int: return self.data[k]
		D=np.zeros(self.size)
		for x in k:
			D+=self.data[x]
		return D
	
	def listChannels(self):
		return self.Peaks

	def show(self, ax, *channels):
		D=self.getChannels(*channels)
		ax.imshow(D)
		if len(channels)>1:
			ax.set_title(','.join([str(z) for z in channels]))
		else:
			ax.set_title(str(channels[0]))
		
