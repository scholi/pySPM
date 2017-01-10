import numpy as np
import struct
import re
import pySPM

class Bruker:
	def __init__(self, path):
		self.path = path
		self.f = open(self.path, 'rb')
		self.Layers=[]
		self.Scanners=[]
		mode=''
		while True:
			l = self.f.readline().rstrip().replace(b'\\',b"")
			if l==b'*Ciao image list':
				self.Layers.append({})
				mode='Image'
			elif l==b'*Scanner list':
				self.Scanners.append({})
				mode='Scanner'
			ll=l.split(b": ")
			if len(ll)>1:
				if mode=='Image':
					self.Layers[-1][ll[0]]=ll[1:]
				elif mode=='Scanner':
					self.Scanners[-1][ll[0]]=ll[1:]
			if l==b"*File list end":
				break

	def getRawLayer(self,i):
			off = int(self.Layers[i][b'Data offset'][0])
			cols = int(self.Layers[i][b'Number of lines'][0])
			rows = int(self.Layers[i][b'Samps/line'][0])
			self.f.seek(off)
			N=cols*rows
			return np.array(struct.unpack("<"+str(N)+"h",self.f.read(N*2)),dtype='float64').reshape((cols,rows))
	
	def loadImage(self, channel, backward=False):	
			for i in range(len(self.Layers)):
				ln=self.Layers[i][b'@2:Image Data'][0].decode('utf8')
				r=re.match(r'([^ ]+) \[([^\]]*)\] "([^"]*)"',ln).groups()
				if r[2]==channel:
					bck=False
					if self.Layers[i][b'Line Direction'][0]==b'Retrace': bck=True
					if bck==backward:
						r = re.match(r'[A-Z]+ \[([^\]]+)\] \([0-9\.]+ [^\)]+\) ([0-9\.]+) V',self.Layers[i][b'@2:Z scale'][0].decode('utf8')).groups()
						Scale=float(r[1])/65536.0
						r=self.Scanners[0][b'@'+r[0].encode('utf8')][0].split()
						Scale*=float(r[1])
						Data=self.getRawLayer(i)*Scale
						
						zscale=r[2].split(b'/')[0]
						s=self.Layers[i][b'Scan Size'][0].split()
						if s[2][0]==126: s[2]=b'u'+s[2][1:]
						Size={'x':float(s[0]),'y':float(s[1]),'unit':s[2].decode('utf8')}
						I=pySPM.SPM_image(channel=channel,backward=backward,BIN=Data,real=Size,_type='Bruker AFM',zscale=zscale.decode('utf8'))
						return I
