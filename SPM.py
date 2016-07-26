import xml.etree.ElementTree as ET
import base64
import numpy as np
import struct
import os
import matplotlib.pyplot as plt
import scipy
import skimage
#import skimage.filters
import copy
from tqdm import tqdm

"""
Library to handle SPM data.
"""

def getSPM(filename, channel, corr=''):
	Fwd = SPM_image(filename, channel)
	Bwd = SPM_image(filename, channel, True)
	Fwd.pixels = (Fwd.pixels + Bwd.pixels)/2
	del Bwd
	if corr.lower() == 'slope':
		Fwd.correctSlope()
	elif corr.lower() == 'lines':
		Fwd.correctLines()
	return Fwd

class SPM_image:
	def __init__(self, filename, channel='Topography', backward=False,corr='none'):
		if not os.path.exists(filename): raise IOError('File Not Found')
		if filename[-4:]!='.xml': raise TypeError("Only xml files are handeled for the moment!")
		self.filename=filename
		tree = ET.parse(filename)
		self.root = tree.getroot()
		self.channel=channel
		direction=['forward','backward'][backward]
		self.direction = direction
		if self.root.tag == "{http://www.nanoscan.ch/SPM}scan":
			namespaces={'spm':"http://www.nanoscan.ch/SPM"}
			self.type = "Nanoscan"
			try:
				RAW = self.root.findall("spm:vector//spm:direction/spm:vector/spm:contents/spm:name[spm:v='%s']/../spm:channel//spm:contents/spm:name[spm:v='%s']/../spm:data/spm:v"%(["forward","backward"][backward], channel), namespaces)[0].text
			except:
				raise 'Channel {0} in {1} scan not found'.format(channel,direction)
				return
			size= [int(z.text) for z in self.root.findall("spm:vector/spm:contents/spm:size/spm:contents//spm:v", namespaces)]
			self.size = {'pixels':{
							'x':size[0],
							'y':size[1]
							},
						 'real':{
							'unit':self.root.findall(".//spm:area/spm:contents/spm:unit/spm:v",namespaces)[0].text,
							'x':float(self.root.findall(".//spm:area//spm:contents/spm:x/spm:v",namespaces)[0].text),
							'y':float(self.root.findall(".//spm:area//spm:contents/spm:y/spm:v",namespaces)[0].text),
							}}
			BIN = base64.b64decode(RAW)
			self.pixels = np.array(struct.unpack("<%if"%(size[0]*size[1]),BIN)).reshape(size)
		elif self.root.tag=="channel_list":# ToF-SIMS data
			self.type = "ToF-SIMS"
			self.channel = "Counts"
			x   = int(self.root.findall("./channel/axis[name='x']/count")[0].text)
			y   = int(self.root.findall("./channel/axis[name='y']/count")[0].text)
			RAW = self.root.findall("./channel/pixels")[0].text
			BIN = base64.b64decode(RAW)
			self.pixels = np.array(struct.unpack("<%if"%(x*y),BIN)).reshape(x,y)
			self.size = {'pixels':{
							'x':x,
							'y':y
							},
						 'real':{
							 'unit':'m',
							 'x':float(self.root.findall("./channel/axis[name='x']/variable/extent")[0].text),
							 'y':float(self.root.findall("./channel/axis[name='y']/variable/extent")[0].text)
							 }}
		if corr.lower() == 'slope':
			self.correctSlope()
		elif corr.lower() == 'lines':
			self.correctLines()
			
	def correctSlope(self):
		s=np.mean(self.pixels,axis=1)
		i=np.arange(len(s))
		fit=np.polyfit(i,s,1)
		self.pixels-=np.tile(np.polyval(fit,i).reshape(len(i),1),len(i))

	def correctLines(self):
		self.pixels-=np.tile(np.mean(self.pixels,axis=1).T,(self.pixels.shape[0],1)).T

	def dist_v2(self, pixel=False):
		if pixel:
			dx=1
			dy=1
		else:
			dx = self.size['real']['x']/self.size['pixels']['x']
			dy = self.size['real']['y']/self.size['pixels']['y']
		x2 = np.arange(self.size['pixels']['x'])
		x2 = (np.minimum(x2, self.size['pixels']['x']-x2) * dx)**2
		y2 = np.arange(self.size['pixels']['y'])
		y2 = (np.minimum(y2, self.size['pixels']['y'] - y2) * dy)**2
		X,Y = np.meshgrid(x2,y2)
		return np.sqrt(X+Y)

	def inv_calc_flat(self,d,l=0.1):
		work_image = self.pixels
		ny,nx = self.pixels.shape
		dx = self.size['real']['x']/self.size['pixels']['x']
		dy = self.size['real']['y']/self.size['pixels']['y']

		k=self.dist_v2()
		k[0,0]=1e-10

		tf = np.exp(-d*k)
		tf[0,0]=np.mean(tf)
		tf/=2
		tf *=1-np.exp(-d * k)
	
		recon_tf = np.ones(tf.shape) / (tf+l*np.ones(tf.shape) / np.conj(tf))
		tf *= recon_tf
		return np.real(np.fft.ifft2(np.fft.fft2(work_image)*recon_tf))

	def show(self, ax=None, sig = None, cmap=None, title=None):
		if ax==None:
			fig, ax = plt.subplots(1,1)
		if title==None:
			title="{0} - {1}".format(self.type,self.channel)
		unit=self.size['real']['unit']
		sunit='afpnum kMGTPE'
		if len(unit)==1: isunit=6
		elif unit[0] in sunit:
			unit=unit[1:]
			isunit=sunit.index(unit[0])
		W = self.size['real']['x']
		H = self.size['real']['y']
		fact=int(round(np.log(W)/np.log(10)/3))
		isunit+=fact
		W,H=W/10**(fact*3),H/10**(fact*3)
		if cmap==None:
			cmap='gray'
			if unit=='m' and self.channel == "Topography":
				cmap='hot'
		extent=(0,W,0,H)
		if sig == None:
			ax.imshow(self.pixels,cmap=cmap,extent=extent)
		else:
			std  = np.std(self.pixels)
			avg  = np.mean(self.pixels)
			vmin  = avg - sig * std
			vmax = avg + sig * std
			ax.imshow(self.pixels,cmap=cmap, vmin=vmin, vmax=vmax, extent=extent)
		if isunit!=6:
			ax.set_xlabel('x [{0}{1}]'.format(sunit[isunit],unit))
			ax.set_ylabel('y [{0}{1}]'.format(sunit[isunit],unit))
		else:
			ax.set_xlabel('x [{0}]'.format(unit))
			ax.set_ylabel('y [{0}]'.format(unit))
		if title != None:
			ax.set_title(title)

	def getProfile(self, x1,y1,x2,y2):
		d=np.sqrt((x2-x1)**2+(y2-y1)**2)
		x,y = np.linspace(x1,x2,int(d)+1),np.linspace(y1,y2,int(d)+1)
		return scipy.ndimage.map_coordinates(self.pixels,np.vstack((y,x)))

	def plotProfile(self, x1,y1,x2,y2, ax=None, col='b-',**kargs):
		if ax==None:
			fig, ax = plt.subplots(1,1)
		d  = np.sqrt((x2-x1)**2+(y2-y1)**2)
		dx = (x2-x1)*self.size['real']['x']/self.size['pixels']['x']
		dy = (y2-y1)*self.size['real']['y']/self.size['pixels']['y']
		rd = np.sqrt(dx**2+dy**2)
		l  = np.linspace(0,rd,int(d)+1)
		x,y = np.linspace(x1,x2,int(d)+1),np.linspace(y1,y2,int(d)+1)
		z=scipy.ndimage.map_coordinates(self.pixels,np.vstack((y,x)))
		l=ax.plot(l,z,col,**kargs)
		ax.set_xlabel("Distance [{0}]".format(self.size['real']['unit']))
		ax.set_ylabel("Height [{0}]".format(self.size['real']['unit']))
		return l

	def getBinThreshold(self, percent, high=True, adaptive=False, binary=True):
		if adaptive:
			if binary:
				return skimage.filters.threshold_adaptive(self.pixels, percent)
			return np.ones(self.pixels.shape)*skimage.filters.threshold_adaptive(self.pixels, percent)
		mi=np.min(self.pixels)
		norm = (self.pixels+mi)/(np.max(self.pixels)-mi)
		if high:
			r=norm > percent
		else:
			r=norm < percent
		if binary: return r
		return np.ones(self.pixels.shape)*r

	def getShadowMask(self, angle, BIN=None, pb=False):
		if type(BIN)!=type(None): 
			BIN=BIN*1.0
		slope = np.tan(np.radians(angle))
		neg   = False
		if slope<0:
			neg   = True
			slope = -slope
			topo  = np.fliplr(self.pixels)
			BIN   = np.fliplr(BIN)
		else:
			topo = self.pixels
		x	= np.linspace(0,self.size['real']['x'],self.pixels.shape[1])
		mask = np.zeros(self.pixels.shape)
		AFM_bin_shadow = np.zeros(self.pixels.shape)
		Y=xrange(self.pixels.shape[0])
		if pb: Y=tqdm(Y)
		for yi in Y:
			for xi in xrange(self.pixels.shape[1]):
				cut   = self.pixels.shape[1]-2
				y_ray = slope*(x-x[xi]) + topo[yi,xi]
				while cut>xi and y_ray[cut]>topo[yi,cut]:
					cut-=1
				if xi==cut:
					AFM_bin_shadow[yi,xi]=BIN[yi,xi]
					continue
				# Cut has been found
				if type(BIN)!=type(None):
					x1 = x[cut]
					x2 = x[cut+1]
					y1 = topo[yi,cut]
					y2 = topo[yi,cut+1]
					x0 = x[xi]
					y0 = topo[yi,xi]
					if y2==y1:
						x_cut = (y1+slope*x0-y0)/slope
						y_cut = y1
					else:
						numerator   = x1/(x2-x1)+(y0-slope*x0-y1)/(y2-y1)
						denominator =  1/(x2-x1)-slope/(y2-y1)
						x_cut = numerator / denominator
						y_cut = slope*(x_cut-x0)+y0
					if x_cut>=x1 and x_cut<=x2:
						y1 = BIN[yi,cut]
						y2 = BIN[yi,cut+1]
						yint = (((y2-y1)/(x2-x1))*(x_cut-x1))+y1
					else: yint = BIN[yi,xi]
					AFM_bin_shadow[yi,xi]=yint
				mask[yi,xi]=1
		if neg:
			mask = np.fliplr(mask)
			AFM_bin_shadow = np.fliplr(AFM_bin_shadow)
		if type(BIN)!=type(None):
			return (mask, AFM_bin_shadow)
		return mask
	
	def adjust_position(self, fixed):
		""" Shift the current pixels to match a fixed image """
		adj = copy.deepcopy(self)
		cor = np.fft.fft2(fixed.pixels)
		cor = np.abs( np.fft.ifft2( np.conj(cor) * np.fft.fft2(self.pixels) ) )
		cor = cor / fixed.pixels.size
		ypeak, xpeak = np.unravel_index(cor.argmax(), cor.shape)
		shift = [-(ypeak-1), -(xpeak-1)]
		adj.pixels = np.roll(self.pixels,shift[0],axis=0)
		adj.pixels = np.roll( adj.pixels,shift[1],axis=1)
		return adj
		

	
def imshow_sig(img,sig=1, ax=None, **kargs):
	if ax==None:
		fig, ax = plt.subplots(1,1)
	std  = np.std(img)
	avg  = np.mean(img)
	vmin  = avg - sig * std
	vmax = avg + sig * std
	ax.imshow(img,vmin=vmin, vmax=vmax, **kargs)
	
def adjust_position(fixed,to_adjust,shift=False):
	""" Shift the current pixels to match a fixed image """
	adj = copy.deepcopy(to_adjust)
	cor = np.fft.fft2(fixed)
	cor = np.abs( np.fft.ifft2( np.conj(cor) * np.fft.fft2(to_adjust) ) )
	cor = cor / to_adjust.size
	ypeak, xpeak = np.unravel_index(cor.argmax(), cor.shape)
	shift = [-(ypeak-1), -(xpeak-1)]
	adj = np.roll(to_adjust,shift[0],axis=0)
	adj = np.roll(adj   ,shift[1],axis=1)
	if shift:	return adj,shift
	return adj

def tukeywin(window_length, alpha=0.5):
	'''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
	that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
	at \alpha = 0 it becomes a Hann window.
	
	We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
	output
	
	Reference
	---------
	http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html
	
	'''
	# Special cases
	if alpha <= 0:
		return np.ones(window_length) #rectangular window
	elif alpha >= 1:
		return np.hanning(window_length)
	
	# Normal case
	x = np.linspace(0, 1, window_length)
	w = np.ones(x.shape)
	
	# first condition 0 <= x < alpha/2
	first_condition = x<alpha/2
	w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))
	
	# second condition already taken care of
	
	# third condition 1 - alpha / 2 <= x <= 1
	third_condition = x>=(1 - alpha/2)
	w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2))) 
	
	return w

if __name__ == "__main__":
	Path = "C:/Users/ols/Dropbox/ToF_SIMS"
	AFM = SPM_image("{0}/CyI5b_0006_ns.xml".format(Path),corr='slope')
	CC2 = SPM_image("{0}/CyI5b_PCBM_CC2.xml".format(Path))
	AFM.show()
