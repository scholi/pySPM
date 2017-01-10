#*-* encoding: utf-8 *-*

import xml.etree.ElementTree as ET
import base64
import numpy as np
import struct
import os
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage
import scipy.optimize
import skimage
import skimage.exposure
import scipy.interpolate
from skimage import transform as tf
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

def funit(v,u=None,iMag=True):
	if u==None:
		u=v['unit']
		v=v['value']
	import math
	shift=int(math.floor(math.log10(v)/3.0))
	mag=u'afpnum1kMGTPE'
	imag=mag.index('1')
	unit=u
	if len(u)>1 and u[0] in mag and u[1:] and iMag:
		imag=mag.index(u[0])
		unit=u[1:]
	value = v/10**(3*shift)
	imag += shift
	if imag<0:
		value *= 10**(imag*3)
		imag=0
	elif imag>=len(mag):
		value *= 10**((imag-len(mag)+1)*3)
		imag = len(mag)-1
	m=mag[imag]
	if m=='1': m=''
	return {'value':value,'unit':u'{mag}{unit}'.format(mag=m,unit=unit)}

def getCurve(filename,channel='Normal Deflection',backward=False):
	"""
	function to retrieve data which are not in the form of images. This is typically used for 1D channel where the normal deflection is recorded while z is swept.
	"""
	tree = ET.parse(filename)
	root = tree.getroot()
	namespace={'spm':'http://www.nanoscan.ch/SPM'}
	RAW=root.findall("spm:vector/spm:contents/spm:direction/spm:vector/spm:contents/spm:name[spm:v='{direction}']/../spm:channel/spm:vector/spm:contents/spm:name[spm:v='{channel}']/../spm:data/spm:v".format(direction=['forward','backward'][backward],channel=channel),namespace)[0].text
	start=float(root.findall("spm:vector/spm:contents/spm:axis/spm:vector/spm:contents/spm:start/spm:vector/spm:v",namespace)[0].text)
	stop=float(root.findall("spm:vector/spm:contents/spm:axis/spm:vector/spm:contents/spm:stop/spm:vector/spm:v",namespace)[0].text)
	unit=root.findall("spm:vector/spm:contents/spm:axis/spm:vector/spm:contents/spm:unit/spm:v",namespace)[0].text
	BIN=base64.b64decode(RAW)
	N=len(BIN)
	vals=np.array(struct.unpack("<"+str(N//4)+"f",BIN))
	x = np.linspace(start,stop,len(vals))
	return x,vals

class SPM_image:
	"""
	Main class to handle SPM images
	"""
	def __init__(self, filename=None, channel='Topography', backward=False,corr='none',BIN=None,real=None,_type=None,zscale='?'):
		if filename is None and not BIN is None:
			self.channel = channel
			self.direction = 'Unknown'
			self.size={'pixels':{'x':BIN.shape[1],'y':BIN.shape[0]}}
			if not real is None:
				self.size['real']=real
			else: self.size['real']={'unit':'pixels','x':BIN.shape[1],'y':BIN.shape[0]}
			self.size['recorded']={'pixels':self.size['pixels'],'real':self.size['real']}
			self.pixels = BIN
			if not _type is None: self.type = _type
			else: self.type = 'RawData'
			self.zscale=zscale
			return
		if not os.path.exists(filename): raise IOError('File "{0}" Not Found'.format(filename))
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
			self.scanSpeed = {'value':float(self.root.findall("spm:vector//spm:direction/spm:vector/spm:contents/spm:name[spm:v='%s']/../spm:point_interval/spm:v"%(["forward","backward"][backward]), namespaces)[0].text) * size[0],
					'unit': self.root.findall("spm:vector//spm:direction/spm:vector/spm:contents/spm:name[spm:v='%s']/../spm:point_interval_unit/spm:v"%(["forward","backward"][backward]), namespaces)[0].text}
			fbPath = "spm:vector/spm:contents/spm:instrumental_parameters/spm:contents/spm:z_control/spm:contents"
			self.feedback = {'channel':self.root.findall('{0}/spm:z_feedback_channel/spm:v'.format(fbPath), namespaces)[0].text}
			self.feedback['P']={'value':float(self.root.findall('{0}/spm:proportional_z_gain/spm:v'.format(fbPath), namespaces)[0].text),
					'unit':self.root.findall('{0}/spm:proportional_z_gain_unit/spm:v'.format(fbPath), namespaces)[0].text}
			self.feedback['I']={'value':float(self.root.findall('{0}/spm:integral_z_time/spm:v'.format(fbPath), namespaces)[0].text),
					'unit':self.root.findall('{0}/spm:integral_z_time_unit/spm:v'.format(fbPath), namespaces)[0].text}
			if self.feedback['channel']=='df':
				self.feedback['channel']=u'Δf'
			uval = float(self.root.findall(".//spm:area//spm:contents/spm:size/spm:contents/spm:fast_axis/spm:v",namespaces)[0].text)
			udispu = self.root.findall(".//spm:area//spm:contents/spm:display_unit/spm:v",namespaces)[0].text
			udisps = float(self.root.findall(".//spm:area/spm:contents/spm:display_scale/spm:v",namespaces)[0].text)
			uname = self.root.findall(".//spm:area/spm:contents/spm:unit/spm:v",namespaces)[0].text
			x=funit(uval*udisps, udispu)
			uval = float(self.root.findall(".//spm:area//spm:contents/spm:size/spm:contents/spm:slow_axis/spm:v",namespaces)[0].text)
			y=funit(uval*udisps,udispu)
			self.size = {'pixels':{
							'x':size[0],
							'y':size[1]
							},
						 'real':{
							'unit':x['unit'],
							'x':x['value'],
							'y':y['value'],
							}}
			BIN = base64.b64decode(RAW)
			recorded_size = len(BIN)/4
			self.size['recorded']={'pixels':{'y':int(recorded_size/size[0]),'x':size[0]}}
			self.size['recorded']['real']={'x':self.size['real']['x'],
				'y':self.size['real']['y']*self.size['recorded']['pixels']['y']/float(self.size['pixels']['y'])}
			self.pixels = np.array(struct.unpack("<%if"%(recorded_size),BIN)).reshape((self.size['recorded']['pixels']['y'],self.size['recorded']['pixels']['x']))
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
							 'y':float(self.root.findall("./channel/axis[name='y']/variable/extent")[0].text)},
						 'recorded':{'real':{
							 'unit':'m',
							 'x':float(self.root.findall("./channel/axis[name='x']/variable/extent")[0].text),
							 'y':float(self.root.findall("./channel/axis[name='y']/variable/extent")[0].text)
							 }}}
		if corr.lower() == 'slope':
			self.correctSlope()
		elif corr.lower() == 'lines':
			self.correctLines()
		elif corr.lower() == 'plane':
			self.correctPlane()

	def Offset(self, profiles, width=1, ax=None, inline=True):
		offset=np.zeros(self.pixels.shape[0])
		counts=np.zeros(self.pixels.shape[0])
		for p in profiles:
			y,D = self.getRowProfile(*p,width=width,ax=ax)
			counts[y]+=1
			offset[y[1:]]+=np.diff(D)
		counts[counts==0] = 1
		offset = offset/counts
		offset = np.cumsum(offset)
		offset = offset.reshape((self.pixels.shape[1],1))
		if inline:
			self.pixels=self.pixels-np.flipud(np.repeat(offset,self.pixels.shape[1],axis=1))
			return True
		else:
			C = copy.deepcopy(self)
			C.pixels = self.pixels - np.flipud(np.repeat(offset,self.pixels.shape[1],axis=1))
			return C

	def getRowProfile(self,x1,y1,x2,y2,width=1,ax=None):
		if y2<y1:
				x1,y1,x2,y2=x2,y2,x1,y1

		if ax!=None: ax.plot((x1,x2),(y1,y2),'w-')
		x = np.arange(self.pixels.shape[1])
		y = np.arange(self.pixels.shape[0])
		I=scipy.interpolate.interp2d(x,y,np.flipud(self.pixels))
		
		Y=np.arange(y1,y2+1)
		V=np.zeros(len(Y))
		for w in np.arange(width):
				xl=np.linspace(x1-(width-1)/2.+w,x2-(width-1)/2.+w,len(Y))
				for i in range(len(Y)):
						Z=I(xl[i],Y[i])
						V[i]+=Z
		return Y,V/width

	def getSummary(self):
		x=funit(self.size['real']['x'],self.size['real']['unit'])
		y=funit(self.size['real']['y'],self.size['real']['unit'])
		P=funit(self.feedback['P'])
		I=funit(self.feedback['I'])
		return u"""Feedback: {feedback[channel]} : P:{P[value]}{P[unit]} : I:{I[value]}{I[unit]}
Size: {size[pixels][x]}×{size[pixels][y]} pixels = {x[value]:.3} {x[unit]}×{y[value]:.3} {y[unit]}
Scan Speed: {scanSpeed[value]}{scanSpeed[unit]}/line""".format(x=x,y=y,P=P,I=I,feedback=self.feedback,size=self.size,scanSpeed=self.scanSpeed)

	def correctMedianDiff(self):
		N=self.pixels
		N2=np.vstack([N[1:,:],N[-1:,:]])-N # Difference of the pixel between two consecutive row
		C=np.cumsum(np.median(N2,axis=1)) # Take the median of the difference and cumsum them
		D=np.tile(C,(N.shape[0],1)).T	 # Extend the vector to a matrix (row copy)
		self.pixels = N-D
					
	def correctSlope(self):
		s=np.mean(self.pixels,axis=1)
		i=np.arange(len(s))
		fit=np.polyfit(i,s,1)
		self.pixels-=np.tile(np.polyval(fit,i).reshape(len(i),1),len(i))

	def correctPlane(self):
		x = np.arange(self.pixels.shape[1])
		y = np.arange(self.pixels.shape[0])
		X,Y = np.meshgrid(x,y)
		Z = self.pixels
		A = np.column_stack((np.ones(Z.ravel().size),X.ravel(),Y.ravel()))
		c, resid, rank, sigma = np.linalg.lstsq(A,Z.ravel())
		self.pixels -= c[0]*np.ones(self.pixels.shape) + c[1] * X + c[2] * Y

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

	def getExtent(self):
		W = self.size['recorded']['real']['x']
		H = self.size['recorded']['real']['y']
		return (0,W,0,H)

	def show(self, ax=None, sig = None, cmap=None, title=None, adaptive=False, dmin=0, dmax=0,pixels=False,flip=False,**kargs):
		if ax==None:
			fig, ax = plt.subplots(1,1)
		if title==None:
			title=u"{0} - {1}".format(self.type,self.channel)
		unit=self.size['real']['unit']
		sunit='afpnum kMGTPE'
		if len(unit)==1 or unit in ['pixels']:
			isunit=6
		elif unit[0] in sunit:
			isunit=sunit.find(unit[0])
			unit=unit[1:]
		else: isunit=6
		W = self.size['recorded']['real']['x']
		H = self.size['recorded']['real']['y']
		fact=int(np.floor(np.log(W)/np.log(10)/3))
		isunit+=fact
		W,H=W/10**(fact*3),H/10**(fact*3)
		if cmap==None:
			cmap='gray'
			if unit=='m' and self.channel == "Topography":
				cmap='hot'
		extent=(0,W,0,H)
		mi,ma = np.min(self.pixels), np.max(self.pixels)
		if adaptive:
			img = np.asarray(256*(self.pixels-mi)/(ma-mi),dtype=np.uint8)
			mi,ma=0,1
			img = skimage.exposure.equalize_adapthist(img, clip_limit=0.03)
		else:
			img=self.pixels
		if sig == None:
			vmin=mi+dmin
			vmax=ma+dmax
		else:
			std  = np.std(img)
			avg  = np.mean(img)
			vmin  = avg - sig * std
			vmax = avg + sig * std
		if not pixels:
			xp = np.linspace(0,self.pixels.shape[1],11)
			xr = np.linspace(0,W,11)
			ax.set_xticks(xp)
			ax.set_xticklabels([str(round(z,2)) for z in xr])
			yp = np.linspace(0,self.pixels.shape[0],11)
			yr = np.linspace(0,H,11)
			ax.set_yticks(yp)
			ax.set_yticklabels([str(round(z,2)) for z in yr])
		if not flip:
			ax.imshow(np.flipud(img),cmap=cmap,vmin=vmin,vmax=vmax,**kargs)
		else:
			ax.imshow(img,cmap=cmap,vmin=vmin,vmax=vmax,**kargs)
			
		if not pixels:
			if isunit!=6:
				u = sunit[isunit]
				if u=='u':
					u='$\\mu$'
				ax.set_xlabel(u'x [{0}{1}]'.format(u,unit))
				ax.set_ylabel(u'y [{0}{1}]'.format(u,unit))
			else:
				ax.set_xlabel(u'x [{0}]'.format(unit))
				ax.set_ylabel(u'y [{0}]'.format(unit))
			if title != None:
				ax.set_title(title)

	def getProfile(self, x1,y1,x2,y2, img=None, imgColor='w-'):
		if not img is None:
			img.plot([x1,x2],[y1,y2],imgColor)
		d=np.sqrt((x2-x1)**2+(y2-y1)**2)
		x,y = np.linspace(x1,x2,int(d)+1),np.linspace(y1,y2,int(d)+1)
		return np.linspace(0,d,int(d)+1),scipy.ndimage.map_coordinates(np.flipud(self.pixels),np.vstack((y,x)))

	def plotProfile(self, x1,y1,x2,y2, ax=None, col='b-', pixels=False,img=None,imgColor='w-',**kargs):
		if ax==None:
			fig, ax = plt.subplots(1,1)
		d  = np.sqrt((x2-x1)**2+(y2-y1)**2)
		dx = (x2-x1)*self.size['real']['x']/self.size['pixels']['x']
		dy = (y2-y1)*self.size['real']['y']/self.size['pixels']['y']
		if not img is None:
			img.plot([x1,x2],[y1,y2],imgColor)
		if pixels:
			rd=d
		else:
			rd = np.sqrt(dx**2+dy**2)
		l  = np.linspace(0,rd,int(d)+1)
		x,y = np.linspace(x1,x2,int(d)+1),np.linspace(y1,y2,int(d)+1)
		z=scipy.ndimage.map_coordinates(np.flipud(self.pixels),np.vstack((y,x)))
		p=ax.plot(l,z,col,**kargs)
		ax.set_xlabel("Distance [{0}]".format(self.size['real']['unit']))
		try:
			ax.set_ylabel("Height [{0}]".format(self.zscale))
		except:
			pass
		return {'plot':p,'l':l,'z':z}

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

	def align(self, tform, cut=True):
		New = copy.deepcopy(self)
		New.pixels = tf.warp(self.pixels, tform, preserve_range=True)
		if not cut:
			return New
		cut=[0,0]+list(self.pixels.shape)
		if tform.translation[0]>=0: cut[2]-=tform.translation[0]
		elif tform.translation[0]<0: cut[0]-=tform.translation[0]
		if tform.translation[1]>=0: cut[1]+=tform.translation[1]
		elif tform.translation[1]<0: cut[3]+=tform.translation[1]
		cut = [int(x) for x in cut]
		New.cut(cut,inplace=True)
		return New, cut

	def getFFT(self):
			return np.fft.fftshift(np.fft.fft2(self.pixels))

	def getRmask(self):
			x = np.linspace(-1,1,self.pixels.shape[1])
			y = np.linspace(-1,1,self.pixels.shape[0])
			X,Y=np.meshgrid(x,y)
			R=np.sqrt(X**2+Y**2)
			return R

	def corrFit2d(self, nx=2,ny=1):
		r,z = fit2d(self.pixels,nx,ny)
		self.pixels -= z

	def filterLowPass(self, p, inline=True):
			F=self.getFFT()
			mask=self.getRmask()<p
			if inline:
				self.pixels=np.real(np.fft.ifft2(np.fft.fftshift(F*mask)))
			else:
				C = copy.deepcopy(self)
				C.pixels=np.real(np.fft.ifft2(np.fft.fftshift(F*mask)))
				return C

	def ResizeInfos(self):
		self.size['real']['x']*=self.pixels.shape[1]/self.size['pixels']['x']
		self.size['real']['y']*=self.pixels.shape[0]/self.size['pixels']['y']
		self.size['recorded']['real']['x']*=self.pixels.shape[1]/self.size['pixels']['x']
		self.size['recorded']['real']['y']*=self.pixels.shape[0]/self.size['pixels']['y']
		self.size['recorded']['pixels']['x']*=self.pixels.shape[1]/self.size['pixels']['x']
		self.size['recorded']['pixels']['y']*=self.pixels.shape[0]/self.size['pixels']['y']
		self.size['real']['y']*=self.pixels.shape[0]/self.size['pixels']['y']
		self.size['pixels']['x']=self.pixels.shape[1]
		self.size['pixels']['y']=self.pixels.shape[0]

	def filterScarsRemoval(self, thresh=.5,inline=True):
		if not inline:
			C = copy.deepcopy(self)
		else:
			C = self
		for y in range(1,self.pixels.shape[0]-1):
			b=self.pixels[y-1,:]
			c=self.pixels[y,:]
			a=self.pixels[y+1,:]
			mask=np.abs(b-a)<.5*(np.abs(c-a))
			C.pixels[y,mask]=b[mask]
		return C
		
	def cut(self, c, inplace=False):
		if not inplace:
			new = copy.deepcopy(self)
			new.pixels = cut(self.pixels, c)
			new.ResizeInfos()
			return new
		else:
			self.pixels = cut(self.pixels, c)
			self.ResizeInfos()

def cut(img, c):
	if c[2]-c[0]==img.shape[1] and c[3]-c[1]==img.shape[0]:
		raise Exception("Reshaping the same array again?")
	return img[c[1]:c[3],c[0]:c[2]]
	
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

def tukeyfy(A, alpha):
	tuky = tukeywin(A.shape[0],alpha)
	tukx = tukeywin(A.shape[1],alpha)
	tuk = np.multiply(tukx[:,None].T,tuky[:,None])
	return A * tuk

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

def Normalize(x):
	return (x-np.min(x))/(np.max(x)-np.min(x))

def overlay(ax,mask,color,**kargs):
	m = ma.masked_array(mask,~mask)
	col = np.array(colors.colorConverter.to_rgba(color))
	I=col[:,None,None].T*m[:,:,None]
	ax.imshow(I,**kargs);

def NormP(x,p, trunk=True):
	thresh_high = np.percentile(x,100-p)
	thresh_low = np.percentile(x,p)
	if thresh_low==thresh_high:
		thresh_high=np.max(x)
		thresh_low=np.min(x)
	if thresh_low==thresh_high:
		thresh_high=thresh_low+1
	r=(x-thresh_low)/(thresh_high-thresh_low)
	if trunk:
		r[r<0]=0
		r[r>1]=1
	return r

def BeamProfile(target, source, mu = 1e-6):
	source = 2*source-1
	tf = np.fft.fft2(source)
	tf /= np.size(tf)
	recon_tf = 1 / ( tf + (mu / np.conj(tf)) )
	return np.fft.fftshift(np.real(np.fft.ifft2( np.fft.fft2(target) * recon_tf )))

def px2real(x,y,size,ext):
    rx=ext[0]+(x/size[1])*(ext[1]-ext[0])
    ry=ext[2]+(y/size[0])*(ext[3]-ext[2])
    return rx,ry

def real2px(x,y,size,ext):
    px = size[1]*(x-ext[0])/(ext[1]-ext[0])
    py = size[0]*(y-ext[2])/(ext[3]-ext[2])
    return px,py

def gaussian(x, mu, sig, A=1):
    return A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def stat(x):
	print("Min: {mi:.3f}, Max: {ma:.3f}, Mean: {mean:.3f}, Std: {std:.3f}".format(mi=np.min(x),ma=np.max(x), mean=np.mean(x), std=np.std(x)))

def fit2d(Z,dx=2,dy=1):
	x = np.arange(Z.shape[1],dtype=np.float)
	y = np.arange(Z.shape[0],dtype=np.float)
	X, Y = np.meshgrid(x,y)
	x2 = X.ravel()
	y2 = Y.ravel()
	A=np.vstack([x2**i for i in range(dx+1)])
	A=np.vstack([A]+[y2**i for i in range(1,dy+1)])
	res = scipy.optimize.lsq_linear(A.T,Z.ravel())
	r=res['x']
	Z2 = r[0]*np.ones(Z.shape)
	for i in range(1,dx+1): Z2+=r[i]*(X**i)
	for i in range(1,dy+1): Z2+=r[dx+i]*(Y**i)
	return r,Z2

def align(img, tform):
	New = tf.warp(img, tform, preserve_range=True)
	Cut=[0,0]+list(img.shape)
	if tform.translation[0]>=0: Cut[2]-=tform.translation[0]
	elif tform.translation[0]<0: Cut[0]-=tform.translation[0]
	if tform.translation[1]>=0: Cut[1]+=tform.translation[1]
	elif tform.translation[1]<0: Cut[3]+=tform.translation[1]
	Cut = [int(x) for x in Cut]
	New = cut(New, Cut)
	return New, Cut

if __name__ == "__main__":
	Path = "C:/Users/ols/Dropbox/ToF_SIMS"
	AFM = SPM_image("{0}/CyI5b_0006_ns.xml".format(Path),corr='slope')
	CC2 = SPM_image("{0}/CyI5b_PCBM_CC2.xml".format(Path))
	AFM.show()
