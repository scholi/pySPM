import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA as PCA1
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy import stats
import matplotlib.pyplot as plt
from pySPM.SPM import SPM_image
import matplotlib as mpl
from matplotlib import cm
from pySPM import collection

class PCA:
	def __init__(self, data):
		self.data = data
		self.pca = None
		self.standX = None
		
	def scatter(self, **kargs):
		pd.tools.plotting.scatter_matrix(self.data, diagonal="kde", **kargs);
		plt.tight_layout();
		
	def corr(self):
		corrmat = self.data.corr()
		return corrmat
	
	def corrShow(self):
		import seaborn as sns
		corrmat = self.corr()
		sns.heatmap(corrmat).xaxis.tick_top()
	
	def hinton(self, max_weight=None, ax=None):
		"""Draw Hinton diagram for visualizing a weight matrix."""
		matrix = self.corr()
		ax = ax if ax is not None else plt.gca()

		if not max_weight:
			max_weight = 2**np.ceil(np.log(np.abs(matrix).max())/np.log(2))

		ax.patch.set_facecolor('lightgray')
		ax.set_aspect('equal', 'box')
		ax.xaxis.set_major_locator(plt.NullLocator())
		ax.yaxis.set_major_locator(plt.NullLocator())

		for (x, y), w in np.ndenumerate(matrix):
			color = 'red' if w > 0 else 'blue'
			size = np.sqrt(np.abs(w))
			rect = plt.Rectangle([x - size / 2, y - size / 2], size, size,
				facecolor=color, edgecolor=color)
			ax.add_patch(rect)

		nticks = matrix.shape[0]
		ax.xaxis.tick_top()
		ax.set_xticks(range(nticks))
		ax.set_xticklabels(list(matrix.columns), rotation=90)
		ax.set_yticks(range(nticks))
		ax.set_yticklabels(matrix.columns)
		ax.grid(False)

		ax.autoscale_view()
		ax.invert_yaxis()
		
	def standardized(self):
		self.standX = pd.DataFrame(scale(self.data), index=self.data.index, columns=self.data.columns)
		return self.standX

	def runPCA(self):
		if self.standX is None:
			self.standardized()
		self.pca = PCA1().fit(self.standX)

	def pca_summary(self):
		if self.pca is None:
			self.runPCA()
		names = ["PC"+str(i) for i in range(1, len(self.pca.explained_variance_ratio_)+1)]
		a = list(np.std(self.pca.transform(self.standX), axis=0))
		b = list(self.pca.explained_variance_ratio_)
		c = [np.sum(self.pca.explained_variance_ratio_[:i]) for i in range(1, len(self.pca.explained_variance_ratio_)+1)]
		columns = pd.MultiIndex.from_tuples([("sdev", "Standard deviation"), ("varprop", "Proportion of Variance"), ("cumprop", "Cumulative Proportion")])
		Z=zip(a, b, c)
		summary = pd.DataFrame(list(Z), index=names, columns=columns)
		return summary
		
	def screeplot(self,ax=None):
		if self.pca is None:
			self.runPCA()
		ax = ax if ax is not None else plt.gca()
		y = np.std(self.pca.transform(self.standX), axis=0)**2
		x = np.arange(len(y)) + 1
		ax.grid(True)
		ax.plot(x,y, "o-")
		ax.set_xticks(x)
		ax.set_xticklabels(["Comp."+str(i) for i in x], rotation=60)
		ax.set_ylabel("Variance")
		
	def pc(self,id=0):
	# find the number of samples in the data set and the number of variables
		if self.pca is None:
			self.runPCA()
		pc = np.matmul(self.standX.as_matrix(),self.pca.components_[id].T)
		return pc
		
	def loadings(self, id=None):
		if self.pca is None:
			self.runPCA()
		if id is not None:
			return pd.DataFrame(self.pca.components_[id,None],columns=self.data.columns)
		return pd.DataFrame(self.pca.components_,columns=self.data.columns,index=["PC{0}".format(i+1) for i in range(len(self.pca.components_))])
	
	def getPCAtransf(self):
		if self.pca is None:
			self.runPCA()
		return self.pca.transform(self.standX)
	
	def showStand(self):
		return pd.DataFrame([self.standX.apply(np.mean),self.standX.apply(np.std)],index=['Mean','Std'])
	
	def pca_scatter(self, classifs=None, light=False):
		import seaborn as sns
		foo = self.getPCAtransf()
		if classifs is None:
			if light:
				plt.scatter(foo[:,0],foo[:,1])
			else:
				bar = pd.DataFrame(list(zip(foo[:, 0], foo[:, 1])), columns=["PC1", "PC2"])
				sns.lmplot("PC1", "PC2", bar, fit_reg=False)
		else:
			if light:
				plt.scatter(foo[:,0],foo[:,1],color=cm.Scalar)
			else:
				bar = pd.DataFrame(list(zip(foo[:, 0], foo[:, 1], classifs)), columns=["PC1", "PC2", "Class"])
				sns.lmplot("PC1", "PC2", bar, hue="Class", fit_reg=False)
		
class ITA_PCA(PCA):
	def __init__(self, c, channels=None):
		"""Init a PCA class from a collection"""
		self.col = c
		if channels is None:
			channels=c.channels.keys()
		mul = self.col.get_multivariate(channels)
		PCA.__init__(self, mul)
	
	def showPCA(self, num=None, ax=None):
		c = self.getPCAcol(num)
		c.show(ax=ax,cmap='hot')
		
	def getPCAcol(self, num=None):
		if num is None:
			num = self.data.shape[1]
		assert num<=self.data.shape[1]
		PCA_col = collection.Collection(cls=self.col,name=self.col.name+"[PCA]")
		for i in range(num):
			PC = self.getPCA(i)
			PCA_col.add(PC,'PC{0}'.format(i+1))
		return PCA_col
			
	def getPCA(self,id=0):
		s = list(self.col.channels.values())[0].shape
		PC = self.pc(id).reshape(s)
		return PC
		
