import sys
from collections import OrderedDict
import numpy as np
from sklearn.cluster import KMeans
from sklearn.utils import resample
from sklearn import metrics
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from scipy import stats
from Bio.Seq import Seq
from .Data import LoadData
from .colors import colors_hex
from .fonts import fonts
from .RunCmdsMP import logger, pool_func

class Cluster:
	def __init__(self, datafile, n_clusters, sg_prefix='SG', sg_assigned={}, 
						replicates=1000, jackknife=80, **kargs):
		data = LoadData(datafile)
		data.load_matrix()
		# normalize
		self.raw_data = data.data
		self.data = self.raw_data.transpose()	# col: kmer, row: chr
		self.data = self.normalize_data(self.data, axis=0)
		
		self.chrs = data.colnames
		self.kmers = data.rownames
		self.d_kmers = data.d_rows
		
		self.n_clusters = n_clusters
		self.sg_prefix = sg_prefix
		self.kargs = kargs
		if sg_assigned:
			logger.info('Skip k-means clustering')
			labels = [sg_assigned[chr] for chr in self.chrs]
			self.n_clusters = len(set(sg_assigned.values()))
			self.d_sg = self.assign_subgenomes(labels=labels)
		else:
			self.kmean = self.fit(self.data, n_clusters, **kargs)
			self.d_sg = self.assign_subgenomes()
		self.d_bs = self.bootstrap(replicates, jackknife)
	def pca(self, outfig, n_components=2, ):
		pca = PCA(n_components=n_components)
		X_pca = pca.fit_transform(self.data)
		percent = pca.explained_variance_ratio_ * 100
		X_pca = self.normalize_data(X_pca, axis=0)
		x = X_pca[:, 0]
		y = X_pca[:, 1]
#		colors = [colors_hex[lab] for lab in self.labels]
		plt.figure(figsize=(7,7), dpi=300, tight_layout=True)
#		plt.scatter(x, y, c=colors, marker='o')
		d_coord = {}
		for _x, _y, _c, _l in zip(x, y, self.chrs, self.labels):
			sg = self.d_sg[_c]
			if sg not in d_coord:
				d_coord[sg] = [[], [], colors_hex[_l]]
			d_coord[sg][0] += [_x]
			d_coord[sg][1] += [_y]
		for sg, (_x, _y, _c) in sorted(d_coord.items()):
			plt.scatter(_x, _y, c=_c, marker='o', label=sg)
		plt.axhline(0, ls='--', c='grey')
		plt.axvline(0, ls='--', c='grey')
		xlabel = 'PC1 ({:.1f}%)'.format(percent[0])
		ylabel = 'PC2 ({:.1f}%)'.format(percent[1])
		plt.xlabel(xlabel, fontsize=fonts['fontsize'])
		plt.ylabel(ylabel, fontsize=fonts['fontsize'], ha='center', va='center')
		plt.legend(fontsize=fonts['fontsize'])
		plt.tick_params(labelsize=fonts['labelsize'])
		plt.savefig(outfig, bbox_inches='tight', dpi=300)
	def normalize_data(self, data, axis=0):
		'''Z normalization: only axis=0 work'''
		mean = data.mean(axis=axis)
		std = data.std(axis=axis)
		return (data-mean)/std

	def bootstrap(self, replicates=1000, jackknife=80):
		logger.info('Performing bootstrap of {} replicates, with each replicate resampling {}% data\
 with replacement'.format(replicates, jackknife))
		jackknife = max(int(jackknife/100 * len(self.kmers)), 100)
		xlabels = []
		scores, measures = [], []
		raw_data = self.data.transpose()
		for i in range(replicates):
			data = resample(raw_data, replace=True, n_samples=replicates)
			data = data.transpose()
#			print(data.shape)
			kmean = self.fit(data, self.n_clusters)
			labels = kmean.labels_
			labels = self.sort_subgenomes(labels)
			xlabels += [labels]
			adjusted_rand_score = metrics.adjusted_rand_score(self.labels, labels)
			scores += [adjusted_rand_score]
			v_measure_score = metrics.v_measure_score(self.labels, labels)
			measures += [v_measure_score]
		xlabels = np.array(xlabels)
#		print(xlabels.shape)
		d_bs = {}
		for i, (label,chr) in enumerate(zip(self.labels, self.chrs)):
			clabels = list(xlabels[:, i])
			bs = clabels.count(label)
			d_bs[chr] = int(100*bs/replicates)
		self.mean_adjusted_rand_score = np.mean(scores)
		self.mean_v_measure_score = np.mean(measures)
		logger.info('Bootstrap: mean Adjusted Rand-Index: {:.4f}; mean V-measure score: {:.4f}'.format(
			self.mean_adjusted_rand_score, self.mean_v_measure_score))
		return d_bs

	def fit(self, data, n_clusters, **kargs):
		'''fit KMeans cluster'''
		kmean = KMeans(n_clusters=n_clusters)
		kmean.fit(data)
		return kmean
	def sort_subgenomes(self, labels):
		assert len(self.chrs) == len(labels)
		d_map = {}
		for label, chr in sorted(zip(labels, self.chrs), key=lambda x:x[1]):
			if label not in d_map:
				try: d_map[label] = max(d_map.values()) + 1
				except ValueError: d_map[label] = 0
		return [d_map[label] for label in labels]

	def assign_subgenomes(self, base=1, labels=None):
		'''labels is the same order as self.chrs'''
		if labels is None:
			labels = self.kmean.labels_
		max_len = len(str(self.n_clusters))
		fmtstr = '{{}}{{:0>{}d}}'.format(max_len)
		d_sg = OrderedDict()
		sg_names = set([])
		self.labels = labels = self.sort_subgenomes(labels)
		assert len(self.chrs) == len(labels)
		for label, chr in zip(labels, self.chrs):
			sg = fmtstr.format(self.sg_prefix, label+base)
			d_sg[chr] = sg
			sg_names.add(sg)
		self.sg_names = sorted(sg_names)
		return d_sg
	def output_subgenomes(self, fout=sys.stdout):
		line = ['#chrom', 'subgenome', 'bootstrap']
		print('\t'.join(line), file=fout)
		for chr, sg in sorted(self.d_sg.items(), key=lambda x:x[1]):
			line = [chr, sg, self.d_bs[chr]]
			line = map(str, line)
			print('\t'.join(line), file=fout)
	def output_kmers(self, fout=sys.stdout, max_pval=0.05, ncpu=4, method='map', test_method='ttest_ind'):
		'''test_method: kruskal, '''
		d_groups = {}
		for i, (chr,sg) in enumerate(self.d_sg.items()):
			try: d_groups[sg] += [i]
			except KeyError: d_groups[sg] = [i]
		d_ksg = {}
		line = ['#kmer', 'subgenome', 'p_value', 'ratios']
		print('\t'.join(line), file=fout)
		test_method = eval('stats.{}'.format(test_method))
		iterable = ((kmer, array, d_groups, test_method) for kmer, array in self.d_kmers.items())
		jobs = pool_func(_output_kmers, iterable, processors=ncpu, method='map', )
		#jobs = list(jobs)
		i = 0
		for kmer, max_sg, pvalue, rc_kmer, mean_vals in jobs:
			i += 1
			if pvalue > max_pval:
				continue
			ratios = ','.join(map(str, mean_vals))
			line = [kmer, max_sg, pvalue, ratios]
			line = map(str, line)
			print('\t'.join(line), file=fout)
#			kmer = tuple(kmer)
			d_ksg[kmer] = max_sg
			d_ksg[rc_kmer] = max_sg
		return d_ksg
#	def _output_kmers(self, args):
def _output_kmers(args):
		kmer, array, d_groups, test_method = args
		grouped = [ [array[i] for i in idx] for sg, idx in sorted(d_groups.items())]
		mean_vals = [np.mean(x) for x in grouped]
		sgs = sorted(d_groups.keys())
		xgrouped = sorted(zip(grouped, sgs), key=lambda x: -sum(x[0])/len(x[0]))
		grouped = [x[0] for x in xgrouped]  # sorted with the same order
		sgs = [x[1] for x in xgrouped]
		max_freqs = grouped[0]
		min_freqs = grouped[1]
		max_sg = sgs[0]
		#ttest = stats.ttest_ind(max_freqs, min_freqs)
		#test = stats.kruskal(max_freqs, min_freqs)
		test = test_method(max_freqs, min_freqs)
		pvalue = test.pvalue
		rc_kmer = str(Seq(kmer).reverse_complement())
		return kmer, max_sg, pvalue, rc_kmer, mean_vals
