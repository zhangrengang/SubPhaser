import sys
from collections import OrderedDict
from sklearn.cluster import KMeans
from scipy import stats
from .Data import LoadData

class Cluster:
	def __init__(self, datafile, n_clusters, **kargs):
		data = LoadData(datafile)
		data.load_matrix()
		self.data = data.data.transpose()
		self.chrs = data.colnames
		self.kmers = data.rownames
		self.d_kmers = data.d_rows
		self.kmean = self.fit(n_clusters, **kargs)
		self.n_clusters = n_clusters
		self.kargs = kargs
		self.d_sg = self.assign_subgenomes()
	def fit(self, n_clusters, **kargs):
		kmean = KMeans(n_clusters=n_clusters)
		kmean.fit(self.data)
		return kmean
	def assign_subgenomes(self, prefix='SG', base=1):
		d_sg = OrderedDict()
		assert len(self.chrs) == len(self.kmean.labels_)
		for label, chr in zip(self.kmean.labels_, self.chrs):
			sg = '{}{}'.format(prefix, label+base)
			d_sg[chr] = sg
		return d_sg
	def output_subgenomes(self, fout=sys.stdout, prefix='SG', base=1):
		line = ['#chrom', 'subgenome']
		print('\t'.join(line), file=fout)
		for chr, sg in sorted(self.d_sg.items(), key=lambda x:x[1]):
			line = [chr, sg]
			print('\t'.join(line), file=fout)
	def output_kmers(self, fout=sys.stdout, min_pval=0.05):
		d_groups = {}
		for i, (chr,sg) in enumerate(self.d_sg.items()):
			try: d_groups[sg] += [i]
			except KeyError: d_groups[sg] = [i]
		d_ksg = {}
		line = ['#kmer', 'subgenome', 'p_value']
		print('\t'.join(line), file=fout)
		for kmer, array in self.d_kmers.items():
			grouped = [ [array[i] for i in idx] for sg, idx in sorted(d_groups.items())]
			sgs = sorted(d_groups.keys())
			xgrouped = sorted(zip(grouped, sgs), key=lambda x: -sum(x[0])/len(x[0]))
			grouped = [x[0] for x in xgrouped]	# sorted with the same order
			sgs = [x[1] for x in xgrouped]
			max_freqs = grouped[0]
			min_freqs = grouped[1]
			max_sg = sgs[0]
			ttest = stats.ttest_ind(max_freqs, min_freqs)
			if ttest.pvalue > min_pval:
				continue
			line = [kmer, max_sg, ttest.pvalue]
			line = map(str, line)
			print('\t'.join(line), file=fout)
			kmer = tuple(kmer)
			d_ksg[kmer] = max_sg
		return d_ksg
