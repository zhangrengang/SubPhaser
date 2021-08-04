import fisher
import copy
import numpy as np

def fisher_test(each, total):
	assert len(each) == len(total)
	pvals = []
	sum_each = sum(each)
	sum_total = sum(total)
	for i in range(len(each)):
		x11 = each[i]
		x12 = sum_each - x11
		x21 = total[i] - x11
		x22 = sum_total-x21 - x12
		pval = fisher.pvalue(x11, x12, x21, x22).right_tail
		pvals += [pval]
	return pvals

def enrich_ltr(fout, *args, **kargs):
	'''Output LTR enrichments'''
	line = ['#ltr', 'SG', 'pval', 'counts']
	fout.write('\t'.join(line)+'\n')
	d_enriched = {}
	for res in enrich(*args, **kargs):
		ltr, *_ = res.rowname
		counts = ','.join(map(str, res.counts))
		line = [ltr, res.key, res.pval, counts]
		line = map(str, line)
		fout.write('\t'.join(line)+'\n')
		d_enriched[ltr] = res.key
	return d_enriched	# significant results

def enrich_bin(fout, *args, **kargs):
	'''Enrich by chromosome bins'''
	line = ['#chrom', 'start', 'end', 'SG', 'pval', 'counts']
	fout.write('\t'.join(line)+'\n')
	last_end = 0
	for res in enrich(*args, **kargs):
		chrom, start, end = res.rowname
		counts = ','.join(map(str, res.counts))
		line = [chrom, start, end, res.key, res.pval, counts]
		line = map(str, line)
		fout.write('\t'.join(line)+'\n')
		if start > last_end:
			line = [chrom, last_end, start, None]
			line = map(str, line)
			fout.write('\t'.join(line)+'\n')
		last_end = end

def enrich(matrix, colnames=None, rownames=None, **kargs):	# kargs: min_pval
	arr = np.array(matrix)
	if colnames is not None and rownames is not None:
		assert arr.shape == (len(rownames), len(colnames))
	total = list(arr.sum(axis=0))	# sum by column
	for row, rowname in zip(matrix, rownames):
		pvals = fisher_test(row, total)
		pvals = Pvalues(pvals, colnames, **kargs)
		_min = pvals.get_enriched()
		if _min is None:
			continue
		_min.counts = row
		_min.pvals = pvals
		_min.rowname = rowname
		yield _min	# min pvalue and max proportion

class Pvalues:
	def __init__(self, pvals, keys, cutoff=100, min_pval=0.05):
		assert len(pvals) == len(keys) > 1
		self.pvals = pvals
		self.keys = keys
		self.cutoff = cutoff
		self.min_pval = min_pval
	def __iter__(self):
		return (Pvalue(pval, key) for pval, key in zip(self.pvals, self.keys))
	def sort(self):
		return sorted(self, key=lambda x:x.pval)
	def get_enriched(self):
		pvals = self.sort()
		_min = pvals[0]
		_submin = pvals[1]
		if _min.pval > self.min_pval:
			return None
		if _min.pval == 0:
			return _min
		if _submin.pval / _min.pval > self.min_pval/_submin.pval:
			return _min
		return None

class Pvalue:
	def __init__(self, pval, key):
		self.pval = pval
		self.key = key
