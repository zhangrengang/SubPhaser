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
		if not res.sig:
			continue
		ltr, *_ = res.rowname
		counts = ','.join(map(str, res.counts))
		line = [ltr, res.key, res.pval, counts]
		line = map(str, line)
		fout.write('\t'.join(line)+'\n')
		d_enriched[ltr] = res.key
	return d_enriched	# significant results

def enrich_bin(fout, *args, **kargs):
	'''Enrich by chromosome bins'''
	line = ['#chrom', 'start', 'end', 'SG', 'pval', 'counts', 'ratios', 'enrich','pvals']
	fout.write('\t'.join(line)+'\n')
#	last_end = 0
	lines = []
	for res in enrich(*args, **kargs):
		chrom, start, end = res.rowname
		key = res.key if res.sig else None
		counts = ','.join(map(str, res.counts))
		enrichs = ','.join(map(str, res.enrich))
		ratios = ','.join(map(str, res.ratios))
		pvals = ','.join(map(str, res.pvals))
		line = [chrom, start, end, key, res.pval, counts, ratios, enrichs, pvals]
		line = list(map(str, line))
		fout.write('\t'.join(line)+'\n')
#		if start > last_end:
#			line = [chrom, last_end, start, None]
#			line = map(str, line)
#			fout.write('\t'.join(line)+'\n')
#		last_end = end
		lines += [line]
	return lines

def enrich(matrix, colnames=None, rownames=None, min_ratio=0.5, **kargs):	# kargs: min_pval
	arr = np.array(matrix)
	if colnames is not None and rownames is not None:
		assert arr.shape == (len(rownames), len(colnames))
	total = list(arr.sum(axis=0))	# sum by column
	for row, rowname in zip(matrix, rownames):
		pvals = fisher_test(row, total)
		pvals = Pvalues(pvals, colnames, **kargs)
		_min = pvals.get_enriched()
		_min.counts = row
		_min.pvals = pvals.pvals
		ratios = np.array(row) / np.array(total)
		_min.ratios = ratios / ratios.sum()
		_min.ratio = _min.ratios[_min.idx]
		if _min.ratio < min_ratio:
			_min.sig = False
		_min.enrich = [0]* (len(colnames)+1)
		if _min.sig:
			_min.enrich[_min.idx] = 1
		else:
			_min.enrich[-1] = 1
		_min.rowname = rowname
		yield _min	# min pvalue and max proportion

class Pvalues:
	def __init__(self, pvals, keys, cutoff=1, min_pval=0.05):
		assert len(pvals) == len(keys) > 1
		self.pvals = pvals
		self.keys = keys
		self.cutoff = cutoff
		self.min_pval = min_pval
	def __iter__(self):
		return (Pvalue(pval, key, i) for i, (pval, key) in enumerate(zip(self.pvals, self.keys)))
	def sort(self):
		return sorted(self, key=lambda x:x.pval)
	def get_enriched(self):
		pvals = self.sort()
		_min = pvals[0]
		_submin = pvals[1]
		_min.sig = True	# significant
		if _min.pval > self.min_pval:
			_min.sig = False
		if _min.pval == 0:	# to avoid divide by zero
			pass
		elif _submin.pval / _min.pval < self.min_pval/_submin.pval*self.cutoff:
			_min.sig = False
		return _min

class Pvalue:
	def __init__(self, pval, key, idx):
		self.pval = pval
		self.key = key
		self.idx = idx
