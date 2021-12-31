import fisher
import re
import copy
import numpy as np
from statsmodels.stats.multitest import multipletests
from .RunCmdsMP import logger, pool_func

MAX_INT = 2147483647 // 10

def correct_pvals(pvals, method='fdr_bh'):
	return multipletests(pvals, method=method)[1]
	
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
		x21 = min(x21, MAX_INT)
		x22 = min(x22, MAX_INT)
		try: pval = fisher.pvalue(x11, x12, x21, x22).right_tail
		except OverflowError as e:
			print(each, total, (x11, x12, x21, x22))
			raise OverflowError(e)
		pvals += [pval]
	return pvals

def enrich_ltr(fout, d_sg, *args, **kargs):
	'''Output LTR enrichments'''
	total, consistent, exchange = 0,0,0
	d_enriched = {}
	lines = []
	pvalues = []
	for res in enrich(*args, **kargs):
		ltr, *_ = res.rowname
		try: chrom = re.compile(r'(\S+?):\d+\-\d+').match(ltr).groups()[0]
		except TypeError: chrom = None
		obs_sg = d_sg.get(chrom)
		sg = res.key if res.sig else None
#			d_enriched[ltr] = 'shared'
#			continue
		potential_exchange = is_exchange(obs_sg, sg)
		counts = ','.join(map(str, res.counts))
		line = [ltr, sg, res.pval, counts, potential_exchange]
		lines += [line]
		pvalues += [res.pval]
		if sg:
			d_enriched[ltr] = sg #res.key
		total += 1
		if potential_exchange == 'yes':
			exchange += 1
		elif potential_exchange == 'no':
			consistent += 1
	if exchange>0 and consistent>0:
		logger.info('Consistent with subgenome assignment: {} ({:.2%}); potential exchange: {} ({:.2%})'.format(
			consistent, consistent/total, exchange, exchange/total))
	# correct p values
	qvals = correct_pvals(pvalues)
	# output
	line = ['#id', 'subgenome', 'p_value', 'counts', 'potential_exchange', 'p_corrected']
	fout.write('\t'.join(line)+'\n')
	for line, qval in zip(lines, qvals):
		line += [qval]
		line = list(map(str, line))
		fout.write('\t'.join(line)+'\n')
	return d_enriched	# significant results

def enrich_bin(fout, d_sg, *args, **kargs):
	'''Enrich by chromosome bins'''
	
	total, consistent, exchange = 0,0,0
	lines = []
	pvalues = []
	for res in enrich(*args, **kargs):
		chrom, start, end = res.rowname
		key = res.key if res.sig else None	# expected SG
		obs_sg = d_sg.get(chrom)
		potential_exchange = is_exchange(obs_sg, key)
		counts = ','.join(map(str, res.counts))
		enrichs = ','.join(map(str, res.enrich))
		ratios = ','.join(map(str, res.ratios))
		pvals = ','.join(map(str, res.pvals))
		line = [chrom, start, end, key, res.pval, counts, ratios, enrichs, pvals, potential_exchange]
		
		lines += [line]
		pvalues += [res.pval]
		total += 1
		if potential_exchange == 'yes':
			exchange += 1
		elif potential_exchange == 'no':
			consistent += 1
	logger.info('Consistent with subgenome assignment: {} ({:.2%}); potential exchange: {} ({:.2%})'.format(
		consistent, consistent/total, exchange, exchange/total))
	# correct p values
	qvals = correct_pvals(pvalues)
	# output
	line = ['#chrom', 'start', 'end', 'subgenome', 'p_value', 'counts', 'ratios', 'enrich','pvals',
				'potential_exchange', 'p_corrected']
	fout.write('\t'.join(line)+'\n')
	for line, qval in zip(lines, qvals):
		line += [qval]
		line = list(map(str, line))
		fout.write('\t'.join(line)+'\n')
	return lines
def is_exchange(obs_sg, exp_sg):
	if not exp_sg or not obs_sg:
		return 'none'
	if obs_sg == exp_sg:
		return 'no'
	return 'yes'
	
def enrich(matrix, colnames=None, rownames=None, ncpu=4, min_ratio=0.5, **kargs):	# kargs: max_pval
	arr = np.array(matrix)
	if colnames is not None and rownames is not None:
		assert arr.shape == (len(rownames), len(colnames)), '{} != {}'.format(
								arr.shape, (len(rownames), len(colnames)))
	total = list(arr.sum(axis=0))	# sum by column
	iterable = ((row, rowname, total, colnames, min_ratio, kargs) for row, rowname in zip(matrix, rownames))
	for _min in pool_func(_enrich, iterable, processors=ncpu, method='map'):
		yield _min
	
def _enrich(args):
	row, rowname, total, colnames, min_ratio, kargs = args
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
	return _min
	
class Pvalues:
	def __init__(self, pvals, keys, cutoff=1, max_pval=0.05):
		assert len(pvals) == len(keys) > 1
		self.pvals = pvals
		self.keys = keys
		self.cutoff = cutoff
		self.max_pval = max_pval
	def __iter__(self):
		return (Pvalue(pval, key, i) for i, (pval, key) in enumerate(zip(self.pvals, self.keys)))
	def sort(self):
		return sorted(self, key=lambda x:x.pval)
	def get_enriched(self):
		pvals = self.sort()
		_min = pvals[0]
		_submin = pvals[1]
		_min.sig = True	# significant
		if _min.pval > self.max_pval:
			_min.sig = False
		if _min.pval == 0:	# to avoid divide by zero
			pass
		elif _submin.pval / _min.pval < self.max_pval/_submin.pval*self.cutoff:
			_min.sig = False
		return _min

class Pvalue:
	def __init__(self, pval, key, idx):
		self.pval = pval
		self.key = key
		self.idx = idx
