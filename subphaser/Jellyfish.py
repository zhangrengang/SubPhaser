import sys,os
from collections import OrderedDict,Counter
import numpy as np
from Bio.Seq import Seq
from scipy import stats
from xopen import xopen as open
from .RunCmdsMP import run_cmd, logger, pool_func
from .small_tools import is_gz, mk_ckp, check_ckp
from .fonts import fonts
from .colors import colors_r

class JellyfishDumpRecord(object):
	def __init__(self, seq, freq):
		self.seq = seq
		self.freq = int(freq)
	@property
	def iseq(self):
		return seq2int(self.seq)
class JellyfishDumpLine:
	def __init__(self, line):
		seq, freq, *plus = line.strip().split()
		self.seq = seq
		self.freq = int(freq)
		if plus:
			self.sample = plus[0]
class JellyfishHistoLine:
	def __init__(self, line):
		depth, freq, *plus = line.strip().split()
		self.depth = int(depth)
		self.freq = int(freq)
class JellyfishHisto:
	def __init__(self, histo,):
		self.histo = histo
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.histo):
			yield JellyfishHistoLine(line)
	def get_total_kmers(self, lower=1, upper=1e12):
		total = 0
		for rc in self:
			if lower <= rc.depth <= upper:
				total += rc.depth * rc.freq
		return total

class JellyfishDump(object):
	def __init__(self, column, ncpu=4, ):
		self.column = column
		self.ncpu = ncpu
	def __iter__(self):
		return self.parse()
	def parse(self):
		return self._parse1()
	def _parse2(self):
		for rc in pool_func(_parse_line, open(self.column), self.ncpu, ):
			yield rc
	def _parse1(self):
		for line in open(self.column):
		#	seq, freq = line.strip().split()
		#	yield JellyfishDumpRecord(seq=seq, freq=freq)
			yield JellyfishDumpLine(line)
	def _parse_line(self, line):
		seq, freq = line.strip().split()
		return JellyfishDumpRecord(seq=seq, freq=freq)
	def get_kmers(self, kmer_set):
		d_kmer = {}
		for rc in self.parse():
			if rc.seq in kmer_set:
				d_kmer[rc.seq] = rc.freq
		return d_kmer
	def to_dict(self):
		d = {}
		for rc in self:
			d[rc.seq] = rc.freq
		return d
class JellyfishDumpLines:
	def __init__(self, lines):
		self.lines = lines
		self.seq = lines[0].seq
	@property
	def freq_dict(self):
		return {line.sample: line.freq for line in self.lines}
	def __len__(self):
		return len(self.lines)
		
def _parse_line(line):
	seq, freq = line.strip().split()
	return JellyfishDumpRecord(seq=seq, freq=freq)

def _to_matrix2(arg):
	dumpfile = arg
	seqs, freqs = [], []
	tot = 0
	for rc in JellyfishDump(dumpfile):
		seqs += [rc.seq]
		freqs += [rc.freq]
		tot += rc.freq
	return seqs, freqs, dumpfile, tot
	
class JellyfishDumpSamples:
	def __init__(self, merged_dumpfile, samples=[], groups={}, 
			ncpu=4, method='imap_unordered', chunksize=10000, **kargs):
		self.merged_dumpfile = merged_dumpfile
		self.samples = samples
		self.nsampels = len(self.samples)
		self.groups = groups
		self.ncpu = ncpu
		self.method = method
		self.chunksize = chunksize
	def __iter__(self):
		return self._parse()
	def _parse(self):
		last_seq = ''
		lines = []
		for rc in JellyfishDump(self.merged_dumpfile):
			if last_seq and rc.seq != last_seq:
				yield JellyfishDumpLines(lines)
				lines = []
			lines += [rc]
			last_seq = rc.seq
	def prefilter(self, min_mean_freq=1, max_missing_rate=0.5, by_group=False):
		for lines in self:
			d_freq = lines.freq_dict
			
			if not by_group:
				missing_rate = 1 - len(lines)/self.nsampels
				if missing_rate > max_missing_rate:
					continue
			else:
				groups = []
				for group, samples in self.groups.items():
					non_missing = len([1 for sample in samples if sample in d_freq])
					missing_rate = 1 - len(non_missing)/len(samples)
					groups += [missing_rate > max_missing_rate]
				if all(groups):	# all with too many missing
					continue
			mean_freq = sum(d_freq.values()) / self.nsampels
			if mean_freq < min_mean_freq:
				continue
			kmer = lines.seq
			freqs = [d_freq.get(sample, 0) for sample in self.samples]
			yield kmer, freqs

	def stats(self, fout=sys.stdout, out_sig=sys.stdout, out_all=sys.stdout, 
			denominators=None, max_pval=0.05, 
			test_method='ttest_ind', filter_args={}):
		samples = self.samples
		d_groups = self.groups
		logger.info('Groups: {}'.format(d_groups.keys()))
		line = ['kmer', 'group', 'p_value', 'ratios', 'fold',]
		print('\t'.join(line), file=fout)
		line = ['kmer'] + samples
		print('\t'.join(line), file=out_all)
		print('\t'.join(line), file=out_sig)
		
		test_method = eval('stats.{}'.format(test_method))
		#iterable = ((kmer, array, denominators, samples, d_groups, test_method ) \
		#			for kmer, array in self.prefilter(**filter_args))
		iterable = ((lines, denominators, samples, d_groups, test_method, filter_args) \
						for lines in self)
		jobs = pool_func(_map_stats, iterable, processors=self.ncpu, method=self.method, chunksize=self.chunksize)
		i,j,k = 0,0,0
		d_ksg = {}
		for result in jobs:
			i += 1
			if i % 10000000 == 0:
				logger.info('Processed {} kmers'.format(i))
			if not result:
				continue
			j += 1
			kmer, freqs, max_group, pvalue, rc_kmer, fold, ratios = result
			line = [kmer] + freqs
			line = list(map(str, line))
			print('\t'.join(line), file=out_all)
			if pvalue > max_pval:
				continue
			k += 1
			print('\t'.join(line), file=out_sig)
			
			line = [kmer, max_group, pvalue, ratios, fold, ]
			line = map(str, line)
			print('\t'.join(line), file=fout)
			
			d_ksg[kmer] = max_group
			d_ksg[rc_kmer] = max_group
		logger.info('Total kmers: {}; After prefilter: {} ({:.2%}); Significant kmers: {} ({:.2%})'.format(
				i, j, j/i, k, k/j))
		return d_ksg
	def heatmap(self, matfile, **kargs):
		_heatmap(matfile, **kargs)
class KmerMatrix:
	def __init__(self, matfile, d_total={},
			ncpu=4, method='imap_unordered', chunksize=10000):
		self.matfile = matfile
		self.d_total = d_total
		self.ncpu = ncpu
		self.method = method
		self.chunksize = chunksize
	def __iter__(self):
		return self._parse()
	def _parse(self):
		i = 0
		for line in open(self.matfile):
			i += 1
			if i == 1:
				samples = line.rstrip().split()[1:]
				continue
			yield line, samples
	def filter(self, fout=sys.stdout, **kargs):
		iterable = ((line, samples, kargs) \
						for line, samples in self)
		jobs = pool_func(_filter_matrix, iterable, 
				processors=self.ncpu, method=self.method, chunksize=self.chunksize)
		i,j,k = 0,0,0
		for result in jobs:
			i += 1
			if i % 10000000 == 0:
				logger.info('Processed {} kmers'.format(i))
			if not result:
				continue
			j += 1
			line = result
			if j == 1:
				line.write_title(fout)
			line.write(fout)
			
		logger.info('Total kmers: {}; After filter: {} ({:.2%});'.format(
				i, j, j/i, ))
		return 
	def to_gemma(self, prefix, genome=None):
		iterable = ((line, samples) \
						for line, samples in self)
		jobs = pool_func(_parse_line, iterable, 
				processors=self.ncpu, method=self.method, chunksize=self.chunksize)
		outgeno = prefix + '.geno'
		outfa = prefix + '.fa'
		outsnp = prefix + '.snps'
		fgeno = open(outgeno, 'w')
		ffa = open(outfa, 'w')
		for line in jobs:
			line.write_gemma(fgeno)
			line.write_fasta(ffa)
		if not genome:
			return
		
		ref= '{}.ref'.format( prefix)
		cmd = '[ ! -f {1}.ok ] && bowtie-build {0} {1} && touch {1}.ok'.format(genome, ref)
		run_cmd(cmd, log=True)
		outsam = prefix + '.sam'
		cmd = '''bowtie {ref} {reads} -f  -p 10 -S --best --sam-nohead 2> {reads}.map.log \
			| awk '{{if ($3=="*"){{$3="chr0";$4=NR*10}}; print $1"\t"$4"\t"$3}}' > {snp}'''.format(
				ref=ref, reads=outfa, snp=outsnp)
		run_cmd(cmd, log=True)

class KmerMatrixLine:
	def __init__(self, line=None, samples=[], freqs=[], seq=None):
		if line is not None:
			if isinstance(line, str):
				self.line = line
				self.seq, *self.freqs = line.strip().split()
				self.freqs = list(map(int, self.freqs))
			else: # list
				self.seq, *self.freqs = line
			
		else:
			self.seq, self.freqs = seq, freqs
		self.samples = samples
		self.recoded = False
	def __len__(self):
		return len(self.freqs)
	def write(self, fout):
		if self.recoded or not hasattr(self, 'line'):
			line = [self.seq] + self.freqs
		else:
			line = self.line

		if not isinstance(line, str):
			line = '\t'.join(map(str,line)) + '\n'
		fout.write(line)
	def pre_write(self, ):
		if self.recoded or not hasattr(self, 'line'):
			line = [self.seq] + self.freqs
			self.line = '\t'.join(map(str,line)) + '\n'
			self.recoded = False
	def write_title(self, fout):
		line = '\t'.join(['kmer']+self.samples) + '\n'
		fout.write(line)
	def write_gemma(self, fout):
		line = [self.seq, 0, 1] + self.freqs
		line = '\t'.join(map(str,line)) + '\n'
		fout.write(line)
	def write_fasta(self, fout):
		fout.write('>{0}\n{0}\n'.format(self.seq))
		
	def subset(self, subset):
		samples = list(subset)
		d_freq = self.freq_dict
		freqs = [d_freq[sample] for sample in samples if sample in d_freq]
		return KmerMatrixLine(samples=samples, freqs=freqs, seq=self.seq)
	@property
	def freq_dict(self):
		return {sample: freq for sample, freq in zip(self.samples, self.freqs)}
	@property
	def missing(self):
		return [sample for sample, freq in zip(self.samples, self.freqs) if freq ==0]
	@property
	def missing_rate(self):
		return len(self.missing) / len(self)
	@property
	def maf(self):
		counts = Counter(self.freqs).values()
		if len(counts) ==1:
			return 0
		return min(counts)/len(self)
	def recode_as_01(self):
		self.freqs = [min(1, freq) for freq in self.freqs]
		self.recoded = True
	def prefilter(self, subset=None, d_groups={}, min_mean_freq=1, max_missing_rate=0.2, by_group=True, **kargs):
		if subset is None:
			line = self
		else:
			line = self.subset(subset)
		d_freq = line.freq_dict
		if not by_group:
			missing_rate = line.missing_rate
			if missing_rate > max_missing_rate:
				return
			mean_freq = sum(d_freq.values()) / nsampels
			if mean_freq < min_mean_freq:
				return
		else:
			groups = []
			for group, _samples in d_groups.items():
				freqs = [d_freq[sample] for sample in _samples if sample in d_freq]
				_line = KmerMatrixLine(samples=_samples, freqs=freqs,)
				missing_rate = _line.missing_rate
				mean_freq = sum(freqs) / len(_samples)
				groups += [missing_rate > max_missing_rate or mean_freq < min_mean_freq]
			if all(groups):	# all with too many missing or two low depth
				return
		return line
	def recode_and_filter_by_maf(self, min_maf=0.05, **kargs):
		self.recode_as_01()
		maf = self.maf
		if maf < min_maf:
			return
		self.pre_write()
		return self
		
def _filter_matrix(args):
	line, samples, kargs = args
	line = KmerMatrixLine(line=line, samples=samples)
	line = line.prefilter(**kargs)
	if line is None:
		return None
	return line.recode_and_filter_by_maf(**kargs)
def _parse_line(args):
	line, samples, = args
	line = KmerMatrixLine(line=line, samples=samples)
	return line
	
def prefilter(lines, samples, d_groups, min_mean_freq=3, max_missing_rate=0.2, min_cv=0.05, by_group=True):
	nsampels = len(samples)
	d_freq = lines.freq_dict
	
	if not by_group:
		missing_rate = 1 - len(lines)/nsampels
		if missing_rate > max_missing_rate:
			return
		mean_freq = sum(d_freq.values()) / nsampels
		if mean_freq < min_mean_freq:
			return
	else:
		groups = []
		for group, _samples in d_groups.items():
			freqs = [d_freq[sample] for sample in _samples if sample in d_freq]
			non_missing = len(freqs)
			missing_rate = 1 - non_missing/len(_samples)
			mean_freq = sum(freqs) / len(_samples)
			groups += [missing_rate > max_missing_rate or mean_freq < min_mean_freq]
		if all(groups):	# all with too many missing or two low depth
			return
	freqs = [d_freq.get(sample, 0) for sample in samples]

	cv = np.std(freqs) / np.mean(freqs)
	if cv < min_cv:
		return
	kmer = lines.seq
	return kmer, freqs

def normalize(freqs, denominators, upper_limit=2):
	'''divide by total kmer or main peak depth'''
	return [min(freq/denominator, upper_limit) for freq, denominator in zip(freqs, denominators)]
def _map_stats(args):
	#kmer, freqs, denominators, samples, d_groups, test_method = args
	lines, denominators, samples, d_groups, test_method, filter_args = args
	results = prefilter(lines, samples, d_groups, **filter_args)
	if not results:
		return 
	kmer, freqs = results
	if denominators is not None:
		freqs = normalize(freqs, denominators)
	d_freq = dict(zip(samples, freqs))
	grouped = []
	groups = []
	for group, _samples in d_groups.items():
		_freqs = [d_freq[sample] for sample in _samples]
		grouped += [_freqs]
		groups += [group]
	mean_vals = [np.mean(x) for x in grouped]
	ratios = ','.join(map(str, mean_vals))
	
	xgrouped = sorted(zip(grouped, groups, mean_vals), key=lambda x: x[2])
	grouped = [x[0] for x in xgrouped]  # sorted with the same order
	groups = [x[1] for x in xgrouped]
	max_freqs = grouped[0]
	min_freqs = grouped[1]
	max_group = groups[0]
	test = test_method(max_freqs, min_freqs)
	pvalue = test.pvalue
	rc_kmer = str(Seq(kmer).reverse_complement())
	mean_vals = sorted(mean_vals, reverse=True)
	fold = mean_vals[1] / mean_vals[0]
	return kmer, freqs, max_group, pvalue, rc_kmer, fold, ratios

class JellyfishDumps:
	def __init__(self, dumpfiles, labels=None, ncpu=4, method='map', chunksize=None, **kargs):
		self.dumpfiles = dumpfiles
		self.labels = labels
		self.ncpu = ncpu
		self.method = method
		self.chunksize = chunksize
	def __len__(self):
		return len(self.dumpfiles)
	def to_matrix(self, array=False):
		ncol = len(self)
		lengths = [0] * ncol
		d_idx = dict(zip(self.dumpfiles, range(ncol)))
		d_mat = {}
		# args = []
		# multiprocessing by chromosomes
		dumps = pool_func(_to_matrix2, self.dumpfiles, self.ncpu, method=self.method)
		for seqs, freqs, dumpfile, tot in dumps:
			i = d_idx[dumpfile]
			lengths[i] = tot
			logger.info('Loading '+ dumpfile)
			for seq, freq in zip(seqs, freqs):
				if seq not in d_mat:
					if array:
						d_mat[seq] = np.array([0] * ncol)
					else:
						d_mat[seq] = [0] * ncol
				d_mat[seq][i] = freq
			del seqs, freqs
		self.lengths = lengths
		return d_mat

	def filter(self, d_mat, lengths, sgs, outfig=None, by_count=False, 
				min_freq=200, max_freq=10000, min_fold=2, baseline=1, 
				min_prop=None, max_prop=None):
		#logger.info([min_freq, max_freq, min_fold])
		#logger.info(sgs)
		tot_lens = sum(self.lengths)
		if min_prop is not None:
			min_freq = min_prop * tot_lens
			logger.info('Adjust `min_freq` to {} according to `min_prop`'.format(min_freq))
		if max_prop is not None:
			max_freq = max_prop * tot_lens
			logger.info('Adjust `max_freq` to {} according to `max_prop`'.format(max_freq))
		if min_freq > max_freq:
			raise ValueError('`min_freq` ({}) should be lower than `max_freq` ({})'.format(min_freq, max_freq))

		i = 0
		for sg in sgs:
			if len(sg) == 1:
				logger.warn('Singleton `{}` is ignored'.format(sg))
				i += 1
		if i == len(sgs):
			raise ValueError('All singletons are not allowed')
		
		d_mat2 = {}
		d_lens = OrderedDict(zip(self.labels, self.lengths))
		args = ((kmer, counts, d_lens, sgs, outfig, by_count, min_freq, max_freq, min_fold, baseline) \
					for kmer, counts in d_mat.items())
		i = 0
		tot_freqs = []
		for kmer, freqs, tot_freq in pool_func(_filter_kmer, args, self.ncpu, 
					method=self.method, chunksize=self.chunksize):
			i += 1
			if i % 10000000 == 0:
				logger.info('Processed {} kmers'.format(i))
			if freqs:
				d_mat2[kmer] = freqs
			if tot_freq:
				tot_freqs += [tot_freq]
		remain, total = len(d_mat2), len(tot_freqs)
		logger.info('After filtering, remained {} ({:.2%}) differential (freq >= {}) and {} ({:.2%}) \
candidate (freq > 0) kmers'.format(remain, remain/i, min_freq, total, total/i))
		# plot
		if outfig is not None:
			if total == 0:
				raise ValueError('0 kmer with fold > {}. Please reset the filter options.'.format(min_fold))
			logger.info('Plot ' + outfig)
			plot_histogram(tot_freqs, outfig, vline=None)
		return d_mat2

	
	def write_matrix(self, d_mat, fout):
		fout.write( '\t'.join(['kmer'] + self.labels) + '\n')
		for kmer, counts in d_mat.items():
			line = [kmer] + counts
			line = map(str, line)
			fout.write( '\t'.join(line) + '\n')
	def heatmap(self, matfile, **kargs):
		_heatmap(matfile, **kargs)
	
def _heatmap(matfile, mapfile=None, kmermapfile=None, figfmt='pdf', color=('green', 'black', 'red'), 
		heatmap_options='scale="row", key=TRUE, density.info="density", trace="none", labRow=F, main="",xlab=""'):
	if len(color) == 3:
		color = 'colorpanel(100, low="{}", mid="{}", high="{}")'.format(*color)
	elif len(color) == 2:
		color = 'colorpanel(100, low="{}", high="{}")'.format(*color)
	else:
		logger.error('Colors must 2 or 3 in size but {} ({}) got'.format(len(color), color))
	outfig = matfile+ '.'+ figfmt
	rsrc = r'''
data = read.table('{matfile}',fill=T,header=T, row.names=1, sep='\t', check.names=F)
map = read.table('{mapfile}', colClasses = "character")
kmermap = read.table('{kmermapfile}', sep='\t', ,colClasses = "character")
colors = {colors}

# subsample
size = 10000
if (nrow(data) > size) {{
	sampels = sample(1:nrow(data), size)
	data = data[sampels,]
}}

# Z scale and transpose
data = as.matrix(data)
z.scale <- function(x) (x - mean(x))/ sqrt(var(x))
data = apply(data, 1, z.scale)	# transpose

# map: sg -> color
cmap = list()
sgs = sort(unique(map[,2]))
for (i in 1:length(sgs)) {{
	color= colors[i]
	sg = sgs[i]
	cmap[[sg]] = color
}}

# map: chrom -> SG -> color
idmap = list()
for (i in 1:nrow(map)) {{
	chrom = as.character(map[i, 1])
	sg = as.character(map[i, 2])
	idmap[[chrom]] = cmap[[sg]] #paste(sg, chrom, sep='|')
}}

# map: kmer -> SG -> color
kmermap = kmermap[kmermap[,1] %in% colnames(data),]	# reduce size
kimap = list()
for (i in 1:nrow(kmermap)) {{
	kmer = as.character(kmermap[i, 1])
	sg = as.character(kmermap[i, 2])
	kimap[[kmer]] = cmap[[sg]]
}}

# map: rownames (chrom) -> color
names = vector()
for (name in rownames(data)){{
	new_name = idmap[[name]]
	names = c(names, new_name)
}}

# map: colnames (kmer) -> color
knames = vector()
for (name in colnames(data)){{
	new_name = kimap[[name]]
			if (is.null(new_name)) {{new_name=NA}}
	knames = c(knames, new_name)
}}

#rownames(data) = names

nr = nrow(data)
cex = min(1.5, 30/nr)

library("gplots")
{dev}('{outfig}')
heatmap.2(data, col={color}, {heatmap_options}, RowSideColors=names, ColSideColors=knames, cexRow = cex)
dev.off()
'''.format(matfile=matfile, mapfile=mapfile, kmermapfile=kmermapfile,
		dev=figfmt, outfig=outfig, colors=colors_r, # colors by SG
		color=color, heatmap_options=heatmap_options, )
	rsrc_file = matfile + '.R'
	with open(rsrc_file, 'w') as f:
		f.write(rsrc)
	cmd = 'Rscript ' + rsrc_file
	run_cmd(cmd, log=True)
	return outfig
		
def _filter_kmer(arg):
	(kmer, counts, d_lens, sgs, outfig, by_count, 
		min_freq, max_freq, min_fold, baseline) = arg
	labels = d_lens.keys()
	lengths = d_lens.values()
	tot = sum(counts)
	if not outfig and (tot < min_freq or tot > max_freq):
		return kmer, False, None
	# normalize
	d_counts = dict(zip(labels, counts))
	includes = []
	for sg in sgs:
		if len(sg) == 1:
			continue
		include = False
		freqs = []
		for chrs in sg:
			if len(chrs) == 1:
				chr = chrs[0]
				count = d_counts[chr]
				lens = d_lens[chr]
				freq = count/lens if not by_count else count
			else:
				count = [d_counts[chr] for chr in chrs]
				lens = [d_lens[chr] for chr in chrs]
				freq = sum(count) / sum(lens) if not by_count else sum(count)
			freqs += [freq]
		freqs = sorted(freqs, reverse=1)
		_max_freq = freqs[0]
		_min_freq = freqs[baseline]
		if 1.0 * _max_freq / (_min_freq+1e-20) >= min_fold:
			include = True
		includes += [include]
	if not all(includes):
		return kmer, False, None
	if outfig and (tot < min_freq or tot > max_freq):
		return kmer, False, tot
	return kmer, [c/l for c,l in zip(counts, lengths)], tot
		
def plot_histogram(data, outfig, step=25, xlim=99, xlabel='Kmer occurrence', ylabel='Count', vline=None):
	from matplotlib import pyplot as plt
	_min, _max = 0, max(data)
	nbins = int((_max-_min)/step)
	plt.figure(figsize=(7,5), dpi=300, tight_layout=True)
	n,bins,patches = plt.hist(data, bins=nbins)
	xlim = np.percentile(data, xlim)
	plt.xlim(_min, xlim)
#	plt.semilogy()
	plt.xlabel(xlabel, fontsize=fonts['fontsize'])
	plt.ylabel(ylabel, fontsize=fonts['fontsize'], ha='center', va='center')
	plt.tick_params(labelsize=fonts['labelsize'])
	plt.ticklabel_format(style='plain')
	if vline is not None:
		plt.axvline(vline, ls='--', c='grey')
	plt.savefig(outfig, bbox_inches='tight', dpi=300)




def run_jellyfish_dumps(seqfiles, ncpu=4, **kargs):
	dumpfiles = []
	args = ((seqfile, kargs) for seqfile in seqfiles)
	for dumpfile in pool_func(_run_jellyfish_dump, args, ncpu):
		dumpfiles += [dumpfile]
	return dumpfiles
def _run_jellyfish_dump(arg):
	'''wrapper for pool_func'''
	seqfile, kargs = arg
	return run_jellyfish_dump(seqfile, **kargs)
def run_jellyfish_dump(seqfile, threads=4, k=17, prefix=None, lower_count=2, method='jellyfish', overwrite=False):
	if isinstance(seqfile, (list, tuple, set)):
		_seqfile = list(seqfile)[0]
		seqfile = ' '.join(seqfile)
	else: # str
		_seqfile = seqfile
		if prefix is None:
			prefix = seqfile
		
	output = '{prefix}_{KMER}.fa'.format(KMER=k, prefix=prefix)
	ckp_file = output + '.ok'
	
	if not overwrite and check_ckp(ckp_file):
		pass
	else:
		xcat = 'zcat' if is_gz(_seqfile) else 'cat'
		cmd = '{XCAT} {seqfile} | jellyfish count -t {NCPU} -m {KMER} -s 100000000  --canonical /dev/stdin -o "{prefix}_{KMER}.jf" && touch "{prefix}_{KMER}.jf.ok" &&\
		jellyfish histo -h 100000 -t {NCPU} -o {prefix}_{KMER}.histo {prefix}_{KMER}.jf &&\
	jellyfish dump -c -o "{prefix}_{KMER}.fa" "{prefix}_{KMER}.jf" -L {lower_count} && touch "{ckp_file}" \
	&& rm "{prefix}_{KMER}.jf"'.format(
			XCAT=xcat, seqfile=seqfile, NCPU=threads, KMER=k, prefix=prefix, 
				ckp_file=ckp_file, lower_count=lower_count)
		run_cmd(cmd, log=True)
	return output
def run_kmc_sort(seqfile, threads=4, k=25, prefix=None, outdir='.', 
		lower_count=3, upper_count=1e6, 
		fmt='fq', overwrite=False):
	upper_count = int(upper_count)
	if isinstance(seqfile, (list, tuple, set)):
		seqfile = ' '.join(seqfile)
	else: # str
		if prefix is None:
			prefix = os.path.basename(seqfile)
	output = '{outdir}/{prefix}_{KMER}.sorted'.format(outdir=outdir, KMER=k, prefix=prefix)
	ckp_file = output + '.ok'
	#print(check_ckp(ckp_file))
	if not overwrite and check_ckp(ckp_file):
		pass
	else:
		cmd = 'realpath {seqfile} > {outdir}/{prefix}.list && cd {outdir} &&\
	mkdir -p {prefix}_{KMER} && \
	kmc -m58 -t{NCPU} -k{KMER} -ci{lower_count} -cs{upper_count} -{format} @{prefix}.list {prefix}_{KMER} {prefix}_{KMER} > {prefix}_{KMER}.stats &&\
        kmc_tools transform {prefix}_{KMER} sort {prefix}_{KMER}.sorted && \
		kmc_tools transform {prefix}_{KMER} histogram {prefix}_{KMER}.histo -cx100000 &&\
		rm {prefix}_{KMER}.kmc_* {prefix}_{KMER} -r && touch {ckp_file}'.format(
		seqfile=seqfile, NCPU=threads, KMER=k, prefix=prefix, outdir=outdir, format=fmt,
				ckp_file=ckp_file, lower_count=lower_count, upper_count=upper_count)
		run_cmd(cmd, log=True)
	return output

def kmc_matrix(dbs, samples=None, outMat=None, bin='matrixer',  overwrite=False):
	ckp_file = outMat + '.ok'
	if not overwrite and check_ckp(ckp_file):
		pass
	else:
		input = outMat + '.input'
		cmd = 'realpath {} > {}'.format(' '.join(dbs), input)
		run_cmd(cmd, log=True)
		
		with open(outMat, 'w') as fout:
			if samples is not None and len(samples) == len(dbs):
				fout.write('\t'.join(['kmer'] + samples)+'\n')
			else:
				fout.write('')
		dbfiles = ['{}.kmc_*'.format(db) for db in dbs]
		cmd = '{bin} {input} | pigz >> {output} && touch {ckp_file} && xrm {dbfiles}'.format(
			bin=bin, input=input, output=outMat, ckp_file=ckp_file, dbfiles=' '.join(dbfiles))
		run_cmd(cmd, log=True)
def run_kmc_dict(d_seqfiles, outdir='.', overwrite=False, **kargs):
	dumpfiles = []
	histofiles = []
	for sample, seqfile in d_seqfiles.items():
		#logger.info('Count kmers of ' + sample)
		dumpfile = run_kmc_sort(seqfile, prefix=sample, outdir=outdir, overwrite=overwrite, **kargs)
		dumpfiles += [dumpfile]
		histofile = os.path.splitext(dumpfile)[0] + '.histo'
		histofiles += [histofile]
		histofiles += [histofile]
	outMat = '{}/kmer.matrix.gz'.format(outdir)
	#logger.info('Pick kmer matrix table')
	kmc_matrix(dumpfiles, samples=list(d_seqfiles.keys()), outMat=outMat, overwrite=overwrite)
	total_kmers = [] #[JellyfishHisto(histofile).get_total_kmers() for histofile in histofiles]
	return outMat, total_kmers
	
def run_jellyfish_dict(d_seqfiles, outdir='.', **kargs):
	dumpfiles = []
	histofiles = []
	for sample, seqfile in d_seqfiles.items():
		prefix = '{}/{}'.format(outdir, sample)
		dumpfile = run_jellyfish_dump(seqfile, prefix=prefix, **kargs)
		histofile = os.path.splitext(dumpfile)[0] + '.histo'
		histofiles += [histofile]
		dumpfile = sort_dump(dumpfile, sample, tmpdir=outdir, **kargs)
		dumpfiles += [dumpfile]
	merged_dumpfile = '{}/merged.dumps'.format(outdir)
	merge_dumps(dumpfiles, merged_dumpfile, tmpdir=outdir, **kargs)
	total_kmers = [JellyfishHisto(histofile).get_total_kmers() for histofile in histofiles]
	return merged_dumpfile, total_kmers
	
def sort_dump(dumpfile, sample, tmpdir='/tmp', threads=10, memory='20G', overwrite=False, **kargs):
	output = dumpfile + '.sorted'
	ckp_file = output + '.ok'
	if not overwrite and check_ckp(ckp_file):
		pass
	else:
		cmd = '''sort {dumpfile} -k 1 --parallel={threads} -S {memory} -T {tmpdir} | \
			awk '{{print $0" {sample}"}}' > {dumpfile}.sorted && touch {ckp_file} && rm {dumpfile}'''.format(
				dumpfile=dumpfile, sample=sample, ckp_file=ckp_file, tmpdir=tmpdir, 
				threads=threads, memory=memory)
		run_cmd(cmd, log=True)
	return output
def merge_dumps(dumpfiles, merged_dumpfile, tmpdir='/tmp', threads=10, memory='20G', overwrite=False, **kargs):
	dumpfiles = ' '.join(dumpfiles)
	ckp_file = merged_dumpfile + '.ok'
	if not overwrite and check_ckp(ckp_file):
		pass
	else:
		cmd = 'sort -k 1 --parallel={threads} -S {memory} -T {tmpdir} -m {dumpfiles} > {merged_dumpfile} \
	&& touch {ckp_file} && rm {merged_dumpfile} '.format(
		dumpfiles=dumpfiles, merged_dumpfile=merged_dumpfile, ckp_file=ckp_file, tmpdir=tmpdir, 
				threads=threads, memory=memory)
		run_cmd(cmd, log=True)
if __name__ == '__main__':
	pass
