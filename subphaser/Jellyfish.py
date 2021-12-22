import sys,os
from collections import OrderedDict
import numpy as np
from Bio.Seq import Seq
from xopen import xopen as open
from .RunCmdsMP import run_cmd, logger, pool_func
from .small_tools import is_gz, mk_ckp, check_ckp
from .fonts import fonts
from .colors import colors_r

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
			seq, freq = line.strip().split()
			yield JellyfishDumpRecord(seq=seq, freq=freq)
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

	def heatmap(self, matfile, mapfile=None, kmermapfile=None, figfmt='pdf', color=('green', 'black', 'red'), 
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
map = read.table('{mapfile}', sep='\t', ,colClasses = "character")
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
cex = min(2, 33/nr)

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
		
def plot_histogram(data, outfig, step=25, xlim=99, xlabel='Kmer Frequency', ylabel='Count', vline=None):
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
	


class JellyfishDumpRecord(object):
	def __init__(self, seq, freq):
		self.seq = seq
		self.freq = int(freq)
	@property
	def iseq(self):
		return seq2int(self.seq)

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
def run_jellyfish_dump(seqfile, threads=4, k=17, prefix=None, lower_count=2, overwrite=False):
	if isinstance(seqfile, (list, tuple, set)):
		seqfile = ' '.join(seqfile)
	else: # str
		if prefix is None:
			prefix = seqfile
	output = '{prefix}_{KMER}.fa'.format(KMER=k, prefix=prefix)
	ckp_file = output + '.ok'
	if not overwrite and check_ckp(ckp_file):
		pass
	else:
		xcat = 'zcat' if is_gz(seqfile) else 'cat'
		cmd = '{XCAT} {seqfile} | jellyfish count -t {NCPU} -m {KMER} -s 100000000  --canonical /dev/stdin -o "{prefix}_{KMER}.jf" && \
	jellyfish dump -c -o "{prefix}_{KMER}.fa" "{prefix}_{KMER}.jf" -L {lower_count} && touch "{ckp_file}" \
	&& rm "{prefix}_{KMER}.jf"'.format(
			XCAT=xcat, seqfile=seqfile, NCPU=threads, KMER=k, prefix=prefix, 
				ckp_file=ckp_file, lower_count=lower_count)
		run_cmd(cmd, log=True)
	return output

if __name__ == '__main__':
	pass
