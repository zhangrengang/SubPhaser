import sys,os
import numpy as np
from Bio.Seq import Seq
from xopen import xopen as open
from .RunCmdsMP import run_cmd, logger
from .small_tools import is_gz

class JellyfishDump(object):
	def __init__(self, column):
		self.column = column
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.column):
			seq, freq = line.strip().split()
			yield JellyfishDumpRecord(seq=seq, freq=freq)
	def get_kmers(self, kmer_set):
		d_kmer = {}
		for rc in self.parse():
			if rc.seq in kmer_set:
				d_kmer[rc.seq] = rc.freq
		return d_kmer
class JellyfishDumps:
	def __init__(self, dumpfiles, labels=None):
		self.dumpfiles = dumpfiles
		self.labels = labels
	def __len__(self):
		return len(self.dumpfiles)
	def to_matrix(self):
		d_mat = {}
		lengths = [0] * len(self)
		for i, dumpfile in enumerate(self.dumpfiles):
			logger.info('Loading '+ dumpfile)
			for j, line in enumerate(JellyfishDump(dumpfile)):
				if line.seq  not in d_mat:
					d_mat[line.seq] = [0] * len(self)
				d_mat[line.seq][i] = line.freq
				lengths[i] += line.freq
		self.lengths = lengths
		return d_mat
	def filter(self, d_mat, lengths, sgs, d_targets={}, 
				min_freq=200, max_freq=10000, min_fold=2,
				min_prop=None, max_prop=None):
		logger.info([min_freq, max_freq, min_fold])
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

		d_mat2 = {}
		d_lens = dict(zip(self.labels, self.lengths))
		i = 0
		for kmer, counts in d_mat.items():
			i += 1
			if i % 1000000 == 0:
				logger.info('Processed {} kmers'.format(i))
#			rv_kmer = str(Seq(kmer).reverse_complement())
#			if rv_kmer in d_mat:
#				logger.warn(kmer+'is not canonical representation')
			tot = sum(counts)
			if tot < min_freq or tot > max_freq:
				continue
			# normalize
			d_counts = dict(zip(self.labels, counts))
			includes = []
			for sg in sgs:
				include = False
				freqs = []
				for chrs in sg:
					chrs = [d_targets.get(chr, chr) for chr in chrs]
					if len(chrs) == 1:
						chr = chrs[0]
						count = d_counts[chr]
						lens = d_lens[chr]
						freq = count/lens
					else:
						count = [d_counts[chr] for chr in chrs] 
						lens = [d_lens[chr] for chr in chrs]
						freq = sum(count) / sum(lens)
					freqs += [freq]
				_min_freq = min(freqs)
				_max_freq = max(freqs)
				if 1.0 * _max_freq / (_min_freq+1e-9) >= min_fold:
					include = True
				includes += [include]
			if not all(includes):
				continue
#			logger.info(tot)
			d_mat2[kmer] = [c/l for c,l in zip(counts, self.lengths)]
		return d_mat2
	def write_matrix(self, d_mat, fout):
		fout.write( '\t'.join(['kmer'] + self.labels) + '\n')
		for kmer, counts in d_mat.items():
			line = [kmer] + counts
			line = map(str, line)
			fout.write( '\t'.join(line) + '\n')

	def heatmap(self, matfile, figfmt='pdf', color=('green', 'black', 'red'), 
			heatmap_options='scale="row", key=TRUE, density.info="density", trace="none", labRow=F, main="",xlab=""'):

		if len(color) == 3:
			color = 'colorpanel(100, low="{}", mid="{}", high="{}")'.format(*color)
		elif len(color) == 2:
			color = 'colorpanel(100, low="{}", high="{}")'.format(*color)
		else:
			logger.error('Colors must 2 or 3 in size but {} ({}) got'.format(len(color), color))
		outfig = matfile+ '.'+ figfmt
		rsrc = '''
data = read.table('{matfile}',fill=T,header=T, row.names=1, sep='\t', check.names=F)
data = as.matrix(data)
#data = log(data+1e-9)
library("gplots")
if (nrow(data) > 10000) {{
	data = data[1:10000, ]
	print('Only use the first 10000 kmers to plot.')
}}
{dev}('{outfig}')
heatmap.2(data, col={color}, {heatmap_options})
dev.off()
'''.format(matfile=matfile, dev=figfmt, outfig=outfig, color=color, heatmap_options=heatmap_options)
		rsrc_file = matfile + '.R'
		with open(rsrc_file, 'w') as f:
			f.write(rsrc)
		cmd = 'Rscript ' + rsrc_file
		run_cmd(cmd, log=True)
		return outfig

class JellyfishDumpRecord(object):
	def __init__(self, seq, freq):
		self.seq = seq
		self.freq = int(freq)

def run_jellyfish_dumps(seqfiles, **kargs):
	dumpfiles = []
	for seqfile in seqfiles:
		dumpfile = run_jellyfish_dump(seqfile, **kargs)
		dumpfiles += [dumpfile]
	return dumpfiles
def run_jellyfish_dump(seqfile, ncpu=4, k=17, prefix=None, lower_count=2, overwrite=False):
	if prefix is None:
		prefix = seqfile
	output = '{prefix}_{KMER}.fa'.format(KMER=k, prefix=prefix)
	ckp_file = output + '.ok'
	if not overwrite and os.path.exists(ckp_file):
		logger.info('Check point {} exists, skipped'.format(ckp_file))
	else:
		xcat = 'cat'
		if is_gz(seqfile):
			return 'zcat'
		cmd = '{XCAT} "{seqfile}" | jellyfish count -t {NCPU} -m {KMER} -s 100000000  --canonical /dev/stdin -o "{prefix}_{KMER}.jf" && \
	jellyfish dump -c -o "{prefix}_{KMER}.fa" "{prefix}_{KMER}.jf" -L {lower_count} && touch "{ckp_file}" \
	&& rm "{prefix}_{KMER}.jf"'.format(
		XCAT=xcat, seqfile=seqfile, NCPU=ncpu, KMER=k, prefix=prefix, ckp_file=ckp_file, lower_count=lower_count)
		run_cmd(cmd, log=True)
	return output

if __name__ == '__main__':
	pass
