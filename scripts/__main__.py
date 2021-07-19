import sys,os
import argparse
import json
from collections import OrderedDict, Counter
from xopen import xopen as open
from Bio import SeqIO
from .split_genomes import split_genomes
from .Jellyfish import run_jellyfish_dumps, JellyfishDumps
from .small_tools import mkdirs, rmdirs
from .RunCmdsMP import logger

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i', "-genomes", dest='genomes', nargs='+', metavar='GENOME', required=True,
					help="Input genome sequences in fasta format [required]")
	parser.add_argument('-c', '-sg_cfgs', dest='sg_cfgs', nargs='+', required=True, metavar='CFGFILE',
					help="Subgenomes config file (one homologous group per line); \
this chromosome set is for identifying differential kmers [required]")
	parser.add_argument('-labels', nargs='+', type=str, metavar='LABEL',
					help="For multiple genomes, provide prefix labels for each genome sequence\
to avoid conficts among chromosome id \
[default: '1- 2- ... n-']")
	parser.add_argument('-no_label', action="store_true", default=False,
					help="Do not use prefix labels for genome sequences as there is \
no conficts among chromosome id [default: %(default)s]")
	parser.add_argument('-pre', '-prefix', default=None, dest='prefix', metavar='STR',
                    help="Prefix for output directories [default=%(default)s]")
	parser.add_argument('-o', '-outdir', default='phase-results', dest='outdir', metavar='DIR',
					help="Output directory [default=%(default)s]")
	parser.add_argument('-tmpdir', default='tmp', type=str, metavar='DIR',
					help="Temporary directory [default=%(default)s]")
	parser.add_argument("-target", default=None, type=str, metavar='FILE',
                    help="Target chromosomes to output; id mapping is allowed; \
this chromosome set is for cluster and phase \
[default: the chromosome set as `-sg_cfgs`]")
	parser.add_argument("-sep", default="|", type=str, metavar='STR',
                    help="Seperator for chromosome ID")
	parser.add_argument('-k', type=int, default=15, metavar='INT',
					 help="Length of kmer [default=%(default)s]")
	parser.add_argument('-ncpu', type=int, default=16, metavar='INT',
					 help="Number of threads [default=%(default)s]")
	parser.add_argument('-lower_count', type=int, default=3, metavar='INT',
					 help="Don't output k-mer with count < lower-count [default=%(default)s]")

	parser.add_argument('-q', '-min_freq', type=int, default=200, metavar='INT', dest='min_freq',
					 help="Minimum total count for each kmer; will not work \
if `-min_prop` is specified [default=%(default)s]")
	parser.add_argument('-min_prop', type=float, default=None, metavar='FLOAT',
					help="Minimum total proportion (< 1) for each kmer [default=%(default)s]")

	parser.add_argument('-max_freq', type=int, default=1e9, metavar='INT',
					 help="Maximum total count for each kmer; will not work \
if `-max_prop` is specified [default=%(default)s]")
	parser.add_argument('-max_prop', type=float, default=None, metavar='FLOAT',
                    help="Maximum total proportion (< 1) for each kmer [default=%(default)s]")

	parser.add_argument('-f', '-min_fold', type=float, default=2, metavar='FLOAT', dest='min_fold',
					 help="Minimum fold [default=%(default)s]")

	parser.add_argument("-figfmt", default='pdf', type=str, metavar='STR',
					choices=['pdf', 'png', 'tiff', 'jpeg', 'bmp'],
					help="Format of figures [default=%(default)s]")
	parser.add_argument('-heatmap_colors', nargs='+', default=('green', 'black', 'red'), metavar='COLOR',
					help="Color panel (2 or 3 colors) for heatmap plot [default=%(default)s]")
	parser.add_argument('-heatmap_options', metavar='STR',
					default="Rowv=T,Colv=T,scale='row',key=T,key.title=NA,density.info='density',trace='none',labRow=F,main=NA,xlab=NA,margins=c(5,5)",
					help='Options for heatmap plot (see more in R shell with `?heatmap.2` of `gplots` package) [default="%(default)s"]')
	parser.add_argument('-clean', action="store_true", default=False,
                    help="Remove the temporary directory [default: %(default)s]")	
	parser.add_argument('-overwrite', action="store_true", default=False,
                    help="Overwrite even if check point files existed [default: %(default)s]")
	args = parser.parse_args()
	if args.prefix is not None:
		args.prefix = args.prefix.replace('/', '_')
		args.outdir = args.prefix + args.outdir
		args.tmpdir = args.prefix + args.tmpdir
	return args

class Pipeline:
	def __init__(self, genomes, sg_cfgs,
				labels=None, **kargs):
		self.genomes = genomes
		self.sg_cfgs = sg_cfgs
		self.__dict__.update(**kargs)
#		for key, value in kargs.items():
#			setattr(self, key, value)
		# labels
		if labels is None:
			if len(genomes) == 1 or self.no_label:
				self.labels = [''] * (len(genomes))
			else:
				self.labels = ['{}-'.format(i, ) for i in range(len(genomes))]
		else:
			self.labels = labels
		# config + label
		if len(self.labels) == len(self.sg_cfgs):
			cfg_labels = self.labels
		else:
			cfg_labels = [None] * len(self.sg_cfgs)
		self.sgs, self.chrs, self.nsg = [], [], 0
		for sgcfg, label in zip(self.sg_cfgs, cfg_labels):
			sgcfg = SGConfig(sgcfg, prefix=label, sep=self.sep)
			self.sgs += sgcfg.sgs
			self.chrs += sgcfg.chrs
			self.nsg += sgcfg.nsg
		# re-labels
		if self.no_label:
			self.labels = [''] * (len(genomes))

		self.kargs = kargs
	def run(self):
		# check 
		check_duplicates(self.genomes)
		check_duplicates(self.labels)
		# mkdir
		mkdirs(self.outdir)
		mkdirs(self.tmpdir)
		self.outdir += '/'
		self.tmpdir += '/'
		if self.prefix is not None:
			self.outdir = self.outdir + self.prefix
			self.tmpdir = self.tmpdir + self.prefix

		# split genome
		logger.info('Target chromosomes: {}'.format(self.chrs))
	#	logger.info('Splitting genomes by chromosom')
	#	ckp_file = '{}/split.ok'.format(self.tmpdir)
	#	if os.path.exists(ckp_file):
	#		logger.info('Check point `{}` exists, skipped'.format(ckp_file))
	#	else:
		d_targets = parse_idmap(self.target)
		chromfiles, labels, d_targets = split_genomes(self.genomes, self.labels, 
					self.chrs, self.tmpdir, d_targets=d_targets, sep=self.sep)
	#		os.mknod(ckp_file)
		logger.info('Split chromosomes {} with {}'.format(chromfiles, labels))
		logger.info('ID map: {}'.format(d_targets))
		if len(chromfiles) == 0:
			raise ValueError('Only 0 chromosome remained after filtering. Please check the inputs.')

		# jellyfish
		logger.info('Counting kmer by jellyfish')
		dumpfiles = run_jellyfish_dumps(chromfiles, k=self.k, ncpu=self.ncpu, lower_count=self.lower_count,
						overwrite=self.overwrite)

		# matrix
		logger.info('Loading kmer matrix of jellyfish')
		dumps = JellyfishDumps(dumpfiles, labels)
		matfile = '{}kmer_k{}_q{}_f{}.mat'.format(self.outdir, self.k, self.min_freq, self.min_fold)
		ckp_file = '{}{}.ok'.format(self.tmpdir, os.path.basename(matfile))
		if self.overwrite or not os.path.exists(ckp_file):
			d_mat = dumps.to_matrix()
			logger.info('{} kmers in total'.format(len(d_mat)))
			lengths = dumps.lengths
			logger.info('Filtering')
			d_mat = dumps.filter(d_mat, lengths, self.sgs, d_targets=d_targets, 
					min_fold=self.min_fold, 
					min_freq=self.min_freq, max_freq=self.max_freq,
					min_prop=self.min_prop, max_prop=self.max_prop)
			logger.info('{} kmers remained after filtering'.format(len(d_mat)))
			if len(d_mat) == 0:
				raise ValueError('Only 0 kmer remained after filtering. Please reset the options.')
			with open(matfile, 'w') as fout:
				dumps.write_matrix(d_mat, fout)
			mk_ckp(ckp_file)
		else:
			logger.info('Check point `{}` exists, skipped'.format(ckp_file))	
		# heatmap
		outfig = dumps.heatmap(matfile, figfmt=self.figfmt, color=self.heatmap_colors, heatmap_options=self.heatmap_options)
		
		# PCA

		# kmeans
	
		# circos
	
		# clean
		if self.clean:
			rmdirs(self.tmpdir)
def mk_ckp(ckgfile, data=None):
	with open(ckgfile, 'w') as f:
		if data is not None:
			json.dump(data, f)
	logger.info('New check point file: `{}`'.format(ckgfile))

def parse_idmap(mapfile=None):
	if not mapfile:
		return None
	d_map = {}
	for line in open(mapfile):
		line = line.strip().split('#')[0]
		if not line:
			continue
		temp = line.split()
		old_id = temp[0]
		try: new_id = temp[1]
		except IndexError: new_id = old_id.split('|')[-1]
		d_map[old_id] = new_id
	return d_map
def check_duplicates(lst):
	count = Counter(lst)
	duplicates = {v: c for v,c in count.items() if c>1}
	if duplicates:
		raise ValueError('Duplicates detected: {}'.format(duplicates))

class SGConfig:
	def __init__(self, sgcfg, **kargs):
		self.sgcfg = sgcfg
		self.kargs = kargs
		self.sgs = list(self)
	def __iter__(self):
		return self._parse()
	def _parse(self):
		self.nsgs = []
		self.nsg = 0
		self.chrs = []
		for line in open(self.sgcfg):
			temp = line.split('#')[0].strip().split()
			if not temp:
				continue
			chrs = [list(map(lambda x: add_prefix(x, **self.kargs), v.split(','))) 
						for v in temp]
			self.nsgs += [len(chrs)]
			if self.nsg == 0:
				self.nsg = len(chrs)
			if len(chrs) != self.nsg:
				logger.warn('Number of subgenome is different in line {}: \
{} in this line but {} in previous line'.format(temp, len(chrs), self.nsg))
			for xchr in chrs:
				for xxchr in xchr:
					self.chrs += [xxchr]
			yield chrs
		self.nsg = max(self.nsgs)

def add_prefix(val, prefix=None, sep='|'):
	if prefix:
		vals = ['{}{}'.format(prefix, v) for v in val.split(sep)]
		return ''.join(val)
	else:
		return val

def main():
	args = makeArgparse()
	logger.info('ARGS: {}'.format(args.__dict__))
	pipeline = Pipeline(**args.__dict__)
	pipeline.run()


if __name__ == '__main__':
	main()

