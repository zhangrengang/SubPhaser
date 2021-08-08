import sys,os
import argparse
import shutil
import json
import multiprocessing
from collections import OrderedDict, Counter
from xopen import xopen as open
from Bio import SeqIO
from . import Seqs
from .Jellyfish import run_jellyfish_dumps, JellyfishDumps
from .Cluster import Cluster
from . import LTR
from .LTR import LTRpipelines
from . import Stats
from . import Circos
from .Circos import circos_plot
from .small_tools import mkdirs, rmdirs, mk_ckp, check_ckp
from .RunCmdsMP import logger

bindir = os.path.dirname(os.path.realpath(__file__))
NCPU = multiprocessing.cpu_count()

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	# input
	group_in = parser.add_argument_group('Input', 'Input genome and config files')
	group_in.add_argument('-i', "-genomes", dest='genomes', nargs='+', metavar='GENOME', required=True,
					help="Input genome sequences in fasta format [required]")
	group_in.add_argument('-c', '-sg_cfgs', dest='sg_cfgs', nargs='+', required=True, metavar='CFGFILE',
					help="Subgenomes config file (one homologous group per line); \
this chromosome set is for identifying differential kmers [required]")
	group_in.add_argument('-labels', nargs='+', type=str, metavar='LABEL',
					help="For multiple genomes, provide prefix labels for each genome sequence \
to avoid conficts among chromosome id \
[default: '1- 2- ... n-']")
	group_in.add_argument('-no_label', action="store_true", default=False,
					help="Do not use default prefix labels for genome sequences as there is \
no conficts among chromosome id [default: %(default)s]")
	group_in.add_argument("-target", default=None, type=str, metavar='FILE',
					help="Target chromosomes to output; id mapping is allowed; \
this chromosome set is for cluster and phase \
[default: the chromosome set as `-sg_cfgs`]")
	group_in.add_argument("-sep", default="|", type=str, metavar='STR',
					help='Seperator for chromosome ID [default="%(default)s"]')
	# output
	group_out = parser.add_argument_group('Output')
	group_out.add_argument('-pre', '-prefix', default=None, dest='prefix', metavar='STR',
					help="Prefix for output [default=%(default)s]")
	group_out.add_argument('-o', '-outdir', default='phase-results', dest='outdir', metavar='DIR',
					help="Output directory [default=%(default)s]")
	group_out.add_argument('-tmpdir', default='tmp', type=str, metavar='DIR',
					help="Temporary directory [default=%(default)s]")
	# kmer
	group_kmer = parser.add_argument_group('Kmer', 'Options to count and filter kmers')
	group_kmer.add_argument('-k', type=int, default=15, metavar='INT',
					 help="Length of kmer [default=%(default)s]")
	group_kmer.add_argument('-f', '-min_fold', type=float, default=2, metavar='FLOAT', dest='min_fold',
					help="Minimum fold [default=%(default)s]")
	group_kmer.add_argument('-q', '-min_freq', type=int, default=200, metavar='INT', dest='min_freq',
					 help="Minimum total count for each kmer; will not work \
if `-min_prop` is specified [default=%(default)s]")
	group_kmer.add_argument('-lower_count', type=int, default=3, metavar='INT',
					 help="Don't output k-mer with count < lower-count [default=%(default)s]")
	group_kmer.add_argument('-min_prop', type=float, default=None, metavar='FLOAT',
					help="Minimum total proportion (< 1) for each kmer [default=%(default)s]")
	group_kmer.add_argument('-max_freq', type=int, default=1e9, metavar='INT',
					help="Maximum total count for each kmer; will not work \
if `-max_prop` is specified [default=%(default)s]")
	group_kmer.add_argument('-max_prop', type=float, default=None, metavar='FLOAT',
					help="Maximum total proportion (< 1) for each kmer [default=%(default)s]")
	group_kmer.add_argument('-low_mem', action="store_true", default=None,
					help="Low MEMory but slower [default: True if genome size > 3G, else False]")
	# cluster
	group_clst = parser.add_argument_group('Cluster', 'Options for cluster to phase')
	group_clst.add_argument('-replicates', type=int, default=1000, metavar='INT',
					help="Number of replicates for bootstrap [default=%(default)s]")
	group_clst.add_argument('-jackknife', type=float, default=80, metavar='FLOAT',
					help="Percent of kmers to resample for bootstrap [default=%(default)s]")
	group_clst.add_argument('-min_pval', type=float, default=0.05, metavar='FLOAT',
					help="Minimum P value for hypothesis test [default=%(default)s]")

	group_clst.add_argument("-figfmt", default='pdf', type=str, metavar='STR',
					choices=['pdf', 'png', 'tiff', 'jpeg', 'bmp'],
					help="Format of figures [default=%(default)s]")
	group_clst.add_argument('-heatmap_colors', nargs='+', default=('green', 'black', 'red'), metavar='COLOR',
					help="Color panel (2 or 3 colors) for heatmap plot [default=%(default)s]")
	group_clst.add_argument('-heatmap_options', metavar='STR',
					default="Rowv=T,Colv=T,scale='col',dendrogram='row',labCol=F,trace='none',\
key=T,key.title=NA,density.info='density',main=NA,xlab=NA,margins=c(5,6)",
					help='Options for heatmap plot (see more in R shell with `?heatmap.2` \
of `gplots` package) [default="%(default)s"]')
	# LTR
	group_ltr = parser.add_argument_group('LTR', 'Options for LTR analyses')
	group_ltr.add_argument('-disable_ltr', action="store_true", default=False,
					help="Disable this step (this step is time-consuming for large genome)\
 [default: %(default)s]")

	group_ltr.add_argument("-ltr_detectors", nargs='+', metavar='PROG', 
					default=['ltr_finder', 'ltr_harvest'], 
					choices=['ltr_finder', 'ltr_harvest'],
					help="Programs to detect LTR-RTs [default=%(default)s]")
	group_ltr.add_argument("-ltr_finder_options", metavar='STR',
					default='-w 2 -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.85',
					help='Options for `ltr_finder` to identify LTR-RTs (see more with \
`ltr_finder -h`)[default="%(default)s"]')
	group_ltr.add_argument("-ltr_harvest_options", metavar='STR',
					default='-similar 85 -vic 10 -seed 20 -seqids yes -minlenltr 100 \
-maxlenltr 7000 -mintsd 4 -maxtsd 6',
					help='Options for `gt ltrharvest` to identify LTR-RTs (see more with \
`gt ltrharvest -help`) [default="%(default)s"]')
	group_ltr.add_argument("-tesorter_options", metavar='STR',
					default='-db rexdb-plant',
					help='Options for `TEsorter` to classify LTR-RTs [default="%(default)s"]')
	group_ltr.add_argument('-mu', metavar='FLOAT', type=float, default=7e-9,
					help='Substitution rate for estimating age of LTR insertion \
[default=%(default)s]')

	# circos
	group_circ = parser.add_argument_group('Circos', 'Options for circos plot')
	group_circ.add_argument('-disable_circos', action="store_true", default=False,
					help="Disable this step [default: %(default)s]")
	group_circ.add_argument('-window_size', type=int, default=1000000, metavar='INT',
					help="Window size for circos plot [default: %(default)s]")

	# others
	group_other = parser.add_argument_group('Other options')
	group_other.add_argument('-p', '-ncpu', type=int, default=NCPU, metavar='INT', dest='ncpu',
					 help="Maximum number of processors to use [default=%(default)s]")
	group_other.add_argument('-cleanup', action="store_true", default=False,
					help="Remove the temporary directory [default: %(default)s]")	
	group_other.add_argument('-overwrite', action="store_true", default=False,
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
		# thread for pre process, i.e. jellyfish
		self.threads = max(1, self.ncpu//len(set(self.chrs)))

	def run(self):
		# check 
		check_duplicates(self.genomes)
		check_duplicates(self.labels)
		# mkdir
		self.outdir = os.path.realpath(self.outdir)
		self.tmpdir = os.path.realpath(self.tmpdir)
		mkdirs(self.outdir)
		mkdirs(self.tmpdir)
		self.outdir += '/'
		self.tmpdir += '/'
		self._outdir = self.outdir
		if self.prefix is not None:
			self.outdir = self.outdir + self.prefix
			self.tmpdir = self.tmpdir + self.prefix

		# split genome
		logger.info('Target chromosomes: {}'.format(self.chrs))
		logger.info('Splitting genomes by chromosome into `{}`'.format(self.tmpdir))
		ckp_file = self.mk_ckpfile('split')
		ckp = check_ckp(ckp_file)
		if isinstance(ckp, list) and len(ckp) ==4 and not self.overwrite:
			chromfiles, labels, d_targets, d_size = data = ckp
		else:
			d_targets = parse_idmap(self.target)
			data = chromfiles, labels, d_targets, d_size = Seqs.split_genomes(self.genomes, self.labels, 
					self.chrs, self.tmpdir, d_targets=d_targets, sep=self.sep)
			mk_ckp(ckp_file, *data)
		self.chromfiles = chromfiles
#`		logger.info('Split chromosomes {} with {}'.format(chromfiles, labels))
		logger.info('ID map: {}'.format(d_targets))
		if len(chromfiles) == 0:
			raise ValueError('Only 0 chromosome remained after filtering. Please check the inputs.')

		# auto set pool method for multiprocessing
		genome_size = sum(d_size.values())
		logger.info('Genome size: {:,}'.format(genome_size))
		self.pool_method = 'map'
		self.chunksize = None
		if self.low_mem is None and genome_size > 3e9:
			self.pool_method = 'imap_unordered'
			self.chunksize = 100
			logger.info('Change pool method to reduce memory')

		logger.info('Step: Kmer Count')
		# jellyfish
		logger.info('Counting kmer by jellyfish')	# multiprocessing by chromfile
		dumpfiles = run_jellyfish_dumps(chromfiles, k=self.k, ncpu=self.ncpu, lower_count=self.lower_count,
						threads=self.threads, overwrite=self.overwrite)	# ckp for each chromfile

		# matrix
		logger.info('Loading kmer matrix from jellyfish')	# multiprocessing by kmer
		dumps = JellyfishDumps(dumpfiles, labels, ncpu=self.ncpu, 
							method=self.pool_method, chunksize=self.chunksize)
		self.para_prefix = '{}k{}_q{}_f{}'.format(self.outdir, self.k, self.min_freq, self.min_fold)
		matfile = self.para_prefix + '.kmer.mat'
		ckp_file = self.mk_ckpfile(matfile)
		if self.overwrite or not check_ckp(ckp_file):
			d_mat = dumps.to_matrix()	# multiprocessing by chrom
			logger.info('{} kmers in total'.format(len(d_mat)))
			lengths = dumps.lengths
			logger.info('Filtering')	# multiprocessing by kmer
			d_mat = dumps.filter(d_mat, lengths, self.sgs, d_targets=d_targets, 
					min_fold=self.min_fold, 
					min_freq=self.min_freq, max_freq=self.max_freq,
					min_prop=self.min_prop, max_prop=self.max_prop)
			logger.info('{} kmers remained after filtering'.format(len(d_mat)))
			if len(d_mat) == 0:
				raise ValueError('Only 0 kmer remained after filtering. Please reset the filter options.')
			with open(matfile, 'w') as fout:
				dumps.write_matrix(d_mat, fout)
			mk_ckp(ckp_file)
		# heatmap	# static mechod
		outfig = dumps.heatmap(matfile, figfmt=self.figfmt, color=self.heatmap_colors, 
					heatmap_options=self.heatmap_options)
		
		# kmeans cluster
		logger.info('Step: Cluster')
		cluster = Cluster(matfile, n_clusters=self.nsg, sg_prefix='SG',
				replicates=self.replicates, jackknife=self.jackknife)
		d_sg = cluster.d_sg	# chrom -> SG
		self.sg_names = cluster.sg_names
		sg_chrs = self.para_prefix + '.subgenomes'
		logger.info('Outputing `chromosome` - `subgenome` assignments to `{}`'.format(sg_chrs))
		with open(sg_chrs, 'w') as fout:
			cluster.output_subgenomes(fout)  # multiprocessing by kmer
		# PCA
		outfig = self.para_prefix + '.pca.' + self.figfmt
		logger.info('Outputing PCA plot to `{}`'.format(outfig))
		cluster.pca(outfig)

		# specific kmers and location
		sg_kmers = self.para_prefix + '.kmers'
		logger.info('Outputing significant differiential `kmer` - `subgenome` maps to `{}`'.format(sg_kmers))
		with open(sg_kmers, 'w') as fout:	# multiprocessing by kmer
			# kmer -> SG
			d_kmers = cluster.output_kmers(fout, min_pval=self.min_pval, ncpu=self.ncpu)
		logger.info('{} significant subgenome-specific kmers'.format(len(d_kmers)//2))
	
		sg_map = self.para_prefix + '.sg.bed'
		ckp_file = self.mk_ckpfile(sg_map)
		if self.overwrite or not check_ckp(ckp_file):	# SG id should be stable
			logger.info('Outputing `coordinate` - `subgenome` maps to `{}`'.format(sg_map))
			chunksize = None if self.pool_method == 'map' else 10
			with open(sg_map, 'w') as fout:	# multiprocessing by chrom chunk
				Seqs.map_kmer3(chromfiles, d_kmers, fout=fout, k=self.k, 
					ncpu=self.ncpu, method=self.pool_method, chunksize=chunksize)
			mk_ckp(ckp_file)
		# enrich by BIN
		bins, counts = Circos.counts2matrix(sg_map, keys=self.sg_names, keycol=3, window_size=self.window_size)
		bin_enrich = self.para_prefix + '.bin.enrich'
		with open(bin_enrich, 'w') as fout:
			sg_lines = Stats.enrich_bin(fout, counts, colnames=self.sg_names, rownames=bins,
					min_pval=self.min_pval)

		# LTR
		ltr_bedlines = self.step_ltr(d_kmers) if not self.disable_ltr else []
		# circos
		if not self.disable_circos:
			self.step_circos(
				bedfile=sg_map, # chrom - coord - SG			circles 3 - n+2
				sg_lines=sg_lines, # SG ratio and enrichment	circles 1 - 2
				ltr_lines=ltr_bedlines, # LTR bed lines			circles n+3
				d_sg = d_sg, # chrom -> SG, for colors
				window_size=self.window_size)

		self.step_final()

	def mk_ckpfile(self, file):
		return '{}{}.ok'.format(self.tmpdir, os.path.basename(file))

	def step_ltr(self, d_kmers):
		# LTR
		logger.info('Step: LTR')
		tmpdir = '{}LTR'.format(self.tmpdir)
		if '-p' not in self.tesorter_options:
			self.tesorter_options += ' -p {}'.format(self.threads)
		kargs = dict(progs=self.ltr_detectors, 
				options={'ltr_finder': self.ltr_finder_options, 'ltr_harvest':self.ltr_harvest_options},
				job_args={'mode': 'local', 'retry': 3, 'cont': 1,  'tc_tasks': self.threads},
				tesorter_options=self.tesorter_options, mu=self.mu,)
		# multiprocessing by chromfile
		ltrs, ltrfiles = LTRpipelines(self.chromfiles, tmpdir=tmpdir, 
					only_ltr=True, ncpu=self.ncpu, **kargs).run()
		ltr_map = self.para_prefix + '.ltr.bed'
		ckp_file = self.mk_ckpfile(ltr_map)
		if self.overwrite or not check_ckp(ckp_file):
#		if True:
			logger.info('Outputing `coordinate` - `LTR` maps to `{}`'.format(ltr_map))
			with open(ltr_map, 'w') as fout:	# multiprocessing by LTR
				Seqs.map_kmer3(ltrfiles, d_kmers, fout=fout, k=self.k, ncpu=self.ncpu, 
								chunk=False, log=False, method=self.pool_method, chunksize=self.chunksize)
			mk_ckp(ckp_file)

		# enrich SG by LTR
		bins, counts = Circos.counts2matrix(ltr_map, keys=self.sg_names, keycol=3, window_size=100000000)
		ltr_enrich = self.para_prefix + '.ltr.enrich'
		with open(ltr_enrich, 'w') as fout:
			d_enriched = Stats.enrich_ltr(fout, counts, colnames=self.sg_names, rownames=bins, 
					min_pval=self.min_pval)
		# plot insert age
		prefix = self.para_prefix + '.ltr.insert'
		LTR.plot_insert_age(ltrs, d_enriched, prefix, mu=self.mu, figfmt=self.figfmt)

		# ltr bed
		ltr_bedlines = [ltr.to_bed() for ltr in ltrs]
		return ltr_bedlines

	def step_circos(self, *args, **kargs):
		# circos
		circos_dir = bindir+'/circos'
		wkdir = self._outdir+'/circos'
		rmdirs(wkdir)
		try:
			logger.info('Copy `{}` to `{}`'.format(circos_dir, self._outdir))
			shutil.copytree(circos_dir, wkdir)
		except FileExistsError:
			pass
		Circos.circos_plot(self.chromfiles, wkdir, *args, **kargs)

	def step_final(self):
		# cleanup
		if self.cleanup:
			rmdirs(self.tmpdir)

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

