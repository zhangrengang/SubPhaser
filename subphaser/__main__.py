import sys,os
import copy
import argparse
import shutil
import json
import math
import multiprocessing
from collections import OrderedDict, Counter
from xopen import xopen as open
from Bio import SeqIO
from . import Seqs
from . import Jellyfish
from .Jellyfish import run_jellyfish_dumps, JellyfishDumps
from .Cluster import Cluster
from . import LTR
from . import Stats
from . import Circos
from . import Blocks
from .small_tools import mkdirs, rmdirs, mk_ckp, check_ckp, test_s
from .RunCmdsMP import logger, available_memory, limit_memory
from .__version__ import version


bindir = os.path.dirname(os.path.realpath(__file__))
NCPU = multiprocessing.cpu_count()
MEM = available_memory()

def makeArgparse():
	parser = argparse.ArgumentParser( 
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='Phase and visualize subgenomes of an allopolyploid or hybrid \
based on the repetitive kmers.',
		)
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
[default: '1-, 2-, ..., n-']")
	group_in.add_argument('-no_label', action="store_true", default=False,
					help="Do not use default prefix labels for genome sequences as there is \
no confict among chromosome id [default=%(default)s]")
	group_in.add_argument("-target", default=None, type=str, metavar='FILE',
					help="Target chromosomes to output; id mapping is allowed; \
this chromosome set is for cluster and phase \
[default: the same chromosome set as `-sg_cfgs`]")
	group_in.add_argument("-sg_assigned", default=None, type=str, metavar='FILE',
					help="Provide subgenome assignments to skip k-means clustering and \
to identify subgenome-specific features \
[default=%(default)s]")
	group_in.add_argument("-sep", default="|", type=str, metavar='STR',
					help='Seperator for chromosome ID [default="%(default)s"]')
	group_in.add_argument('-custom_features', nargs='+', metavar='FASTA', default=None, 
					help="Custom features in fasta format to enrich subgenome-specific kmers, \
such as TE and gene [default: %(default)s]")
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
	group_kmer.add_argument('-baseline', type=int, default=1, 
					 help="Use sub-maximum (1) or minimum (-1) as the baseline of fold \
[default=%(default)s]")

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
	group_kmer.add_argument('-by_count', action="store_true", default=False,
					help="Calculate fold by count instead of by proportion [default=%(default)s]")
	group_kmer.add_argument('-re_filter', action="store_true", default=False,
					help="Re-filter with subset of chromosomes (subgenome assignments are expected to change)\
 [default=%(default)s]")
 
	# cluster
	group_clst = parser.add_argument_group('Cluster', 'Options for clustering to phase')
	group_clst.add_argument('-nsg', type=int, default=None, metavar='INT',
					help="Number of subgenomes (>1) [default: auto]")
	group_clst.add_argument('-replicates', type=int, default=1000, metavar='INT',
					help="Number of replicates for bootstrap [default=%(default)s]")
	group_clst.add_argument('-jackknife', type=float, default=50, metavar='FLOAT',
					help="Percent of kmers to resample for each bootstrap [default=%(default)s]")
	group_clst.add_argument('-max_pval', type=float, default=0.05, metavar='FLOAT',
					help="Maximum P value for all hypothesis tests [default=%(default)s]")
	group_clst.add_argument("-test_method", default='ttest_ind', 
					choices=['ttest_ind', 'kruskal', 'wilcoxon', 'mannwhitneyu'],
					help="The test method to identify differiential kmers[default=%(default)s]")
					
	group_clst.add_argument("-figfmt", default='pdf', type=str, 
					choices=['pdf', 'png', ], # 'svg','tiff', 'jpeg', 'bmp'],
					help="Format of figures [default=%(default)s]")
	group_clst.add_argument('-heatmap_colors', nargs='+', default=('green', 'black', 'red'), metavar='COLOR',
					help="Color panel (2 or 3 colors) for heatmap plot [default: %(default)s]")
	group_clst.add_argument('-heatmap_options', metavar='STR',
					default="Rowv=T,Colv=T,scale='col',dendrogram='row',labCol=F,trace='none',\
key=T,key.title=NA,density.info='density',main=NA,xlab='Differential kmers',margins=c(2.5,12)",
					help='Options for heatmap plot (see more in R shell with `?heatmap.2` \
of `gplots` package) [default="%(default)s"]')
	group_clst.add_argument('-just_core', action="store_true", default=False,
					help="Exit after the core phasing module [default=%(default)s]")

	# LTR
	group_ltr = parser.add_argument_group('LTR', 'Options for LTR analyses')
	group_ltr.add_argument('-disable_ltr', action="store_true", default=False,
					help="Disable this step (this step is time-consuming for large genome)\
 [default=%(default)s]")

	group_ltr.add_argument("-ltr_detectors", nargs='+',  
					default=['ltr_harvest'], 
					choices=['ltr_finder', 'ltr_harvest'],
					help="Programs to detect LTR-RTs [default=%(default)s]")
	group_ltr.add_argument("-ltr_finder_options", metavar='STR',
				#	default='-w 2 -D 20000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.8',
					default='-w 2 -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.8',
					help='Options for `ltr_finder` to identify LTR-RTs (see more with \
`ltr_finder -h`) [default="%(default)s"]')
	group_ltr.add_argument("-ltr_harvest_options", metavar='STR',
				#	default='-seqids yes -similar 80 -vic 10 -seed 20 -minlenltr 100 \
#-maxlenltr 7000 -maxdistltr 20000 -mindistltr 1000 -mintsd 4 -maxtsd 20',
					default='-seqids yes -similar 80 -vic 10 -seed 20 -minlenltr 100 \
-maxlenltr 7000 -mintsd 4 -maxtsd 6',
					help='Options for `gt ltrharvest` to identify LTR-RTs (see more with \
`gt ltrharvest -help`) [default="%(default)s"]')
	group_ltr.add_argument("-tesorter_options", metavar='STR',
					default='-db rexdb -dp2',
					help='Options for `TEsorter` to classify LTR-RTs (see more with `TEsorter -h`) \
[default="%(default)s"]')
					
	group_ltr.add_argument('-all_ltr', action="store_true", default=False,
                    help="Use all LTR-RTs identified by `-ltr_detectors` (more LTR-RTs but slower) \
[default: only use LTR as classified by `TEsorter`]")
	group_ltr.add_argument('-intact_ltr', action="store_true", default=False,
                    help="Use completed LTR-RTs classified by `TEsorter` (less LTR-RTs but faster) \
[default: the same as `-all_ltr`]")
	group_ltr.add_argument('-exclude_exchanges', action="store_true", default=False,
                    help="Exclude potential exchanged LTRs for insertion age estimation and phylogenetic trees \
[default=%(default)s]")
	group_ltr.add_argument('-shared_ltr', action="store_true", default=False,
                    help="Identify shared LTR-RTs among subgenomes (experimental) [default=%(default)s]")
	group_ltr.add_argument('-mu', metavar='FLOAT', type=float, default=13e-9,
					help='Substitution rate per year in the intergenic region, \
for estimating age of LTR insertion \
[default=%(default)s]')
	group_ltr.add_argument('-disable_ltrtree', action="store_true", default=False,
					help="Disable subgenome-specific LTR tree (this step is time-consuming when \
subgenome-specific LTR-RTs are too many, so `-subsample` is enabled by defualt)\
 [default=%(default)s]")
	group_ltr.add_argument('-subsample', type=int, default=1000, metavar='INT',
					help="Subsample LTR-RTs to avoid too many to construct a tree \
[default=%(default)s] (0 to disable)")
	group_ltr.add_argument("-ltr_domains", nargs='+', 
					default=['INT', 'RT', 'RH'], 
					choices=['GAG', 'PROT', 'INT', 'RT', 'RH', 'AP', 'RNaseH'],
					help="Domains for LTR tree (Note:  for domains identified by `TEsorter`, \
PROT (rexdb) = AP (gydb), RH (rexdb) = RNaseH (gydb)) [default: %(default)s]")
	group_ltr.add_argument("-trimal_options", metavar='STR',
					default='-automated1',
					help='Options for `trimal` to trim alignment (see more with `trimal -h`) \
[default="%(default)s"]')
	group_ltr.add_argument("-tree_method",  
					default='FastTree', 
					choices=['iqtree', 'FastTree'],
					help="Programs to construct phylogenetic trees [default=%(default)s]")

	group_ltr.add_argument("-tree_options", metavar='STR',
					default='',
					help='Options for `-tree_method` to construct phylogenetic trees \
(see more with `iqtree -h` or `FastTree -expert`) [default="%(default)s"]')
	group_ltr.add_argument("-ggtree_options", metavar='STR',
					default="branch.length='none', layout='circular'",
					help='Options for `ggtree` to show phylogenetic trees \
(see more from `https://yulab-smu.top/treedata-book`) [default="%(default)s"]')

	# circos
	group_circ = parser.add_argument_group('Circos', 'Options for circos plot')
	group_circ.add_argument('-disable_circos', action="store_true", default=False,
					help="Disable this step [default=%(default)s]")
	group_circ.add_argument('-window_size', type=int, default=1000000, metavar='INT',
					help="Window size (bp) for circos plot [default=%(default)s]")
	group_circ.add_argument('-disable_blocks', action="store_true", default=False,
					help="Disable to plot homologous blocks [default=%(default)s]")
	group_circ.add_argument("-aligner", metavar='PROG', 
					default='minimap2', 
					choices=['minimap2', 'unimap'],
					help="Programs to identify homologous blocks [default=%(default)s]")
	group_circ.add_argument("-aligner_options", metavar='STR',
					default='-x asm20 -n 10',
					help='Options for `-aligner` to align chromosome sequences [default="%(default)s"]')
	group_circ.add_argument('-min_block', type=int, default=100000, metavar='INT',
					help="Minimum block size (bp) to show [default=%(default)s]")
	group_circ.add_argument('-alt_cfgs', nargs='+', metavar='CFGFILE', default=None,
					help="An alternative config file for identifying homologous blocks \
[default=%(default)s]")
	group_circ.add_argument("-chr_ordered", default=None, type=str, metavar='FILE',
					help="Provide a chromosome order to plot circos \
[default=%(default)s]")

	# others
	group_other = parser.add_argument_group('Other options')
	group_other.add_argument('-p', '-ncpu', type=int, default=NCPU, metavar='INT', dest='ncpu',
					 help="Maximum number of processors to use [default=%(default)s]")
	group_other.add_argument('-max_memory', type=str, default=MEM, metavar='MEM', 
					 help="Maximum memory to use where limiting can be enabled. [default=%(default)s]")
	group_other.add_argument('-cleanup', action="store_true", default=False,
					help="Remove the temporary directory [default=%(default)s]")	
	group_other.add_argument('-overwrite', action="store_true", default=False,
					help="Overwrite even if check point files existed [default=%(default)s]")
	group_other.add_argument('-v', '-version', action='version', version=version)
	
	args = parser.parse_args()
	if args.prefix is not None:
		args.prefix = args.prefix.replace('/', '_')
		args.outdir = args.prefix + args.outdir
		args.tmpdir = args.prefix + args.tmpdir
	args.ltr_detectors = sorted(set(args.ltr_detectors))
	return args

class Pipeline:
	def __init__(self, genomes, sg_cfgs,
				labels=None, **kargs):
		self.genomes = genomes
		self.sg_cfgs = sg_cfgs
		self.__dict__.update(**kargs)
		
		# check 
		check_duplicates(genomes)
		check_duplicates(labels)
		
		# labels
		if labels is None:
			if len(genomes) == 1 or self.no_label:
				self.labels = [''] * (len(genomes))
			else:
				self.labels = ['{}-'.format(i+1, ) for i in range(len(genomes))]
		else:
			self.labels = labels
		# config + label
		if len(self.labels) == len(self.sg_cfgs):
			cfg_labels = self.labels
		else:
			cfg_labels = [None] * len(self.sg_cfgs)
		self.sgs, self.chrs, _nsg = [], [], 0
		for sgcfg, label in zip(self.sg_cfgs, cfg_labels):
			sgcfg = SGConfig(sgcfg, prefix=label, sep=self.sep)
			self.sgs += sgcfg.sgs
			self.chrs += sgcfg.chrs
			_nsg += sgcfg.nsg
		if not self.alt_cfgs:
			self.alt_sgs = self.sgs
		else:
			self.alt_sgs = []
			for sgcfg in self.alt_cfgs:
				sgcfg = SGConfig(sgcfg, sep=self.sep)
				self.alt_sgs += sgcfg.sgs
		if not self.nsg or self.nsg < 2:
			self.nsg = _nsg
		
		
		# re-labels
		if self.no_label:
			self.labels = [''] * (len(genomes))
		self.kargs = kargs
		# thread for pre process, i.e. jellyfish
		self.threads = max(1, self.ncpu//len(set(self.chrs)))

	def update_sgs(self, sgs, d_targets):
		sgs = copy.deepcopy(sgs)
		for sg in sgs:
			for i, chrs in enumerate(sg):
				chrs = [d_targets.get(chr, chr) for chr in chrs]
				sg[i] = chrs
		return sgs
	def parse_assigned(self, d_targets):
		d_assigned = {}
		if not self.sg_assigned:
			return d_assigned
		for line in open(self.sg_assigned):
			if line.startswith('#'):
				continue
			chr, sg = line.strip().split()[:2]
			chr = d_targets.get(chr, chr)
			d_assigned[chr] = sg
		return d_assigned
	def parse_ordered(self, d_targets):
		chroms = []
		if not self.chr_ordered:
			return chroms
		for line in open(self.chr_ordered):
			if line.startswith('#'):
				continue
			chr = line.strip().split()[0]
			chr = d_targets.get(chr, chr)
			chroms += [chr]
		return chroms
		
	def run(self):
		
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
		split = False
		ckp_file = self.mk_ckpfile('split')
		ckp = check_ckp(ckp_file)
		if isinstance(ckp, list) and len(ckp) ==4 and not self.overwrite:
			chromfiles, labels, d_targets, d_size = data = ckp
			if set(d_targets) != set(self.chrs):
				self.re_filter = True
				#if set(self.chrs)-set(d_targets):	# new chrs, have to re-split
				#	logger.info('Chromosome set is not equal to the checkpoint, re-splitting..')
				split = True
		else:
			split = True
		
		if split:
			d_targets = parse_idmap(self.target)
			outdir = '{}chromosomes/'.format(self.tmpdir)
			mkdirs(outdir)
			data = chromfiles, labels, d_targets, d_size = Seqs.split_genomes(self.genomes, self.labels, 
					self.chrs, outdir, d_targets=d_targets, sep=self.sep)
			mk_ckp(ckp_file, *data)
		labels, chromfiles = self.sort_labels(d_targets.values(), labels, chromfiles)
		logger.info('Chromosomes: {}'.format(labels))
		logger.info('Chromosome Number: {}'.format(len(labels)))
		
		self.chromfiles = chromfiles	# ordered as the genome files
		# update
		self.labels = labels
		self.sgs = self.update_sgs(self.sgs, d_targets)
		self.alt_sgs = self.update_sgs(self.alt_sgs, d_targets)
		self.sg_assigned = self.parse_assigned(d_targets)
		logger.info('CONFIG: {}'.format(self.sgs))
		self.d_chromfiles = OrderedDict(zip(labels, chromfiles))
		self.chr_ordered = self.parse_ordered(d_targets)
		# if self.chr_ordered:
			# labels = self.chr_ordered
			# chromfiles = [self.d_chromfiles[label] for chr in self.chr_ordered]
			# self.d_chromfiles = OrderedDict(zip(labels, chromfiles))
		self.d_size =  d_size
#		logger.info('Split chromosomes {} with {}'.format(chromfiles, labels))
#		logger.info('ID map: {}'.format(d_targets))
		if len(chromfiles) == 0:
			raise ValueError('0 chromosome remained after filtering. Please check the inputs.')

		# auto set pool method for multiprocessing
		genome_size = sum(d_size.values())
		logger.info('Genome size: {:,} bp'.format(genome_size))
		self.pool_method = 'map'
		self.chunksize = None
		if self.low_mem is None and genome_size > 3e9:
			self.pool_method = 'imap_unordered'
			self.chunksize = 100
			logger.info('Change pool method to reduce memory')

		logger.info('###Step: Kmer Count')
		# jellyfish
		logger.info('Counting kmer by jellyfish')	# multiprocessing by chromfile
		dumpfiles = run_jellyfish_dumps(chromfiles, k=self.k, ncpu=self.ncpu, lower_count=self.lower_count,
						threads=self.threads, overwrite=self.overwrite)	# ckp for each chromfile

		# matrix
		logger.info('Loading kmer matrix from jellyfish')	# multiprocessing by kmer
		chunksize = None if self.pool_method == 'map' else 20000
		dumps = JellyfishDumps(dumpfiles, labels, ncpu=self.ncpu, 
							method=self.pool_method, chunksize=chunksize)
		self.basename = 'k{}_q{}_f{}'.format(self.k, self.min_freq, self.min_fold)
		self.para_prefix = '{}{}'.format(self.outdir, self.basename)
		matfile = self.para_prefix + '.kmer.mat'
		ckp_file = self.mk_ckpfile(matfile)
		if self.overwrite or self.re_filter or not check_ckp(ckp_file) or not test_s(matfile):
			d_mat = dumps.to_matrix()	# multiprocessing by chrom
			kmer_count = len(d_mat)
			logger.info('{} kmers in total'.format(kmer_count))
			lengths = dumps.lengths
			logger.info('Filtering differential kmers')	# multiprocessing by kmer
			histfig = self.para_prefix + '.kmer_freq.' + self.figfmt
			d_mat = dumps.filter(d_mat, lengths, self.sgs, outfig=histfig, #d_targets=d_targets, 
					min_fold=self.min_fold, baseline=self.baseline,
					min_freq=self.min_freq, max_freq=self.max_freq,
					min_prop=self.min_prop, max_prop=self.max_prop)
			kmer_count2 = len(d_mat)
	#		logger.info('{} ({:.2%}) kmers remained after filtering'.format(kmer_count2, kmer_count2/kmer_count))
			if len(d_mat) == 0:
				raise ValueError('0 kmer remained after filtering. Please reset the filter options.')
			with open(matfile, 'w') as fout:
				dumps.write_matrix(d_mat, fout)
			mk_ckp(ckp_file)
		
		
		# kmeans cluster
		logger.info('###Step: Cluster')
		cluster = Cluster(matfile, n_clusters=self.nsg, sg_prefix='SG', sg_assigned=self.sg_assigned,
				replicates=self.replicates, jackknife=self.jackknife)
		self.d_sg = d_sg = cluster.d_sg	# chrom -> SG
		logger.info('Subgenome assignments: {}'.format(d_sg))
		self.sg_names = cluster.sg_names
		sg_chrs = self.para_prefix + '.chrom-subgenome.tsv'
		logger.info('Outputing `chromosome` - `subgenome` assignments to `{}`'.format(sg_chrs))
		with open(sg_chrs, 'w') as fout:
			cluster.output_subgenomes(fout)  # multiprocessing by kmer
			

		# specific kmers and location
		sg_kmers = self.para_prefix + '.sig.kmer-subgenome.tsv'
		logger.info('Outputing significant differiential `kmer` - `subgenome` maps to `{}`'.format(sg_kmers))
		with open(sg_kmers, 'w') as fout:	# multiprocessing by kmer
			# kmer -> SG
			d_kmers = cluster.output_kmers(fout, max_pval=self.max_pval, ncpu=self.ncpu,
							test_method=self.test_method)
		logger.info('{} significant subgenome-specific kmers'.format(len(d_kmers)//2))
		for sg, count in sorted(Counter(d_kmers.values()).items()):
			logger.info('\t{} {}-specific kmers'.format(count//2, sg))
		
		
		# heatmap	# static method
		outfig = dumps.heatmap(matfile, mapfile=sg_chrs, kmermapfile=sg_kmers,
					figfmt=self.figfmt, color=self.heatmap_colors, 
					heatmap_options=self.heatmap_options)
		
		# PCA
		outfig = self.para_prefix + '.kmer_pca.' + self.figfmt
		logger.info('Outputing PCA plot to `{}`'.format(outfig))
		cluster.pca(outfig, n_components=self.nsg)


		if self.just_core:
			self.step_final()
			logger.info('Pipeline completed early')
			return
		
		sg_map = self.para_prefix + '.subgenome.bin.count'
		ckp_file = self.mk_ckpfile(sg_map)
		if self.overwrite or self.re_filter or not check_ckp(ckp_file) or not test_s(sg_map):	
		# SG id should be stable
			logger.info('Outputing `coordinate` - `subgenome` maps to `{}`'.format(sg_map))
			chunksize = None if self.pool_method == 'map' else 10
			with open(sg_map, 'w') as fout:	# multiprocessing by chrom chunk
				Seqs.map_kmer3(chromfiles, d_kmers, fout=fout, k=self.k, 
					bin_size=10000, sg_names=self.sg_names,
					ncpu=self.ncpu, method=self.pool_method, chunksize=chunksize)
			mk_ckp(ckp_file)
		# enrich by BIN
		logger.info('Enriching subgenome by chromosome window (size: {})'.format(self.window_size))
	#	bins, counts = Circos.counts2matrix(sg_map, keys=self.sg_names, keycol=3, window_size=self.window_size)
		bins, counts = Circos.stack_matrix(sg_map, window_size=self.window_size)
	#	logger.info('Matrix loaded')
		bin_enrich = self.para_prefix + '.bin.enrich'
		bin_exchange = self.para_prefix + '.bin.group'
		with open(bin_enrich, 'w') as fout:	# multiprocessing by chrom bin
			with open(bin_exchange, 'w') as fout2:	# group exchanges
				sg_lines = Stats.enrich_bin(fout, fout2, self.d_sg, counts, colnames=self.sg_names, rownames=bins,
					max_pval=self.max_pval, ncpu=self.ncpu)
		logger.info('Output: {}'.format(bin_enrich))

		# custom
		if self.custom_features is not None:
	#		for i, feature in enumerate(self.custom_features):
			feat_map = self.para_prefix + '.custom.bin.count'
			logger.info('Mapping subgenome-specific kmers to custom features: {}'.format(self.custom_features))
			pool_method = self.pool_method
			chunksize = None if pool_method == 'map' else 10000
			with open(feat_map, 'w') as fout:	# multiprocessing by LTR
				Seqs.map_kmer3(self.custom_features, d_kmers, fout=fout, k=self.k, ncpu=self.ncpu, 
								bin_size=10000000, sg_names=self.sg_names,
								chunk=False, log=False, method=pool_method, chunksize=chunksize)
			# enrich SG by feature
			logger.info('Enriching subgenome-specific features')
			bins, counts = Circos.stack_matrix(feat_map, window_size=100000000)
			feat_enrich = self.para_prefix + '.custom.enrich'
			with open(feat_enrich, 'w') as fout:
				d_enriched, *_ = Stats.enrich_ltr(fout, self.d_sg, counts, colnames=self.sg_names, rownames=bins, 
						max_pval=self.max_pval, ncpu=self.ncpu)
			logger.info('Output: {}'.format(feat_enrich))
			
			logger.info('{} significant subgenome-specific features'.format(len(d_enriched)))
			for sg, count in sorted(Counter(d_enriched.values()).items()):
				suffix = {'shared':''}.get(sg, '-specific')
				logger.info('\t{} {}{} features'.format(count, sg, suffix))
		
		# LTR
		ltr_bedlines, enrich_ltr_bedlines = self.step_ltr(d_kmers) if not self.disable_ltr else ([],[])

		# circos
		if self.chr_ordered:
			self.chromfiles = [self.d_chromfiles[chr] for chr in self.chr_ordered]
		if not self.disable_circos:
			self.step_circos(
				bedfile=sg_map, # chrom - coord - SG			circles 3 - n+2
				sg_lines=sg_lines, # SG ratio and enrichment	circles 1 - 2
				ltr_lines=ltr_bedlines, # LTR bed lines			circles n+3
				enrich_ltr_bedlines=enrich_ltr_bedlines, 	#	circles n+3
				d_sg = d_sg, # chrom -> SG, for colors
				prefix=self.para_prefix + '.circos',
				figfmt=self.figfmt,
				window_size=self.window_size)

		self.step_final()
		logger.info('Pipeline completed')

	def mk_ckpfile(self, file):
		return '{}{}.ok'.format(self.tmpdir, os.path.basename(file))

	def step_ltr(self, d_kmers):
		# LTR
		logger.info('###Step: LTR')
		tmpdir = '{}LTR'.format(self.tmpdir)
		if ' -p ' not in self.tesorter_options:
			self.tesorter_options += ' -p {}'.format(self.ncpu)
		cont = not self.overwrite
		job_args = {'mode': 'local', 'retry': 3, 'cont': cont,  'tc_tasks': self.ncpu}
		kargs = dict(progs=self.ltr_detectors, 
				options={'ltr_finder': self.ltr_finder_options, 'ltr_harvest':self.ltr_harvest_options},
				job_args=job_args, tesorter_options=self.tesorter_options, mu=self.mu,)
		# multiprocessing by chromfile
		logger.info('Identifying LTR-RTs by {}'.format(self.ltr_detectors))
		pipeline = LTR.LTRpipeline(self.chromfiles, tmpdir=tmpdir, 
					all_ltr=self.all_ltr, intact_ltr=self.intact_ltr,
					ncpu=self.ncpu, **kargs)
		ltrs, ltrfile = pipeline.run()
		
		ltr_map = self.para_prefix + '.ltr.bin.count'
		ckp_file = self.mk_ckpfile(ltr_map)
		if self.overwrite or self.re_filter or not check_ckp(ckp_file) or not test_s(ltr_map):
			logger.info('Outputing `coordinate` - `LTR` maps to `{}`'.format(ltr_map))
			pool_method = self.pool_method
			chunksize = None if pool_method == 'map' else 2000
			with open(ltr_map, 'w') as fout:	# multiprocessing by LTR
				Seqs.map_kmer3([ltrfile], d_kmers, fout=fout, k=self.k, ncpu=self.ncpu, 
								bin_size=10000000, sg_names=self.sg_names,
								chunk=False, log=False, method=pool_method, chunksize=chunksize)
			# only mapped LTRs are output; the others are ignored
			mk_ckp(ckp_file)

		# enrich SG by LTR
		logger.info('Enriching subgenome-specific LTR-RTs')
#		bins, counts = Circos.counts2matrix(ltr_map, keys=self.sg_names, keycol=3, window_size=100000000)
		bins, counts = Circos.stack_matrix(ltr_map, window_size=100000000)
		ltr_enrich = self.para_prefix + '.ltr.enrich'
		with open(ltr_enrich, 'w') as fout:
			# ltr_id -> SG
			d_enriched, d_exchange = Stats.enrich_ltr(fout, self.d_sg, counts, 
					colnames=self.sg_names, rownames=bins, 
					max_pval=self.max_pval, ncpu=self.ncpu)
		logger.info('Output: {}'.format(ltr_enrich))
		
		logger.info('{} significant subgenome-specific LTR-RTs'.format(len(d_enriched)))
		for sg, count in sorted(Counter(d_enriched.values()).items()):
			suffix = {'shared':''}.get(sg, '-specific')
			logger.info('\t{} {}{} LTR-RTs'.format(count, sg, suffix))
		
		# shared
		if self.shared_ltr:
			ckp_file = tmpdir + '.' + self.basename + '.shared.ok'
			ckp = check_ckp(ckp_file)
			if self.overwrite or self.re_filter or not ckp:
				d_shared = self.identify_shared_ltrs(self.d_sg, self.d_chromfiles, ltrfile, d_enriched)
				mk_ckp(ckp_file, d_shared)
			else:
				d_shared, *_ = ckp
		else:
			d_shared = {}
		
		# plot insert age
		prefix = self.para_prefix + '.ltr.insert'
		enrich_ltrs = LTR.plot_insert_age(ltrs, d_enriched, prefix, shared=d_shared, 
						exclude_exchanges=self.exclude_exchanges, d_exchange=d_exchange, 
						mu=self.mu, figfmt=self.figfmt)
		
		# ltr tree
		if not self.disable_ltrtree:
			domfile = pipeline.int_seqs + '.cls.pep'
			overwrite = (self.overwrite or self.re_filter)
			prefix = tmpdir + '.' + self.basename
			tree = LTR.LTRtree(enrich_ltrs, domains=self.ltr_domains, domfile=domfile, prefix=prefix,
					trimal_options=self.trimal_options, 
					tree_method=self.tree_method, tree_options=self.tree_options, 
					subsample=self.subsample, exclude_exchanges=self.exclude_exchanges,
					ncpu=self.ncpu, overwrite=overwrite)
			job_args['cont'] = not (self.overwrite or self.re_filter)
			d_files = tree.build(job_args=job_args)
			for key, (treefile, mapfile) in d_files.items():
				key = tuple([v for v in key if v])
				outfig = '{}.{}.tree.pdf'.format(self.para_prefix, '_'.join(key))
				tree.visualize_treefile(treefile, mapfile, outfig, 
						ggtree_options=self.ggtree_options)
		# ltr bed for circos
		ltr_bedlines = [ltr.to_bed() for ltr in ltrs]
		# bin ltrs by sg
		d_enrich_ltr_bedlines = {}
		for ltr in enrich_ltrs:
			try: d_enrich_ltr_bedlines[ltr.sg] += [ltr.to_bed()]
			except KeyError: d_enrich_ltr_bedlines[ltr.sg] = [ltr.to_bed()]
		#enrich_ltr_bedlines = [ltr.to_bed() for ltr in enrich_ltrs]
		enrich_ltr_bedlines = [v for k,v in sorted(d_enrich_ltr_bedlines.items()) if v]
		return ltr_bedlines, enrich_ltr_bedlines
	def identify_shared_ltrs(self, d_sg, d_chromfiles, ltrfile, d_enriched):
		d_chromfiles_sg = {}
		for chrom, sg in d_sg.items():
			chromfile = d_chromfiles[chrom]
			try: d_chromfiles_sg[sg] += [chromfile]
			except KeyError: d_chromfiles_sg[sg] = [chromfile]
		sgs = sorted(d_chromfiles_sg.keys())
		dumpfiles = []
		for sg, chromfiles in sorted(d_chromfiles_sg.items()):
			prefix = '{}{}.{}'.format(self.tmpdir, self.basename, sg)
			dumpfile = Jellyfish.run_jellyfish_dump(chromfiles, threads=self.ncpu, k=self.k, 
				prefix=prefix, lower_count=self.lower_count*2, overwrite=self.overwrite)
			dumpfiles += [dumpfile]
		chunksize = None if self.pool_method == 'map' else 10000
		dumps = Jellyfish.JellyfishDumps(dumpfiles, sgs, ncpu=self.ncpu,
                            method=self.pool_method, chunksize=chunksize)
		d_mat = dumps.to_matrix(array=False)
		
		# count kmer freq
		logger.info('Identifying shared LTR-RTs')
		pool_method = 'map' #self.pool_method
		chunksize = 2000 #None if pool_method == 'map' else 2000
		out = self.tmpdir + 'LTR.' + self.basename + '.shared'
		with open(out, 'w') as fout:
			d_shared = Seqs.count_kmer([ltrfile], d_mat, lengths=dumps.lengths, k=self.k, ncpu=self.ncpu,
                                exclude=d_enriched, min_prob=0.75, min_count=10, max_fold=1.05, fout=fout,
                                method=pool_method, chunksize=chunksize)
		logger.info('{} shared LTR-RTs'.format(len(d_shared)))
		return d_shared
		
	def step_circos(self, *args, **kargs):
		logger.info('###Step: Circos')
		# blocks
		if not self.disable_blocks:
			pafs, paf_offsets = self.step_blocks()
			kargs['pafs'] = pafs
			kargs['paf_offsets'] = paf_offsets
			kargs['min_block'] = self.min_block

		# circos
		circos_dir = bindir+'/circos'
		wkdir = self.para_prefix + '.circos' #self._outdir+'/circos'
		rmdirs(wkdir)
		try:
			logger.info('Copy `{}` to `{}`'.format(circos_dir, self._outdir))
			shutil.copytree(circos_dir, wkdir)
		except FileExistsError:
			pass
		Circos.circos_plot(self.chromfiles, wkdir, *args, **kargs)
	
	def step_blocks(self):
		outdir = '{}Blocks/'.format(self.tmpdir)
		multiple = {'minimap2': 20, 'unimap': 40}	# relative with repeat amount
		max_size = max(self.d_size.values())
		mem_per_cmd = max_size*math.log(max_size, 10) * multiple[self.aligner]
		ncpu = min(self.ncpu, limit_memory(mem_per_cmd, self.max_memory))
		logger.info('Using {} processes to align chromosome sequences'.format(ncpu))
		thread = int(self.ncpu // ncpu)
		
		mkdirs(outdir)
		pafs, paf_offsets = Blocks.run_align(self.alt_sgs, self.d_chromfiles, outdir, aligner=self.aligner,
						ncpu=ncpu, thread=thread, d_size=self.d_size, overlap=self.min_block*5,
						opts=self.aligner_options, overwrite=self.overwrite)
		#print(paf_offsets)
		return pafs, paf_offsets
	
	def step_final(self):
		# cleanup
		if self.cleanup:
			logger.info('Cleaning {}'.format(self.tmpdir))
			rmdirs(self.tmpdir)
	def sort_labels(self, order, labels, chromfiles):
		d = dict(zip(labels, chromfiles))
		labels, chromfiles = [], []
		for lab in order:
			if not lab in d:
				continue
			chromfile = d[lab]
			labels += [lab]
			chromfiles += [chromfile]
		return labels, chromfiles

def parse_idmap(mapfile=None):
	'''idmap: old_id new_id'''
	if not mapfile:
		return None
	d_map = OrderedDict()
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
			chrs = [list(map(lambda x: add_prefix(x, **self.kargs), v.strip(',').split(','))) 
						for v in temp]
			self.nsgs += [len(chrs)]
			if self.nsg == 0:
				self.nsg = len(chrs)
			if len(chrs) != self.nsg:
				logger.warn('Number of column is different in line {}: \
{} in this line but {} in previous line'.format(temp, len(chrs), self.nsg))
			for xchr in chrs:
				for xxchr in xchr:
					self.chrs += [xxchr]
			yield chrs
		self.nsg = max(self.nsgs)
		for chr, count in Counter(self.chrs).items():
			if count > 1:
				logger.warn('Chromsome id {} repeat {} times'.format(chr, count))

def add_prefix(val, prefix=None, sep='|'):
	if prefix:
		vals = ['{}{}'.format(prefix, v) for v in val.split(sep) if v]
		return ''.join(vals)
	else:
		return val

def main():
	args = makeArgparse()
	logger.info('Command: {}'.format(' '.join(sys.argv)))
	logger.info('Version: {}'.format(version))
	logger.info('Arguments: {}'.format(args.__dict__))
	pipeline = Pipeline(**args.__dict__)
	pipeline.run()


if __name__ == '__main__':
	main()

