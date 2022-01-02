import sys, re, os
import math
import itertools
import random
from collections import OrderedDict
import numpy as np
from Bio import SeqIO
try:	# TEsorter need to be updated
	from TEsorter.app import CommonClassifications
except ImportError:
	from .api.TEsorter.app import CommonClassifications
try:
	from TEsorter.modules.concatenate_domains import concat_domains
except :
	from .api.TEsorter.modules.concatenate_domains import concat_domains
	
from .split_records import bin_split_fastx_by_chunk_num
from .RunCmdsMP import run_job, run_cmd, pool_func, logger
from .small_tools import mkdirs, rmdirs, mk_ckp, check_ckp
from .colors import colors_r
from .fonts import fonts, fonts_r

job_args = {
	'tc_tasks': 40,
	'mode': 'grid',
	'grid_opts': '-tc {tc} -l h_vmem={mem} -pe mpi {cpu}',
	'retry': 3,
	'cont': 0,
	'cpu': 1,
	'mem': '2g',
	'template': 'if [ $SGE_TASK_ID -eq {id} ]; then\n{cmd}\nfi',
	}
templates = {
	'ltr_finder': 'ltr_finder {options} {input} > {output}',
	'ltr_harvest': 'gt suffixerator -db {input} -indexname {input} -tis -suf -lcp -des -ssp -sds -dna \
&& gt ltrharvest -index {input} {options} > {output}'
	}
options = {
	'ltr_finder': '-w 2 -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.85',
	'ltr_harvest': '-similar 85 -vic 10 -seed 20 -seqids yes -minlenltr 100 \
-maxlenltr 7000 -mintsd 4 -maxtsd 6 ', # -motif TGCA -motifmis 1
}

def detect_ltr_pp(chunk_files, prefix, 
		progs=['ltr_finder', 'ltr_harvest'], job_args=job_args, options=options ):
	# pp runing
	cmds = []
	d_outputs = {}
	for prog in progs:
		outputs = []
		for chunk_file in chunk_files:
			output = '{}.{}.scn'.format(chunk_file, prog)
			opts = options[prog]
			cmd = templates[prog].format(input=chunk_file, output=output, options=opts)
			cmds.append(cmd)
			outputs.append(output)
		d_outputs[prog] = outputs

	cmd_file = '{}ltr_denovo.sh'.format(prefix)
	run_job(cmd_file, cmds, **job_args)
	return d_outputs

def chunk_fastas(inSeqs, d_len, per_bin=20e6, tmpdir='/io/tmp/share',
			unique=True, window_size=10e6, window_ovl=1e5):
	if unique:	# unique dir
		tmpdir = '{}/ltr_{}'.format(tmpdir, os.getpid())
	mkdirs(tmpdir)
	# binning
	seq_len = sum(d_len.values())
	nbins = int(seq_len/per_bin + 1)
		# prefix = '{}/{}'.format(tmpdir, os.path.basename(inSeq))
	prefix = '{}/chunks'.format(tmpdir, )
	*_, chunk_files = bin_split_fastx_by_chunk_num(
			inSeqs, prefix=prefix, chunk_num=nbins, window_size=window_size, window_ovl=window_ovl, 
			seqfmt='fasta', suffix='')
		# xchunk_files += chunk_files
		# d_chunks[inSeq] = chunk_files
	return chunk_files

def detect_ltr(inSeqs, harvest_out, d_len, progs=['ltr_finder', 'ltr_harvest'],  
		job_args=job_args, options=options, 
		tmpdir='/io/tmp/share', per_bin=20e6, 
		unique=True, window_size=10e6, window_ovl=1e5, **kargs):
	chunk_files = chunk_fastas(inSeqs, d_len, per_bin=per_bin, 
			tmpdir=tmpdir, unique=unique, window_size=window_size, window_ovl=window_ovl)
	prefix = tmpdir +'/'
	d_outputs = detect_ltr_pp(chunk_files, prefix, progs=progs, job_args=job_args, options=options)

	d_idmap = {raw_id: i for i, raw_id in enumerate(d_len.keys())}
	lines = []
	parsers = {
		'ltr_finder': LTRFinder,
		'ltr_harvest': LTRHarvest,
	}
	d_ltrs = OrderedDict((raw_id, []) for raw_id in d_len.keys())
	for prog in progs:
		for input, output in zip(chunk_files, d_outputs[prog]):
			args = {'ltr_finder': {'finder_out': output}, 
					'ltr_harvest': {'harvest_out': output, 'harvest_in': input}
					}
			for rc in parsers[prog](**args[prog]):
				raw_id, start, end = re.compile(r'^(\S+):(\d+)\-(\d+)$').match(rc.seq_id).groups()
				start = int(start)
				rc.seq_id = raw_id
				rc.start, rc.end = rc.start+start-1, rc.end+start-1
				rc.seq_nr = d_idmap[raw_id]
				#	print(inSeq, prefix,chunk_files, d_outputs[prog], e, d_idmap)
				rc.source = [prog]
				d_ltrs[raw_id] += [rc]

	# merge, remove overlaps and output
	print('''# LTR_pp
# Note: overlap between two LTRs is resolved by removing the partial one and the shorter one
#start end len lLTR_str lLTR_end lLTR_len rLTR_str rLTR_end rLTR_len similarity seqid chr direction TSD lTSD rTSD motif superfamily family age(ya)
# LTR_FINDER args=$finder_para
# predictions are reported in the following way
# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr
# where:
# s = starting position
# e = ending position
# l = length
# ret = LTR-retrotransposon
# lLTR = left LTR
# rLTR = right LTR
# sim = similarity
# seq-nr = sequence order''', file=harvest_out)

	all_ltrs = []
	for chrom in d_len.keys():
		ltrs = d_ltrs[chrom]
		if len(progs) > 1:
			ltrs = resolve_overlaps(ltrs, max_ovl=95)	#relaxed
		all_ltrs += ltrs
		for rc in ltrs:
			line = tuple(rc.harvest_output())
			source = rc.source
			source = '#' + ','.join(source)
			line = list(map(str, line)) + [source]
			print(' '.join(line), file=harvest_out)

	#rmdirs(tmpdir)
	return all_ltrs

class LTRtree:
	tree_template = {
		'iqtree': 'iqtree -s {alnfile}.trimal -pre {alnfile} {tree_options} -nt {ncpu} -redo',
		'FastTree': 'FastTree {tree_options} {alnfile}.trimal  > {alnfile}.treefile'}
	tree_template['fasttree'] = tree_template['FastTree']
	def __init__(self, ltrs, domains, domfile, prefix='ltrtree', 
			overwrite=False, ncpu=10, subsample=None, 
			categories=[('LTR', 'Copia', None), ('LTR', 'Gypsy', None)],
			trimal_options='-automated1', 
			tree_method='iqtree', tree_options='-mset JTT -nt AUTO'):
		'''ltrs: enrich_ltrs (list)'''
		self.ltrs = ltrs
		self.domains = domains
		self.domfile = domfile
		self.prefix = prefix
		self.categories = categories
		self.trimal_options = trimal_options
		self.tree_method = tree_method
		self.tree_options = tree_options
		self.overwrite = overwrite
		self.ncpu = ncpu
		self.subsample = subsample
	
	def build(self, job_args):
		d_ltrs = {ltr.id: ltr for ltr in self.ltrs}
#		ltr_ids = {ltr.id for ltr in ltrs}
		subsample = self.subsample
		nseqs = []
		alnfiles = []
		for key in self.categories:
			order, superfamily, clade = key
			key = tuple([v for v in key if v])
			logger.info('Extracting and aligning protein domain sequences of {}'.format('/'.join(key)))
			prefix = '{}.{}'.format(self.prefix, '_'.join(key))
			alnfile = prefix + '.aln'
			mapfile = prefix + '.map'
			ckp_file = alnfile + '.ok'
			ckp = check_ckp(ckp_file)
			if ckp and isinstance(ckp, list) and not self.overwrite:
				nseq, *_ = ckp
			else:
				with open(alnfile, 'w') as fout:
					_, d_idmap = concat_domains(self.domfile, self.domains, outSeq=fout, 
						targets=d_ltrs, unique=True, prefix=prefix, raw=True, format_id=True,
						order=order, superfamily=superfamily, clade=clade, subsample=subsample)
				with open(mapfile, 'w') as f:
					line = ['label', 'Clade', 'Subgenome']
					f.write('\t'.join(line)+'\n')
					for raw_id, new_id in d_idmap.items():
						ltr = d_ltrs[raw_id]
						clade, sg = ltr.clade, ltr.sg
						line = [new_id, clade, sg]
						f.write('\t'.join(line)+'\n')
				nseq = len(d_idmap)
				mk_ckp(ckp_file, nseq)
			nseqs += [nseq]
			alnfiles += [(key, alnfile, mapfile, nseq)]
			
		# assign CPU
		prop = [v**2 for v in nseqs]
		tprop = sum(prop)
		ncpus = [max(1, int(self.ncpu*v/tprop)) for v in prop]
		#print(ncpus)
		
		# tree
		d_files = {}
		cmds = []
		for ncpu, (key, alnfile, mapfile, nseq) in zip(ncpus, alnfiles):
			if nseq < 4:
				continue
			treefile = alnfile + '.rooted.tre'
			template = 'trimal -in {alnfile} {trimal_options} > {alnfile}.trimal && ' + \
					self.tree_template[self.tree_method] + \
					'&& nw_reroot {alnfile}.treefile > {treefile}'
			cmd = template.format(
				alnfile=alnfile, treefile=treefile, 
				trimal_options=self.trimal_options,
				tree_options=self.tree_options, ncpu=ncpu)
			
			cmds +=[cmd]
			d_files[key] = (treefile, mapfile)
			
		cmd_file = '{}.LTRtree.sh'.format(self.prefix)
		run_job(cmd_file, cmds, **job_args)
		return d_files

	def visualize_treefile(self, treefile, mapfile, outfig, 
			ggtree_options="branch.length='none', layout='circular'"):
		rsrc_file = os.path.splitext(outfig)[0] + '.R'
		rsrc = '''treefile = "{treefile}"
mapfile = "{mapfile}"
branch_color = 'Subgenome'
library(ggplot2)
library(ggtree)
library(treeio)

map = read.table(mapfile, head=T, fill=T)
tree <- read.tree(file = treefile)

if (branch_color == 'Clade') {{
	clades = sort(unique(map$Clade))

	grp = list()
	for (clade in clades){{
			labels = map[which(map$Clade==clade), ]
			labels = labels$label
			grp[[clade]] = labels
	}}
	tree3 = groupOTU(tree, grp, 'Clade')
	map3 = data.frame(label=map$label, Subgenome=map$Subgenome)
	p = ggtree(tree3 , aes(color=Clade) , {ggtree_options} ) %<+% map3 +
	  theme(legend.position="right")  + 
	  geom_tippoint(aes(fill=Subgenome), pch=21, stroke=0, size=1.1, color='#00000000') +
	  scale_fill_manual(values={colors}) + scale_colour_discrete(limits=clades, labels=clades) +
	  scale_fill_hue(l=35) +
	  guides(colour=guide_legend(order = 1), fill=guide_legend(order = 2))

}} else {{	# branch_color == 'Subgenome'
	subgenomes = sort(unique(map$Subgenome))

	grp = list()
	for (subgenome in subgenomes){{
			labels = map[which(map$Subgenome==subgenome), ]
			labels = labels$label
			grp[[subgenome]] = labels
	}}
	tree3 = groupOTU(tree, grp, 'Subgenome')
	map3 = data.frame(label=map$label, Clade=map$Clade)
	p = ggtree(tree3 , aes(color=Subgenome) , {ggtree_options} ) %<+% map3 +
	  theme(legend.position="right")  + 
	  scale_colour_manual(values={colors},limits=subgenomes, labels=subgenomes) +
	  geom_tippoint(aes(fill=Clade), pch=21, stroke=0, size=1.2, color='#00000000') +
	  scale_fill_hue(l=35) +
	  guides(colour=guide_legend(order = 1), fill=guide_legend(order = 2))

}}
p = p + theme(plot.margin=margin(0,0,0,0)) +
	theme(legend.position=c(1.145,0.9), legend.justification=c(1.145,0.9)) +
	theme(legend.background=element_blank(), legend.key=element_blank()) +
	theme(legend.text=element_text(size={fontsize}), legend.title=element_text(size={title_fontsize}))

#	theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())

ggsave("{outfig}", p, width=10.2, height=8.4, dpi=300, units="in")
'''.format(treefile=treefile, mapfile=mapfile, outfig=outfig, colors=colors_r, 
			ggtree_options=ggtree_options, **fonts_r)
		with open(rsrc_file, 'w') as f:
			f.write(rsrc)
		cmd = 'Rscript ' + rsrc_file
		run_cmd(cmd, log=True)



class LTRpipeline:
	def __init__(self, genomes, tmpdir='./tmp', mu=7e-9, tesorter_options='', 
			all_ltr=False, intact=False, only_ltr=True, overwrite=False, **kargs):
		'''all_ltr: use all LTR identified by LTR detectors
only_ltr: use LTR as classified by TEsorter
intact: only use completed LTR as classified by TEsorter'''
		self.genomes = genomes
		self.tmpdir = tmpdir
		self.tesorter_options = tesorter_options
		self.all_ltr = all_ltr
		self.intact = False if only_ltr else intact 
		self.only_ltr = only_ltr	# only LTR classified by tesorter
		self.mu = mu
		self.overwrite = overwrite
		self.kargs = kargs
	def run(self):
		mkdirs(self.tmpdir)
		self.d_seqs = OrderedDict((rc.id, rc.seq) for genome in self.genomes \
									 for rc in SeqIO.parse(genome, 'fasta'))
		self.d_len = OrderedDict((id, len(seq)) for id, seq in self.d_seqs.items())
		self.prefix = '{}'.format(self.tmpdir, )
		ltr_out = self.prefix +'.scn'
		ckp = ltr_out + '.ok'
		if check_ckp(ckp) and not self.overwrite:
			ltrs = list(LTRHarvest(ltr_out))
		else:
			with open(ltr_out, 'w') as fout:
				ltrs = self.identify(fout)
			mk_ckp(ckp)
		ltr_count = len(ltrs)
		logger.info('{} LTRs identified'.format(ltr_count))
		
		d_class = self.classfify(ltrs)

			
		filtered_ltrs = []
		i, j = 0,0
		for ltr in ltrs:
			cls = d_class.get(ltr.id, )
			if cls:
				ltr.__dict__.update(**cls.__dict__)
			else:
				ltr.__dict__.update(order=None, superfamily=None,
					clade=None, completed=None, strand=None, domains=None)
			order = getattr(cls, 'order', None)
			completed = getattr(cls, 'completed', None)
			if order == 'LTR':
				i += 1
			if completed == 'yes':
				j += 1
			if self.all_ltr:	# no filter
				pass
			elif self.only_ltr and order != 'LTR':
				continue
			elif self.intact and completed != 'yes':
				continue
			filtered_ltrs += [ltr]
		logger.info('By TEsorter, {} ({:.1%}) are classified as LTRs, of which {} ({:.1%}) are intact \
with complete protein domains'.format(i, i/ltr_count, j, j/i))
		
		i, j = len(ltrs), len(filtered_ltrs)
	#	logger.info('After filter, {} / {} ({:.1f}%) LTRs retained'.format(j, i, 1e2*j/i))
		
		ltrs = group_resolve_overlaps(filtered_ltrs)
		j = len(ltrs)
		logger.info('After filtering, {} / {} ({:.1%}) LTRs retained'.format(j, i, j/i))
		ltr_seqs = self.prefix + '.filtered.LTR.fa' 
		with open(ltr_seqs, 'w') as fout:
			self.get_full_seqs(ltrs, fout)	# for kmer enrich
		return ltrs, ltr_seqs


	def identify(self, fout):
		return detect_ltr(self.genomes, fout, self.d_len, tmpdir=self.prefix, unique=False, **self.kargs)

	def classfify(self, ltrs):
		inseq = self.int_seqs = int_seqs = self.prefix + '.inner.fa'
		ckp = self.prefix + '.tesort.ok'
		if check_ckp(ckp) and not self.overwrite:
			pass
		else:
			logger.info('Extracting inner sequences of LTRs to classify by `TEsorter`')
			with open(int_seqs, 'w') as fout:
				self.get_int_seqs(ltrs, fout)
				
			cmd = 'TEsorter {seqfile} {options} -pre {seqfile} -tmp {tmpdir} &> \
{seqfile}.tesort.log'.format(
					seqfile=inseq, options=self.tesorter_options, tmpdir=self.prefix)
			if self.intact and '-dp2' not in cmd and "--disable-pass2" not in cmd:
				cmd += ' -dp2'
			run_cmd(cmd, log=True)
			mk_ckp(ckp)
		clsfile = '{seqfile}.cls.tsv'.format(seqfile=inseq)
		d_class = {}
		for classification in CommonClassifications(clsfile):
			d_class[classification.id] = classification
		return d_class

	def get_seqs(self, ltrs, fout, method='get_full_seq'):
		d_seqs = self.d_seqs
		for ltr in ltrs:
			seq = getattr(ltr, method)(d_seqs=d_seqs)
			if seq is None:
				continue
			fout.write('>{}\n{}\n'.format(ltr.id, seq))
	def get_int_seqs(self, ltrs, fout):
		return self.get_seqs(ltrs, fout, method='get_int_seq')
	def get_full_seqs(self, ltrs, fout):
		return self.get_seqs(ltrs, fout, method='get_full_seq')

def group_resolve_overlaps(ltrs):
	'''assume multiple chromsomes'''
	resolved_ltrs = []
	for chrom, items in itertools.groupby(ltrs, key=lambda x:x.seq_id):
		resolved_ltrs += resolve_overlaps(list(items))
	return resolved_ltrs

def resolve_overlaps(ltrs, max_ovl=10):
	'''assume only one chromsome'''
	last_ltr = None
	discards = []
	ie, io = 0, 0
	for ltr in sorted(ltrs, key=lambda x:x.start):
		discard = None
	#	print(last_ltr, ltr)
		if last_ltr:
			both_completed = is_completed(ltr) and is_completed(last_ltr)
			both_uncompleted = not (is_completed(ltr) or is_completed(last_ltr))
			 
			if ltr == last_ltr:	# equal
				ie += 1
				ltr_pair = [last_ltr, ltr]	# retain, discard
			elif both_completed or both_uncompleted:
				if ltr.overlap(last_ltr) > max_ovl: # overlaps
					io += 1
					if ltr.element_len > last_ltr.element_len:
						ltr_pair = [ltr, last_ltr]
					else:
						ltr_pair = [last_ltr, ltr]
				else:	# no overlap or too short overlap
					last_ltr = ltr
					continue
			else:
				if ltr.overlap(last_ltr) > max_ovl:
					io += 1
					if is_completed(ltr):	# completed is prior
						ltr_pair = [ltr, last_ltr]
					else: # is_completed(last_ltr)
						ltr_pair = [last_ltr, ltr]
				else:	# no overlap or too short overlap
					last_ltr = ltr
					continue
			
			retain, discard = ltr_pair
			try:
				retain.source += discard.source
			except AttributeError:
				pass
			discards += [discard]

		if not last_ltr or discard != ltr:
			last_ltr = ltr
#	logger.info('Discard {} equal and {} overlapped LTRs; {} in total'.format(ie, io, ie+io))
	return sorted(set(ltrs) - set(discards), key=lambda x:x.start)

def is_completed(ltr):
	completed = getattr(ltr, 'completed', None)
	return True if completed == 'yes' else False

def plot_insert_age(ltrs, d_enriched, prefix, mu=7e-9, shared={}, figfmt='pdf'):
	datfile = prefix + '.data'
	fout = open(datfile, 'w')
	line = ['ltr', 'sg', 'age']
	fout.write('\t'.join(line) + '\n')
	d_data = {}
	enriched_ltrs = []
	for ltr in ltrs:
		age = ltr.estimate_age(mu=mu)
		if ltr.id in d_enriched:
			sg = d_enriched[ltr.id]
			enriched_ltrs += [ltr]
		elif ltr.id in shared:
			sg = 'shared'
		else:
			continue
#		try: sg = d_enriched[ltr.id]
#		except KeyError: continue
		ltr.sg = sg
		age = age/1e6
#		enriched_ltrs += [ltr]
		line = [ltr.id, sg, age]
		line = map(str, line)
		fout.write('\t'.join(line) + '\n')
		try: d_data[sg] += [age]
		except KeyError: d_data[sg] = [age]
	fout.close()
	# summary
	sumfile = prefix + '.summary'
	with open(sumfile, 'w') as fout:
		d_info = summary_ltr_time(d_data, fout)
	text = 'Summary: median (95% CI)\\n'
	for sg, info in sorted(d_info.items()):
		text += '{}: {}\\n'.format(sg, info)
		
	rsrc_file = prefix + '.R'
	outfig = prefix + '.' + figfmt
	rsrc = '''library(ggplot2)
data = read.table('{datfile}',fill=T,header=T, sep='\\t')
p <- ggplot(data, aes(x = age, color=sg)) + geom_line(stat="density", size=1.5) + 
	xlab('LTR insertion age (million years)') + ylab('Density') + 
	scale_colour_manual(values={colors}) + labs(color='Subgenome') + 
	annotate('text',x=Inf, y=Inf, label="{annotate}", hjust=1.1, vjust=1.1)+
	theme_bw() + 
	theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
	theme(legend.position=c(0.95,0.8), legend.justification=c(0.95,0.8)) +
	theme(legend.background=element_blank(), legend.key=element_blank()) +
	theme(legend.text=element_text(size={fontsize}), axis.title=element_text(size={fontsize}),
		axis.text=element_text(size={tick_fontsize})) +
	guides(color=guide_legend(title=NULL))

ggsave('{outfig}', p, width=7, height=7, dpi=300, units="in") 
'''.format(datfile=datfile, outfig=outfig, colors=colors_r, annotate=text, **fonts_r)
	with open(rsrc_file, 'w') as f:
		f.write(rsrc)
	cmd = 'Rscript ' + rsrc_file
	run_cmd(cmd, log=True)
	
	return enriched_ltrs

def summary_ltr_time(d_data, fout):
	fout.write('# Summary of LTR insertion age (million years)\n')
	line = ['#subgenome', 'mean', 'median', 'standard_deviation', '95%-CI', '75%-CI']
	fout.write('\t'.join(line) + '\n')
	d_info = {}
	xages = []
	_medians = []
	_tile2_5s = []
	_tile97_5s = []
	for sg, ages in sorted(d_data.items()):
		xages += ages
		ages = np.array(ages)
		_mean = ages.mean()
		_median = np.median(ages)
		_medians += [_median]
		_std = np.std(ages)
		_mean = '{:.3f}'.format(_mean)
		_median = '{:.3f}'.format(_median)
		_std = '{:.3f}'.format(_std)
		_tile2_5 = np.percentile(ages, 2.5)
		_tile2_5s += [_tile2_5]
		_tile97_5 = np.percentile(ages, 97.5)
		_tile97_5s += [_tile97_5]
		_tile12_5 = np.percentile(ages, 12.5)
		_tile87_5 = np.percentile(ages, 87.5)
		_ci95 = '{:.3f}-{:.3f}'.format(abs(_tile2_5), _tile97_5)
		_ci75 = '{:.3f}-{:.3f}'.format(_tile12_5, _tile87_5)
		line = [sg, _mean, _median, _std, _ci95, _ci75]
		fout.write('\t'.join(line) + '\n')
		d_info[sg] = '{} ({})'.format(_median, _ci95)
	logger.info('Summary of overall LTR insertion age (million years):')
	logger.info('\tmedian: {:.3f}\t95% CI (percentile-based): {:.3f}-{:.3f}'.format(
		np.median(xages), abs(np.percentile(xages, 2.5)), np.percentile(xages, 97.5)))
	logger.info('A rough estimation of the divergence–hybridization period: {:.3f}-{:.3f} ({:.3f})'.format(
		np.mean(_tile2_5s), np.mean(_tile97_5s), np.mean(_medians)))
	return d_info


class LTRHarvest():
	def __init__(self, harvest_out, harvest_in=None):
		self.harvest_out = harvest_out
		self.harvest_in = harvest_in
		self.idmap = self._parse_idmap()
	def __iter__(self):
		return self.parse()
	def _parse_idmap(self):
		des = str(self.harvest_in) + '.des'
		d = {}
		if os.path.exists(des):
			for i, line in enumerate(open(des, 'rb')):
				try: line = line.decode('UTF-8')
				except: continue
				d[i] = line.split()[0]
		elif not self.harvest_in or not os.path.exists(self.harvest_in):
			pass
		else:
			for i, rc in enumerate(SeqIO.parse(self.harvest_in, 'fasta')):
				d[i] = rc.id
		return d

	def parse(self):
		for line in open(self.harvest_out):
			if not line.startswith('#'):
				yield LTRHarvestRecord(line, self.idmap)
class LTRHarvestRecord(object):
	def __init__(self, line, idmap=None):
		self.line = line.strip().split()
		self.start, self.end, self.element_len, self.start, \
			self.lltr_e0, self.lltr, self.rltr_s0, self.end, self.rltr, \
			self.similarity, self.seq_nr = self.line[:11]
		self.seq_nr = int(self.seq_nr)
		if len(self.line) == 11: # raw format
			self.seq_id = idmap[self.seq_nr]
		elif len(self.line) >= 12:	# modified format output by pipeline
			self.seq_id = self.line[11]
		else:
			raise ValueError('Unrecognized LTRHarvest format: {}'.format(self.line))
		self.start, self.end = int(self.start), int(self.end)
		self.lltr, self.rltr =  int(self.lltr), int(self.rltr)
		self.similarity = float(self.similarity)
		self.element_len = int(self.element_len)
		self.lltr_e0, self.rltr_s0 = int(self.lltr_e0), int(self.rltr_s0)
	def harvest_output(self, fout=None):
		line = [self.start, self.end, self.element_len, self.start,
				self.lltr_e, self.lltr, self.rltr_s, self.end, self.rltr,
				self.similarity, self.seq_nr, self.seq_id]
		if fout is not None:
			line = list(map(str, line))
			print(' '.join(line), file=fout)
		else:
			return line
	@property
	def key(self):
		return (self.seq_id, self.start, self.end, self.lltr_e, self.rltr_s)
	@property
	def id(self):
		return '{seq_id}:{start}-{end}:{lltr_e}-{rltr_s}'.format(
			lltr_e=self.lltr_e, rltr_s=self.rltr_s, **self.__dict__)
	def to_bed(self):
		return [self.seq_id, self.start, self.end, self.id]
	def __str__(self):
		return self.id
	def __hash__(self):
		return hash(self.key)
	def __eq__(self, other):
		return self.key == other.key
	def overlap(self, other):
		ovl = max(0, min(self.end, other.end) - max(self.start, other.start))
		return 100*ovl/(min(self.element_len, other.element_len))
	def estimate_age(self, mu=7e-9, method='JC69'):
		div = 1- self.similarity/100
		if div >= 0.75:
			dist = div
		else:
			dist = -3/4* math.log(1-4*div/3)
		return dist / (mu *2)
	@property
	def lltr_e(self):	# compitable with ltrfinder
		return self.start + self.lltr - 1
	@property
	def rltr_s(self):
		return self.end - self.rltr + 1
	def get_seq(self, seq=None, start=None, end=None, d_seqs=None):
		if seq is None and d_seqs is not None:
			try: seq = d_seqs[self.seq_id]
			except KeyError: return
		return seq[start:end]
	def get_int_seq(self, **kargs):
		return self.get_seq(start=self.lltr_e, end=self.rltr_s, **kargs)
	def get_full_seq(self, **kargs):
		return self.get_seq(start=self.start, end=self.end, **kargs)

class LTRFinder():
	def __init__(self, finder_out):
		self.finder_out = finder_out
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.finder_out):
			if line.startswith('['):
				yield LTRFinderRecord(line)
class LTRFinderRecord(LTRHarvestRecord):
	def __init__(self, line):
		self.line = line.strip().split('\t')
		self.index, self.seq_id, self.location, self.ltr_len, self.element_len, \
			self.TSR, self.PBS, self.PPT, self.RT, self.IN_core, self.IN_term, \
			self.RH, self.strand, self.score, self.sharpness, self.similarity = self.line
		self.similarity = round(100 * float(self.similarity),1 )
		self.start, self.end = self.location.split('-')
		self.start, self.end = int(self.start), int(self.end)
		self.lltr, self.rltr = self.ltr_len.split(',')
		self.lltr, self.rltr =  int(self.lltr), int(self.rltr)
		self.element_len = int(self.element_len)
		self.seq_nr = None

if __name__ == '__main__':
	inSeq = sys.argv[1]
	harvest_out = sys.stdout
	main(inSeq, harvest_out)
