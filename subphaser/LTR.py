import sys, re, os
import math
from collections import OrderedDict
from Bio import SeqIO
from TEsorter.app import CommonClassifications
from .split_records import bin_split_fastx_by_chunk_num
from .RunCmdsMP import run_job, run_cmd, pool_func, logger
from .small_tools import mkdirs, rmdirs, mk_ckp,check_ckp

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

def detect_ltr(inSeq, harvest_out, progs=['ltr_finder', 'ltr_harvest'],  per_bin=20e6, 
			tmpdir='/io/tmp/share', unique=True, window_size=10e6, window_ovl=1e5,
			job_args=job_args, options=options):
	if unique:	# unique dir
		tmpdir = '{}/ltr_{}'.format(tmpdir, os.getpid())
	mkdirs(tmpdir)
	d_len = OrderedDict((rc.id, len(rc.seq)) for rc in SeqIO.parse(inSeq, 'fasta'))
	# binning
	seq_len = sum(d_len.values())
	nbins = int(seq_len/per_bin + 1)
	prefix = '{}/{}'.format(tmpdir, os.path.basename(inSeq))
	*_, chunk_files = bin_split_fastx_by_chunk_num(
		inSeq, prefix=prefix, chunk_num=nbins, window_size=window_size, window_ovl=window_ovl, seqfmt='fasta', suffix='')

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

	cmd_file = '{}.ltr_denovo.sh'.format(prefix)
	run_job(cmd_file, cmds, **job_args)

	print('''# LTR_pp {}
# Note: overlap between two LTRs is resolved by removing the shorter one
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
# seq-nr = sequence order'''.format(inSeq), file=harvest_out)
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
	all_ltrs = []
	for chrom in d_len.keys():
		ltrs = d_ltrs[chrom]
		if len(progs) > 1:
			ltrs = resolve_overlaps(ltrs)
		all_ltrs += ltrs
		for rc in ltrs:
			line = tuple(rc.harvest_output())
			source = rc.source
			source = '#' + ','.join(source)
			line = list(map(str, line)) + [source]
			print(' '.join(line), file=harvest_out)

	#rmdirs(tmpdir)
	return all_ltrs

def resolve_overlaps(ltrs, max_ovl=10):
	last_ltr = None
	discards = []
	ie, io = 0, 0
	for ltr in sorted(ltrs, key=lambda x:x.start):
		discard = None
	#	print(last_ltr, ltr)
		if last_ltr:
			if ltr == last_ltr:	# equal
				ie += 1
				discard = ltr
				last_ltr.source += ltr.source
			elif ltr.overlap(last_ltr) > max_ovl: # overlaps
				io += 1
				if ltr.element_len > last_ltr.element_len:
					discard = last_ltr
					ltr.source += last_ltr.source
				else:
					last_ltr.source += ltr.source
					discard = ltr
			else:	# no overlap or too short overlap
				last_ltr = ltr
				continue
			discards += [discard]
#			if last_ltr.source != ltr.source:
#				print(last_ltr, ltr, ltr.overlap(last_ltr))
#		else:
#			last_ltr = ltr
#		if not discard or discard == last_ltr:
		if not last_ltr or discard != ltr:
			last_ltr = ltr
	logger.info('Discard {} equal and {} overlapped; {} in total'.format(ie, io, ie+io))
	return sorted(set(ltrs) - set(discards), key=lambda x:x.start)

class LTRpipelines:
	def __init__(self, genomes, ncpu=8, tmpdir='tmp', **kargs):
		self.genomes = genomes
		self.ncpu = ncpu
		self.tmpdir = tmpdir
		self.kargs = kargs
	def run(self):
		ltrs = []
		seqs = []
		for _ltrs, _seqs in pool_func(self._run, self.genomes, self.ncpu):
			ltrs += _ltrs
			seqs += [_seqs]
		return ltrs, seqs
	def _run(self, genome):
#		tmpdir = '{}/{}'.format(self.tmpdir, os.path.basename(genome))
		return LTRpipeline(genome, tmpdir=self.tmpdir, **self.kargs).run()

def plot_insert_age(ltrs, d_enriched, prefix, mu=7e-9, figfmt='pdf'):
	datfile = prefix + '.data'
	fout = open(datfile, 'w')
	line = ['ltr', 'sg', 'age']
	fout.write('\t'.join(line) + '\n')
	enriched_ltrs = []
	for ltr in ltrs:
		age = ltr.estimate_age(mu=mu)
		try: sg = d_enriched[ltr.id]
		except KeyError: continue
		enriched_ltrs += [ltr]
		line = [ltr.id, sg, age]
		line = map(str, line)
		fout.write('\t'.join(line) + '\n')
	fout.close()
	rsrc_file = prefix + '.R'
	outfig = prefix + '.' + figfmt
	rsrc = '''library(ggplot2)
data = read.table('{datfile}',fill=T,header=T, sep='\t')
p <- ggplot(data, aes(x = age, color=sg)) + geom_line(stat="density", size=1.15) + \
xlab('LTR insertion age (years)') + scale_colour_hue(l=45)
ggsave('{outfig}', p, width=12, height=7)
'''.format(datfile=datfile, outfig=outfig)
	with open(rsrc_file, 'w') as f:
		f.write(rsrc)
	cmd = 'Rscript ' + rsrc_file
	run_cmd(cmd, log=True)
	return enriched_ltrs

class LTRpipeline:
	def __init__(self, genome, tmpdir='./tmp', mu=7e-9, tesorter_options='', 
			all_ltr=False, intact=False, only_ltr=True, **kargs):
		'''all_ltr: use all LTR identified by LTR detectors
only_ltr: use LTR as classified by TEsorter
intact: only use completed LTR as classified by TEsorter'''
		self.genome = genome
		self.tmpdir = tmpdir
		self.tesorter_options = tesorter_options
		self.all_ltr = all_ltr
		self.intact = False if only_ltr else intact 
		self.only_ltr = only_ltr	# only LTR classified by tesorter
		self.mu = mu
		self.kargs = kargs
	def run(self):
		mkdirs(self.tmpdir)
		self.prefix = '{}/{}'.format(self.tmpdir, os.path.basename(self.genome))
		ltr_out = self.prefix +'.scn'
		ckp = ltr_out + '.ok'
		if check_ckp(ckp):
			ltrs = list(LTRHarvest(ltr_out))
		else:
			with open(ltr_out, 'w') as fout:
				ltrs = self.identify(fout)
			mk_ckp(ckp)
		int_seqs = self.prefix + '.INT.fa'
		with open(int_seqs, 'w') as fout:
			self.get_int_seqs(ltrs, fout)
		d_class = self.classfify(int_seqs)

			
		filtered_ltrs = []
		for ltr in ltrs:
			if self.all_ltr:	# no filter
				ltr.age = ltr.estimate_age(mu=self.mu)
				
				filtered_ltrs += [ltr]
				continue
			if ltr.id not in d_class:	# filter by classification
				continue
			cls = d_class[ltr.id]
			if self.only_ltr and cls.order != 'LTR':
				continue
			elif self.intact and cls.completed != 'yes':
				continue
			ltr.superfamily = cls.superfamily
			ltr.age = ltr.estimate_age(mu=self.mu)
			filtered_ltrs += [ltr]
		ltrs = filtered_ltrs
		ltr_seqs = self.prefix + '.LTR.fa' 
		with open(ltr_seqs, 'w') as fout:
			self.get_full_seqs(ltrs, fout)
		return ltrs, ltr_seqs

	def identify(self, fout):
		return detect_ltr(self.genome, fout, tmpdir=self.prefix, unique=False, **self.kargs)
	def classfify(self, inseq):
		ckp = self.prefix + '.tesort.ok'
		if check_ckp(ckp):
			pass
		else:
			cmd = 'TEsorter {seqfile} {options} -pre {seqfile} -nocln -tmp {tmpdir} > \
{seqfile}.tesort.log'.format(
					seqfile=inseq, options=self.tesorter_options, tmpdir=self.prefix)
			if self.intact:
				cmd += ' -dp2'
			run_cmd(cmd, log=True)
			mk_ckp(ckp)
		clsfile = '{seqfile}.cls.tsv'.format(seqfile=inseq)
		d_class = {}
		for classification in CommonClassifications(clsfile):
			d_class[classification.id] = classification
		return d_class

	def get_seqs(self, ltrs, fout, method='get_full_seq'):
		d_seqs = {rc.id: rc.seq for rc in SeqIO.parse(self.genome, 'fasta')}
		for ltr in ltrs:
			seq = getattr(ltr, method)(d_seqs=d_seqs)
			fout.write('>{}\n{}\n'.format(ltr.id, seq))
	def get_int_seqs(self, ltrs, fout):
		return self.get_seqs(ltrs, fout, method='get_int_seq')
	def get_full_seqs(self, ltrs, fout):
		return self.get_seqs(ltrs, fout, method='get_full_seq')

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
			dist = -3/4* math.log(1-4*div/3, 10)
		return dist / (mu *2)
	@property
	def lltr_e(self):	# compitable with ltrfinder
		return self.start + self.lltr - 1
	@property
	def rltr_s(self):
		return self.end - self.rltr + 1
	def get_seq(self, seq=None, start=None, end=None, d_seqs=None):
		if seq is None and d_seqs is not None:
			seq = d_seqs[self.seq_id]
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
