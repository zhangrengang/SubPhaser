import sys
import copy
from Bio import SeqIO
from Bio.Seq import Seq
from xopen import xopen as open
from . import Kmer
from .Kmer import get_subgenomes
from .RunCmdsMP import logger, pp_func, pool_func

def split_genomes(genomes, prefixes, targets, outdir, d_targets=None, sep='|'):
	# allow renaming id split by `|`
	d_targets2 = {}
	if not d_targets:
		d_targets = {}
		for t in targets:
			temp = t.split(sep, 1)
			id, new_id = temp[-1], temp[0]
			d_targets[id] = new_id
			d_targets2[t] = new_id
	elif set(targets) - set(d_targets):
		for t in set(targets) - set(d_targets):
			temp = t.split(sep, 1)
			id, new_id = temp[-1], temp[0]
			d_targets[id] = new_id
			d_targets2[t] = new_id
	else:
		d_targets2 = copy.deepcopy(d_targets)

	outfas, labels = [], []
	d_size = {}
	for genome, prefix in zip(genomes, prefixes):
		for rc in SeqIO.parse(open(genome), 'fasta'):
			rc.id = '{}{}'.format(prefix, rc.id)
			if d_targets and rc.id not in d_targets:
				continue
			rc.id = d_targets[rc.id]
			outfa = '{}{}.fasta'.format(outdir, rc.id)
			with open(outfa, 'w') as fout:
				SeqIO.write(rc, fout, 'fasta')
			outfas += [outfa]
			labels += [rc.id]
			d_size[rc.id] = len(rc.seq)
	return outfas, labels, d_targets2, d_size

def map_kmer_each(args):
	rc, d_kmers, k = args
	seq = list(rc.seq.upper())
	lines = []
	for s,e,sg in get_subgenomes(d_kmers, seq, k):
		line = [rc.id, s, e, sg]
		line = list(map(str, line))
		lines += [line]
	return rc.id, lines
def map_kmer_each2(rc, d_kmers, k, ncpu=4):
	seq = list(rc.seq.upper())
#	lines = []
#	iterable = ((seq, i, k) for i in range(len(seq)))
#	jobs = pool_func(Kmer.get_2kmer, iterable, processors=ncpu)
#	for kmers in jobs:
#		for kmer in kmers:
#			try: sg = d_kmers[kmer]
#			except KeyError: continue
#			yield rc.id, s, e, sg
	iterable = ((seq, i, k, d_kmers) for i in range(len(seq)))
	jobs = pool_func(_get_2kmer, iterable, processors=ncpu)
	for kmers in jobs:
		for kmer in kmers:
			if kmer is None:
				continue
			yield (rc.id,) + kmer
def _get_2kmer(args):
	seq, i, k, d_kmers = args
	kmers = Kmer.get_2kmer(seq, i, k)
	for kmer in kmers:
		try: sg = d_kmers[kmer]
		except KeyError: yield None
		yield s, e, sg
def map_kmer(chromfiles, d_kmers, fout=sys.stdout, k=None, ncpu='autodetect', method='map'):
	if k is None:
		k = len(list(d_kmers.keys())[0])
	iterable = ((rc, d_kmers, k) for chromfile in chromfiles for rc in SeqIO.parse(open(chromfile), 'fasta'))
#	logger.info('here?')
#	jobs = pool_func(map_kmer_each, iterable, processors=ncpu, method=method)
#	logger.info('or here?')
	xlines = []	
#	for job in jobs:
#		id, lines = job #()
	j = 0
	for (rc, d_kmers, k) in iterable:
		line = map_kmer_each2(rc, d_kmers, k, ncpu=ncpu)
		line = list(map(str, line))
#		j = 0
#		for line in lines:
		print('\t'.join(line), file=fout)
#			j += 1
#		xlines += lines		# ImportError: cannot import name 'map_kmer_each'
		xlines += [line]
		logger.info('Mapped {} kmers to {}'.format(j, id))
	return xlines


