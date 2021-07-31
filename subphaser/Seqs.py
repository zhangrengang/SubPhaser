import sys
import copy
from Bio import SeqIO
from Bio.Seq import Seq
from xopen import xopen as open
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
	return outfas, labels, d_targets2

def map_kmer_each(args):
	rc, d_kmers, k = args
	seq = list(rc.seq.upper())
	lines = []
	for s,e,sg in get_subgenomes(d_kmers, seq, k):
		line = [rc.id, s, e, sg]
		line = map(str, line)
		lines += [line]
	return rc.id, lines

def map_kmer(chromfiles, d_kmers, fout=sys.stdout, k=None, ncpu='autodetect'):
	if k is None:
		k = len(list(d_kmers.keys())[0])
	iterable = ((rc, d_kmers, k) for chromfile in chromfiles for rc in SeqIO.parse(open(chromfile), 'fasta'))
	jobs = pool_func(map_kmer_each, iterable, processors=ncpu)

	xlines = []	
	for job in jobs:
		id, lines = job #()
		j = 0
		for line in lines:
			print('\t'.join(line), file=fout)
			j += 1
		xlines += lines		# ImportError: cannot import name 'map_kmer_each'
		logger.info('Mapped {} kmers to {}'.format(j, id))
	return xlines


