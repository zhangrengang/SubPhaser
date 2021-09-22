import sys
import copy
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from xopen import xopen as open
from .RunCmdsMP import logger, pp_func, pool_func

def split_genomes(genomes, prefixes, targets, outdir, d_targets=None, sep='|'):
	# allow renaming id split by `|`
	d_targets2 = OrderedDict()
	if not d_targets:
		d_targets = OrderedDict()
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
	got_ids = set([])
	for genome, prefix in zip(genomes, prefixes):
		for rc in SeqIO.parse(open(genome), 'fasta'):
			old_id, new_id = rc.id, '{}{}'.format(prefix, rc.id)
#			rc.id = '{}{}'.format(prefix, rc.id)
			if d_targets:
				if new_id in d_targets:
					rc.id = new_id
				elif old_id in d_targets:
					pass
				else:
					continue
			got_ids.add(rc.id)
			rc.id = d_targets[rc.id]
			outfa = '{}{}.fasta'.format(outdir, rc.id)
			with open(outfa, 'w') as fout:
				SeqIO.write(rc, fout, 'fasta')
			outfas += [outfa]
			labels += [rc.id]	# new id
			d_size[rc.id] = len(rc.seq)
	ungot_ids = set(d_targets) - got_ids
	if ungot_ids:
		logger.error('Chromosomes {} are not found in sequences files'.format(ungot_ids))
	return outfas, labels, d_targets2, d_size


def map_kmer3(chromfiles, d_kmers, fout=sys.stdout, k=None, window_size=10e6, 
		bin_size=10000, sg_names=[], 
		ncpu='autodetect', method='map', log=True, chunk=True, chunksize=None):
	if k is None:
		for key in d_kmers.keys():
			k = len(key)
			break
	overlap = k-1
	if chunk:
		chunks = chunk_chromfiles(chromfiles, window_size=window_size, overlap=overlap)
	else:
		chunks = unchunk_chromfiles(chromfiles)

	iterable = ((id, start, seq, k, d_kmers, bin_size, sg_names) for id, start, seq in chunks)
	last_id = ''
	i, j = 0, 0
	mapped_num = 0	# number of mappped kmers
	mapped_seqs = 0	# number of mappped sequences
	mapped_cat = set([])
	#method='map'
	line = ['#chrom', 'start', 'end'] + sg_names
	fout.write('\t'.join(line)+'\n')
	for id, c, mapped_kmers, lines in pool_func(map_kmer_each4, iterable, 
					processors=ncpu, method=method, chunksize=chunksize):
		if last_id and id != last_id:
			if log:
				logger.info('Mapped {} kmers to chromsome {}'.format(j, last_id))
			j = 0
		fout.write(lines)	# not lines when no kmer mapped
		j += c
		i += 1
		if i % 100000 == 0:
			logger.info('Processed {} sequences'.format(i))
		del lines
		last_id = id
		mapped_cat |= mapped_kmers
		mapped_num += c
		if c > 0:
			mapped_seqs += 1
	logger.info('Processed {} sequences'.format(i))
	logger.info('{} ({:.2%}) sequences contain subgenome-specific kmers'.format(mapped_seqs, mapped_seqs/i))
	mapped_cat, total = len(mapped_cat), len(d_kmers)
	logger.info('{:.2%} of {} subgenome-specific kmers are mapped'.format(mapped_cat/total, total//2))

def chunk_chromfiles(chromfiles, window_size=10e6, overlap=0):
	'''too large genome need to chunk'''
	window_size = int(window_size)
	j = 0
	for chromfile in chromfiles:
#		logger.info('Loading ' + chromfile)
		for rc in SeqIO.parse(open(chromfile), 'fasta'):
			logger.info('Chunking chromsome {}: {:,} bp'.format(rc.id, len(rc.seq)))
			rc_seq = str(rc.seq).upper()
			x = 0
			for i in range(0, len(rc_seq), window_size):
				j += 1
				x += 1
				start = max(0, i-overlap)
				end = i+window_size
				seq = rc_seq[start:end]
				yield rc.id, start, seq	# offset
		#	logger.info('Chunk of {}: {}'.format(rc.id, x))
	logger.info('Chunk {} by window size {}'.format(j, window_size))
def unchunk_chromfiles(chromfiles, exclude=None, rv=False):
	j = 0
	for chromfile in chromfiles:
		for rc in SeqIO.parse(open(chromfile), 'fasta'):
			if exclude and rc.id in exclude:
				continue
			j += 1
			seq = str(rc.seq).upper()	#.replace('U', 'T')
			if rv:
				rc_seq = str(rc.seq.reverse_complement()).upper()
				yield rc.id, 0, (seq, rc_seq)
			else:
				yield rc.id, 0, seq
	#logger.info('{} sequences in total'.format(j))

def count_kmer(chromfiles, d_kmers, lengths, fout=sys.stdout, k=None, 
		sg_names=[], exclude=None, 
		min_prob=0.5, min_count=5, max_fold=1.2,
		ncpu='autodetect', method='map', log=True, chunksize=None):
	if k is None:
		for key in d_kmers.keys():
			k = len(key)
			break
	chunks = unchunk_chromfiles(chromfiles, exclude=exclude, rv=True)
	iterable = ((id, seqs, k, d_kmers,lengths, min_prob, min_count, max_fold) for id, _, seqs in chunks)
	i,j = 0,0
	d_counts = {}
	for id, counts in pool_func(map_kmer_sum, iterable, 
					processors=ncpu, method=method, chunksize=chunksize):
	#	logger.info( (id, counts))
		i += 1
		if i % 100000 == 0:
			logger.info('Processed {} sequences'.format(i))
		if counts is None:
			continue
		j += 1
		d_counts[id] = counts
		print(id, counts, counts/lengths, file=fout)
	logger.info('Processed {} sequences'.format(i))
	logger.info('Identified {} ({:.2%}) shared sequences'.format(j, j/i))
	return d_counts
	
def map_kmer_sum(args):
	id, seqs, k, d_kmers,lengths, min_prob, min_count, max_fold = args
	i = 0
	counts = []
	for seq in seqs:
		seq_len = len(seq)
	#	logger.info( (id, seq_len))
		for s, kmer in _get_kmer(seq, k):
			try: _counts = d_kmers[kmer]
			except KeyError: continue
			i += 1
			counts += [_counts]
			# if i == 1:
				# counts = _counts
			# else:
				# counts =[v1+v2 for v1, v2 in zip(counts, _counts)] #+= _counts
	counts = np.array(counts).sum(axis=0)
	#logger.info( (id, counts, counts/lengths))
	if i/seq_len < min_prob:	# exclude too low copy
		return id, None
	if min(counts/i) < min_count:	# exclude too low copy
		return id, None
	ratios = sorted(counts/lengths)
	if ratios[-1] / ratios[0] > max_fold:	# exclude too different
		return id, None
	return id, counts
		
def map_kmer_each4(args):
	'''bin counts'''
	id, offset, seq, k, d_kmers, bin_size, sg_names = args
	size = len(seq) + offset
	mapped_kmers = set([])
	c = 0	# mapped kmers
	d_bin = OrderedDict()
#	logger.info( (id, size))
	for s, kmer in _get_kmer(seq, k):
		try: sg = d_kmers[kmer]
		except KeyError: continue
		c += 1	# count
		s += offset
		bin = s // bin_size
		if bin not in d_bin:
			d_bin[bin] = OrderedDict((v, 0) for v in sg_names)
		d_bin[bin][sg] += 1	# count by sg
		mapped_kmers.add(kmer)
	
	lines = []
	for bin, d_counts in d_bin.items():
		s = bin * bin_size
		e = min(s + bin_size, size)
		counts = d_counts.values()
		line = [id, s,e] + list(counts)
		line = map(str, line)
		line = '\t'.join(line) + '\n'
		lines += [line]
	return id, c, mapped_kmers, ''.join(lines)
	

def _get_kmer(seq, k):
	for i in range(len(seq)):
		kmer = seq[i:i+k]
#		yield i, ''.join(kmer)
		yield i, kmer
