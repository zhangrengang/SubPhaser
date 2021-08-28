import sys
import copy
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from xopen import xopen as open
# from . import Kmer
# from .Kmer import get_subgenomes
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
			
			rc.id = d_targets[rc.id]
			outfa = '{}{}.fasta'.format(outdir, rc.id)
			with open(outfa, 'w') as fout:
				SeqIO.write(rc, fout, 'fasta')
			outfas += [outfa]
			labels += [rc.id]	# new id
			d_size[rc.id] = len(rc.seq)
	return outfas, labels, d_targets2, d_size

# def map_kmer_each(args):
	# rc, d_kmers, k = args
	# seq = list(rc.seq.upper())
	# lines = []
	# for s,e,sg in get_subgenomes(d_kmers, seq, k):
# #		line = [rc.id, s, e, sg]
		# line = [s, sg]
	# #	line = list(map(str, line))
		# lines += [line]
	# return rc.id, lines
# def map_kmer_each2(rc, d_kmers, k, ncpu=4):
	# seq = list(rc.seq.upper())
	# iterable = ((seq, i, k, d_kmers) for i in range(len(seq)))
	# jobs = pool_func(_get_2kmer, iterable, processors=ncpu, )
	# for kmers in jobs:
		# for kmer in kmers:
			# line = (rc.id,) + kmer
			# line = list(map(str, line))
			# yield line
# def _get_2kmer(args):
	# seq, i, k, d_kmers = args
	# s,e = i, i+k
	# kmers = Kmer.get_2kmer(seq, i, k)
	# res = []
	# for kmer in kmers:
		# try: sg = d_kmers[kmer]
		# except KeyError: continue
		# res += [( s, e, sg)]
	# return res
# def map_kmer(chromfiles, d_kmers, fout=sys.stdout, k=None, ncpu='autodetect', method='map'):
	# if k is None:
		# k = len(list(d_kmers.keys())[0])
	# iterable = ((rc, d_kmers, k) for chromfile in chromfiles for rc in SeqIO.parse(open(chromfile), 'fasta'))
	# xlines = []	
	# for (rc, d_kmers, k) in iterable:
		# j = 0
		# for line in map_kmer_each2(rc, d_kmers, k, ncpu=ncpu):
			# print('\t'.join(line), file=fout)
			# j += 1
			# xlines += [line]
		# logger.info('Mapped {} kmers to {}'.format(j, rc.id))
	# return xlines
# def map_kmer2(chromfiles, d_kmers, fout=sys.stdout, k=None, ncpu='autodetect', method='map'):
	# if k is None:
		# k = len(list(d_kmers.keys())[0])
	# iterable = ((rc, d_kmers, k) for chromfile in chromfiles for rc in SeqIO.parse(open(chromfile), 'fasta'))
	# jobs = pool_func(map_kmer_each, iterable, processors=ncpu, method=method)
	# xlines = []
	# for job in jobs:
		# id, lines = job
		# j = 0
# #		for line in lines:
		# for s, sg in lines:
			# line = [id, s, s+k, sg]
			# line = list(map(str, line))
			# print('\t'.join(line), file=fout)
			# j += 1
			# xlines += [line]
	# return xlines

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
	#method='map'
	line = ['#chrom', 'start', 'end'] + sg_names
	fout.write('\t'.join(line)+'\n')
	for id, c, lines in pool_func(map_kmer_each4, iterable, 
					processors=ncpu, method=method, chunksize=chunksize):
		if last_id and id != last_id:
			if log:
				logger.info('Mapped {} kmers to chromsome {}'.format(j, last_id))
			j = 0
		fout.write(lines)	# not lines when no kmer mapped
		j += c
		i += 1
		if i % 10000 == 0:
			logger.info('Processed {} sequences'.format(i))
		del lines
		last_id = id
	# iterable = ((id, start, seq, k, d_kmers) for id, start, seq in chunks)
	# last_id = ''
	# i, j = 0, 0
	# #method='map'
	# for id, c, lines in pool_func(map_kmer_each3, iterable, 
					# processors=ncpu, method=method, chunksize=chunksize):
		# if last_id and id != last_id:
			# if log:
				# logger.info('Mapped {} kmers to chromsome {}'.format(j, last_id))
			# j = 0
		# fout.write(lines)	# not lines when no kmer mapped
		# j += c
		# i += 1
		# if i % 10000 == 0:
			# logger.info('Processed {} sequences'.format(i))
		# del lines
		# last_id = id

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
def unchunk_chromfiles(chromfiles):
	for chromfile in chromfiles:
		for rc in SeqIO.parse(open(chromfile), 'fasta'):
			seq = str(rc.seq).upper()
			yield rc.id, 0, seq

def map_kmer_each4(args):
	'''bin counts'''
	id, offset, seq, k, d_kmers, bin_size, sg_names = args
	size = len(seq) + offset
	c = 0
	d_bin = OrderedDict()
	for s, kmer in _get_kmer(seq, k):
		try: sg = d_kmers[kmer]
		except KeyError: continue
		c += 1	# count
		s += offset
		bin = s // bin_size
		if bin not in d_bin:
			d_bin[bin] = OrderedDict((v, 0) for v in sg_names)
		d_bin[bin][sg] += 1	# count by sg
	
	lines = []
	for bin, d_counts in d_bin.items():
		s = bin * bin_size
		e = s + bin_size
		if e > size:
			e = size
		counts = d_counts.values()
		line = [id, s,e] + list(counts)
		line = map(str, line)
		line = '\t'.join(line) + '\n'
		lines += [line]
	return id, c, ''.join(lines)
	
# def map_kmer_each3(args):
	# id, offset, seq, k, d_kmers = args
	# lines = []
	# c = 0
	# for s, kmer in _get_kmer(seq, k):
		# try: sg = d_kmers[kmer]
		# except KeyError: continue
		# c += 1	# count
		# s += offset
		# e = s + k
		# line = [id, s,e, sg]
		# line = map(str, line)
		# line = '\t'.join(line) + '\n'
		# lines += [line]
	# return id, c, ''.join(lines)

def _get_kmer(seq, k):
	for i in range(len(seq)):
		kmer = seq[i:i+k]
#		yield i, ''.join(kmer)
		yield i, kmer
