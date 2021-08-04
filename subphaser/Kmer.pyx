#cimport cython
#import numpy as np
# cython: initializedcheck=False
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

cdef dict ambiguous_dna_complement = {
	"A": "T",
	"C": "G",
	"G": "C",
	"T": "A",
	"M": "K",
	"R": "Y",
	"W": "W",
	"S": "S",
	"Y": "R",
	"K": "M",
	"V": "B",
	"H": "D",
	"D": "H",
	"B": "V",
	"X": "X",
	"N": "N",
	}
#cdef dict ambiguous_dna_complement = {
#	65: 84,
#	67: 71,
#	71: 67,
#	84: 65,
#}
#for k, v in list(ambiguous_dna_complement.items()):
#	ambiguous_dna_complement[k.lower()] = v.lower()

cdef str blank = ''

cdef tuple reverse_complement(list seq):
#	cdef char cseq[k]
#	cdef char rcseq[k]
#	cdef char[::1] cseq = np.empty(k, dtype=np.str)
#	cdef char[::1] rcseq = cseq
	cdef list cseq, rcseq
	cdef tuple kmer
	cdef str v
	cdef int i
#	for i in range(k):
#		cseq[i] = ambiguous_dna_complement[seq[i]]
	cseq = [ambiguous_dna_complement[v] for v in seq]
#	for i in range(k):
#		rcseq[i] = cseq[k-i]
	rcseq = cseq[::-1]
	kmer = tuple(rcseq) # blank.join(rcseq)
	return kmer

cdef tuple get_kmer(list seq, int s, int k):
	#DEF N=k
	#cdef char kseq[N]
#	cdef str[::1] kseq = np.empty(k, dtype=np.str)
	cdef list kseq
	cdef tuple kmer, rc_kmer
	cdef int i
	kseq = seq[s:s+k]
#	for i in range(k):
		#print(kseq[i])
		#print(seq[s+i])
#		kseq[i] = seq[s+i]
	kmer = tuple(kseq) #blank.join(kseq)
	rc_kmer = reverse_complement(kseq)
	return kmer, rc_kmer

def _reverse_complement(seq):
	cseq = seq #[ambiguous_dna_complement[v] for v in seq]
	rcseq = cseq[::-1]
	kmer = blank.join(rcseq)
def _get_kmer(seq, si, ei):
	blank = ''
	kseq = seq[si:ei]
	kmer = blank.join(kseq)
	rc_kmer = _reverse_complement(kseq)
# %time x=[ _get_kmer(seq, i, i+k) for i in range(size)]

cdef tuple get_kmer2(list seq, int s, int k):
#	cdef int s,e
#	s,e = i, i+k
	return get_kmer(seq, s, k)

def get_2kmer(list seq, int s, int k):
	return get_kmer(seq, s, k)

def get_kmers(list seq, int size, int k):
	cdef int i
	cdef list xkmers
	for i in range(size):
		yield get_kmer2(seq, i, k)

def get_subgenomes(dict d_kmers, list seq, int k):
	cdef int s,e,i,size
	cdef str sg
	cdef tuple kmers, kmer
	size = len(seq)
#	cdef str sg
#	cdef list lines = []
	for i in range(size):
		s,e = i, i+k
		kmers = get_kmer2(seq, i, k)
#	for s, e, kmers in get_kmers(seq, size, k):
		for kmer in kmers:
#			kmer = kmer.upper()
			try: sg = d_kmers[kmer]
			except KeyError: continue
			yield s, e, sg
#			lines += [(s, e, sg)]
#	return lines

#def get_subgenomes(dict d_kmers, list seq, int size, int k):
#	return _get_subgenomes(d_kmers, seq, size, k)
