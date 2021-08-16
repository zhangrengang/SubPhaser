import sys, io
import re
import numpy as np
from xopen import xopen as open
def best_hit(inPaf=sys.stdin, outPaf=sys.stdout):
	last_id = ''
	for line in inPaf:
		temp = line.rstrip().split('\t')
		qid = temp[0]
		if qid == last_id:
			continue
		outPaf.write(line)
		last_id = qid

class PafRecord():
	'''define each reacord of PAF format
	line = ST-E00276:351:HHVCHCCXY:1:1101:22993:3876	   150	 43	  77	  +	   NC_034747.1	 133730  133522  133556	34	  34	  3	   NM:i:0  ms:i:68 AS:i:68 nn:i:0  tp:A:P  cm:i:2  s1:i:29 s2:i:0  cg:Z:34M
	'''
	def __init__(self, line, max_hang=500, paired=False):
		self.line = line.rstrip()
		temp = line.split('\t')
		qid, qlen, qstart, qend, strand, \
		  tid, tlen, tstart, tend, \
		  match, alen, ascore = temp[:12]
		ilist = [qlen, qstart, qend, qlen, tstart, tend, tlen, match, alen, ascore]
		[qlen, qstart, qend, qlen, tstart, tend, tlen, match, alen, ascore] = list(map(int, ilist))
		self.qid = qid
		self.qlen = qlen
		self.qstart = qstart
		self.qend = qend
		self.strand = strand
		self.tid = tid
		self.tlen = tlen
		self.tstart = tstart
		self.tend = tend
		self.match = match
		self.alen = alen
		self.ascore = ascore
		self._parse_sam(temp[12:])
		self.max_hang = max_hang
		if paired:
			self._deal_paired()
	def _deal_paired(self):
		if self.qid.endswith('/1'):
			self.qid = self.qid[:-2]
		elif self.qid.endswith('/2'):
			self.qid = self.qid[:-2]
			self.strand = '-' if self.strand == '+' else '+'
		else:
			pass
	def _parse_sam(self, fields):
		d_type = {'i':int, 'f':float, 'A':str, 'Z':str, 'H':str, 'B':str}
		self.pairs = {}
		for field in fields:
			temp = field.split(':', 2)
			tag = temp[0]
			_type = temp[1]
		#	value = ':'.join(temp[2:])
			value = temp[2]
			func = d_type.get(_type, str)
			value = func(value)
			self.pairs[tag] = value
		#	exec 'self.%s = value' % (tag)
			setattr(self, tag, value)
	@property
	def q5hang(self):
		return self.qstart
	@property
	def is_primary(self):
		return self.pairs['tp'] == 'P'
	@property
	def q3hang(self):
		return self.qlen - self.qend
	@property
	def t5hang(self):
		return self.tstart
	@property
	def t3hang(self):
		return self.tlen - self.tend
	@property
	def qmatch(self):
		return self.qend-self.qstart
	@property
	def tmatch(self):
		return self.tend-self.tstart
	@property
	def qcov(self):
		return 1.0* self.qmatch / self.qlen
	@property
	def tcov(self):
		return 1.0* self.tmatch / self.tlen
	@property
	def div(self):
		try: return self.dv
		except: return self.de
	@property
	def t_contains_q(self):
		# t ---------
		# q   <--->
		if self.q5hang < self.max_hang and self.q3hang < self.max_hang:
			return True
		return False
	@property
	def q_contains_t(self):
		# t	---
		# q <------->
		if self.t5hang < self.max_hang and self.t3hang < self.max_hang:
			return True
		return False
	@property
	def contain(self):
		return self.t_contains_q or self.q_contains_t
	@property
	def overlap(self):
		if self.strand == '+':
			if   self.t5hang < self.max_hang and self.q3hang < self.max_hang:
				# t   ---> 5
				# q --->   3
				return (5,3)
			elif self.t3hang < self.max_hang and self.q5hang < self.max_hang:
				# t --->   3
				# q   ---> 5
				return (3,5)
		elif self.strand == '-':
			if   self.t5hang < self.max_hang and self.q5hang < self.max_hang:
				# t   ---> 5
				# q <---   5
				return (5,5)
			elif self.t3hang < self.max_hang and self.q3hang < self.max_hang:
				# t --->   3
				# q  <---  3
				return (3,3)
		return False
	def coverage(self):
		if self.t_contains_q:
			return self.qcov
		elif self.q_contains_t:
			return self.tcov
		overlap = self.overlap
		if not overlap:
			return min(self.qcov, self.tcov)
		if   overlap == (5,3):
			tcov = 1.0*(self.tmatch) / self.tend
			qcov = 1.0*(self.qmatch) / (self.qlen-self.qstart)
		elif overlap == (3,5):
			tcov = 1.0*(self.tmatch) / (self.tlen-self.tstart)
			qcov = 1.0*(self.qmatch) / self.qend
		elif overlap == (5,5):
			tcov = 1.0*(self.tmatch) / self.tend
			qcov = 1.0*(self.qmatch) / self.qend
		elif overlap == (3,3):
			tcov = 1.0*(self.tmatch) / (self.tlen-self.tstart)
			qcov = 1.0*(self.qmatch) / (self.qlen-self.qstart)
		return (tcov+qcov) / 2
	def parse_cs(self):
		cs = self.cs
		toffset = self.tstart	# 1~ t, 2~ q
		strand = self.strand
		if strand == '-':
			qoffset = self.qend
		else:
			qoffset = self.qstart
		return CsBlocks(cs, toffset, qoffset, strand)
	
def test_cs(strand='+'):
	line = 'q\t27\t0\t27\t'+strand+'\tr\t27\t0\t27\t24\t30\t60\tcs:Z::6-ata:10+gtc:4*at:3'
	rc = PafRecord(line)
	blocks = rc.parse_cs()
	for block in blocks:
		print(block.tstart, block.tend, block.qstart, block.qend)
class CsBlocks:
	def __init__(self, cs, toffset=0, qoffset=0, strand='+'):
		self.cs = cs
		self.toffset = toffset
		self.qoffset = qoffset
		self.strand = strand
	def __iter__(self):
		return self._parse()
	def _parse(self):
		expression = '(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)'
		i = 0
		toffset = self.toffset
		qoffset = self.qoffset
		for string in re.compile(expression).findall(self.cs):
			i += 1
			block = CsBlock(string, toffset=toffset, qoffset=qoffset, strand=self.strand)
			yield block
			toffset = block.tend
			if self.strand == '-':
				qoffset = block.qstart
			else:
				qoffset = block.qend
	def call(self, chrom, variant=False):	# per block
		from Vcf import VcfRecord
		for block in self:
			if variant and block.type == 'identical':	# only variant sites
				continue
#			chrom = block.tid
			qual, fliter, info, format = '.', '.', '.', 'GT'
			if block.type == 'identical':
				alt = '.'
				for start, base in zip(list(range(block.tstart, block.tend)), block.tseq):
					ref = base
					sample = '0'
					line = [chrom, start+1, '.', ref, alt, qual, fliter, info, format, sample]
					yield VcfRecord(line)
			else:
				ref = block.tseq.upper()
				alt = block.qseq.upper()
				#if ref == 'N' or alt =='N':	# gap
				if 'N' in set(ref) or 'N' in set(alt):
					continue
				sample = '1'
				start = block.tstart
				if block.type == 'ins':
					ref = '-'
				elif block.type == 'del':
					alt = '-'
				line = [chrom, start+1, '.', ref, alt, qual, fliter, info, format, sample]
				yield VcfRecord(line)
class CsBlock:
	def __init__(self, string, toffset=0, qoffset=0, strand='+'):
		self.key = string[0]
		self.value = string[1:]
		self.tstart = self.toffset = toffset
		self.qoffset = qoffset
		if strand =='-':
			self.qend = qoffset
		else:
			self.qstart = qoffset
		self.strand = strand
		self._parse()
	def _parse(self):
		self.match = 0
		if self.key == ':':		# identical, --cs or --cs=short
			self.tlen = self.qlen = self.match = int(self.value)	# 1=ref, 2=qry
			self.tseq = self.qseq = ''
			self.type = 'identical' # dentical block
		elif self.key == '=':	# identical, --cs=long
			self.tlen = self.qlen = self.match = len(self.value)
			self.tseq = self.qseq = self.value
			self.type = 'identical'
		elif self.key == '-':	# delete
			self.tlen, self.qlen = len(self.value), 0
			self.tseq, self.qseq = self.value, ''
			self.type = 'del' # deletion
		elif self.key == '+':	# insert
			self.tlen, self.qlen = 0, len(self.value)
			self.tseq, self.qseq = '', self.value
			self.type = 'ins' # insertion 
		elif self.key == '*':	# snp
			self.tlen = self.qlen = self.match = 1
			self.tseq, self.qseq = self.value[0], self.value[1]
			self.type = 'snp'
		else:
			raise ValueError('unknown mark :{}; must be in [:=+-*].'.format(self.key))
		self.tend = self.tstart + self.tlen  # 0-based
		if self.strand =='-':
			self.qstart = self.qend - self.qlen
		else:
			self.qend = self.qstart + self.qlen
	def is_adj(self, other):
		return self.tstart == other.tend
			
class PafParser():
	'''parse PAF file from minimap2'''
	def __init__(self, inPaf, **record_kargs):
		self.inPaf = inPaf
		#if type(self.inPaf) == file:
		if isinstance(self.inPaf, io.TextIOWrapper):
			self.inPaf = self.inPaf
		else:
			self.inPaf = open(self.inPaf)
		self.record_kargs = record_kargs
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in self.inPaf:
			if line.startswith('#'):
				continue
			yield PafRecord(line, **self.record_kargs)
	def call_vcf(self, fout=sys.stdout, sample="sample", primary=True, min_alen=5000, **kargs):
		print('##fileformat=VCFv4.1', file=fout)
		#lines = []
		d_contigs = {}
		d_lines = {}
		for rc in self:
			if primary and not rc.is_primary:
				continue
			if rc.alen < min_alen:
				continue
			d_contigs[rc.tid] = rc.tlen
			blocks = rc.parse_cs()
			for line in blocks.call(rc.tid, **kargs):
				#line.write(fout)
				#lines += [line]
				try: d_lines[rc.tid] += [line]
				except KeyError: d_lines[rc.tid] = [line]
		for tid, tlen in sorted(d_contigs.items()):
			print('##contig=<ID=' + tid + ',length=' + str(tlen) + '>', file=fout)
		print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=fout)
		print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample, file=fout)
		for tid, lines in sorted(d_lines.items()): #, key=lambda x: (x.CHROM, x.POS)):
			sites = set([])
			lines = sorted(lines, key=lambda x: x.POS)
			new_lines = []
			for line in sorted(lines, key=lambda x: (x.POS, x.REF)):
				site = (line.POS, ) #line.REF, line.ALT)
				if site in sites:	# remove duplicates
					continue
				if line.REF == '-': # ins
					last_line = new_lines.pop()
					try: assert last_line.POS == line.POS - 1
					except AssertionError:	# unknown
						last_line.write(sys.stderr)
						line.write(sys.stderr)
						continue
					line.REF = last_line.REF
					if last_line.ALT[0] == '.':
						line.ALT = (last_line.REF + line.ALT[0], )
					else:
						line.ALT = (last_line.ALT[0] + line.ALT[0], )
					line.POS = line.POS - 1
				elif line.ALT[0] == '-': # del
					last_line = new_lines.pop()
					try: assert last_line.POS == line.POS - 1
					except AssertionError: 
						last_line.write(sys.stderr)
						line.write(sys.stderr)
						continue
					line.REF = last_line.REF + line.REF
					if last_line.ALT[0] == '.':
						line.ALT = (last_line.REF, )
					else:
						line.ALT = (last_line.ALT[0], )
					line.POS = line.POS - 1

				new_lines += [line]
				site = (line.POS, )
				sites.add(site)
			for line in new_lines:
				line.write(fout)
class LinkPafParser(object):
	'''yield records by query'''
	def __init__(self, inPaf, primary=True, sortby='q', **parser_kargs):
		self.inPaf = inPaf
		self.primary = primary
		self.sortby = sortby
		self.parser_kargs = parser_kargs
	def __iter__(self):
		return self._parse()
	def _parse(self):
		last_key = ''
		records = []
		for rc in PafParser(self.inPaf, **self.parser_kargs):
			if self.primary and not rc.is_primary:
				continue
			key = self.by_key(rc)
			if last_key and key != last_key:
				yield LinkPafRecords(records, sortby=self.sortby)
				records = []
			records += [rc]
			last_key = key
		yield LinkPafRecords(records, sortby=self.sortby)
	def by_key(self, rc):
		return rc.qid
class LinkPafParser2(LinkPafParser):
	def __init__(self, inPaf, primary=True, sortby='t', **parser_kargs):
		super(LinkPafParser2, self).__init__(inPaf, primary=primary, sortby=sortby, **parser_kargs)
	def by_key(self, rc):
		return (rc.qid, rc.tid, rc.strand)
	
def linked_reads(inPaf, fout=sys.stdout):
	for records in LinkPafParser(inPaf):
		if not records.has_good_path():
			continue
		print(records.qid, records.path, file=fout)
			
class LinkPafRecords():
	def __init__(self, records, primary=True, sortby=None):
#		if primary:
#			records = [rc for rc in records if rc.is_primary]
		if sortby == 'q':
			self.records = sorted(records, key=lambda x:x.qstart)
		elif sortby == 't':
			self.records = sorted(records, key=lambda x:x.tstart)
		else:
			self.records = records
	def __iter__(self):
		return iter(self.records)
	def __len__(self):
		return len(self.records)
	def __getitem__(self, index):
		if isinstance(index, int):
			return self.records[index]
		else:
			return self.__class__(self.records[index])
	@property
	def qid(self):
		return self.records[0].qid
	@property
	def tid(self):
		return self.records[0].tid
	@property
	def strand(self):
		return self.records[0].strand
	@property
	def tlen(self):
		return self.records[0].tlen
	@property
	def qlen(self):
		return self.records[0].qlen

	@property
	def tstart(self):
		return self.records[0].tstart
	@property
	def tend(self):
		return self.records[-1].tend
	@property
	def qstart(self):
		return min(self.records[0].qstart, self.records[-1].qstart)
	@property
	def qend(self):
		return max(self.records[0].qend, self.records[-1].qend)
	@property
	def key(self):
		return (self.qid, self.qstart, self.qend, self.strand, self.tid, self.tstart, self.tend)
	def has_good_path(self):
		if len(self.records) < 2:
			return False
		ends = self.records[0], self.records[-1]
		internals = self.records[1:-1]
		end_overlap = [record.overlap for record in ends]
		internal_contian = [record.contain for record in internals]
		if all(end_overlap+internal_contian):
			return True
		return False
	@property
	def path(self):
		pth = [(rec.strand, rec.tid) for rec in self.records]
#		return [int('{}{}'.format(rec.strand, rec.tid)) for rec in self.records]
		try: return [int('{}{}'.format(*x)) for x in pth]
		except ValueError: return pth
	@property
	def score(self):
		'''total chain score'''
		return sum([rec.s1 for rec in self.records])
	@property
	def match(self):
		'''total match'''
		return sum([rec.match for rec in self.records])
	def sum_property(self, key)	:
		return sum([getattr(self, key) for rec in self.records])
	@property
	def alen(self):
		return self.sum_property('alen')
	@property
	def dv(self):
		if len(self.dvs) == 0:
			return 0
		return np.mean(self.dvs)
	@property
	def dvs(self):
		try: return [rec.dv for rec in self.records]
		except AttributeError: 
			try: return [rec.de for rec in self.records]
			except AttributeError: return [0 for rec in self.records]
	
	@property
	def similarity(self):
		alens = [rec.alen for rec in self.records]
		dvs = self.dvs
		return 1 - sum([alen*dv for alen, dv in zip(alens, dvs)]) / sum(alens)
	@property
	def cov(self):
		return np.mean([rec.coverage for rec in self.records])
	@property
	def idt(self):
		return 1 - self.dv

def paf2ovl(inPaf, outOvl, max_hang=5000, min_match=500):
	for rc in PafParser(inPaf):
		line = \
		(qry_name, qry_strand, qry_length, qry_beg, qry_end,
		 sbj_name, sbj_strand, sbj_length, sbj_beg, sbj_end,
		 score, idenity, n_mat, n_mis, n_ins, n_del, cigar) = \
		(rc.qid,   '+',	rc.qlen, rc.qstart-1, rc.qend,
		 rc.tid,   rc.strand,  rc.tlen, rc.tstart-1, rc.tend,
		 rc.match, 1-rc.dv, rc.match, 0,  0,  0,   '0M')
		qhang = min(qry_beg, qry_length-qry_end)
		thang = min(sbj_beg, sbj_length-sbj_end)
		if not (qhang <= max_hang or thang <= max_hang): # "or" for wtclp
			continue
		if n_mat < min_match:
			continue
	#	line = (qry_name, qry_strand, qry_length, qry_beg, qry_end,
	#			sbj_name, sbj_strand, sbj_length, sbj_beg, sbj_end,
	#			score, idenity, n_mat, n_mis, n_ins, n_del, cigar)
		line = list(map(str, line))
		print('\t'.join(line), file=outOvl)
	
def filter(inPaf=sys.stdin, outPaf=sys.stdout, max_dv=1, min_qcov=0, min_match=0, min_qlen=0,
		min_tcov=0, is_primary=True):
	for rc in PafParser(inPaf):
		if is_primary and not rc.is_primary:
			continue
		if rc.dv > max_dv or rc.qcov < min_qcov or rc.match < min_match or \
		   rc.qlen < min_qlen or rc.tcov < min_tcov:
			continue
		print(rc.line, file=outPaf)

def filter_alignment(inPaf, outPaf, max_dv=0.025, max_hang=500, min_match=500):
	ctq, cqt, o35, o53, o55, o33 = 0,0,0,0,0,0
	for rc in PafParser(inPaf):
		if rc.dv > max_dv or rc.match < min_match:
			continue
#		if (rc.qstart < max_hang or rc.qlen - rc.qend < max_hang) or \
#			(rc.tstart < max_hang or rc.tlen - rc.tend < max_hang):
#			print >> outPaf, rc.line
		out = False
		q_start_hang, q_end_hang = rc.qstart, rc.qlen - rc.qend
		t_start_hang, t_end_hang = rc.tstart, rc.tlen - rc.tend
		if q_start_hang < max_hang and q_end_hang < max_hang:  # t contains q
				# q   --->
				# t --------
			out = True
			ctq += 1
		elif t_start_hang < max_hang and t_end_hang < max_hang:  # q contains t
				# q --------
				# t   --->
			out = True
			cqt += 1
		elif rc.strand == '+':
			if q_end_hang < max_hang and t_start_hang < max_hang:
				# q --->
				# t   --->
				out = True
				o35 += 1
			elif q_start_hang < max_hang and t_end_hang < max_hang:
				# q   --->
				# t --->
				out = True
				o53 += 1
		elif rc.strand == '-':
			if q_start_hang < max_hang and t_start_hang < max_hang:
				# q <---
				# t	--->
				out = True
				o55 += 1
			elif q_end_hang < max_hang and t_end_hang < max_hang:
				# q --->
				# t   <---
				out = True
				o33 += 1
		if out:
			print(rc.line, file=outPaf)
	print('{} t contains q, {} q contains t, {} 35, {} 53, {} 55, {} 33'.format(ctq, cqt, o35, o53, o55, o33), file=sys.stderr)
def check_pafs(inPafs):
	for inPaf in inPafs:
		print('\tchecking', inPaf, file=sys.stderr)
		try:
			i = 0
			for rc in PafParser(inPaf):
				i += 1
				continue
		except Exception as e:
			print('ERROR:', inPaf, '--', 'LINE', i, '--', e, file=sys.stdout)

def main():
	import sys
	subcmd = sys.argv[1]
	if subcmd == 'paf2ovl':
		try: inPaf = sys.argv[2]
		except IndexError: inPaf = sys.stdin
		paf2ovl(inPaf=inPaf, outOvl=sys.stdout)
	elif subcmd == 'filter_organ':
		try: max_dv = float(sys.argv[2])
		except IndexError: max_dv = 0.025
		filter_alignment(inPaf=sys.stdin, outPaf=sys.stdout, max_dv=max_dv)
	elif subcmd == 'filter':
		sys.argv.pop(1)
		args = makeArgs()
		list(filter(**args.__dict__))
	elif subcmd == 'check':
		inPafs = sys.argv[2:]
		check_pafs(inPafs)
	elif subcmd == 'best':
		best_hit()
	elif subcmd == 'test':
		try: strand = sys.argv[2]
		except IndexError: strand= '+'
		test_cs(strand)
	elif subcmd == 'linked_reads':
		inPaf = sys.argv[2]
		linked_reads(inPaf)
	elif subcmd == 'call_vcf':
		try: inPaf = sys.argv[2]
		except IndexError: inPaf = sys.stdin
		PafParser(inPaf).call_vcf()
	else:
		raise ValueError('unknown command: {}'.format(subcmd))
def makeArgs():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", action="store",type=str,
					dest="inPaf", default=sys.stdin,
					help="input [default=%(default)s]")
	parser.add_argument("-o", action="store",type=str,
					dest="outPaf", default=sys.stdout,
					help="output [default=%(default)s]")
	parser.add_argument("-max_dv", action="store",type=float,
					default=1,
					help="max dv [default=%(default)s]")
	parser.add_argument("-min_qcov", action="store",type=float,
					default=0,
					help="min qcov [default=%(default)s]")
	parser.add_argument("-min_tcov", action="store",type=float,
                    default=0,
                    help="min tcov [default=%(default)s]")

	parser.add_argument("-min_match", action="store",type=int,
					default=0,
					help="min match [default=%(default)s]")
	parser.add_argument("-min_qlen", action="store",type=int,
					default=0,
					help="min qlen [default=%(default)s]")
	parser.add_argument("-is_primary", action="store_true", default=False,
					help="is primary alignment [default=%(default)s]")
	options = parser.parse_args()
	return options
if __name__ == '__main__':
	main()
