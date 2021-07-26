import sys, os
import argparse
#from xopen import xopen as open
from collections import OrderedDict
import numpy as np
from .small_tools import open_file as open
__version__='0.1'
__LastModified__='20190115'
__Example__=''

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter,\
		epilog="Version: {}\nLast Modification Date: {}".format(__version__,__LastModified__),\
		version="Version: {}".format(__version__),\
		description="Example: {}".format(__Example__))
	parser.add_argument('--genome-fasta', metavar='genome fasta file', type=str, default=None, help="")
	parser.add_argument('--chr-mark', metavar='chr string', type=str, default='chr', help="")
	parser.add_argument('--gene-gff3', metavar='gene gff3 file', type=str, default=None, help="")
	parser.add_argument('--gene-features', nargs='+', metavar='gene featurs to parse', type=str, default=['gene'], help="")
	parser.add_argument('--out-dir', metavar='out dir', type=str, default='.', help="")
	parser.add_argument('--window-size', metavar='window-size', type=int, default=1000000, help="")
	parser.add_argument('--repeat-out', metavar='RepeatMasker out file', type=str, default=None, help="")
	parser.add_argument('--repeat-class', nargs='+', metavar='repeat class [ClassI, ClassII]', type=str, default=['ClassI', 'ClassII'], help="")
	parser.add_argument('--repeat-order', nargs='+', metavar='repeat class [LTR, LINE, DNA ...]', type=str, default=None, help="")
	parser.add_argument('--repeat-superfamily', nargs='+', metavar='repeat superfamily [Copia, Gypsy...]', type=str, default=None, help="")
	parser.add_argument('--by-sites', action='store_true', default=False, help="by sites instead by elements for gene or TE")
	parser.add_argument('--no_trim', action='store_true', default=False, help="not trim abnormal values")
	parser.add_argument('--variant-vcf', metavar='variant vcf file', type=str, default=None, help="")
	parser.add_argument('--variant-is-het', action='store_true', default=None, help="")
	parser.add_argument('--mcscan-collinearity', metavar='mcscan collinearity file', type=str, default=None, help="")
	parser.add_argument('--mcscan-gff', metavar='mcscan gff file', type=str, default=None, help="")
	parser.add_argument('--mcscan-chrmap', metavar='mcscan chr map file', type=str, default=None, help="")
	parser.add_argument('--mcscan-kaks', metavar='mcscan kaks file', type=str, default=None, help="")
	parser.add_argument('--max-ks', metavar='max ks; block with median_ks > --max-ks will be removed', type=str, default=None, help="")
	parser.add_argument('--min_genes', metavar='min genes', type=int, default=10, help="min number of genes in a block")
	parser.add_argument('--sp1', metavar='species 1', type=str, default=None, help="")
	parser.add_argument('--sp2', metavar='species 2', type=str, default=None, help="")
	parser.add_argument('--sp', metavar='species to target --mcscan-chrmap', type=str, default=None, help="")
	parser.add_argument('--rechr', action='store_true', default=False, help="rename chromosome id according to --mcscan-chrmap")
	parser.add_argument('--bed', metavar='bed file', type=str, default=None, help="")
	args = parser.parse_args()
	return args

global d_rechr
d_rechr = {}
global no_trim
no_trim = False

def main(args):
	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)
	if args.rechr:
		#if args.mcscan_chrmap is None or args.sp is None:
		#	raise ValueError('--mcscan-chrmap or --sp is not specified for --rechr')
		d_rechr = parse_chrmap(args.mcscan_chrmap, args.sp)
	if args.genome_fasta is not None:
		outfile1 = '%s/genome_karyotype.txt' % (args.out_dir, )
		outfile2 = '%s/genome_gc.txt' % (args.out_dir, )
		genome_base(args.genome_fasta, out_karyotype=outfile1, out_gc=outfile2, \
				window_size=args.window_size, chr_mark=args.chr_mark, rechrom=args.rechr)
	if args.gene_gff3 is not None:
		outfile = '%s/gene_density.txt' % (args.out_dir, )
		featurs = set(args.gene_features)
		gene_density(args.gene_gff3, outfile, window_size=args.window_size, featurs=featurs, 
			by_sites=args.by_sites, )
	if args.bed is not None:
		outfile = '%s/%s_density.txt' % (args.out_dir, os.path.basename(args.bed))
		bed_density(args.bed, outfile, window_size=args.window_size, by_sites=args.by_sites,
			)
	if args.repeat_out is not None:
		if args.repeat_class is not None:
			outfile = '%s/repeat_%s_density.txt' % (args.out_dir, '_'.join(args.repeat_class))
			repeat_density(args.repeat_out, outfile, window_size=args.window_size, \
					repeat_class=args.repeat_class, repeat_order=None, repeat_superfamily=None, 
					by_sites=args.by_sites, )
		if args.repeat_order is not None:
			outfile = '%s/repeat_%s_density.txt' % (args.out_dir, '_'.join(args.repeat_order))
			repeat_density(args.repeat_out, outfile, window_size=args.window_size, \
					repeat_class=None, repeat_order=args.repeat_order, repeat_superfamily=None, 
					by_sites=args.by_sites, )
		if args.repeat_superfamily is not None:
			outfile = '%s/repeat_%s_density.txt' % (args.out_dir, '_'.join(args.repeat_superfamily))
			repeat_density(args.repeat_out, outfile, window_size=args.window_size, \
					repeat_class=None, repeat_order=None, repeat_superfamily=args.repeat_superfamily, 
					by_sites=args.by_sites, )
	if args.variant_vcf is not None:
		outfile = '%s/variant_density.txt' % (args.out_dir, )
		variant_density(args.variant_vcf, outfile, window_size=args.window_size, 
			variant_is_het=args.variant_is_het,)
	if (args.gene_gff3 is not None or args.sp1 is not None) and args.mcscan_collinearity is not None:
		outfile = '%s/segment_link.txt' % (args.out_dir, )
		collinearity(args.mcscan_collinearity, outfile, gene_gff3=args.gene_gff3, 
				mcscan_gff=args.mcscan_gff, mcscan_chrmap=args.mcscan_chrmap, 
				mcscan_kaks=args.mcscan_kaks, max_ks=args.max_ks,
				min_genes=args.min_genes, rechrom=args.rechr, 
				sp1=args.sp1, sp2=args.sp2)
def parse_chrmap(mcscan_chrmap, sp):
	d = {}
	for line in open(mcscan_chrmap):
		temp = line.strip().split()
		new_chr, raw_chr, spn, ng = temp[:4]
		raw_chr = raw_chr.split('|')[-1]
		if spn != sp:
			continue
		d[raw_chr] = new_chr
	return d
def rechr(chr):
	return d_rechr.get(chr, chr)
def collinearity(mcscan_collinearity, outfile, gene_gff3=None, min_genes=20, rechrom=False, 
			mcscan_kaks=None, max_ks=None,
			mcscan_gff=None, mcscan_chrmap=None, sp1=None, sp2=None, d_rechr={}):
	'''gene_gff3 or mcscan_gff + chrmap + sp1 [+ sp2]'''
	if sp1 is not None and sp2 is None:
		sp2 = sp1
	from mcscan import Collinearity
	if gene_gff3 is not None:
		d_genes = gff2dict(gene_gff3)
	else:
		d_genes = None
	if rechrom:
		blocks = Collinearity(mcscan_collinearity, gff=mcscan_gff, kaks=mcscan_kaks)
	else:
		blocks = Collinearity(mcscan_collinearity, gff=mcscan_gff, chrmap=mcscan_chrmap, kaks=mcscan_kaks)
	f = open(outfile, 'w')
	lines = []
	for rc in blocks.parse():
		if rc.N < min_genes:
			continue
		if max_ks is not None and rc.median_ks > max_ks:
			continue
		if d_genes is not None:
			from_g1 = rc.head1.split('|')[1]
			from_g2 = rc.head2.split('|')[1]
			to_g1   = rc.tail1.split('|')[1]
			to_g2   = rc.tail2.split('|')[1]
		if d_genes is not None and not (from_g1 in d_genes and from_g2 in d_genes):
			continue
		if sp1 is not None and not sorted([rc.species1, rc.species2]) == sorted([sp1,sp2]):
			continue
		if d_genes is not None:
			chr1, start11, end11 = d_genes[from_g1]
			chr1, start12, end12 = d_genes[to_g1]
			chr2, start21, end21 = d_genes[from_g2]
			chr2, start22, end22 = d_genes[to_g2]
		elif mcscan_gff is not None:
			chr1, start11, end12 = rc.chr1, rc.start1, rc.end1
			chr2, start21, end22 = rc.chr2, rc.start2, rc.end2
			chr1 = chr1.split('|')[-1]
			chr2 = chr2.split('|')[-1]
		if rc.strand == '-':
			start21, end22 = end22, start21
		line = [rechr(chr1), start11, end12, rechr(chr2), start21, end22, rc.N]
		lines += [line]
	lines = sorted(lines, key=lambda x: x[-1])
	for line in lines:
		line = line[:-1]
		line = list(map(str, line))
		print(' '.join(line), file=f)
	f.close()
def gff2dict(gene_gff3):
	d = {}
	for line in open(gene_gff3):
		if line.startswith('#'):
			continue
		temp = line.strip().split('\t')
		CHR = temp[0]
		TYPE = temp[2]
		if not TYPE == 'gene':
			continue
		START = temp[3]
		END = temp[4]
		INFO = temp[8]
		d_info = dict([v.split('=') for v in INFO.split(';')])
		ID = d_info['ID']
		d[ID] = [CHR, START, END]
	return d
def genome_base(genome_fasta, out_karyotype, out_gc, window_size=None, chr_mark=None, rechrom=False):
	from Bio import SeqIO
	from Bio.SeqUtils import GC
	colors = ['chr'+ str(chr) for chr in list(range(1, 23))+['X', 'Y']]
	f1 = open(out_karyotype, 'w')
	f2 = open(out_gc, 'w')
	i = 0
	for rc in SeqIO.parse(open(genome_fasta), 'fasta'):
		seq_len = len(rc.seq)
		for start in range(0, seq_len, window_size):
			end = start + window_size
			if end > seq_len:
				end = seq_len
			gc = GC(rc.seq[start:end])
			line2 = [rc.id, start, end, gc]
			line2 = list(map(str, line2))
			print(' '.join(line2), file=f2)
		if chr_mark is not None and not rc.id.lower().startswith(chr_mark.lower()):
			continue
		try:
			color = colors[i]
		except IndexError:
			color = colors[i-25]
		if rechrom and rc.id not in d_rechr:
			continue
		line1 = ['chr', '-', rechr(rc.id), rechr(rc.id), 0, seq_len, color]
		line1 = list(map(str, line1))
		print(' '.join(line1), file=f1)
		i += 1
	f1.close()
	f2.close()
def genomes_base(genome_fastas, out_karyotype):
	from Bio import SeqIO
	colors = ['chr'+ str(chr) for chr in list(range(1, 23))+['X', 'Y']]
	f1 = open(out_karyotype, 'w')
	i = 0
	for genome_fasta in genome_fastas:
		for rc in SeqIO.parse(open(genome_fasta), 'fasta'):
			seq_len = len(rc.seq)
			color = colors[i%len(colors)]
			line1 = ['chr', '-', rechr(rc.id), rechr(rc.id), 0, seq_len, color]
			line1 = map(str, line1)
			print(' '.join(line1), file=f1)
			i += 1
	f1.close()
def circos_plot(genomes, bedfiles, wddir='circos', window_size=100000):
	from vualize_akchr import colors_rgb
	from .RunCmdsMP import run_cmd
	from .small_tools import mkdirs
	datadir = '{}/data'.format(wddir)
	mkdirs(datadir)
	karyotype_file = '{}/genome_karyotype.txt'.format(datadir)
	genomes_base(genomes, karyotype_file)
	outpre = '{}/subgenome'.format(datadir)
	d_outfiles = bed_density_by_col(bedfiles, outpre, window_size=window_size)
	conf_file = datadir = '{}/histogram.conf'.format(wddir)
	fout = open(conf_file, 'w')
	fout.write('<plots>\n\n')
	start = 0.99
	step = 0.5/len(d_outfiles)
	for i, (key, datafile) in enumerate(sorted(d_outfiles.items())):
		r1, r0 = start, start-step
		color = colors_rgb[i]
		circle = '''<plot>
type        = line
file        = {datafile}
r1          = {r1}r
r0          = {r0}r
min         = 0
extend_bin  = no
color = no
fill_color  = {color}
thickness   = 1p
orientation = out

<axes>
<axis>
color     = black
thickness = 2
position  = 0
</axis>
</axes>

</plot>\n\n'''.format(datafile=datafile, r1=r1, r0=r0, color=color)
		fout.write(circle)
		start = r0
	fout.write('</plots>\n\n')
	fout.close()
	cmd = 'cd {} && circos -conf ./circos.conf'.format(wddir)
	run_cmd(cmd, log=True)
	annofile = '{}/../circos_legend.txt'.format(wddir)
	with open(annofile, 'w') as fout:
		fout.write('Rings from outer to inner:\n\t1. Karyotypes\n')
		for i, sg in enumerate(sorted(d_outfiles.keys())):
			fout.write('\t{}. Density of {}-specific kmers\n'.format(i+2, sg))
		fout.write('Window size: {} bp\n'.format(window_size))

def _bed_density(inBed, window_size=None, chr_col=0, start_col=1, end_col=2, based=0, 
		by_sites=False, split_by=None):
	if split_by is not None:	# split by one col
		d_counts = {}
	d_count = OrderedDict()
	for line in open(inBed):
		temp = line.strip().split()
		try:
			CHR, START = temp[chr_col], temp[start_col]
			START = int(START) - based
			if end_col is not None:
				END = temp[end_col]
				END = int(END)
		except IndexError: continue
		except ValueError: continue
		if by_sites and end_col is not None:
			for POS in range(START, END):
				BIN = POS // window_size
				add_pos(d_count, CHR, BIN, POS)
			continue
		if split_by is not None:
			key = temp[split_by]
			try: d_count = d_counts[key]
			except KeyError: d_counts[key] = d_count = OrderedDict()
		BIN = START // window_size
		try: d_count[CHR][BIN] += 1
		except KeyError:
			try: d_count[CHR][BIN] = 1
			except KeyError:
				d_count[CHR] = {BIN: 1}
	if by_sites:
		count_pos(d_count)
	if split_by:
		return d_counts
	return d_count
def bed_density_by_col(inBed, outpre, keycol=3, window_size=100000, **kargs):
	d_counts = _bed_density(inBed, window_size=window_size, split_by=keycol, **kargs)
	d_outfiles = {}
	for key, d_count in d_counts.items():
		outfile = '{}.{}.txt'.format(outpre, key)
		write_density(d_count, outfile, window_size)
		d_outfiles[key] = outfile
	return d_outfiles

def bed_density(inBed, outfile, window_size=100000, **kargs):
	d_count = _bed_density(inBed, window_size=window_size, **kargs)
	write_density(d_count, outfile, window_size)

def variant_density(variant_vcf, outfile, window_size=None, variant_is_het=None):
	d_count = OrderedDict()
	for line in open(variant_vcf):
		if line.startswith('#'):
			continue
		temp = line.split('\t')
		CHR, POS = temp[:2]
		POS = int(POS)
		try:
			sample = temp[9]
			if variant_is_het is not None and is_het(sample) != variant_is_het:
				continue
		except IndexError:
			pass
		BIN = POS / window_size
		try: d_count[CHR][BIN] += 1
		except KeyError:
			try: d_count[CHR][BIN] = 1
			except KeyError:
				d_count[CHR] = {BIN: 1}
	write_density(d_count, outfile, window_size)

def is_het(sample):
	gt = sample.split(':')[0]
	if '/' in set(gt):
		sep = '/'
	elif '|' in set(gt):
		sep = '|'
	else:
		return 
	if len(set(gt.split('/'))) > 1:
		return True
def repeat_density(repeat_out, outfile, window_size=None, 
		repeat_class=None, repeat_order=None, repeat_superfamily=None, by_sites=False):
	import re
	d_class = dict([
			('LTR', 'ClassI'),
			('LINE', 'ClassI'),
			('SINE', 'ClassI'),
#			('', 'ClassI'),
			('DNA', 'ClassII'),
			('RC', 'ClassII'),
			('Unknown', 'Others'),
			('Satellite', 'Others'),
			('Simple_repeat', 'Others'),
			('Low_complexity', 'Others'),
			('rRNA', 'Others'),
			('tRNA', 'Others'),
			('snRNA', 'Others'),
			('Other', 'Others'),
			('Unspecified', 'Others')
			])
	d_count = OrderedDict()
	for line in open(repeat_out):
		line = line.strip()
		if not re.compile(r'^\d').match(line):
			continue
		temp = line.split()
		CHR = temp[4]
		START = int(temp[5])
		TYPE = temp[10]
		ORDER = TYPE.split('/')[0]
		ORDER = ORDER.rstrip('?')
		try: SUPERFAMILY = TYPE.split('/')[1]
		except IndexError: SUPERFAMILY = ''
		if ORDER not in d_class:
			print('unknown TE order:', ORDER, file=sys.stderr)
			continue
		if repeat_class is not None and not d_class[ORDER] in set(repeat_class):
			continue
		if repeat_order is not None and not ORDER in set(repeat_order):
			continue
		if repeat_superfamily is not None and not SUPERFAMILY in set(repeat_superfamily):
			continue
		if by_sites:
			END = int(temp[6])
			for POS in range(START-1, END):
				BIN = POS // window_size
				add_pos(d_count, CHR, BIN, POS)
			continue
		BIN = START // window_size
		try: d_count[CHR][BIN] += 1
		except KeyError:
			try: d_count[CHR][BIN] = 1
			except KeyError:
				d_count[CHR] = {BIN: 1}
	if by_sites:
		count_pos(d_count)
	write_density(d_count, outfile, window_size)	
def write_density(d_count, outfile, window_size):
	f = open(outfile, 'w')
	counts = []
	for CHR, d_bin in list(d_count.items()):
		last_BIN = -1
		for BIN, count in sorted(d_bin.items()):
			counts += [count]
			if BIN - last_BIN > 1:
				for bin in range(last_BIN+1, BIN):
					count = 0
					d_count[CHR][bin] = count
					counts += [count]
			last_BIN = BIN
	if not no_trim:
		upper,lower = abnormal(counts)
		print('using cutoff: upper {} and lower {}'.format(upper,lower), file=sys.stderr)
	for CHR, d_bin in list(d_count.items()):
		for BIN, count in sorted(d_bin.items()):
			if not no_trim and count > upper:
				count = min(upper*1.1, count)
			elif not no_trim and count < lower:
				count = max(lower*0.9, count)
			start = int(BIN * window_size)
			end = start + window_size
			line = [rechr(CHR), start, end, count]
			line = list(map(str, line))
			print(' '.join(line), file=f)
	f.close()
def abnormal(data, k=1.5):
	q1 = np.percentile(data, 25)
	q3 = np.percentile(data, 75)
	iqr = q3 - q1
	upper = q3 + iqr*k
	upper = np.percentile(data, 99)
	lower = np.percentile(data, 1)
	return upper, lower
def gene_density(gene_gff3, outfile, window_size=None, featurs=None, by_sites=False):
	d_count = OrderedDict()
	for line in open(gene_gff3):
		if line.startswith('#'):
			continue
		temp = line.strip().split('\t')
		CHR = temp[0]
		TYPE = temp[2]
		START = int(temp[3])
		if not TYPE in featurs:
			continue
		if by_sites:
			END = int(temp[4])
			for POS in range(START-1, END):
				BIN = POS / window_size
				add_pos(d_count, CHR, BIN, POS)
			continue
		BIN = START / window_size
		try: d_count[CHR][BIN] += 1
		except KeyError:
			try: d_count[CHR][BIN] = 1
			except KeyError:
				d_count[CHR] = {BIN: 1}
	if by_sites:
		count_pos(d_count)
	write_density(d_count, outfile, window_size)
def add_pos(d_count, CHR, BIN, POS):
	try: d_count[CHR][BIN].add(POS)
	except KeyError:
		try: d_count[CHR][BIN] = {POS}
		except KeyError:
			d_count[CHR] = {BIN: {POS}}
def count_pos(d_count):
	for CHR, d_bin in list(d_count.items()):
		for BIN, POSs in list(d_bin.items()):
			d_count[CHR][BIN] = len(POSs)
if __name__ == '__main__':
	main(makeArgparse())
