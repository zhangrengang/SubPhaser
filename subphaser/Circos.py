import sys, os, re, shutil, math
import argparse
#from xopen import xopen as open
import collections
from collections import OrderedDict
import numpy as np
#from .small_tools import open_file as open
from xopen import xopen as open
from .small_tools import lazy_open


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
# defualt colors
defualt_colors = ['chr'+ str(chr) for chr in list(range(1, 23))+['X', 'Y']]

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
def genome_base(genome_fasta, out_karyotype, out_gc, window_size=None, chr_mark=None, rechrom=False, color=None):
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
		if chr_mark is not None and not re.compile(chr_mark, re.I).match(rc.id):
			continue
		try:
			chr_color = colors[i]
		except IndexError:
			chr_color = colors[i%len(colors)]
		if rechrom and rc.id not in d_rechr:
			continue
		line1 = ['chr', '-', rechr(rc.id), rechr(rc.id), 0, seq_len, chr_color]
		line1 = list(map(str, line1))
		print(' '.join(line1), file=f1)
		i += 1
	if i == 0:
		raise ValueError('No chromosomes are recognized. Please reset `-chr_prefix`')
	f1.close()
	f2.close()
def genomes_base(genome_fastas, out_karyotype, d_cmap={}):
	from Bio import SeqIO
	colors = ['chr'+ str(chr) for chr in list(range(1, 23))+['X', 'Y']]
	f1 = open(out_karyotype, 'w')
	i = 0
	for genome_fasta in genome_fastas:
		for rc in SeqIO.parse(open(genome_fasta), 'fasta'):
			seq_len = len(rc.seq)
			color = colors[i%len(colors)]
			color = d_cmap.get(rc.id, color)
			line1 = ['chr', '-', rechr(rc.id), rechr(rc.id), 0, seq_len, color]
			line1 = map(str, line1)
			print(' '.join(line1), file=f1)
			i += 1
	f1.close()


		
CIRCLE = '''<plot>
type        = {type}
file        = {datafile}
r1          = {r1}r
r0          = {r0}r
min         = 0
extend_bin  = no
color = no
fill_color  = {color}
thickness   = 0
orientation = out
stroke_color     = undef
stroke_thickness = 0

# <axes>
# <axis>
# color     = black
# thickness = 2
# position  = 0
# </axis>
# </axes>

</plot>\n\n'''
class CircleLegend:
	def __init__(self, labels, colors, title=None, **kargs):
		self.labels = labels
		self.colors = colors
		self.title = title
		self.kargs = kargs
		self.ncols = math.ceil(len(labels) / 6)
	def plot(self, ax):
		for clade, color in zip(self.labels, self.colors):
			ax.barh(0, 0, height=0, color=color, left=0, align='center', label=clade)
		ncols = len(self.labels)//5 + 1
		ax.legend(loc='upper left',fancybox=False, frameon=False, ncol=self.ncols)
		ax.xaxis.set_tick_params(length=0)
		ax.spines['right'].set_color('none')
		ax.spines['top'].set_color('none')
		ax.spines['left'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.xaxis.set_ticks([])
		ax.yaxis.set_ticks([])
		ax.set_title(self.title, loc='left')
class CirclesLegend:
	def __init__(self, legends):
		self.legends = legends
	def plot(self, outfig):
		import matplotlib.pyplot as plt
		plt.switch_backend('agg')
		n = len(self.legends)
		plt.figure()

		for i, legend in enumerate(self.legends):
			ax = plt.subplot(n, 1, i+1)
			legend.plot(ax)
		plt.savefig(outfig, bbox_inches='tight', transparent=True)
		plt.close()
		
def centomics_plot(genome, wddir='circos', tr_bed='', tr_labels='',
		hic_bed='', chip_bed='', 
		prefix='circos', figfmt='pdf',
		chr_prefix='chr', window_size=100000):
	from .RunCmdsMP import run_cmd
	from .colors import Colors
	from .small_tools import mkdirs
	
	datadir = '{}/data'.format(wddir)
	mkdirs(datadir)
	
	legends = []
	# karyotype
	karyotype_file = '{}/genome_karyotype.txt'.format(datadir)
	gc_file = '{}/genome_gc.txt'.format(datadir, )
	genome_base(genome, karyotype_file, gc_file, chr_mark=chr_prefix, window_size=window_size)
	
	# tr
	tr_file = '{}/tr_density.txt'.format(datadir, )
	_n = n = stack_bed(tr_bed, tr_file, window_size=window_size)
	
	
	# histogram
	conf_file = '{}/histogram.conf'.format(wddir)
	fout = open(conf_file, 'w')
	fout.write('<plots>\n\n')
	
	n_circles = 1
	if hic_bed:
		n_circles += 2-1
		n += 1
	if chip_bed:
		n_circles += 1
		n += 1
	start = 0.99
	step = round(0.55/n_circles, 2)
	
	# circle 1
	colors = ['c{}'.format(i) for i in range(1, n+1)]
	colors_mat = Colors(n)
	colors_hex = colors_mat.to_hex()
	colors_rgb = colors_mat.to_rgb()
	_colors, d_cmap = fmt_color(colors, colors_rgb)
	create_color_conf(wddir, d_cmap)
	r1, r0 = start, start-step
	color = ','.join(colors[:_n])
	circle = CIRCLE.format(type='histogram', datafile=tr_file, r1=r1, r0=r0, color=color)
	fout.write(circle)
	start = r0-0.01
	legend = CircleLegend(labels=tr_labels, colors=colors_hex, title='Tandem Repeats')
	legends += [legend]
	
	xcolors = list(zip(colors, colors_hex))
	colors = xcolors[_n:]
	# hic
	if hic_bed:
		# inter
		hic_file = '{}/hic_inter.density.txt'.format(datadir, )
		stack_bed(hic_bed[0], hic_file, window_size=window_size)
		r1, r0 = start, start-step
		color, color_hex = colors[0]
		circle = CIRCLE.format(type='histogram', datafile=hic_file, r1=r1, r0=r0, color=color)
		fout.write(circle)
		start = r0-0.01
		legend = CircleLegend(labels=['inter-chromosomal interaction'], colors=[color_hex], title='Hi-C')
		legends += [legend]
		# intra
		# hic_file = '{}/hic_intra.density.txt'.format(datadir, )
		# stack_bed(hic_bed[1], hic_file, window_size=window_size)
		# r1, r0 = start, start-step
		# color = 'lred'
		# circle = CIRCLE.format(type='histogram', datafile=hic_file, r1=r1, r0=r0, color=color)
		# fout.write(circle)
		# start = r0-0.01
		colors = colors[1:]
	
	# chip
	if chip_bed:
		chip_file = '{}/chip_density.txt'.format(datadir, )
		stack_bed(chip_bed, chip_file, window_size=window_size)
		r1, r0 = start, start-step
		color, color_hex = colors[0] #'green'
		circle = CIRCLE.format(type='histogram', datafile=chip_file, r1=r1, r0=r0, color=color)
		fout.write(circle)
		start = r0-0.01
		legend = CircleLegend(labels=['coverage depth'], colors=[color_hex], title='ChIP-Seq')
		legends += [legend]
	
	fout.write('</plots>\n\n')
	fout.close()
	
	# plot
	cmd = 'cd {} && circos -conf ./circos.conf'.format(wddir)
	exit = run_cmd(cmd, log=True)
	
	svgfile = '{}/circos.svg'.format(wddir)
	fmt_svg(svgfile)
	figfile = '{}/circos.{}'.format(wddir, figfmt)
	if figfmt == 'pdf':
		svg2pdf(svgfile, figfile)
	figfmts = set([figfmt, 'png'])
	for figfmt in figfmts:
		# link file
		figfile = '{}/circos.{}'.format(wddir, figfmt)
		dstfig = '{}.{}'.format(prefix, figfmt)
		try: os.remove(dstfig)
		except FileNotFoundError: pass
		try: os.link(figfile, dstfig)
		except PermissionError: shutil.move(figfile, dstfig)
	
	# legend
	lefig = '{}_legend.pdf'.format(prefix)
	CirclesLegend(legends).plot(lefig)
	
	annofile = '{}_legend.txt'.format(prefix)
	fout = open(annofile, 'w')
	fout.write('Rings from outer to inner:\n\t1. Karyotypes\n')
	ring = 2
	legend = 'Density of tandem repeats'
	fout.write('\t{}. {}\n'.format(ring, legend))
	if hic_bed:
		ring += 1
		legend = 'Density of Hi-C links'
		fout.write('\t{}. {}\n'.format(ring, legend))
	if chip_bed:
		ring += 1
		legend = 'Density of ChIP reads'
		fout.write('\t{}. {}\n'.format(ring, legend))
	fout.close()
	
	
def stack_bed(bedfile, outfile, window_size=100000, trim=True, high_tile=99, low_tile=0):
	bins, counts = stack_matrix(bedfile, window_size=window_size)
	# limit
	if trim:
		counts = np.array(counts)
		uppers = []
	#	print(counts.shape, counts[:10])
		for i in range(counts.shape[1]):
			_count = counts[:,i]
			upper, lower = abnormal(_count, high_tile=high_tile, low_tile=low_tile)
			uppers += [upper]

	fout = open(outfile, 'w')
	for bin, count in zip(bins, counts):
		n = len(count)
		count = limit_upper(count, uppers)
		count = ','.join(map(str, count))
		line = list(bin) + [count]
		line = map(str, line)
		print('\t'.join(line), file=fout)
	fout.close()
	return n
def limit_upper(count, uppers):
	return [min(c,u) for c,u in zip(count, uppers)]
	
# subphaser
def circos_plot(genomes, wddir='circos', bedfile='', sg_color=None,
		sg_lines=[], d_sg={}, prefix='circos', figfmt='pdf',
		ltr_lines=[], enrich_ltr_bedlines=[[]],	# list
		pafs=[], paf_offsets={}, min_block=10000, # blocks
		window_size=100000):
#	from .colors import colors_rgb
	from .RunCmdsMP import run_cmd
	from .small_tools import mkdirs
	colors_rgb = sg_color.colors_rgb
	sgs = sorted(set(d_sg.values()))
	if len(sgs) > len(colors_rgb):
		colors_rgb = colors_rgb * (len(sgs)//len(colors_rgb) + 1)
	colors, d_cmap = fmt_color(sgs, colors_rgb)
	datadir = '{}/data'.format(wddir)
	mkdirs(datadir)
	
	# karyotype
	karyotype_file = '{}/genome_karyotype.txt'.format(datadir)
	genomes_base(genomes, karyotype_file, d_cmap=d_sg)

	# bedfile, split by keys
	outpre = '{}/subgenome'.format(datadir)
	d_outfiles = stack_bed_density(bedfile, outpre, colnames=sgs, window_size=window_size)
	
	# colors
	create_color_conf(wddir, d_cmap)
	
	# histogram
	conf_file = '{}/histogram.conf'.format(wddir)
	
	
	fout = open(conf_file, 'w')
	fout.write('<plots>\n\n')

	n_circles = len(d_outfiles)
	if sg_lines:
		n_circles += 2
	if ltr_lines :
		n_circles += 1
	start = 0.99
	step = round(0.55/n_circles, 2)
	if sg_lines:
		# ratio and enrich
		ratio_file, enrich_file = out_sg_lines(sg_lines, datadir)
		# circle 1
		r1, r0 = start, start-step
		color = ','.join(colors)
		circle = CIRCLE.format(type='histogram', datafile=enrich_file, r1=r1, r0=r0, color=color)
		fout.write(circle)
		start = r0-0.01
		# circle 2
		r1, r0 = start, start-step
		color = ','.join(colors[:-1])
		circle = CIRCLE.format(type='histogram', datafile=ratio_file, r1=r1, r0=r0, color=color)
		fout.write(circle)
		start = r0-0.01
	# circle 3 - n+2
	for i, (key, datafile) in enumerate(sorted(d_outfiles.items())):
		r1, r0 = start, start-step
		color = colors_rgb[i]
		circle = CIRCLE.format(type='histogram', datafile=datafile, r1=r1, r0=r0, color=color)
		fout.write(circle)
		start = r0-0.01

	# circle LTR
	if ltr_lines:
		ltr_file = '{}/ltr_density.txt'.format(datadir)
	#	bed_density(ltr_lines, ltr_file, window_size=window_size)
		bed_density_minus(ltr_lines,enrich_ltr_bedlines, ltr_file, window_size=window_size)
		start = start-0.01
		r1, r0 = start, start-step
		color = ','.join(colors[:-1]) + ',grey' #colors_rgb[i+1]
		circle = CIRCLE.format(type='histogram', datafile=ltr_file, r1=r1, r0=r0, color=color)
		fout.write(circle)
		start = r0-0.01
	
	fout.write('</plots>\n\n')
	fout.close()

	# Links
	conf_file = '{}/link.conf'.format(wddir)
	fout = open(conf_file, 'w')
	fout.write('''
radius    = {}r
crest     = 0.25
ribbon    = yes
twist     = yes
flat      = yes
thickness = 1
color     = 255,0,0,0.5
bezier_radius        = 0r
bezier_radius_purity = 1\n'''.format(start))

	if pafs:
		linkfile = '{}/block_link.txt'.format(datadir)
		paf2blocks(pafs, linkfile, paf_offsets, min_block=min_block)
		fout.write('''<link>
file       = {}
</link>
'''.format(linkfile))
	
	#fout.write('</links>')
	fout.close()
	
	# plot
	cmd = 'cd {} && circos -conf ./circos.conf'.format(wddir)
	exit = run_cmd(cmd, log=True)
	#print(exit)
	
	svgfile = '{}/circos.svg'.format(wddir)
	fmt_svg(svgfile)
	figfile = '{}/circos.{}'.format(wddir, figfmt)
	if figfmt == 'pdf':
		svg2pdf(svgfile, figfile)
	figfmts = set([figfmt, 'png'])
	for figfmt in figfmts:
		# link file
		figfile = '{}/circos.{}'.format(wddir, figfmt)
		dstfig = '{}.{}'.format(prefix, figfmt)
		try: os.remove(dstfig)
		except FileNotFoundError: pass
		try: os.link(figfile, dstfig)
		except PermissionError: shutil.move(figfile, dstfig)
	
	# legend
	annofile = '{}/../circos_legend.txt'.format(wddir)
	with open(annofile, 'w') as fout:
		fout.write('Rings from outer to inner:\n\t1. Karyotypes\n')
		ring = 1
		if sg_lines:
			for legend in ['Enriched subgenome', 'Normalized proportion of each subgenome']:
				ring += 1
				fout.write('\t{}. {}\n'.format(ring, legend))
		for i, sg in enumerate(sorted(d_outfiles.keys())):
			ring += 1
			fout.write('\t{}. Density of {}-specific kmers\n'.format(ring, sg))
		if ltr_lines:
			ring += 1
			legend = 'Density of LTR-RTs'
			fout.write('\t{}. {}\n'.format(ring, legend))
		if pafs:
			ring += 1
			legend = 'Homologous blocks'
			fout.write('\t{}. {}\n'.format(ring, legend))
		if pafs:
			fout.write('Window size: {} bp\n'.format(window_size))
			

def create_color_conf(wddir, d_cmap):
	# colors
	conf_file = '{}/colors.conf'.format(wddir)
	with open(conf_file, 'w') as f:
		f.write('<colors>\n\n')
		for name, color in sorted(d_cmap.items()):
			f.write('{} = {}\n'.format(name, color))
		f.write('</colors>\n')
def fmt_color(sgs, color_set, white='white'):
	if len(sgs) > len(color_set):
		color_set = color_set * (len(sgs)//len(color_set) + 1)
	colors, d_cmap = [], {}
	for i, sg in enumerate(sgs):
		cname = str(sg)
		d_cmap[cname] = color_set[i]
		colors += [cname]
	colors += [white]
#	d_cmap[white] = white
	return colors, d_cmap

def out_sg_lines(sg_lines, datadir, ratio_col=6, enrich_col=7):
	ratio_file = '{}/sg_ratio.txt'.format(datadir)
	enrich_file = '{}/sg_enrich.txt'.format(datadir)
	fr = open(ratio_file, 'w')
	fe = open(enrich_file, 'w')
	for line in sg_lines:
		chrom, start, end = line[:3]
		ratio = line[ratio_col]
		enrich = line[enrich_col]
		for f, dat in zip((fr, fe), (ratio, enrich)):
			line = [chrom, start, end, dat]
			line = map(str, line)
			f.write('\t'.join(line) + '\n')
	fr.close()
	fe.close()
	return ratio_file, enrich_file
	
def fmt_svg(svgfile):
	import re
	from .small_tools import backup_file
	bksvgfile, svgfile = backup_file(svgfile)
	svg = open(bksvgfile).read()
	subsvg = re.compile(r'<tspan.*?>(\S+)<\/tspan>').sub(r'_\1', svg)
	with open(svgfile, 'w') as fout:
		fout.write(subsvg)
	os.remove(bksvgfile)
	return svgfile
	
def svg2pdf(svgfile, pdfile):
	from svglib.svglib import svg2rlg
	from reportlab.graphics import renderPDF
	drawing = svg2rlg(svgfile)
	renderPDF.drawToFile(drawing, pdfile)
	
	
def paf2blocks(paf_groups, linkfile, paf_offsets={}, min_block=10000, colors=None):
	from .Paf import PafParser
	#print(vars(), file=sys.stderr)
	if colors is None:
		colors = defualt_colors
	f = open(linkfile, 'w')
	for i, paf_group in enumerate(paf_groups):
		blocks = []
		for paf in paf_group:
			for rc in PafParser(paf):
				if not rc.is_primary:
					continue
				rc.size = rc.alen
				if rc.size < min_block:
					continue
				blocks += [rc]
		color = colors[i%len(colors)]
		for rc in sorted(blocks, key=lambda x:x.size):
			options = 'color={}'.format(color)
			if rc.qid in paf_offsets:
				rc.qid, offset = paf_offsets[rc.qid]
				rc.qstart, rc.qend = rc.qstart+offset, rc.qend+offset
			if rc.tid in paf_offsets:
				rc.tid, offset = paf_offsets[rc.tid]
				rc.tstart, rc.tend = rc.tstart+offset, rc.tend+offset
			line = [rc.qid, rc.qstart, rc.qend, rc.tid, rc.tstart, rc.tend, options]
			line = map(str, line)
			f.write(' '.join(line) + '\n')
	f.close()

def bed_density_minus(totBed, setBeds, outfile, window_size=10000, **kargs):
	d_count = _bed_density_minus(totBed, setBeds, window_size=window_size, **kargs)
	write_density(d_count, outfile, window_size, trim=True)
	
def _bed_density_minus(totBed, setBeds, **kargs):
	import copy
	d_count0 = _bed_density(totBed, **kargs)
	d_count_out = copy.deepcopy(d_count0)
	for setBed in setBeds:
		d_count = _bed_density(setBed, **kargs)
		for CHR, d_bin in list(d_count_out.items()):
			for BIN, count in list(d_bin.items()):
				try: value = d_count[CHR][BIN]
				except KeyError: value = 0
				try: d_count_out[CHR][BIN].append(value)
				except: d_count_out[CHR][BIN] = [value]
	for CHR, d_bin in list(d_count_out.items()):
		for BIN, count in list(d_bin.items()):
			total = d_count0[CHR][BIN]
			values = d_count_out[CHR][BIN]
			last = total - sum(values)
			values.append(last)
			d_count_out[CHR][BIN] = ','.join(map(str, values))
	return d_count_out

def _bed_density(inBed, window_size=None, chr_col=0, start_col=1, end_col=2, based=0, 
		split_by=None, by_sites=False, stack=False):
	'''
split_by: split and count by one key column
by_sites: count by site instead of line
stack: stack short bins into long window
	'''
	if split_by is not None:	# split by one col
		d_counts = OrderedDict()
	d_count = OrderedDict()
	for line in lazy_open(inBed):
		if isinstance(line, str):
			if line.startswith('#'):
				continue
			temp = line.strip().split()
		elif isinstance(line, collections.abc.Iterable):
			temp = list(line)
		try:
			CHR, START = temp[chr_col], temp[start_col]
			START = int(START) - based
			if end_col is not None:
				END = temp[end_col]
				END = int(END)
		except IndexError: continue
		except ValueError: continue
		if stack:
			values = list(map(tr_numeric, temp[3:]))	# assume bed format
			values = np.array(values)
			BIN = int(START // window_size)
			if CHR not in d_count:
				d_count[CHR] = OrderedDict()
			try: d_count[CHR][BIN] += values
			except KeyError: d_count[CHR][BIN] = values
			continue
		if by_sites and end_col is not None:
			for POS in range(START, END):
				BIN = int(POS // window_size)
				add_pos(d_count, CHR, BIN, POS)
			continue
		if split_by is not None:
			key = temp[split_by]
			try: d_count = d_counts[key]
			except KeyError: d_counts[key] = d_count = OrderedDict()
		BIN = int(START // window_size)
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
def tr_numeric(val):
    try: return int(val)
    except:
        try: return float(val)
        except: return val

def bed_density_by_col(inBed, outpre, keycol=3, window_size=100000, **kargs):
	d_counts = _bed_density(inBed, window_size=window_size, split_by=keycol, **kargs)
	d_outfiles = {}
	for key, d_count in d_counts.items():
		outfile = '{}.{}.txt'.format(outpre, key)
		write_density(d_count, outfile, window_size)
		d_outfiles[key] = outfile
	return d_outfiles
def stack_bed_density(inBedCount, outpre, colnames, window_size=100000, trim=True):
	coords, counts = stack_matrix(inBedCount, window_size=window_size)
	d_outfiles = {key: '{}.{}.txt'.format(outpre, key) for key in colnames}
	d_handles = {key: open(d_outfiles[key], 'w') for key in colnames}
	
	# limit
	if trim:
		d_count = {key: [] for key in colnames}
		for coord, _count in zip(coords, counts):
			for key, count in zip(colnames, _count):
				d_count[key] += [count]
		d_upper = {}
		for key, _count in d_count.items():
			upper, lower = abnormal(_count)
			d_upper[key] = upper
			print('using cutoff: upper {} for {}'.format(upper,key), file=sys.stderr)
	# 
	for coord, _count in zip(coords, counts):
		coord = list(coord)
		assert len(colnames) == len(_count), '{} != {}'.format(len(colnames), len(_count))
		for key, count in zip(colnames, _count):
			if trim:
				count = min(count, d_upper[key])
			line = coord + [count]
			line = map(str, line)
			fout = d_handles[key]
			print(' '.join(line), file=fout)
	for f in d_handles.values():
		f.close()
	return d_outfiles

def counts2matrix(inBed, keys=None, keycol=3, window_size=100000, **kargs):
	d_counts = _bed_density(inBed, window_size=window_size, split_by=keycol, **kargs)
	_keys = set([])
	d_bincounts = {}
	# convert data
	for key, d_count in d_counts.items():
		_keys.add(key)
		for CHR, d_bin in list(d_count.items()):
			for BIN, count in sorted(d_bin.items()):
				_key = (CHR, BIN)
				try: d_bincounts[_key][key] = count
				except KeyError: d_bincounts[_key] = {key: count}
	if keys is None:
		keys = sorted(_keys)

	coords, counts = [], []
	for (CHR, BIN), d_count in d_bincounts.items():
		start = int(BIN * window_size)
		end = int(start+window_size)
		count = [d_count.get(key, 0) for key in keys]
		coords += [(CHR, start, end)]
		counts += [count]
	return coords, counts
def stack_matrix(inBedCount, window_size=100000):
	'''stack short bins'''
	d_bincounts = _bed_density(inBedCount, window_size=window_size, stack=True)
	coords, counts = [], []
	for CHR, d_count in d_bincounts.items():
		for BIN, count in d_count.items():
			start = int(BIN * window_size)
			end = int(start+window_size)
			coords += [(CHR, start, end)]
			count = list(count)	# [ordered count]
			counts += [count]
	return coords, counts

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
def write_density(d_count, outfile, window_size, trim=None):
	if trim is None:
		_no_trim = no_trim
	else:
		_no_trim = not trim
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
	if not _no_trim:
		try: 
			upper,lower = abnormal(counts)
		#	print('using cutoff: upper {} and lower {}'.format(upper,lower), file=sys.stderr)
		except Exception as e: 
			_no_trim = True
		#	print('Exception: {} for output `{}`'.format(e, outfile), file=sys.stderr)
	for CHR, d_bin in list(d_count.items()):
		for BIN, count in sorted(d_bin.items()):
			if not _no_trim and count > upper:
				count = min(upper*1.1, count)
			elif not _no_trim and count < lower:
				count = max(lower*0.9, count)
			start = int(BIN * window_size)
			end = start + window_size
			line = [rechr(CHR), start, end, count]
			line = list(map(str, line))
			print(' '.join(line), file=f)
	f.close()
def abnormal(data, k=1.5, high_tile=99, low_tile=1):
	q1 = np.percentile(data, 25)
	q3 = np.percentile(data, 75)
	iqr = q3 - q1
	upper = q3 + iqr*k
	upper = np.percentile(data, high_tile)
	lower = np.percentile(data, low_tile)
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
