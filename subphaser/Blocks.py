import sys, os
import itertools
from .RunCmdsMP import run_cmd, logger, pool_run
from .small_tools import mk_ckp, check_ckp
from .split_records import cut_seqs

def run_align(sgs, d_chromfiles, outdir, ncpu=8, thread=2, overwrite=False,
		aligner='unimap', opts='-x asm20', d_size={}, max_size=10e6, overlap=100e3):
	opts += ' -t {}'.format(thread)
	cmds = []
	paf_groups = []
	d_offsets = {}
	d_ckp = {}
	for sg in sgs:
		paf_group = []
		for sg1, sg2 in itertools.combinations(sg, 2):
			for chr1, chr2 in itertools.product(sg1, sg2):
				outpaf = '{}{}-{}.paf'.format(outdir, chr1, chr2)
				ckp_file = outpaf + '.ok'
				ckp = check_ckp(ckp_file)
				if ckp and not overwrite:
					d_offset, = ckp
					d_offsets.update(d_offset)
					continue
				fa1, fa2 = d_chromfiles[chr1], d_chromfiles[chr2]
				d_offset = {}
				size = d_size.get(chr2, 0)
				if size > max_size:
					fa2_cut = fa2+'.cut'
					with open(fa2_cut, 'w') as fout:
						d_offset = cut_seqs(fa2, fout, window_size=max_size, window_ovl=overlap)
					fa2 = fa2_cut
					logger.info('Split {} ({} Mb) into {} chunks'.format(chr2, int(size/1e6), len(d_offset)))
				
				cmd = '{} {} {} {} > {} && touch {}'.format(
						  aligner, fa1, fa2, opts, outpaf, ckp_file)
				cmds += [cmd]
				paf_group += [outpaf]
				d_offsets.update(d_offset)
				d_ckp[ckp_file] = d_offset
				
		paf_groups += [paf_group]
	pool_run(cmds, ncpu, log=True)
	
	# write to reload
	for ckp_file, d_offset in d_ckp.items():
		if check_ckp(ckp_file, log=False):
			mk_ckp(ckp_file, d_offset, log=False)
	return paf_groups, d_offsets
