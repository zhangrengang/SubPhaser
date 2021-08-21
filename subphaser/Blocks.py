import sys
import itertools
from .RunCmdsMP import run_cmd, logger, pool_run
from .small_tools import mk_ckp, check_ckp

def run_align(sgs, d_chromfiles, outdir, ncpu=8, thread=2, overwrite=False,
		aligner='unimap', opts='-x asm20'):
	opts += ' -t {}'.format(thread)
	cmds = []
	paf_groups = []
	for sg in sgs:
		paf_group = []
		for sg1, sg2 in itertools.combinations(sg, 2):
			for chr1, chr2 in itertools.product(sg1, sg2):
				fa1, fa2 = d_chromfiles[chr1], d_chromfiles[chr2]
				outpaf = '{}{}-{}.paf'.format(outdir, chr1, chr2)
				ckp_file = outpaf + '.ok'
				if not check_ckp(ckp_file) or overwrite:
					cmd = '{} {} {} {} > {} && touch {}'.format(
						  aligner, fa1, fa2, opts, outpaf, ckp_file)
					cmds += [cmd]
#				run_cmd(cmd, log=True)
				paf_group += [outpaf]
		paf_groups += [paf_group]
	pool_run(cmds, ncpu, log=True)
#	for stdout, stderr, status in 	pool_run(cmds, ncpu):
#		print(stdout, stderr, status)
	return paf_groups
