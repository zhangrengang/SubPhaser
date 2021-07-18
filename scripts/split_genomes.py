import sys
from Bio import SeqIO
from xopen import xopen as open

def split_genomes(genomes, prefixes, targets, outdir, d_targets=None):
	# allow renaming id split by `|`
	if not d_targets:
		d_targets = {}
		for t in targets:
			temp = t.split('|')
			id, new_id = temp[-1], temp[0]
			d_targets[id] = new_id
	elif set(targets) - set(d_targets):
		for t in set(targets) - set(d_targets):
			temp = t.split('|')
			id, new_id = temp[-1], temp[0]
			d_targets[id] = new_id

	outfas, labels = [], []
	for genome, prefix in zip(genomes, prefixes):
		for rc in SeqIO.parse(open(genome), 'fasta'):
			rc.id = '{}{}'.format(prefix, rc.id)
			if d_targets and rc.id not in d_targets:
				continue
			rc.id = d_targets[rc.id]
			outfa = '{}/{}.fasta'.format(outdir, rc.id)
			with open(outfa, 'w') as fout:
				SeqIO.write(rc, fout, 'fasta')
			outfas += [outfa]
			labels += [rc.id]
	return outfas, labels

