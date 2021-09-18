import sys
import numpy as np

def main(inTsv=sys.argv[1], outStat=sys.stdout):
	d_count = {}
	ids, sgs = set([]), set([])
	for line in open(inTsv):
		if line.startswith('#'):
			continue
		temp = line.strip().split()
		id, subgenome, p_value, counts = temp
		ann = id.split('-')[0]
		counts = list(map(int, counts.split(',')))
		counts = np.array(counts)
		key = (ann, subgenome)
		if key not in d_count:
			d_count[key] = [1, counts]
		else:
			d_count[key][0] += 1
			d_count[key][1] += counts
		ids.add(key[0])
		sgs.add(key[1])
	for ann in sorted(ids):
		num = []
		for i, sg in enumerate(sorted(sgs)):
			key = (ann, sg)
			if key in d_count:
				_num, _count = d_count[key]
			else:
				_num, _count = 0, np.array([0]*len(sgs))
			num += [_num]
			if i == 0:
				count = _count
			else:
				count += _count
		line = [ann] + num + list(count)
		outStat.write('\t'.join(map(str, line))+'\n')

if __name__ == '__main__':
	main()
	
