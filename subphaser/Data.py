import numpy as np

class LoadData:
	def __init__(self, datafile):
		self.datafile = datafile
	def load_matrix(self):
		self.data = []
		self.colnames = []
		self.rownames = []
		self.d_rows = {}
		for i, line in enumerate(open(self.datafile)):
			temp = line.strip().split()
			if i == 0:
				self.colnames = temp[1:]
				continue
			rkey = temp[0]
			self.rownames += [rkey]
			line = list(map(float, temp[1:]))
			self.d_rows[rkey] = line
			self.data += [line]
		self.data = np.array(self.data)

