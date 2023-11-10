from reportlab.lib.units import cm
from reportlab.lib import colors
import numpy as np
from matplotlib import cm
import matplotlib

COLORS_HEX = ['#f9c00c', '#00b9f1', '#7200da', '#f9320c',
	'#980000','#00ffff','#0000ff','#ff0000','#4a86e8','#ff9900','#ffff00',
	'#00ff00','#9900ff','#ff00ff','#20124d','#274e13','#000000','#cccccc',
	'#7f6000','#a64d79','#6aa84f','#fff2cc','#47a952','#3ea6b6','#a5b805','#8f9276','#ca8d7c']
white = colors.white

class HexColors:
	def __init__(self, colors_hex=None):
		if colors_hex is None:
			colors_hex = COLORS_HEX
		elif isinstance(colors_hex, str):
			colors_hex = colors_hex.split(',')
		self.colors_hex = colors_hex
	@property
	def colors_rgb(self):
		colors_lst = [colors.HexColor(v) for v in self.colors_hex]
		return [','.join(map(str, v.bitmap_rgb())) for v in colors_lst]
	@property 
	def colors_r(self):
		return 'c({})'.format(','.join(map(repr, self.olors_hex)))

class Colors:
	def __init__(self, n):
		self.n = n
		self.colors = self.get_colors()
		
	def get_colors(self):
		if self.n <=10:
			colors = cm.get_cmap('tab10')
		elif self.n <= 20:
			colors = cm.get_cmap('tab20')
		else:
			colors = cm.get_cmap('viridis', self.n)
		colors = colors(np.linspace(0, 1, self.n))
		return colors
	def to_hex(self):
		colors_hex = [matplotlib.colors.to_hex(c) for c in self.colors]
		return colors_hex
	def to_rgb(self):
		colors_hex = self.to_hex()
		colors_lst = [colors.HexColor(v) for v in colors_hex]
		colors_rgb = [','.join(map(str, v.bitmap_rgb())) for v in colors_lst]
		return colors_rgb


