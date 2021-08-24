from reportlab.lib.units import cm
from reportlab.lib import colors

colors_hex = ['#f9c00c', '#00b9f1', '#7200da', '#f9320c',
	'#980000','#00ffff','#0000ff','#ff0000','#4a86e8','#ff9900','#ffff00',
	'#00ff00','#9900ff','#ff00ff','#20124d','#274e13','#000000','#cccccc',
	'#7f6000','#a64d79','#6aa84f','#fff2cc','#47a952','#3ea6b6','#a5b805','#8f9276','#ca8d7c']

colors_lst = [colors.HexColor(v) for v in colors_hex]
colors_rgb = [','.join(map(str, v.bitmap_rgb())) for v in colors_lst]
white = colors.white
colors_r = 'c({})'.format(','.join(map(repr, colors_hex)))
