import numpy as np
from matplotlib import cm
from matplotlib import colors

COLORS_HEX = ['#f9c00c', '#00b9f1', '#7200da', '#f9320c', '#00b8a9',
              "#F4A460", '#009999', '#00C02E',
              '#980000', '#00ffff', '#0000ff', '#ff0000', '#4a86e8', '#ff9900', '#ffff00',
              '#00ff00', '#9900ff', '#ff00ff', '#20124d', '#274e13', '#000000', '#cccccc',
              '#7f6000', '#a64d79', '#6aa84f', '#fff2cc', '#47a952', '#3ea6b6', '#a5b805',
              '#8f9276', '#ca8d7c']
white = '#FFFFFF'


class HexColors:
    '''get default colors'''

    def __init__(self, colors_hex=None):
        if colors_hex is None:	# HEX colors
            colors_hex = COLORS_HEX
        elif isinstance(colors_hex, str):
            colors_hex = colors_hex.split(',')
        self.colors_hex = colors_hex  # HEX colors

    @property
    def colors_rgb(self):
        '''RGB colors'''
        return [hex2rgb(v) for v in self.colors_hex]

    @property
    def colors_r(self):
        '''colors for R'''
        return 'c({})'.format(','.join(map(repr, self.colors_hex)))

    def __str__(self):
        return str(self.colors_hex)

    def __repr__(self):
        return str(self.colors_hex)


def hex2rgb(Hex):
    rgb = colors.to_rgb(Hex)
    rgb255 = ','.join(map(str, [int(v*255) for v in rgb]))
    return rgb255


class Colors:
    '''get n colors'''

    def __init__(self, n):
        self.n = n
        self.colors = self.get_colors()

    def get_colors(self):
        if self.n <= 10:
            _colors = cm.get_cmap('tab10')
        elif self.n <= 20:
            _colors = cm.get_cmap('tab20')
        else:
            _colors = cm.get_cmap('viridis', self.n)
        _colors = _colors(np.linspace(0, 1, self.n))
        return _colors

    def to_hex(self):
        '''HEX colors'''
        colors_hex = [colors.to_hex(c) for c in self.colors]
        return colors_hex

    def to_rgb(self):
        '''RGB colors'''
        colors_rgb = [hex2rgb(v) for v in self.colors]
        return colors_rgb
