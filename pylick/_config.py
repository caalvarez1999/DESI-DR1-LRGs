import os
import inspect

# directories
__ROOT__ = './'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])

dir_lib  = './aux/'

__table__ = dir_lib + 'tableall.dat'
__style__ = dir_lib + 'plotbelli.style'

__dirRes__  = '../outputs/'
__dirPlot__ = '../outputs/plots/'
