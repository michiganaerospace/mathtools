'''Let's try to make plots look nice, shall we?
To preview the colors defined below, check out:
    http://flatuicolors.com/
'''
# Import some common plotting commands, but not all of pylab.
from pylab import plot, xlabel, ylabel, grid, title, ion, figure, close, xlim,\
        ylim, matplotlib, hold, imshow, subplot, axis, text, legend, savefig, \
        clf


# Define a reasonable, consistent color scheme.
turquoise       = '#1abc9c'
green_sea       = '#16a085'
sunflower       = '#f1c40f'
orange          = '#f39c12'

emerald         = '#2ecc71'
nephritis       = '#27ae60'
carrot          = '#e67e22'
pumpkin         = '#d35400'

peter_river     = '#3498db'
belize_hole     = '#2980b9'
alizarin        = '#e74c3c'
pomegranate     = '#c0392b'

amethyst        = '#9b59b6'
wisteria        = '#8e44ad'
clouds          = '#ecf0f1'
silver          = '#bdc3c7'

wet_asphalt     = '#34495e'
midnight_blue   = '#2c3e50'
concrete        = '#95a5a6'
asbestos        = '#7f8c8d'


def setup_plotting():
    '''Set up plotting. Make things look nice.'''
    # Turn on interactive plotting.
    ion()
    close('all')

    # Set font to something reasonable.
    font = {'family': 'Arial', 'weight': 'light', 'size': 16}
    matplotlib.rc('font', **font)
