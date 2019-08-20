#set user libraries
import ipywidgets as widgets
from IPython.display import display, HTML, clear_output, Image
from ipywidgets import HBox, Label
import subprocess
# import ipyleaflet

import pysmm,os
from pysmm.derive_SM import get_map
import sys
import ee
import itertools
from IPython.core.display import HTML 

ee.Initialize()
path = os.path.dirname(pysmm.__file__)
user = path.split("/")[2]

def intdates(arg_index):
    return [int(s) for s in sys.argv[arg_index].split(', ')]

year1 = intdates(1)
month1 = intdates(2)
day1 = intdates(3)

# year1= sys.argv[1].split(', ')
# month1= sys.argv[2].split(', ')
# day1= sys.argv[3].split(', ')

minlon= float(sys.argv[4])
minlat= float(sys.argv[5])
maxlon= float(sys.argv[6])
maxlat= float(sys.argv[7])



# year1= [2018]
# month1= [6]
# day1= [1]

for a, b, c in itertools.product( year1, month1, day1):

    print a
    print b
    print c

    get_map(minlon , minlat ,maxlon, maxlat,
            '/data/home/' + user + '/',
            sampling=100,
        year=a, month=b, day=c,
       # year=int(sys.argv[1]), month=int(sys.argv[2]), day=int(sys.argv[3]),
            tracknr=None,
            tempfilter=True,
            mask='Globcover',
            masksnow=False,
            overwrite=True)
    print('123')
print('4567')