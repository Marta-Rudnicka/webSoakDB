from bokeh.embed import components
import numpy as np
import scipy.special

from bokeh.layouts import gridplot
from bokeh.plotting import figure, output_file, show

from API.models import LibraryPlate

def make_plot(title, hist, edges):
    p = figure(title=title, tools='', background_fill_color="#fafafa")
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
           fill_color="navy", line_color="white", alpha=0.5)

    p.y_range.start = 0
    p.xaxis.axis_label = 'Value'
    p.yaxis.axis_label = 'Compounds'
    p.grid.grid_line_color="white"
    p.height = 300
    p.width = 300
    return p

from datetime import datetime
def graph():
	plate = LibraryPlate.objects.get(pk=2)
	
	t1 = datetime.now()
	x = [c.compound.mol_wt for c in plate.compounds.all()]
	t2 = datetime.now()
	print('list comprehension: ', t2 - t1)
	#y = [c.compound.ring_count for c in plate.compounds.all()]
	
	t3 = datetime.now()
	unique = set(x)

	bins = len(unique)
	if bins > 20:
		bins = 20
	t4 = datetime.now()
	
	print('getting bins: ', t4 - t3)
	
	t5 = datetime.now()
	hist, edges = np.histogram(x, density=False, bins=bins)
	p = make_plot(plate.barcode, hist, edges)
	t6 = datetime.now()
	print('getting plot: ', t6 - t5)
	
	t7 = datetime.now()
	script, div = components(p)
	x = '<script src="https://cdn.bokeh.org/bokeh/release/bokeh-2.3.1.min.js"></script>'
	script = x + script
	
	y = '<div style="height: 100px, width: auto">' + script + div + '</div>'
	t8 = datetime.now()
	print('making component: ', t8 - t7)
	
	return y
	#return x
