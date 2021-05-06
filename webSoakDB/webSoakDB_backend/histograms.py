from API.models import Library, LibraryPlate

from bokeh.embed import components
import numpy as np
from bokeh.plotting import figure

from datetime import datetime

def get_histogram(obj, obj_type, attr):
	compounds = get_compounds(obj, obj_type)
	data = dataset(compounds, attr)
	bins = len(set(data))
	if bins > 20:
		bins = 20
	plot = make_plot(data, bins, properties_dict[attr])
	script, div = components(plot)
	script = '<script src="https://cdn.bokeh.org/bokeh/release/bokeh-2.3.1.min.js"></script>' + script
	
	return '<div style="height: 100px, width: auto">' + script + div + '</div>'


def get_compounds(obj, obj_type):
	if obj_type=="library":
		return get_current_lib_compounds(obj)
	if obj_type=="preset":
		return get_preset_compounds(obj)
	if obj_type=="selection":
		return get_selection_compounds(obj)
	print('get_compounds: Unrecognised type of compounds set')
	return False

def get_selection_compounds(selection):
	#not implemented yet
	pass
	


#the bottleneck
def get_current_lib_compounds(library):
	current_plates = LibraryPlate.objects.filter(library = library, current=True)
	source_wells = []
	for p in current_plates:
		source_wells += p.compounds.all().filter(active=True).prefetch_related("compound")
	return [sw.compound for sw in source_wells] #this is the slow part

#alternative bottleneck
def get_preset_compounds(preset):
	c = []
	for s in preset.subsets.all():
		c += s.compounds.all()
	return c

def dataset(compounds, attr):
	return [getattr(c, attr) for c in compounds]

'''
#I also tried the this:

def get_current_lib_compounds(library):
	current_plates = LibraryPlate.objects.filter(library = library, current=True)
	source_wells = []
	for p in current_plates:
		source_wells += p.compounds.all().filter(active=True)
	return source_wells

def dataset(compounds, attr):
	return [getattr(c.compound, attr) for c in compounds]


The effect was that dataset() became the bottleneck.
I also tried using different loops to get my list of numbers, but I didn't get anything 
significantly better than a list comprehension.

'''

def make_plot(data, bins, title):
	hist, edges = np.histogram(data, density=False, bins=bins)
	p = figure(title=title, tools='', background_fill_color="#fafafa")
	p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
		fill_color="navy", line_color="white", alpha=0.5)

	p.y_range.start = 0
	p.xaxis.axis_label = 'Value'
	p.yaxis.axis_label = 'Compounds'
	p.grid.grid_line_color="white"
	p.height = 280	
	p.width = 280
	return p


properties_dict = {
	"mol_wt" : "Molecular weight" ,
	"tpsa" : "TPSA" ,
	"log_p" : "LogP" ,
	"num_val_electrons" : "Number of valence electrons" ,
	"num_h_acceptors" : "Number of hydrogen bond acceptors" ,
	"num_h_donors" : "Number of hydrogen bond donors" ,
	"num_het_atoms" : "Number of heteroatoms" ,
	"num_rot_bonds" : "Number of rotable bonds" ,
	"ring_count" : "Number of rings" ,
	"heavy_atom_count" : "Number of heavy atoms" ,
	"heavy_atom_mol_wt" : "Average mol. weight of heavy atoms" ,
	"nhoh_count" : "Number of NH and OH" ,
	"no_count" : "Number of nitrogens and oxygens" ,
}


#for performance testing purposes
def get_all_histograms(obj, obj_type):
	compounds = get_compounds(obj, obj_type)
	source = '<div style="display: flex; flex-wrap: wrap"><script src="https://cdn.bokeh.org/bokeh/release/bokeh-2.3.1.min.js"></script>'
	for attr in properties_dict:
		data = dataset(compounds, attr)
		bins = len(set(data))
		if bins > 20:
			bins = 20
		plot = make_plot(data, bins, properties_dict[attr])
		script, div = components(plot)
		source = source + '<div>' + script + div + '</div>'
	
	source = source + '</div>'
	return source
