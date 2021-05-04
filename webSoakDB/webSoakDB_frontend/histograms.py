#THESE ARE ONLY NOTES FOR FUTURE CODE
#In future this file will contain functions needed to create histograms
#of compound properties in a library plate
from API.models import Library, LibraryPlate
#making graphs
import matplotlib #for generating images
import seaborn as sns
import os
'''
mwts = [] #to suppress errors
g = sns.displot(mwts)# create a histogram from the list mwts
g = sns.displot(mwts, bins=30)


matplotlib.pyplot.show() #opens a new window and displays the last generated graph as an image
g.savefig(name) #saves grapg g as an image with the name name (recognises file extention)
g.savefig("pic.svg")

with sns.plotting_context("notebook", font_scale=10):
     g=sns.displot(mol_wt)

#make a graph with a larger font


#checking age of the file
import os
import datetime

s = os.stat('9_mol_wt.png')#get file metadata
s.st_mtime #date of the last modification
value = datetime.datetime.fromtimestamp(s.st_mtime)#turn it into a timestamp
now = datetime.datetime.now()#get current timestamp
age = now - value #how long ago the file was last modified - produces timedelta object
str(age) #human-readable timedelta
d = datetime.timedelta(days=7) #make a timedelta object that means the duration of 7 days
age > d #check if the file has been modified in the last 7 days


for filename in os.listdir('.'):
     print(filename)
#print all files in current directory
'''
#######################################
libraries = Library.objects.filter(public=True)
img_path = "images/graphs/"

properties = ["mol_wt", 
			"log_p", 
			"heavy_atom_count", 
			"heavy_atom_mol_wt", 
			"nhoh_count", 
			"no_count", 
			"num_h_acceptors", 
			"num_h_donors", 
			"num_het_atoms",
			"num_rot_bonds",
			"num_val_electrons",
			"ring_count",
			"tpsa"
			]
		

def get_compounds(obj, obj_type):
	if obj_type=="library":
		return get_current_lib_compounds(obj)
	if obj_type=="preset":
		return get_preset_compounds(obj)
	if obj_type=="selection":
		return get_selection_compounds(obj)
	print('get_compounds: Unrecognised type of compounds set')
	return False
	
def get_preset_compounds(preset):
	c = []
	for s in preset.subsets.all():
		print('subset: ', s.name)
		print('size: ', s.compounds.all().count())
		c += s.compounds.all()
	print('preset size: ', len(c))
	return c

def get_current_lib_compounds(library):
	print('library: ', library)
#	for plate in library.plates:
#		print(plate.name + "current :" + plate.current)
	current_plates = LibraryPlate.objects.filter(library = library, current=True)
#	print('current_plates: ', current_plates)
	source_wells = []
	for p in current_plates:
		source_wells += p.compounds.all().filter(active=True)
	
#	print('source_wells: ', source_wells)
	
	return [sw.compound for sw in source_wells]

def get_selection_compounds(selection):
	return []


def dataset(compounds, attr):
	return [getattr(c, attr) for c in compounds]


def make_images(obj, obj_type):
	compounds = get_compounds(obj, obj_type)
	
	for p in properties:
		name = img_path + p + "/" + str(obj.id) + ".svg"
		
		try:
			os.remove(name)
		except OSError:
			pass
			
		data = dataset(compounds, p)
		with sns.plotting_context("notebook", font_scale=1.3):
			if not p in ["log_p", "mol_wt", "heavy_atom_mol_wt", "tpsa"]:
				graph = sns.displot(data, discrete=True)
			else:
				graph = sns.displot(data)
			
		graph.savefig(name, dpi=100)


def set_image_fields(lib):
	for p in properties:
		attr = p + "_graph"
		path = img_path + p + "/" + str(lib.id) + ".svg"
		setattr(lib, attr, path)
	
	lib.save()

