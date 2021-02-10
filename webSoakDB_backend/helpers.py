from API.models import Library, LibraryPlate, Compounds, SourceWell, Proposals, LibrarySubset
import django.db
import csv
import django.core.exceptions


#prepate list of libraries formatted as froms.ChoiceField choices
def create_lib_selection():
	output = []
	for library in Library.objects.filter(public=True):
		lib = (library.id, library.name)
		output.append(lib)
	return output
	
#prepare a list of library plates in compound selection
def get_selection_details(proposal):
	plate_dict = {}
	for library in proposal.libraries.all():
		if library.public: 
			origin = "in-house library"
		else:
			origin = "user-submitted library"
		
		for plate in LibraryPlate.objects.filter(library=library, current=True):
			plate_dict[plate] = origin
	
	return plate_dict

#import data from a csv file
def upload_plate(file_name, plate):
	with open(file_name, newline='') as csvfile:
		compound_reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		i = 0
		for row in compound_reader:
			try:
				compound = Compound.objects.get(code = row[0])	
			except django.core.exceptions.ObjectDoesNotExist:
				compound = Compound.objects.create(smiles = row[2], code = row[0])
				
			s_well = SourceWell.objects.create(compound = compound, library_plate = plate, well = row[1])
			
			try:
				s_well.concentration = int(row[3])
				s_well.save()
			except (IndexError, ValueError):
				pass
			i += 1

def create_subset(file_name, library, name, origin):
	new_subset = LibrarySubset.objects.create(library=library, name=name, origin=origin)
	
	with open(file_name, newline='') as csvfile:
		compound_reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		for row in compound_reader:
			try:
				compound = Compound.objects.get(code = row[0])
				new_subset.compounds.add(compound)
				new_subset.save()
			except django.core.exceptions.ObjectDoesNotExist:
				print('No such compound found')

'''
def copy_compounds_to_experiment(proposal):
	selection = CompoundSelection.objects.get(proposal=proposal)
	
	for plate in selection.plates:
		for compound in plate.compounds:
			library = compound.library_plate.library
			plate = compound.library_plate.name
			well = compound.well
			code = compound.compound.code
			smiles = compound.compound.smiles
			concentration = compound.concentration
			print('data:', library, plate, well, code, smiles, concentration)
'''
