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

#import data from a csv file - should be used after validation with data_is_valid
def upload_plate(file_name, plate):
	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		for row in compound_reader:
			try:
				compound = Compounds.objects.get(code = row[0].strip(), smiles = row[2].strip())	
			except django.core.exceptions.ObjectDoesNotExist:
				print('Compound not found, creating new one.')
				compound = Compounds.objects.create(smiles = row[2].strip(), code = row[0].strip())
				
			s_well = SourceWell.objects.create(compound = compound, library_plate = plate, well = row[1].strip())
			
			try:
				s_well.concentration = int(row[3])
				s_well.save()
			except (IndexError, ValueError):
				pass

def upload_subset(file_name, library_id, name, origin):
	library = Library.objects.get(id=library_id)
	new_subset = LibrarySubset.objects.create(library=library, name=name, origin=origin)
	
	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		for row in compound_reader:
			try:
				compound = Compounds.objects.get(code = row[0], smiles=row[1])
				new_subset.compounds.add(compound)
				new_subset.save()
			except django.core.exceptions.ObjectDoesNotExist:
				print('No such compound found')
			except django.core.exceptions.MultipleObjectsReturned:
				duplicates = Compounds.objects.filter(code = row[0], smiles=row[1])
				print('duplicates: ', duplicates)
				for compound in duplicates:
					print('smiles: ', compound.smiles)
					print('locations: ')
					for l in compound.locations.all():
						print(l.library_plate.library.name)

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
