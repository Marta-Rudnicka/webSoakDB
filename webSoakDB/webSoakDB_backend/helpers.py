from API.models import Library, LibraryPlate, Compounds, SourceWell, Proposals, LibrarySubset
import django.db
import csv
import re
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
				print('Error: No such compound found.')
				print('Error: upload_subset running on unvalidated data. Use validators.selection_is_valid() on the input data first!')
			except django.core.exceptions.MultipleObjectsReturned:
				print('Error: There are duplicate compounds in the database: ', row[0], ':', row[1])
				print('Error: upload_subset running on unvalidated data. Use validators.selection_is_valid() on the input data first!')
	
	return new_subset
	

def import_full_libraries(proposal):
	proposal = Proposals.objects.get(name=proposal)
	#proposal = Proposals.objects.get(proposal=proposal)
	full_plates = []
	
	for l in proposal.libraries.all():
		current_plates = LibraryPlate.objects.filter(library=l, current=True)
		[full_plates.append(plate) for plate in current_plates]
	
	#print('full_plates in import_full_libraries: ', full_plates)
	
	s_wells = []
	for p in full_plates:
		#print('p.compounds.all()', p.compounds.all())
		s_wells = s_wells + [c for c in p.compounds.all()]
	return s_wells

def find_source_wells(subset, plate_id):
	plate = LibraryPlate.objects.get(pk=plate_id)
	sw = []
	for compound in subset.compounds.all():
		try:
			c = plate.compounds.get(compound = compound)
			sw.append(c)
		except django.core.exceptions.ObjectDoesNotExist:
			pass
	
	return sw
		

def import_library_parts(proposal, data):
	proposal = Proposals.objects.get(name=proposal)
	sw = []
	
	for s in proposal.subsets.all():
		if data.get(str(s.library.id), False):
			sw = sw + find_source_wells(s, data.get(str(s.library.id), False))
		else:							#multiple current plates
			lib = Library.objects.get(pk=s.library.id)
			plates_count = lib.plates.filter(current=True).count()
			for i in range(1, plates_count + 1):
				option_name = str(s.library.id) + '-' + str(i)
				plate_id = data.get(option_name, False)
				sw = sw + find_source_wells(s, plate_id)
	return sw

def export_form_is_valid(post_data):
	proposal = None
	subset_lib_ids = []
	for key in post_data:
		if key=='csrfmiddlewaretoken':
			pass
		elif key=='proposal':
			try:
				proposal = Proposals.objects.get(name=post_data.get(key))
				#proposal = Proposals.objects.get(proposal=post_data.get(key))
				subset_lib_ids = [s.library.id for s in proposal.subsets.all()]
			except django.core.exceptions.ObjectDoesNotExist:
				print('No proposal found')
				return False
		else:
			if re.fullmatch('[0-9]+', key):
				value = post_data.get(str(key), False)
			elif re.fullmatch('([0-9]+)(\-[0-9]+)', key):
				value = post_data.get(str(key), False)
				old_key = re.fullmatch('([0-9]+)(\-[0-9]+)', key)
				key = old_key.group(1)
			else:
				print('key not matching any regex:', key)
				return False
				 
			if not int(key) in subset_lib_ids:
				print('key not found in subset libraries: ', key)
				return False
				
			try:
				plate = LibraryPlate.objects.get(pk=value)
			except django.core.exceptions.ObjectDoesNotExist:
				print('plate does not exist: ', value)
				return False
			
			if plate.library.id != int(key):
				print('plate library not matching subset library')
				return False
	return True
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
