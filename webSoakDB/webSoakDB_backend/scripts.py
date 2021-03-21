from API.models import Library, LibraryPlate, Compounds, SourceWell, Proposals, SoakDBCompound, Crystal, CrystalPlate
import django.core.exceptions

def import_compounds(proposal_name):
	'''create SoakDBCompound objects based on compound selection for the proposal
	and the current state of the inventory. For selected libraries, make SoakDBCompound 
	objects based on all the compounds in the current plates of the library. For 
	selected subsets, find all desired compounds in the current plates, and make
	a list of all the compounds currently unavailable. Returns the list of missing 
	compounds'''
	
	proposal = Proposals.objects.get(name=proposal_name)
	print('proposal: ', proposal)
	
	full_plates = []
	compounds = []
	
	#get current plates of selected full libraries
	for lib in proposal.libraries.all():
		#full_plates.append(lib.plates.all().filter(current=True))
		full_plates = (*full_plates, *lib.plates.all().filter(current=True))
	
	print('full_plates: ', full_plates)
	
	#create SoakDBCompound objects from all compounds in current plate	
	for plate in full_plates:
		for compound in plate.compounds.all():
			print('Creating compound: ')
			print('proposal: ', proposal, ', library_name: ', plate.library.name, ', library_plate: ', plate.name, 
				'well: ', compound.well, ' code:', compound.compound.code, ' smiles: ', compound.compound.smiles, ', crystal: ', None)
			#compound = SoakDBCompound.create(proposal=proposal, library_name=plate.library.name, library_plate=plate.name, 
			#								well=compound.well, code=compound.compound.code, smiles=compound.compound.smiles,
			#								crystal=None)
			#compound.save()
			#compounds.append(compound)
	
	subset_dict = ensure_unique_compounds(proposal)
	
	missing_compounds = []
	
	#get current plates of libraries for which subsets are selected
	for lib in subset_dict:
		plates =  LibraryPlate.objects.filter(library = lib, current=True)
		
		#find desired compounds in the plate, make list of missing items
		for compound in subset_dict[lib]:
			available = False
			c = Compounds.objects.get(code=compound.code, smiles=compound.smiles)
			for plate in plates:
				try:
					found = plate.compounds.get(compound = c)
					print('Creating compound: ')
					print('proposal: ', proposal, ', library_name: ', plate.library.name, ', library_plate: ', plate.name, 
							'well: ', found.well, ' code:', found.compound.code, ' smiles: ', found.compound.smiles, ', crystal: ', None)
					
					#new = SoakDBCompound.create(proposal=proposal, library_name=plate.library.name, library_plate=plate.name, 
					#						well=found.well, code=found.compound.code, smiles=found.compound.smiles,
					#						crystal=None)
					#new.save()
					#compounds.append(new)
					available = True
					break
				except django.core.exceptions.ObjectDoesNotExist:
					pass
			
			if available == False:
				missing_compounds.append(compound)
	
	return missing_compounds


def ensure_unique_compounds(proposal):
	subset_dict = {}
	for subset in proposal.subsets.all():
		if not subset.library in proposal.libraries.all():
			if subset.library not in subset_dict:
				subset_dict[subset.library] = set(subset.compounds.all())
			else:
				subset_dict[subset.library] = subset_dict[subset.library].update(subset.compounds.all())
	
	return subset_dict


def import_crystals(visit, file_name, plate_name, plate_type, drop_volume):
	plate = CrystalPlate.objects.create(name=plate_name, plate_type=plate_type, drop_volume=drop_volume)
	
	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		crystal_reader = csv.reader(csvfile, dialect)
		for row in crystal_reader:
			try:
				crystal = Crystal.objects.create(visit=visit, crystal_plate=plate, well=row[0], echo_x=row[1], echo_y=row[2], score=row[3])
				crystal.save()
			except (IndexError): 					#no crystal found
				pass
