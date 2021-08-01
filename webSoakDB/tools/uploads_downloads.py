'''helper functions for uploading and downloading files'''
import django.core.exceptions
from API.models import Compounds, SourceWell, Library, LibrarySubset, LibraryPlate, Project
from .data_storage_classes import BasicTemporaryCompound
import csv
import re
from datetime import date 
from django.http import HttpResponse

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Draw
from rdkit.Chem.SaltRemover import SaltRemover

def standardize_smiles(raw_string):
	'''Canonicalize SMILES and remove salts - use only on validated SMILES strings;
	For parsing CSV files'''

	#remove whitespace and convert to canonical form
	smiles = Chem.CanonSmiles(raw_string.strip())

	#convert to RDKit mol object 
	mol = Chem.MolFromSmiles(smiles)

	#remove salt
	remover = SaltRemover()
	desalted_mol = remover.StripMol(mol, dontRemoveEverything=True)

	#get SMILES string of the desalted mol and return it
	return Chem.MolToSmiles(desalted_mol)




#LOCATING COMPOUND LIST IN INVENTORY / EXPORTING LOCATION DATA (INVENTORY)
def upload_temporary_subset(file_name):
	compounds = set()

	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		dialect.delimiter = ','
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		for row in compound_reader:
			compound = BasicTemporaryCompound(row[0].strip(), standardize_smiles(row[1]))
			compounds.add(compound)
			
	return compounds

def parse_compound_list(string):
	compounds = set()
	for item in string.split(','):
		data = [x for x in item.split(':')]
		try:
			compounds.add(BasicTemporaryCompound(data[0], data[1]))
		except IndexError:
			pass
	
	return compounds

def parse_id_list(string):
	return [int(i) for i in re.findall('[0-9]+', string)]

#IMPORTING A PLATE MAP (FRONTEND & INVENTORY)

#import data from a csv file - should be used after validation with data_is_valid
def upload_plate(file_name, plate):
	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		dialect.delimiter = ','
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		for row in compound_reader:
			smiles = standardize_smiles(row[2])
			try:
				compound = Compounds.objects.get(code = row[0].strip(), smiles = smiles)	
			except django.core.exceptions.ObjectDoesNotExist:
				compound = Compounds.objects.create(smiles = smiles, code = row[0].strip())
				if compound.smiles:
					add_molecular_properties(compound)
					
			s_well = SourceWell.objects.create(compound = compound, library_plate = plate, well = row[1].strip())
			
			try:
				s_well.concentration = int(row[3])
				s_well.save()
			except (IndexError, ValueError): #ignore if concentration is not provided
				pass

def add_molecular_properties(compound):
	sanitized_mol = Chem.MolFromSmiles(compound.smiles)
	compound.log_p = Crippen.MolLogP(sanitized_mol)
	compound.mol_wt = float(Chem.rdMolDescriptors.CalcExactMolWt(sanitized_mol))
	compound.heavy_atom_count = Chem.Lipinski.HeavyAtomCount(sanitized_mol)
	compound.heavy_atom_mol_wt = float(Descriptors.HeavyAtomMolWt(sanitized_mol))
	compound.nhoh_count = Chem.Lipinski.NHOHCount(sanitized_mol)		
	compound.no_count = Chem.Lipinski.NOCount(sanitized_mol)
	compound.num_h_acceptors = Chem.Lipinski.NumHAcceptors(sanitized_mol)
	compound.num_h_donors = Chem.Lipinski.NumHDonors(sanitized_mol)
	compound.num_het_atoms = Chem.Lipinski.NumHeteroatoms(sanitized_mol)
	compound.num_rot_bonds = Chem.Lipinski.NumRotatableBonds(sanitized_mol)
	compound.num_val_electrons = Descriptors.NumValenceElectrons(sanitized_mol)
	compound.ring_count = Chem.Lipinski.RingCount(sanitized_mol)
	compound.tpsa = Chem.rdMolDescriptors.CalcTPSA(sanitized_mol)
	compound.save()

def upload_subset(file_name, library_id, name, origin):
	'''assumes list have been validated before'''
	library = Library.objects.get(id=library_id)
	new_subset = LibrarySubset.objects.create(library=library, name=name, origin=origin)
	
	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		dialect.delimiter = ','
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		for row in compound_reader:
			try:
				#compound = Compounds.objects.get(code = row[0].strip(), smiles=row[1].strip())

				#find all compounds with matching SMILES string
				matching_smiles = Compounds.objects.filter(smiles=standardize_smiles(row[0]))

				#if you find only one, use it
				if matching_smiles.count() == 1:
					compound = matching_smiles.all()[0]
				
				#if there are more, take the first one you find that belongs to <library>
				else:
					for c in matching_smiles.all():
						if library.id in [sw.library_plate.library.id for sw in c.locations.all()]:
							compound = c
							break

				new_subset.compounds.add(compound)
				new_subset.save()
			except django.core.exceptions.ObjectDoesNotExist:
				print('Error: No such compound found.')
				print('Error: upload_subset running on unvalidated data. Use validators.selection_is_valid() on the input data first!')
			except django.core.exceptions.MultipleObjectsReturned:
				print('Error: There are duplicate compounds in the database: ', row[0], ':', row[1])
				print('Error: upload_subset running on unvalidated data. Use validators.selection_is_valid() on the input data first!')
	
	return new_subset

#EXPORTING SOAK-DB COMPATIBLE COMPOUND LISTS (AS CSV)
def import_full_libraries(proposal):
	print("import_full_libraries", proposal)
	proposal = Project.objects.get(id=proposal)
	full_plates = []
	
	for l in proposal.libraries.all():
		current_plates = LibraryPlate.objects.filter(library=l, current=True)
		[full_plates.append(plate) for plate in current_plates]
	
	s_wells = []
	for p in full_plates:
		s_wells = s_wells + [c for c in p.compounds.filter(active=True)]
	return s_wells

def find_source_wells(subset, plate_ids):
	plates = []
	for id in plate_ids:
		plates.append(LibraryPlate.objects.get(pk=id))

	sw = []
	if type(subset) == set:
		compounds = subset
	else:
		compounds = subset.compounds.all()
	
	for compound in compounds:
		for plate in plates:
			try:
				c = plate.compounds.get(compound__code = compound.code, compound__smiles = compound.smiles, active=True)
				sw.append(c)
				break 		#to avoid duplicates in plate combinations
			except django.core.exceptions.ObjectDoesNotExist:
				pass
	sw.sort(key=lambda x: x.library_plate.id)

	return sw
		
def import_library_parts(proposal, data):
	proposal = Project.objects.get(id=proposal)
	sw = []
	
	for s in proposal.subsets.all():
		if data.get(str(s.library.id), False):
			sw = sw + find_source_wells(s, [ data.get(str(s.library.id))])
		else:							#multiple current plates
			lib = Library.objects.get(pk=s.library.id)
			plates_count = lib.plates.filter(current=True).count()
			for i in range(1, plates_count + 1):
				option_name = str(s.library.id) + '-' + str(i)
				plate_id = data.get(option_name, False)
				sw = sw + find_source_wells(s, [plate_id])
	return sw

def source_wells_to_csv(source_wells, file_path, filename_prefix):
	with open(file_path, 'r+') as f:
		f.truncate(0)
		for c in source_wells:
			line = c.library_plate.barcode + ',' + c.well + ',' + c.library_plate.library.name + ',' + c.compound.smiles + ',' + c.compound.code + "\n"
			f.write(line)
		f.close() #to ensure the whole file will get served
	
	with open(file_path, 'r+') as f:
		filename = filename_prefix + '-soakDB-source-export-' + str(date.today()) + '.csv'
		response = HttpResponse(f, content_type='text/csv')
		response['Content-Disposition'] = "attachment; filename=%s" % filename
	
	return response

#CREATE CSV FILES WITH COMPOUND LISTS (for plates, subsets or preset)
def basic_csv_headers():
	return "Compound Code, Well, SMILES, Concentration"

def extended_csv_headers():
	return ", Molecular weight (Da), \
TPSA, \
LogP, \
Number of valence electrons, \
Number of hydrogen bond acceptors, \
Number of hydrogen bond donors, \
Number of heteroatoms, \
Number of rotable bonds, \
Number of rings , \
Number of heavy atoms , \
Mol. weight of heavy atoms (Da), \
Number of NH and OH , \
Number of nitrogens and oxygens"

def append_mol_properties_to_csv(compound):
	return ', ' + str(compound.mol_wt) + ', ' \
		+ str(compound.tpsa) + ', ' \
		+ str(compound.log_p) + ', '\
		+ str(compound.num_val_electrons) + ', '\
		+ str(compound.num_h_acceptors) + ', '\
		+ str(compound.num_h_donors) + ', '\
		+ str(compound.num_het_atoms) + ', '\
		+ str(compound.num_rot_bonds) + ', '\
		+ str(compound.ring_count) + ', '\
		+ str(compound.heavy_atom_count) + ', '\
		+ str(compound.heavy_atom_mol_wt) + ', '\
		+ str(compound.nhoh_count) + ', '\
		+ str(compound.no_count)

def serve_csv_compound_list(header, collection, filename, include_details):
	with open('files/plate-map.csv', 'r+') as f:
		f.truncate(0)
		f.write(header + "\n")

		#add appropriate data depending on the type of compound collection
		if "LibraryPlate" in str(type(collection)):
			f = add_plate_map_data_to_file(collection, f, include_details)
		if "LibrarySubset" in str(type(collection)):
			f = add_subset_data_to_file(collection, f, include_details, False)
		if "Preset" in str(type(collection)):
			for subset in collection.subsets.all():
				f = add_subset_data_to_file(subset, f, include_details, True)
		f.close()
	
	with open('files/plate-map.csv', 'r+') as f:
		response = HttpResponse(f, content_type='text/csv')
		response['Content-Disposition'] = "attachment; filename=%s" % filename
		return response

def add_plate_map_data_to_file(plate, file, include_details):
	for compound in plate.compounds.filter(active=True):
		line = compound.compound.code + ',' + compound.well + ',' + compound.compound.smiles + ',' + str(compound.concentration)
		if include_details:
			line += append_mol_properties_to_csv(compound.compound)
		line += "\n"
		file.write(line)
	return file

def add_subset_data_to_file(subset, file, include_details, is_preset):
	if is_preset:
		lib_name = subset.library.name + ','
	else:
		lib_name = ""

	for compound in subset.compounds.all():
		line = lib_name + compound.code + ',' + compound.smiles
		if include_details:
			line = line + append_mol_properties_to_csv(compound)
		line += "\n"
		file.write(line)
	return file