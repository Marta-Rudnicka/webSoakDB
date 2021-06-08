from API.models import Library, LibraryPlate, Compounds, SourceWell, Proposals, LibrarySubset
import django.db
import csv
import re
import django.core.exceptions
from datetime import date
from django.http import HttpResponse

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import Draw

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
				compound = Compounds.objects.create(smiles = row[2].strip(), code = row[0].strip())
				if compound.smiles:
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
	proposal = Proposals.objects.get(proposal=proposal)
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
		print(compound, compound.code, compound.smiles)
		for plate in plates:
			try:
				c = plate.compounds.get(compound__code = compound.code, compound__smiles = compound.smiles, active=True)
				sw.append(c)
				break 		#to avoid duplicates in plate combinations
			except django.core.exceptions.ObjectDoesNotExist:
				print('not found')
				pass
	sw.sort(key=lambda x: x.library_plate.id)

	return sw
		
def import_library_parts(proposal, data):
	proposal = Proposals.objects.get(proposal=proposal)
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

def export_form_is_valid(post_data):
	proposal = None
	subset_lib_ids = []
	for key in post_data:
		if key=='csrfmiddlewaretoken':
			pass
		elif key=='proposal':
			try:
				proposal = Proposals.objects.get(proposal=post_data.get(key))
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