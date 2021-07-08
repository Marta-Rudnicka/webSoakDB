from ast import iter_fields
from datetime import datetime
from API.models import Compounds, SourceWell
from django.core.exceptions import ObjectDoesNotExist
import itertools
'''
Classes used by backend views to store model data so they can be modified without affecting the database or associated with additional information
'''

class SourceWellCopy:
	def __init__(self, compound):
		self.well = compound.well         
		self.code = compound.compound.code
		self.active = compound.active
		self.deactivation_date = compound.deactivation_date
		self.changes = [{"date" : ch.date, "activation" : ch.activation} for ch in compound.status_changes.all()]

class BasicTemporaryCompound:
	def __init__(self, code, smiles):
		self.code = code
		self.smiles = smiles
		try:
			self.id = Compounds.objects.get(code=code, smiles=smiles).id
		except ObjectDoesNotExist:
			self.id = 0

class PresetCopy:
	def __init__(self, preset):
		self.id = preset.id         
		self.name = preset.name
		self.description = preset.description
		self.subsets = []
	def __str__(self):
		return "Preset copy: " + self.name
	
class SubsetCopy:
	def __init__(self, subset):
		self.id = subset.id         
		self.name = subset.name
		self.library_name = subset.library.name
		self.library_id = subset.library.id
		self.compounds = []
		self.unavailable_count = 0
		
	def __str__(self):
		return "Subset copy: " + self.name + ":" + self.library_name

class CompoundCopy:
	def __init__(self, compound):
		self.id = compound.id         
		self.code = compound.code
		self.smiles = compound.smiles
		self.available = True
	
	def __str__(self):
		return "Compound copy: " + str(self.id) + ":" + self.code

class PlateCopy:
	def __init__(self, plate, missing_compounds):
		self.id = plate.id
		self.library = plate.library
		self.barcode = plate.barcode
		self.current = plate.current
		self.missing_compounds = missing_compounds
	def __str__(self):
		return "Plate copy: " + self.library.name + " " + self.library.barcode

class SubsetCopyWithAvailability:
	''' a library subset with added information about 10 library plates that have the largest number of the
	desired compounds available; if a library is large and takes up n plates, the availability data include combinations
	of up to n plates'''

	def __init__(self, *args):
		#for LbrarySubset objects saved in the database (which have library attribute available)
		if len(args) == 1:
			self.id = args[0].id
			self.name = args[0].name
			self.library = args[0].library
			self.library_id = args[0].library.id
			self.library_name = args[0].library.name
			self.compounds = set(args[0].compounds.all())
			self.availability = self.get_compound_availability(self.compounds, self.library)

		#for unsaved cherry-picking list uploaded from a file (where library data has to be manually entered)
		elif len(args) == 2:
			self.id = 0
			self.name = "temporary selection"
			self.library = args[1]
			self.library_id = args[1].id
			self.library_name = args[1].name
			self.compounds = args[0]
			self.availability = self.get_compound_availability(self.compounds, self.library)

	def get_compound_availability(self, desired_compounds, library):
		t1 = datetime.now()
		plate_count = self.plates_in_library()
		t2 = datetime.now()
		print(t2-t1, 'counted plates')
		single_plates = self.rank_single_plates(desired_compounds, library)
		t3 = datetime.now()
		print(t3-t2, 'got single plates')
		if plate_count in [0, 1]:
			return single_plates[:9]
		else:
			if single_plates[0].availability == 100:
				return single_plates[:9]
			
			else:
				return self.rank_combinations(single_plates, desired_compounds, plate_count)

	def plates_in_library(self):
		'''find out how many plates are needed to store one instance of the library'''
		return self.library.plates.filter(current=True).count()

	def rank_single_plates(self, desired_compounds, library):
		plate_copies = []

		for plate in library.plates.all():
			plate_copy = self.get_plate_copy(desired_compounds, plate)

			if plate_copy:
				plate_copies.append(plate_copy)

		plate_copies = self.sort_by_availability(plate_copies)
		return plate_copies

	def get_plate_copy(self, desired_compounds, plate):
		missing_compounds = self.find_missing_compounds(desired_compounds, plate)
		availability = ((len(desired_compounds) - len(missing_compounds)) / len(desired_compounds)) * 100
		if availability > 0:
			plate_copy = PlateCopy(plate, missing_compounds)
			plate_copy.availability = round(availability, 2)
			return plate_copy
		else:				#don't include plates that have none of the desired compounds
			return None	

	def find_missing_compounds(self, desired_compounds, plate):
		missing_in_plate = set()
		for c in desired_compounds:
			try:
				plate.compounds.get(compound__code = c.code, compound__smiles = c.smiles, active=True)
			except ObjectDoesNotExist:
				missing_in_plate.add(c)
			except AttributeError: #combined plate
				found = False
				for q in plate.compounds:
					try: 
						q.get(compound__code = c.code, compound__smiles = c.smiles, active=True)
						found = True
						continue
					except ObjectDoesNotExist:
						pass
				if not found:
					missing_in_plate.add(c)	
		return missing_in_plate

	def rank_combinations(self, single_plates, desired_compounds, count):
		'''Rank combinations of plates by compound availability. Start with combinations
		of 2 plates, and then 3, etc. until you reach the number of plates into which the whole 
		library fits. If 100% of desired compounds coverage is found, don't check for higher 
		number of plates (e.g. if a combination of 2 plates has all the desired compounds, don't check
		3-plate combinations at all)'''

		best_10 = single_plates[:]
		found_100 = False

		for i in range(2, count + 1):

			if found_100:
				break
			
			combinations = itertools.combinations(single_plates, i)
			for comb in combinations:
				combined_plate = CombinedPlate(comb)
				plate_copy = self.get_plate_copy(desired_compounds, combined_plate)
				if plate_copy.availability > best_10[-1].availability:
					best_10.append(plate_copy)
					best_10 = sorted(best_10, key=lambda x: x.availability, reverse=True)
					best_10.pop()
				if plate_copy.availability == 100:
					found_100 = True
			i += 1
		
		best_10 = self.sort_by_availability(best_10)

		return best_10

	def sort_by_availability(self, plate_copies):
		'''Sort PlateCopy objects by availability from highest to lowest; if there are multiple plates with the best availability, 
		avoid placing the current plate at the beginning of the list '''
		plate_copies.sort(key=lambda x: x.availability, reverse=True)

		if plate_copies[0].current:
			#try to find an old plate with the same availability and swap them
			for i in range(1, len(plate_copies)):
				if plate_copies[0].availability == plate_copies[i].availability and not plate_copies[i].current:
					plate_copies[0], plate_copies[i] = plate_copies[i], plate_copies[0]
					break
				else:
					break

		return plate_copies

class CombinedPlate:
	def __init__(self, plate_copies):
		compounds = []
		barcodes = ''
		ids = []
		library = plate_copies[0].library
		current = False

		for plate in plate_copies:
			#create a queryset with all relevant source wells (so they can be later processed with QuerySet API)
			compounds.append(SourceWell.objects.filter(library_plate__id = plate.id))
			barcodes = barcodes + plate.barcode + ',' 
			ids.append(plate.id)
			if plate.current:
				current = True
	
		self.id = ids
		self.library = library
		self.barcode = barcodes
		self.current = current
		self.compounds = compounds
