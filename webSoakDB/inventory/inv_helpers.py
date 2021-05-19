from API.models import Library, LibraryPlate, SourceWell, Compounds, Preset
from .forms import LibraryForm, LibraryPlateForm, PlateUpdateForm
from webSoakDB_backend.validators import data_is_valid
from webSoakDB_backend.helpers import upload_plate
from datetime import date, datetime
from django.core.files.storage import FileSystemStorage
import string

def fake_compounds_copy(queryset):
	'''Produces regular Python objects (not django model objects) to store a copy of the data from
	SourceWell objects in quesryset. This copy is used to manipulate the data without chaging 
	anything in the database'''
	class Copy:
		def __init__(self, compound):
			self.well = compound.well         
			self.code = compound.compound.code
			self.active = compound.active
			self.deactivation_date = compound.deactivation_date
			self.changes = [{"date" : ch.date, "activation" : ch.activation} for ch in compound.status_changes.all()]
	
	return [Copy(c) for c in queryset]

def fake_preset_copy(preset, missing_compounds):
	'''Produces a regular Python object (not django model object) to store a copy of the data from
	the Preset object (preset). Then, it adds information about the availability of the compounds 
	in the current LibraryPlates, and alternative locations of the compound if the compound is missing
	in the current plate, based on missing_compounds'''
	
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
	
	preset_copy = PresetCopy(preset)
	
	for subset in preset.subsets.all():
		subset_copy = SubsetCopy(subset)
		preset_copy.subsets.append(subset_copy)
		current = LibraryPlate.objects.filter(library=subset.library.id, current=True)
		for plate in current:
			for compound in subset.compounds.all():
				c = CompoundCopy(compound)
				if (plate.id, c.id) in missing_compounds:
					c.available = False
					subset_copy.unavailable_count += 1
					c.alternatives = [( w.library_plate.library.name + ' : ' + w.library_plate.barcode ) for w in compound.locations.filter(library_plate__library__public = True) if w.active ]
				subset_copy.compounds.append(c)
		
	return preset_copy

def get_plate_size(queryset):
	rows = [char for char in string.ascii_uppercase] + ['AA', 'AB', 'AC', 'AD', 'AE', 'AF']
	columns = ['0' + str(i) for i in range(1, 10)] + [str(i) for i in range(10, 49)]

	large = False
	for i in columns[24:]:
		if queryset.filter(well__contains=str(i)).count() > 0:
			large = True
			
	for i in rows[16:]:
		if queryset.filter(well__contains=str(i)).count() > 0:
			large = True
	
	if not large:
		rows = rows[0:16]
		columns = columns[0:24]
	
	return {'rows' : rows, 'columns' : columns}
		
def get_change_dates(queryset):
		
	changes = []
	for compound in queryset:
		changes += compound.changes
		
	change_dates = set([ change["date"] for change in changes])
	return sorted(change_dates, reverse=True)

def set_status(compound, date):
	compound.active = True
	
	for change in sorted(compound.changes, key=lambda change: change["date"], reverse=True):
		#inspect changes in reverse chronological order
		if change["date"] <= date:
			compound.active = change["activation"]
			break
			
	return compound

def get_usage_stats(queryset):
	
	count = len(queryset)
	if count == 0:
		return {"count" : 0, "active" : 0, "inactive" : 0, "availability" : 0 }
	
	active, inactive = 0, 0
	
	for compound in queryset:
		if compound.active:
			active += 1
		else:
			inactive += 1
	
	availability = round(active / count * 100)
	
	return {"count" : count, "active" : active, "inactive" : inactive, "availability" : availability }
	
def current_library_selection(boolean):
	if boolean == True:
		first_tuple = [("", "Select library...")]
	else:
		first_tuple = []
		
	return  first_tuple + [(library.id, library.name) for library in Library.objects.filter(public=True)]
