from django.core.exceptions import ObjectDoesNotExist
from API.models import Library, LibraryPlate
from tools.data_storage_classes import SourceWellCopy, PresetCopy, SubsetCopy, CompoundCopy, BasicTemporaryCompound, SubsetCopyWithAvailability
import string
import csv
import re

def sw_copies(queryset):
	return [SourceWellCopy(c) for c in queryset]

def fake_preset_copy(preset, missing_compounds):
	'''Copy of a preset with alternative locations for compounds that are missing from the current plate(s) in a library'''

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
					c.alternatives = [( 
						w.library_plate.library.name + 
						' : ' + 
						w.library_plate.barcode ) 
						for w in 
						compound.locations.filter(library_plate__library__public = True) 
						if w.active 
						]
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

def get_subsets_with_availability(*args):
	subset_copies = []
	if len(args) == 1:
		for s in args[0]:
			subset_copies.append(SubsetCopyWithAvailability(s))
	elif len(args) == 2:
		subset_copies.append(SubsetCopyWithAvailability(args[0], args[1]))
	
	return subset_copies

def upload_temporary_subset(file_name):
	compounds = set()

	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		for row in compound_reader:
			compound = BasicTemporaryCompound(row[0], row[1])
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
