from datetime import datetime
from django.core.exceptions import ObjectDoesNotExist
from API.models import Library, LibraryPlate
from tools.data_storage_classes import SourceWellCopy, PresetCopy, SubsetCopyWithAvailability
import string

def sw_copies(queryset):
	return [SourceWellCopy(c) for c in queryset]

def fake_preset_copy(preset):
	preset_copy = PresetCopy(preset)
	for subset in preset.subsets.all():
		subset_copy = SubsetCopyWithAvailability(subset)
		preset_copy.subsets.append(subset_copy)
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
	
	return []
	#return  first_tuple + [(library.id, library.name) for library in Library.objects.filter(public=True)]

def get_subsets_with_availability(*args):
	subset_copies = []
	if len(args) == 1:
		for s in args[0]:
			subset_copies.append(SubsetCopyWithAvailability(s))
	elif len(args) == 2:
		subset_copies.append(SubsetCopyWithAvailability(args[0], args[1]))
	
	return subset_copies

def parse_fedids(string):
	return [s.strip() for s in string.split(',')]
	#TODO: validate fedids
	