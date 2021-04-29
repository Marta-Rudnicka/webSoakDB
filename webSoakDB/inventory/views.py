from django.shortcuts import render
from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponseRedirect, HttpRequest, HttpResponse
from django.core.files.storage import FileSystemStorage
from django.contrib.admin.views.decorators import staff_member_required
from datetime import date, datetime
import string
import re
from .inv_helpers import fake_compounds_copy, get_plate_size, get_change_dates, get_usage_stats, fake_preset_copy, current_library_selection
from .dt import get_well_dictionary
from .forms import LibraryForm, LibraryPlateForm, PlateUpdateForm, NewPresetForm, EditPresetForm, DTMapForm
from webSoakDB_backend.validators import data_is_valid, selection_is_valid
from webSoakDB_backend.helpers import upload_plate, upload_subset
from API.models import Library, LibraryPlate, SourceWell, Compounds, Preset



fake_well_dictionary = {} #

#VIEWS DISPLAYING PAGES (HANDLING GET REQUESTS)
@staff_member_required
def index(request):
	return render(request, "inventory/inventory-index.html")

@staff_member_required
def dispense_testing(request):
	rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	columns = [ '0' + str(i) for i in range(1, 10)] + ['10', '11', '12']
	
	
	l = [('A', 'B', 'C', 'D'), ('E', 'F', 'G', 'H')]
	
	return render(request, "inventory/dispense-testing.html", {
		'columns': columns,
		'rows': rows,
		'l': l,
	})

@staff_member_required
def libraries(request):
	libraries = Library.objects.filter(public=True).order_by("name")
	library_form = LibraryForm()
	form_dict = {}
	for lib in libraries:
		form_dict[lib] = LibraryForm(initial={"name": lib.name, "for_industry": lib.for_industry})
		
	return render(request, "inventory/libraries.html", {"libraries": libraries, "library_form": library_form, "form_dict" : form_dict})
	
@staff_member_required
def plates(request):
	libraries = Library.objects.filter(public=True).order_by("name")
	plate_form = LibraryPlateForm(libs=current_library_selection(True))
	return render(request, "inventory/plates.html", {"libraries": libraries, "plate_form" : plate_form})

@staff_member_required
def presets(request):
	presets = Preset.objects.all().order_by("name")
	all_libs = current_library_selection(False)
	new_preset_form = NewPresetForm(libs=([("", "Select library...")] + all_libs))
	
	#produce a set of all compounds missing from the current library plates (only for libraries used in any of the presets)
	missing_compounds = set()
	
	for preset in presets:
		for subset in preset.subsets.all():
			current = LibraryPlate.objects.filter(library=subset.library.id, current=True)			
			for plate in current:
				m = plate.compounds.filter(active=False).all()
				missing_compounds = {(int(plate.id), int(c.compound.id)) for c in m} | missing_compounds
				
	# tip: the info from which plate the compound is missing is needed for libraries that share the same compounds,
	# which happens when there are different versons of the same library for different solvents. Otherwise, e.g., 
	# a compound that is missing from the current plate DSI DMSO would be recorded as missing also in the current
	# plates of DSI EG, even if it's still available in the EG plate
	
	#make copy of presets data, and add information about availability of the compounds in the current library plates
	presets_copy = [fake_preset_copy(preset, missing_compounds) for preset in presets]
	
	#produce a dictionary that matches each presets with appropriate form to edit this preset
	form_dict = {}
	
	for p in presets_copy:
		#generate lists of valid library choices for the preset
		old_libs_set = {(s.library_id, s.library_name) for s in p.subsets}
		new_libs = [("", "Select library...")] + list(set(all_libs) - old_libs_set)
		old_libs = [("", "Select library...")] + list(old_libs_set)
		
		form_dict[p] = EditPresetForm(old_libs=old_libs, new_libs=new_libs, initial={
			"name": p.name, 
			"description": p.description, 
			"new_library": "", 
			"new_compound_list": None, 
			"edited_library": "",
			"updated_compound_list": None,
			"delete_library": "",
			})
				
	return render(request, "inventory/presets.html", {
		"presets": presets_copy, 
		"new_preset_form": new_preset_form,
		"form_dict" : form_dict,
		})

@staff_member_required
def update_plate(request, pk):
	'''show information about compounds in plate, provides forms to edit plate data, delete a plate, and deactivate compounds'''
	plate = LibraryPlate.objects.get(pk=pk)
	compounds_count = plate.compounds.all().count()
	
	#for compounds that are missing, find in which library plates they are still available
	alternatives = {}
	
	for compound in plate.compounds.all():
		if not compound.active:
			c = Compounds.objects.get(code=compound.compound.code, smiles = compound.compound.smiles)
			alts = [( w.library_plate.library.name + ' : ' + w.library_plate.barcode ) for w in c.locations.filter(library_plate__library__public = True) if w.active ]
			alternatives[compound.id] = alts
	
	#generate info about availablity of compounds
	if compounds_count > 0:
		active_count = plate.compounds.filter(active=True).count()
		inactive_count = plate.compounds.filter(active=False).count()
		availability = round(active_count / compounds_count * 100)
	else:
		active_count, inactive_count, availability = 0, 0, 0
	
	#generate form to update the basic information about the plate
	libs = current_library_selection(True)
	plate_form = PlateUpdateForm(libs=libs, initial={"library" : plate.library.id, "barcode" : plate.barcode, "current" : plate.current})
	dt_map_form = DTMapForm()
		
	return render(request, "inventory/update_plate.html", {
		"plate": plate, 
		"compounds" : plate.compounds.all().order_by("deactivation_date"),
		"plate_form" : plate_form, 
		"active_count" : active_count, 
		"inactive_count" : inactive_count,
		"availability" : availability,
		"alternatives" : alternatives,
		"dt_map_form" : dt_map_form,
		})

@staff_member_required
def track_usage(request, pk, date, mode):
	'''Shows information about availability of compounds in the plate on a given date in three ways: availability statistics,
	a list of compounds on a given date with availability info, and a table graphically representing a physical plate with the
	locations of the active and inactive compounds. There are two modes: general view, which shows everything, and graphical view, 
	in which some of the page elements are hidden to bring the graphical representation of a plate to the top of the page. The 
	mode argument is passed to the template and used as id of the <main> HTML element. The difference between the page layouts 
	in different modes is determined by the CSS used for the children of the <main> element '''
	
	plate = LibraryPlate.objects.get(pk=pk)
	compounds = fake_compounds_copy(plate.compounds.all()) #make a copy of the data to edit it without touching the db
	timestamp = datetime.strptime(date, "%Y-%m-%d").date() #make a datetime object based on the url
	change_dates = get_change_dates(compounds, plate) #produce the list of all dates on which any compounds were deactivated
	
	#generate strings needed to switch between the general and the graphic view
	if mode == "general-view":
		switch_view = "graphic view"
		switch_view_url = "graphic-view"
	else:
		switch_view = "general view"
		switch_view_url = "general-view"
	
	#activate all the compounds that were still active on the inspected date
	for c in compounds:
		if c.deactivation_date and c.deactivation_date > timestamp:
			c.active = True
	
	#decide what size the table graphically representing the library plate should be
	size = get_plate_size(plate.compounds.all())
	rows = size['rows']
	columns = size['columns']
	
	#produce basic usage statistics
	usage_stats = get_usage_stats(compounds)
		
	return render(request, "inventory/track_usage.html", {
		"change_dates": change_dates,
		"compounds" : compounds,
		"plate": plate, 
		"date" : date,
		"timestamp" : timestamp,
		"active_count" : usage_stats['active'], 
		"inactive_count" : usage_stats['inactive'],
		"availability" : usage_stats['availability'],
		"columns" : columns,
		"rows": rows,
		"main_id" : mode,
		"switch_view": switch_view,
		"switch_view_url": switch_view_url,
		})

#FORM ACTION VIEWS (HANDLING POST REQUESTS)
@staff_member_required
def add_library(request):
	form = LibraryForm(request.POST)
	
	if request.method == "POST":	
		if form.is_valid():
			name = form.cleaned_data['name']
			for_industry = form.cleaned_data['for_industry']
			print('name: ', name, ", for_industry: ", for_industry)
			Library.objects.create(name=name, for_industry=for_industry, public=True)

			return HttpResponseRedirect('../libraries')

@staff_member_required
def add_plate(request):
	today = str(date.today())
	
	if request.method == "POST":
		form = LibraryPlateForm(data=request.POST, files=request.FILES, libs=current_library_selection(True))
		
		if form.is_valid():
			log = []
			fs = FileSystemStorage()
			source = request.FILES["plate_map"]
			filename = fs.save(source.name, source)
			if data_is_valid(filename, log):
				
				barcode = form.cleaned_data['barcode']
				library_id = form.cleaned_data['library']
				current = form.cleaned_data['current']
				library = Library.objects.get(pk=library_id)
				today = str(date.today())
				
				plate = LibraryPlate.objects.create(library = library, barcode = barcode, current = current, last_tested = today)
				
				upload_plate(filename, plate)
				fs.delete(filename)
				
				return HttpResponseRedirect('../plates')
			else:
				
				fs.delete(filename)
				
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})
		else:
			print('invalid form')
			print(form.errors)
			return HttpResponseRedirect('../plates')

@staff_member_required
def add_preset(request):
	if request.method == "POST":
		form = NewPresetForm(data=request.POST, files=request.FILES, libs=current_library_selection(True))
		if form.is_valid():
			print('FILES: ', request.FILES)
			log = []
			fs = FileSystemStorage()
			source = request.FILES["new_compound_list"]
			library_id = form.cleaned_data['new_library']
			filename = fs.save(source.name, source)
			print('library_id: ', library_id, 'description: ', form.cleaned_data['description'], 'preset_name:', form.cleaned_data['name'])
			if selection_is_valid(filename, log, library_id):
				library = Library.objects.get(pk=library_id)
				description = form.cleaned_data['description']
				preset_name = form.cleaned_data['name']
			
				#create new Subset object and upload data to it
				subset_name = library.name + " selection"
				origin = "preset: " + preset_name
				first_subset = upload_subset(filename, library_id, subset_name, origin)
				print('create subset: ', library_id, subset_name, origin)
				fs.delete(filename)
								
				#create preset object
				preset = Preset.objects.create(name=preset_name, description=description)
				preset.subsets.add(first_subset)
				preset.save()
				print('create preset: ', preset_name, description)
				
				return HttpResponseRedirect('../presets/')
			else:
				fs.delete(filename)
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})
		else:
			print(form.errors)
			return render(request, 'inventory/presets.html', {'form_errors' : form.errors})

@staff_member_required
def edit_preset(request):
	if request.method == "POST":
	
		#create lists of valid library choices
		preset = Preset.objects.get(pk=request.POST["id"])
		all_libs = current_library_selection(False)
		old_libs_set = {(s.library.id, s.library.name) for s in preset.subsets.all()}
		new_libs = [("", "Select library...")] + list(set(all_libs) - old_libs_set)
		old_libs = [("", "Select library...")] + list(old_libs_set)
		
		#create form object
		form = EditPresetForm(data=request.POST, files=request.FILES, old_libs=old_libs, new_libs=new_libs)
		#https://django-gotchas.readthedocs.io/en/latest/forms.html <-- explanation for the weird kwargs; 
		
		if form.is_valid():
			
			#check if upload files are valid
			new_library = form.cleaned_data["new_library"]
			new_compound_list = request.FILES.get("new_compound_list")
			edited_library = form.cleaned_data["edited_library"]
			deleted_library = form.cleaned_data["deleted_library"]
			updated_compound_list = request.FILES.get("updated_compound_list")
			fs_new, fs_edit = None, None
			
			if new_compound_list:
				fs_new =  FileSystemStorage()
				filename_new = fs_new.save(new_compound_list.name, new_compound_list) 
			if updated_compound_list:
				fs_edit =  FileSystemStorage()
				filename_edit = fs_edit.save(updated_compound_list.name, updated_compound_list)
				
			log = []
			if ( new_compound_list and not selection_is_valid(filename_new, log, new_library)) or (updated_compound_list and not selection_is_valid(filename_edit, log, edited_library)):
				if fs_new:
					fs_new.delete(filename_new)
				if fs_edit:
					fs_edit.delete(filename_edit)
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})
			
			#update submitted data
			preset.name = form.cleaned_data["name"]
			preset.description = form.cleaned_data["description"]
			
			if new_compound_list:
				
				#create new Subset object and upload data to it
				library = Library.objects.get(pk=new_library)
				subset_name = library.name + " selection"
				origin = "preset: " + preset.name
				new_subset = upload_subset(filename_new, library.id, subset_name, origin)
				fs_new.delete(filename_new)
				
				#add new subset to preset
				preset.subsets.add(new_subset)
			
			if updated_compound_list:
				
				#find and remove old list
				library = Library.objects.get(pk=edited_library)
				old_subset = preset.subsets.get(library = library)
				preset.subsets.remove(old_subset)	
				
				#create new Subset object and upload data to it
				library = Library.objects.get(pk=edited_library)
				subset_name = library.name + " selection"
				origin = "preset: " + preset.name
				new_subset = upload_subset(filename_edit, library.id, subset_name, origin)
				fs_edit.delete(filename_edit)
				
				#add new subset to preset
				preset.subsets.add(new_subset)
				
			if deleted_library:
				library = Library.objects.get(pk=deleted_library)
				deleted_subset = preset.subsets.get(library = library)
				preset.subsets.remove(deleted_subset)	
			
			preset.save()
			return HttpResponseRedirect('../presets/')
		else:
			return render(request, "webSoakDB_backend/error_log.html", {"form_errors": form.errors, "non_field_errors": form.non_field_errors})			

@staff_member_required
def edit_library(request):
	form = LibraryForm(request.POST)
	
	if request.method == "POST":	
		if form.is_valid():
			pk = request.POST.get('id')
			name = form.cleaned_data['name']
			for_industry = form.cleaned_data['for_industry']
			lib = Library.objects.get(pk=pk)
			lib.name=name
			lib.for_industry=for_industry
			lib.save()

			return HttpResponseRedirect('../libraries')
	else:
		return render(request, "webSoakDB_backend/error_log.html", {"form_errors": form.errors, "non_field_errors": form.non_field_errors})	

@staff_member_required
def edit_plate(request):
	form = PlateUpdateForm(data=request.POST, libs=current_library_selection(True))
	
	if request.method == "POST":	
		if form.is_valid():
			pk = request.POST.get('id')
			library_id = form.cleaned_data['library']
			library = Library.objects.get(pk=library_id)
			barcode = form.cleaned_data['barcode']
			current = form.cleaned_data['current']
			print("pk: ", pk, ", library: ", library, ", barcode : ", barcode, ", current: ", current)
			plate = LibraryPlate.objects.get(pk=pk)
			plate.library = library
			plate.barcode = barcode
			plate.current = current
			plate.save()
			
			redirect_url = '/inventory/update-plate/' + str(plate.id) + '/'
			return HttpResponseRedirect(redirect_url)
		else:
			return render(request, "webSoakDB_backend/error_log.html", {"form_errors": form.errors, "non_field_errors": form.non_field_errors})	

def update_preset(request):
	pass

@staff_member_required
def delete_library(request):
	if request.method == "POST":	
		pk = request.POST.get('id')
		lib = Library.objects.get(pk=pk)
		if len(lib.plates.all()) == 0:
			print('Good to delete')
			lib.delete()
			return HttpResponseRedirect('../libraries/')
		else:
			print('Cannot be deleted')
			return HttpResponseRedirect('../library-deletion-error/')

@staff_member_required
def delete_plate(request):
	if request.method == "POST":	
		pk = request.POST.get('id')
		plate = LibraryPlate.objects.get(pk=pk)
		plate.delete()
		return HttpResponseRedirect('../plates/')

@staff_member_required
def delete_preset(request):
	if request.method == "POST":	
		pk = request.POST.get('id')
		preset = Preset.objects.get(pk=pk)
		preset.delete()
		return HttpResponseRedirect('../presets/')

def library_deletion_error(request):
	return render(request, "inventory/library_deletion_error.html");

def dummy(request):
	return render(request, "webSoakDB_backend/dummy.html");

@staff_member_required
def deactivate_compounds(request):
	if request.method == "POST":
		today = str(date.today())
		plate = LibraryPlate.objects.get(pk=request.POST.get('plate_id'))
		redirect_url = '/inventory/update-plate/' + str(plate.id) + '/'
		plate.last_tested = today
		plate.save()
		for key in request.POST:
			if key != 'csrfmiddlewaretoken' and key != "plate_id" and key != "":
				compound = SourceWell.objects.get(pk=key)
				print(compound)
				compound.active = False
				compound.deactivation_date = today
				compound.save()
			
		return HttpResponseRedirect(redirect_url)

@staff_member_required
def dispense_testing_map(request):
	rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	columns = [ '0' + str(i) for i in range(1, 10)] + ['10', '11', '12']
	errors = []
	
	form = DTMapForm(request.POST, request.FILES)
	
	if request.method == "POST":
		if form.is_valid():
			source = request.FILES["dt_map"]
			if not re.fullmatch('(.*)+\.csv$', source.name):
				return render(request, "inventory/dispense-testing.html", {"errors": ["Invalid file: the well map must be a csv file."], "filename" : source.name})
			
			pk = request.POST.get('id')
			fs = FileSystemStorage()
			source = request.FILES["dt_map"]
			filename = fs.save(source.name, source)
			dt_dict = get_well_dictionary(filename)
			fs.delete(filename)
			plate=LibraryPlate.objects.get(pk=pk)
			
			if not isinstance(dt_dict, dict):
				return render(request, "inventory/dispense-testing.html", {"errors": dt_dict, "plate": plate, "filename" : source.name})
				
			compounds = plate.compounds.all()
			data = []
			for key in dt_dict:
				try:
					dw = key
					sw = dt_dict[key]
					c =  compounds.get(well=sw)
					active = c.active
					code = c.compound.code
					_id = c.id
					data.append((dw, sw, active, code, _id))
				except(ObjectDoesNotExist):
					errors.append("Invalid mapping: this source plate has no compound recorded in well " + sw)
			
			if errors != []:
				data = None
				rows = None
				columns = None
					
			return render(request, "inventory/dispense-testing.html", {"errors": errors, "data" : data, "rows": rows, "columns": columns, "plate": plate, "filename" : filename})
		else:
			return render(request, "webSoakDB_backend/error_log.html", {"form_errors": form.errors, "non_field_errors": form.non_field_errors})	
