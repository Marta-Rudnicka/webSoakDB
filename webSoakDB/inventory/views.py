from tools.data_storage_classes import SubsetCopyWithAvailability
from tools.compounds import standardize_smiles
from django.shortcuts import redirect, render
from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponseRedirect
from django.core.files.storage import FileSystemStorage
from django.contrib.admin.views.decorators import staff_member_required
from webSoakDB_stack.settings import MEDIA_ROOT
from datetime import date, datetime
import re
import shutil
from .inv_helpers import (
	sw_copies, 
	get_plate_size, 
	get_change_dates, 
	get_subsets_with_availability, 
	get_usage_stats, 
	fake_preset_copy, 
	current_library_selection,
	set_status,
	get_project_by_proposal,
	get_plate_stats
)
from tools.uploads_downloads import upload_temporary_subset, parse_compound_list, parse_id_list
from .dt import get_well_dictionary, manage_sw_status_change
from .forms import (
	AddVisitForm,
	FindCompoundForm,
	LibraryForm, 
	LibraryPlateForm, 
	PlateUpdateForm, 
	NewPresetForm, 
	EditPresetForm, 
	DTMapForm, 
	PlateOpeningForm, 
	OldProjectForm,
	NewProjectForm,
	UnsavedSubsetForm,
	FindPlateForm,
	ComparePlatesForm,
)
from tools.validators import data_is_valid, selection_is_valid
from tools.histograms import update_histograms
from tools.uploads_downloads import upload_plate, upload_subset, find_source_wells, source_wells_to_csv 
from API.models import Library, LibraryPlate, LibrarySubset, SourceWell, Compounds, Preset, PlateOpening, Project, IspybAuthorization
from django.contrib.auth.models import User

#VIEWS HANDLING GET REQUESTS
@staff_member_required
def index(request):
	return render(request, "inventory/inventory-index.html")

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
def projects(request):
	old_project_form = OldProjectForm()
	new_project_form = NewProjectForm()
	return render(request, "inventory/projects.html", {
		"old_project_form" : old_project_form, 
		"new_project_form" : new_project_form
		})

@staff_member_required
def browse_data(request):
	return render(request, "inventory/browse-data.html")

@staff_member_required
def presets(request):

	presets = Preset.objects.all().order_by("name")
	all_libs = current_library_selection(False)
	new_preset_form = NewPresetForm(libs=([("", "Select library...")] + all_libs))

	#make copy of presets data, and add information about availability of the compounds
	presets_copy = presets 
	form_dict = {}
	
	for p in presets:
		#generate lists of valid library choices for the preset
		old_libs_set = {(s.library.id, s.library.name) for s in p.subsets.all()}
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
		"presets": presets, 
		"new_preset_form": new_preset_form,
		"form_dict" : form_dict,
		})

@staff_member_required
def update_plate(request, pk):
	plate = LibraryPlate.objects.get(pk=pk)
	compounds_count = plate.compounds.all().count()
	
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
	
	plate_opening_form = PlateOpeningForm(initial={"date" : datetime.today()})
		
	return render(request, "inventory/update_plate.html", {
	#return render(request, "inventory/background.html", {
		"plate": plate, 
		"compounds" : plate.compounds.all().order_by("deactivation_date"),
		"plate_form" : plate_form,
		"plate_opening_form" : plate_opening_form, 
		"active_count" : active_count, 
		"inactive_count" : inactive_count,
		"availability" : availability,
		"dt_map_form" : dt_map_form,
		})

@staff_member_required
def track_usage(request, pk, date, mode):
		
	plate = LibraryPlate.objects.get(pk=pk)
	compounds = sw_copies(plate.compounds.all()) #make a copy of the data to edit it without touching the db
	modified_compounds = [compound for compound in compounds if len(compound.changes) > 0 ]
	timestamp = datetime.strptime(date, "%Y-%m-%d").date() #make a datetime object based on the url
	change_dates = get_change_dates(modified_compounds) #produce the list of all dates on which any compounds were (de)activated
	
	#generate strings needed to switch between the general and the graphic view
	if mode == "general-view":
		switch_view = "graphic view"
		switch_view_url = "graphic-view"
	else:
		switch_view = "general view"
		switch_view_url = "general-view"
	
	#activate all the compounds that were still active on the inspected date
	for c in modified_compounds:
		c = set_status(c, timestamp)
	
	#decide what size the table graphically representing the library plate should be
	size = get_plate_size(plate.compounds.all())
	rows = size['rows']
	columns = size['columns']
	
	#produce basic usage statistics
	usage_stats = get_usage_stats(compounds)
	opened = len([o for o in plate.opened.all() if o.date <= timestamp])
		
	return render(request, "inventory/track_usage.html", {
		"change_dates": change_dates,
		"compounds" : compounds,
		"plate": plate, 
		"date" : date,
		"timestamp" : timestamp,
		"active_count" : usage_stats['active'], 
		"inactive_count" : usage_stats['inactive'],
		"availability" : usage_stats['availability'],
		"opened" : opened,
		"columns" : columns,
		"rows": rows,
		"main_id" : mode,
		"switch_view": switch_view,
		"switch_view_url": switch_view_url,
		})

@staff_member_required
def proposal(request):
	form = OldProjectForm(request.POST)
	if request.method == "POST":	
		if form.is_valid():
			try:
				pr = form.cleaned_data['proposal']
				project = get_project_by_proposal(pr) #Project.objects.get(proposal=pr)
				subsets = get_subsets_with_availability(set(project.subsets.all()))
				
				return render(request, "inventory/proposal.html", {'proposal' : project, 'subsets': subsets, 'visit_form' : AddVisitForm()})
			except(ObjectDoesNotExist):
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': ['Proposal not found']})
				
		else:
			return render(request, "webSoakDB_backend/error_log.html", {'error_log': [form.errors, form.non_field_errors]})

#FORM ACTION VIEWS (HANDLING POST REQUESTS)
@staff_member_required
def add_library(request):
	form = LibraryForm(request.POST)
	
	if request.method == "POST":	
		if form.is_valid():
			name = form.cleaned_data['name']
			for_industry = form.cleaned_data['for_industry']
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
			filename = MEDIA_ROOT + '/' + fs.save(source.name, source)
			if data_is_valid(filename, log):
				
				barcode = form.cleaned_data['barcode']
				library_id = form.cleaned_data['library']
				current = form.cleaned_data['current']
				library = Library.objects.get(pk=library_id)
				today = str(date.today())
				
				plate = LibraryPlate.objects.create(library = library, barcode = barcode, current = current, last_tested = today)
				
				upload_plate(filename, plate)
				fs.delete(filename)
				if plate.current:
					update_histograms(plate.library, "library")
				
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
			log = []
			fs = FileSystemStorage()
			source = request.FILES["new_preset_compound_list"]
			library_id = form.cleaned_data['new_preset_library']
			filename = MEDIA_ROOT + '/' + fs.save(source.name, source)
			if selection_is_valid(filename, log, library_id):
				library = Library.objects.get(pk=library_id)
				description = form.cleaned_data['description']
				preset_name = form.cleaned_data['new_preset_name']
			
				#create new Subset object and upload data to it
				subset_name = library.name + " selection"
				origin = "preset: " + preset_name
				first_subset = upload_subset(filename, library_id, subset_name, origin)
				fs.delete(filename)
								
				#create preset object
				preset = Preset.objects.create(name=preset_name, description=description)
				preset.subsets.add(first_subset)
				preset.save()
				
				update_histograms(preset, "preset")
				
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
				filename_new = MEDIA_ROOT + '/' + fs_new.save(new_compound_list.name, new_compound_list) 
			if updated_compound_list:
				fs_edit =  FileSystemStorage()
				filename_edit = MEDIA_ROOT + '/' +  fs_edit.save(updated_compound_list.name, updated_compound_list)
				
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
			
			if new_library or edited_library or deleted_library:
				update_histograms(preset, "preset")
			
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
			plate = LibraryPlate.objects.get(pk=pk)
			old_current = plate.current
			old_library = plate.library
			plate.library = library
			plate.barcode = barcode
			plate.current = current
			plate.save()
			
			if old_current != current or old_library != library:
				update_histograms(library, "library")
			
			if old_library != library:
				update_histograms(old_library, "library")
			
			redirect_url = '/inventory/update-plate/' + str(plate.id) + '/'
			return HttpResponseRedirect(redirect_url)
		else:
			return render(request, "webSoakDB_backend/error_log.html", {"form_errors": form.errors, "non_field_errors": form.non_field_errors})	

@staff_member_required
def delete_library(request):
	if request.method == "POST":	
		pk = request.POST.get('id')
		lib = Library.objects.get(pk=pk)
		if len(lib.plates.all()) == 0:
			
			#remove all cached histograms
			path = MEDIA_ROOT + '/html_graphs/library/' + str(lib.id) + '/'
			try:
				shutil.rmtree(path)
			except FileNotFoundError:
				pass #in rare cases there are no graphs and no action needs to be taken

			lib.delete()
			return HttpResponseRedirect('../libraries/')
		else:
			return HttpResponseRedirect('../library-deletion-error/')

@staff_member_required
def delete_plate(request):
	if request.method == "POST":	
		pk = request.POST.get('id')
		plate = LibraryPlate.objects.get(pk=pk)
		plate.delete()
		
		if plate.current:
			update_histograms(plate.library, "library")
		return HttpResponseRedirect('../plates/')

@staff_member_required
def delete_preset(request):
	if request.method == "POST":	
		pk = request.POST.get('id')
		preset = Preset.objects.get(pk=pk)
		
		path = MEDIA_ROOT + '/html_graphs/preset/' + str(preset.id) + '/'

		try:
			shutil.rmtree(path)
		except FileNotFoundError:
			pass #in cases there are no graphs

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
		
		not_dispensed = set()
		#deactivate compounds that were not dispensed (unless already inactive)
		for key in request.POST:
			if key not in ['csrfmiddlewaretoken', "plate_id", "", "already_inactive"]:
				compound = SourceWell.objects.get(pk=key)
				not_dispensed.add(int(key))
				if compound.active:
					compound.active = False
					compound.save()
					
					manage_sw_status_change(compound, today, False)
		
		#activate inactive compounds that were dispensed anyway
		raw_string = request.POST.get('already_inactive')
		already_inactive = set([ int(n) for n in raw_string.split() ])
		activated = already_inactive - not_dispensed
		for sw_id in activated:
			compound = SourceWell.objects.get(pk=sw_id)
			compound.active = True
			compound.save()
			manage_sw_status_change(compound, today, True)
			
		PlateOpening.objects.create(plate=plate, date=today, reason="dispense test")
		if plate.current:
			update_histograms(plate.library, "library")
		
		return HttpResponseRedirect(redirect_url)

@staff_member_required
def deactivate_compounds_manually(request):
	if request.method == "POST":
		today = str(date.today())
		plate = LibraryPlate.objects.get(pk=request.POST.get('plate_id'))
		redirect_url = redirect_url = '/inventory/update-plate/' + str(plate.id) + '/'

		for key in request.POST:
			if key not in ['csrfmiddlewaretoken', "plate_id", ""]:
				compound = SourceWell.objects.get(pk=key)
				if compound.active:
					compound.active = False
					compound.save()
					
					manage_sw_status_change(compound, today, False)
		
		if plate.current:
			update_histograms(plate.library, "library")
		
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
			filename = fs.save(source.name, source)
			dt_dict = get_well_dictionary(MEDIA_ROOT + '/' + filename)
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

@staff_member_required
def open_plate(request):
	if request.method == "POST":
		try:
			plate_id = int(request.POST["id"])
			plate = LibraryPlate.objects.get(pk=plate_id)
			date_input = request.POST["date"]
			date = datetime.strptime(date_input, "%Y-%m-%d").date()
			reason = str(request.POST["reason"])
		except(ValueError):
			return render(request, "webSoakDB_backend/error_log.html", {"error_log": "<p>The values submitted in the form were invalid. Please try again.</p>"})	
			
		PlateOpening.objects.create(plate=plate, date=date, reason=reason)
		redirect_url = '/inventory/update-plate/' + str(plate_id) + '/'
		return HttpResponseRedirect(redirect_url)

@staff_member_required
def delete_multiple_plates(request):
	if request.method == "POST":
		deleted = []
		for key in request.POST:
			if key != "csrfmiddlewaretoken":
				deleted += [ int(n) for n in request.POST.getlist(key) ]
		
		for plate_id in deleted:
			plate = LibraryPlate.objects.get(pk=plate_id)
			plate.delete()
			if plate.current:
				update_histograms(plate.library, "library")

		return HttpResponseRedirect('/inventory/plates/')

@staff_member_required
def add_project(request):
	if request.method == "POST":
		form = NewProjectForm(request.POST)
		if form.is_valid():
			proposal = request.POST["proposal"]
			title = request.POST["title"]
			if request.POST.get("industry_user", False):
				industry_user = True
			else:
				industry_user = False

			first_visit = proposal + '-1'

			new_project = Project.objects.create(industry_user=industry_user)
			
			new_auth = IspybAuthorization.objects.create(project=title, proposal_visit=first_visit)			
			new_project.auth.add(new_auth)
			new_auth.save()
			
			return HttpResponseRedirect('../projects')
		else:
			return render(request, "webSoakDB_backend/error_log.html")	

def get_subset_map(request):
	if request.method == "POST":
		plate_id = request.POST["plate_id"]
		try:
			plate_ids = [ int(plate_id) ]
		except ValueError:
			plate_ids = parse_id_list(plate_id)
		subset_id = int(request.POST["subset_id"])
		
		if int(subset_id) == 0:
			subset = parse_compound_list(request.POST["compound_list"])
		else:
			subset = LibrarySubset.objects.get(pk=subset_id)
			
		source_wells = find_source_wells(subset, plate_ids)
		try:
			filename_prefix = subset.library.name + '-selection'
		except AttributeError:
			filename_prefix = "selection"
		response = source_wells_to_csv(source_wells, "files/soakdb-export.csv", filename_prefix)
		return response

def locate_compounds(request):
	libraries = [("", "Select library...")] + [(library.id, library.name) for library in Library.objects.filter(public=True)]
	
	if request.method == "GET":
		form = UnsavedSubsetForm(libs=libraries)
		return render(request, "inventory/locate_compounds.html", {"form": form, "subsets" : []})
	
	if request.method == "POST":
		form = UnsavedSubsetForm(data=request.POST, files=request.FILES, libs=libraries)
		if form.is_valid():
			log = []
			fs = FileSystemStorage()
			source = request.FILES["compound_list"]
			library_id = form.cleaned_data['library']
			filename = MEDIA_ROOT + '/' + fs.save(source.name, source)
			if selection_is_valid(filename, log, library_id):
				compounds = [c for c in upload_temporary_subset(filename, library_id)]
				lib = Library.objects.get(pk=library_id)
				subset = get_subsets_with_availability(compounds, lib)						
				fs.delete(filename)
				return render(request, "inventory/locate_compounds.html", {
					"form": form, 
					"subsets" : subset,
					"compound_list": compounds,
					})
			else:
				fs.delete(filename)
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})
		else:
			return render(request, "inventory/locate_compounds.html", {
				"form": form, 
				"errors" : [form.errors, form.non_field_errors ]
				})

def find_single_compound(request):
	if request.method == "GET":
		form =FindCompoundForm()
		return render(request, "inventory/find_single_compound.html", {'form' : form})
	
	if request.method == "POST":
		form =FindCompoundForm(data=request.POST)
		if form.is_valid():
			smiles = standardize_smiles(form.cleaned_data['smiles'])
			code = form.cleaned_data['code']

			compounds = []
			if smiles:
				compounds = Compounds.objects.filter(smiles=smiles)
			if code:
				compounds = Compounds.objects.filter(code=code)

			return render(request, "inventory/find_single_compound.html", {
				'form' : form, 
				'compounds': compounds,
				'smiles': smiles,
				'code': code})

def find_plate(request):
	if request.method == "GET":
		form =FindPlateForm()
		return render(request, "inventory/find_plate.html", {'form' : form})
	
	if request.method == "POST":
		form =FindPlateForm(data=request.POST)
		if form.is_valid():
			barcode = form.cleaned_data['barcode']
			print(barcode)

			plates = LibraryPlate.objects.filter(barcode=barcode)
			print(plates)
			return render(request, "inventory/find_plate.html", {
				'form' : form, 
				'plates': plates,
				'barcode': barcode})

def compare_plates(request):
	if request.method == "GET":
		form =ComparePlatesForm()
		return render(request, "inventory/compare_plates.html", {'form' : form})
	
	if request.method == "POST":
		form =ComparePlatesForm(data=request.POST)
		if form.is_valid():
			results = {}
			p1_data = {}
			p2_data = {}
			plate1_id = form.cleaned_data['plate1']
			plate2_id = form.cleaned_data['plate2']

			plate1 = LibraryPlate.objects.get(pk=plate1_id)
			plate2 = LibraryPlate.objects.get(pk=plate2_id)
			
			p1_compounds = set([sw.compound for sw in plate1.compounds.all() ])
			p2_compounds = set([sw.compound for sw in plate2.compounds.all() ])
			p1_smiles = set([sw.compound.smiles for sw in plate1.compounds.all() ])
			p2_smiles = set([sw.compound.smiles for sw in plate2.compounds.all() ])
			p1_available = set([sw.compound.smiles for sw in plate1.compounds.all() if sw.active])
			p2_available = set([sw.compound.smiles for sw in plate2.compounds.all() if sw.active])
			p1_unavailable = set([sw.compound.smiles for sw in plate1.compounds.all() if not sw.active])
			p2_unavailable = set([sw.compound.smiles for sw in plate2.compounds.all() if not sw.active])

			common_compounds = p1_compounds.intersection(p2_compounds)
			common_smiles = p1_smiles.intersection(p2_smiles)
			common_available = p1_available.intersection(p2_available)
			common_unavailable = p1_unavailable.intersection(p2_unavailable)
			smiles_with_different_codes = [s for s in common_smiles if s not in [c.smiles for c in common_compounds]]
			
			results["Compounds in common"] = len(common_smiles)
			results["Compounds in common (including the same code)"] = len(common_compounds)
			results["Compounds that have a different code in another plate"] = len(common_unavailable)
			results["Available compounds in common"] = len(common_available)
			results["Unvailable compounds in common"] = len(common_unavailable)
			results["Smiles with different codes"] = len(smiles_with_different_codes)
			p1_data = get_plate_stats(plate1, common_smiles, smiles_with_different_codes)
			p2_data = get_plate_stats(plate2, common_smiles, smiles_with_different_codes)
			
			return render(request, "inventory/compare_plates.html", {
				'form' : form, 
				'results' : results,
				'p1_data' : p1_data,
				'p2_data' : p2_data,
				'plate1' : plate1.barcode + ': ' + plate1.library.name,
				'plate2' : plate2.barcode + ': ' + plate2.library.name,
				})
		else:
			print('invalid form')


def compound_lookup(request, pk):
	compound = Compounds.objects.get(pk=pk)
	return render(request, "compound_lookup.html", {'compound': compound})


def preset_availability(request, pk):
	preset = Preset.objects.get(pk=pk)
	preset = fake_preset_copy(preset)
	return render(request, "inventory/preset-availability.html", {'preset': preset})

def subset_availability(request, pk):
	subset = LibrarySubset.objects.get(pk=pk)
	subset_copy = SubsetCopyWithAvailability(subset)
	return render(request, "inventory/subset_availability.html", {'subset': subset_copy})