from django.shortcuts import render
from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponseRedirect
from django.core.files.storage import FileSystemStorage
from django.contrib.admin.views.decorators import staff_member_required
from webSoakDB_stack.settings import MEDIA_ROOT
from datetime import date, datetime
import re
from . import inv_helpers as h
from .dt import get_well_dictionary, manage_sw_status_change
from .. import forms as forms
from tools.validators import data_is_valid
from tools.histograms import update_histograms
from tools.uploads_downloads import upload_plate, current_library_selection
from API.models import Library, LibraryPlate, SourceWell, PlateOpening
	
@staff_member_required
def plates(request):
	libraries = Library.objects.filter(public=True).order_by("name")
	plate_form = forms.LibraryPlateForm(libs=current_library_selection(True))
	return render(request, "inventory/plates.html", {"libraries": libraries, "plate_form" : plate_form})

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
	plate_form = forms.PlateUpdateForm(libs=libs, initial={
		"library" : plate.library.id, 
		"barcode" : plate.barcode, 
		"name" : plate.name, 
		"current" : plate.current
		})
	dt_map_form = forms.DTMapForm()
	
	#plate_opening_form = forms.PlateOpeningForm(initial={"date" : datetime.today()})
		
	return render(request, "inventory/update_plate.html", {
	#return render(request, "inventory/background.html", {
		"plate": plate, 
		"compounds" : plate.compounds.all().order_by("deactivation_date"),
		"plate_form" : plate_form,
		#"plate_opening_form" : plate_opening_form, 
		"active_count" : active_count, 
		"inactive_count" : inactive_count,
		"availability" : availability,
		"dt_map_form" : dt_map_form,
		})

@staff_member_required
def track_usage(request, pk, date, mode):
		
	plate = LibraryPlate.objects.get(pk=pk)
	compounds = h.sw_copies(plate.compounds.all()) #make a copy of the data to edit it without touching the db
	#print('compounds[0]: ', compounds[0], compounds[0].changes)
	modified_compounds = [compound for compound in compounds if len(compound.changes) > 0 ]
	#print("modified_compounds: ", modified_compounds)
	timestamp = datetime.strptime(date, "%Y-%m-%d").date() #make a datetime object based on the url
	change_dates = h.get_change_dates(modified_compounds) #produce the list of all dates on which any compounds were (de)activated
	
	#print('change_dates: ', change_dates)
	#generate strings needed to switch between the general and the graphic view
	if mode == "general-view":
		switch_view = "graphic view"
		switch_view_url = "graphic-view"
	else:
		switch_view = "general view"
		switch_view_url = "general-view"
	
	#activate all the compounds that were still active on the inspected date
	for c in modified_compounds:
		c = h.set_status(c, timestamp)
	
	#decide what size the table graphically representing the library plate should be
	size = h.get_plate_size(plate.compounds.all())
	rows = size['rows']
	columns = size['columns']
	
	#produce basic usage statistics
	usage_stats = h.get_usage_stats(compounds)
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
def add_plate(request):
	today = str(date.today())
	
	if request.method == "POST":
		form = forms.LibraryPlateForm(data=request.POST, files=request.FILES, libs=current_library_selection(True))
		
		if form.is_valid():
			log = []
			fs = FileSystemStorage()
			source = request.FILES["plate_map"]
			filename = MEDIA_ROOT + '/' + fs.save(source.name, source)
			if data_is_valid(filename, log):
				
				barcode = form.cleaned_data['barcode']
				name = form.cleaned_data['name']
				library_id = form.cleaned_data['library']
				current = form.cleaned_data['current']
				library = Library.objects.get(pk=library_id)
				today = str(date.today())
				
				plate = LibraryPlate.objects.create(library = library, name=name, barcode = barcode, current = current, last_tested = today)
				
				upload_plate(filename, plate)
				fs.delete(filename)
				if plate.current:
					update_histograms(plate.library, "library")
				
				return HttpResponseRedirect('../plates')
			else:
				
				fs.delete(filename)
				
				return render(request, "inventory/errors.html", {'error_log': log})
		else:
			print('invalid form')
			print(form.errors)
			return HttpResponseRedirect('../plates')

@staff_member_required
def edit_plate(request):
	form = forms.PlateUpdateForm(data=request.POST, libs=current_library_selection(True))
	
	if request.method == "POST":	
		if form.is_valid():
			pk = request.POST.get('id')
			library_id = form.cleaned_data['library']
			library = Library.objects.get(pk=library_id)
			barcode = form.cleaned_data['barcode']
			name = form.cleaned_data['name']
			current = form.cleaned_data['current']
			plate = LibraryPlate.objects.get(pk=pk)
			old_current = plate.current
			old_library = plate.library
			plate.library = library
			plate.barcode = barcode
			plate.name = name
			plate.current = current
			plate.save()
			
			if old_current != current or old_library != library:
				update_histograms(library, "library")
			
			if old_library != library:
				update_histograms(old_library, "library")
			
			redirect_url = '/inventory/update-plate/' + str(plate.id) + '/'
			return HttpResponseRedirect(redirect_url)
		else:
			return render(request, "inventory/errors.html", {"form_errors": form.errors, "non_field_errors": form.non_field_errors})	

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
	print('fired deactivate_compounds_manually')
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
def activate_single_compound(request):
	print('fired activate_compounds_manually')
	if request.method == "POST":
		print(request.POST)
		today = str(date.today())
		plate = LibraryPlate.objects.get(pk=request.POST.get('plate_id'))
		print(plate)
		sw = SourceWell.objects.get(pk=request.POST.get('c_id'))
		print(sw)
		sw.active = True
		sw.save()
		manage_sw_status_change(sw, today, True)

		if plate.current:
			update_histograms(plate.library, "library")
		redirect_url = '/inventory/update-plate/' + str(plate.id) + '/'
		return HttpResponseRedirect(redirect_url)

@staff_member_required
def dispense_testing_map(request):
	rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	columns = [ '0' + str(i) for i in range(1, 10)] + ['10', '11', '12']
	errors = []
	
	form = forms.DTMapForm(request.POST, request.FILES)
	
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
			return render(request, "inventory/errors.html", {"form_errors": form.errors, "non_field_errors": form.non_field_errors})	

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
			return render(request, "inventory/errors.html", {"error_log": "<p>The values submitted in the form were invalid. Please try again.</p>"})	
			
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
