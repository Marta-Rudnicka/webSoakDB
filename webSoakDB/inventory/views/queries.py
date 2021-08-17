from tools.compounds import standardize_smiles
from tools.validators import parse_smiles
from django.shortcuts import redirect, render
from django.core.files.storage import FileSystemStorage
from django.contrib.admin.views.decorators import staff_member_required
from webSoakDB_stack.settings import MEDIA_ROOT
from . import inv_helpers as h
from tools.uploads_downloads import upload_temporary_subset
from .. import forms as forms
from tools.validators import selection_is_valid
from API.models import Library, LibraryPlate, Compounds, LibrarySubset
from tools.uploads_downloads import upload_temporary_subset, parse_compound_list, parse_id_list, find_source_wells, source_wells_to_csv
from tools.data_storage_classes import SubsetCopyWithAvailability
#VIEWS HANDLING GET REQUESTS

@staff_member_required
def browse_data(request):
	return render(request, "inventory/browse-data.html")

def locate_compounds(request):
	libraries = [("", "Select library...")] + [(library.id, library.name) for library in Library.objects.filter(public=True)]
	
	if request.method == "GET":
		form = forms.UnsavedSubsetForm(libs=libraries)
		return render(request, "inventory/locate_compounds.html", {"form": form, "subsets" : []})
	
	if request.method == "POST":
		form = forms.UnsavedSubsetForm(data=request.POST, files=request.FILES, libs=libraries)
		if form.is_valid():
			log = []
			fs = FileSystemStorage()
			source = request.FILES["compound_list"]
			library_id = form.cleaned_data['library']
			filename = MEDIA_ROOT + '/' + fs.save(source.name, source)
			if selection_is_valid(filename, log, library_id):
				compounds = [c for c in upload_temporary_subset(filename, library_id)]
				lib = Library.objects.get(pk=library_id)
				subset = h.get_subsets_with_availability(compounds, lib)						
				fs.delete(filename)
				return render(request, "inventory/locate_compounds.html", {
					"form": form, 
					"subsets" : subset,
					"compound_list": compounds,
					})
			else:
				fs.delete(filename)
				return render(request, "inventory/errors.html", {'error_log': log})
		else:
			return render(request, "inventory/locate_compounds.html", {
				"form": form, 
				"errors" : [form.errors, form.non_field_errors ]
				})

def find_single_compound(request):
	if request.method == "GET":
		form = forms.FindCompoundForm()
		return render(request, "inventory/find_single_compound.html", {'form' : form})
	
	if request.method == "POST":
		form = forms.FindCompoundForm(data=request.POST)
		if form.is_valid():
			string = form.cleaned_data['string']
			if parse_smiles(string): 	#with a valid SMILES string, it returns None
				compounds = Compounds.objects.filter(code=string)
			else:
				smiles = standardize_smiles(string)
				compounds = Compounds.objects.filter(smiles=smiles)
			
			if compounds.count() == 0:
				not_found = 'Compound not found'
			else:
				not_found = ''
			
			return render(request, "inventory/find_single_compound.html", {
				'form' : form, 
				'compounds': compounds,
				'string': string,
				'not_found' : not_found,
				'results': True
				})

def find_plate(request):
	if request.method == "GET":
		form = forms.FindPlateForm()
		return render(request, "inventory/find_plate.html", {'form' : form})
	
	if request.method == "POST":
		form = forms.FindPlateForm(data=request.POST)
		if form.is_valid():
			barcode = form.cleaned_data['barcode']
			print(barcode)

			plates = LibraryPlate.objects.filter(barcode=barcode)

			return render(request, "inventory/find_plate.html", {
				'form' : form, 
				'plates': plates,
				'barcode': barcode,
				})

def compare_plates(request):
	if request.method == "GET":
		form = forms.ComparePlatesForm()
		return render(request, "inventory/compare_plates.html", {'form' : form})
	
	if request.method == "POST":
		form = forms.ComparePlatesForm(data=request.POST)
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
			results["Available compounds in common"] = len(common_available)
			results["Unvailable compounds in common"] = len(common_unavailable)
			results["Smiles with different codes"] = len(smiles_with_different_codes)
			p1_data = h.get_plate_stats(plate1, common_smiles, smiles_with_different_codes)
			p2_data = h.get_plate_stats(plate2, common_smiles, smiles_with_different_codes)

			plate1 = h.make_plate_name(plate1)
			plate2 = h.make_plate_name(plate2)
			
			return render(request, "inventory/compare_plates.html", {
				'form' : form, 
				'results' : results,
				'p1_data' : p1_data,
				'p2_data' : p2_data,
				'plate1' : plate1,
				'plate2' : plate2
				})
		else:
			print('invalid form')


def compound_lookup(request, pk):
	compound = Compounds.objects.get(pk=pk)
	return render(request, "compound_lookup.html", {'compound': compound})


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


def subset_availability(request, pk):
	subset = LibrarySubset.objects.get(pk=pk)
	subset_copy = SubsetCopyWithAvailability(subset)
	return render(request, "inventory/subset_availability.html", {'subset': subset_copy})