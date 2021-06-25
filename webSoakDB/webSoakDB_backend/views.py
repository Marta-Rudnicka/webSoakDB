from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.core.files.storage import FileSystemStorage
from rest_framework.response import Response

from .forms import ExternalLibraryForm, SubsetForm
from tools.histograms import get_histogram, get_all_histograms, get_selection_histogram
from tools.validators import data_is_valid, selection_is_valid, export_form_is_valid
from tools.uploads_downloads import source_wells_to_csv, upload_plate, upload_subset, import_full_libraries, import_library_parts
from API.models import Library, LibraryPlate, LibrarySubset, Proposals, Compounds, Preset
from webSoakDB_stack.settings import MEDIA_ROOT

from slugify import slugify
from datetime import date, datetime
from rdkit import Chem
from rdkit.Chem import Draw

def upload_user_library(request):
	today = str(date.today())
	
	if request.method == "POST":
		form = ExternalLibraryForm(request.POST, request.FILES)
		if form.is_valid():
			print('form valid')
			log = []
			fs = FileSystemStorage()
			source = request.FILES["data_file"]
			filename = MEDIA_ROOT + '/' + fs.save(source.name, source)
			print(filename)
			if data_is_valid(filename, log):
				print('data is valid')
				#data to be submitted
				submitted_name = form.cleaned_data['name']
				proposal_name = form.cleaned_data['proposal']
				name = submitted_name + '(' + proposal_name + ')'
				today = str(date.today())
				
				#create new Library and LibraryPlate objects
				user_lib = Library.objects.create(name=name, public=False, for_industry=True)
				user_plate = LibraryPlate.objects.create(library = user_lib, barcode = name, current = True, last_tested = today)
				
				#upload the compound data for the new library plate
				upload_plate(filename, user_plate)
				fs.delete(filename)
				
				#add new library to user's proposal
				proposal = Proposals.objects.get(proposal=proposal_name)
				proposal.libraries.add(user_lib)
				proposal.save()
			
				fs.delete(filename)
				return render(request, "webSoakDB_backend/upload_success.html")
			else:
				fs.delete(filename)
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})

def upload_user_subset(request):
	if request.method == "POST":
		form = SubsetForm(request.POST, request.FILES)
		
		if form.is_valid():
			log = []
			fs = FileSystemStorage()
			source = request.FILES["data_file"]
			library_id = form.cleaned_data['lib_id']
			filename = MEDIA_ROOT + '/' + fs.save(source.name, source)
			if selection_is_valid(filename, log, library_id):
								
				#data to be submitted
				name = form.cleaned_data['name']
				proposal_name = form.cleaned_data['proposal']
				origin = "User selection for proposal " + proposal_name
				
				#create new LibrarySubset object and upload data to it
				subset = upload_subset(filename, library_id, name, origin)
				
				#add new subset to user's proposal
				proposal = Proposals.objects.get(proposal=proposal_name)
				proposal.subsets.add(subset)
				proposal.save()
			
				fs.delete(filename)
				return render(request, "webSoakDB_backend/upload_success.html")
			else:
				fs.delete(filename)
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})

def download_current_plate_map(request, pk):
	
	plate = LibraryPlate.objects.get(pk=pk)
	compounds = plate.compounds.filter(active=True)
	
	with open('files/plate-map.csv', 'r+') as f:
		f.truncate(0)
		for compound in compounds:
			line = compound.compound.code + ',' + compound.well + ',' + compound.compound.smiles + ','
			if compound.concentration:
				line += str(compound.concentration)
			line += "\n"
			f.write(line)
		f.close()
	
	
	with open('files/plate-map.csv', 'r+') as f:
		filename = slugify(plate.library.name) + '-' + slugify(plate.barcode) + '-map-' + str(date.today()) + '.csv'
		response = HttpResponse(f, content_type='text/csv')
		response['Content-Disposition'] = "attachment; filename=%s" % filename
		return response

def export_selection_for_soakdb(request):
	if request.method == "POST":
		if not export_form_is_valid(request.POST):
			error_str = ["Application error:  cannot generate download file. Please report to developers."]
			return render(request, "webSoakDB_backend/error_log.html", {'error_log': error_str})
		
		prop = request.POST.get('proposal', False)
		source_wells = import_full_libraries(prop) + import_library_parts(prop, request.POST)
		response = source_wells_to_csv(source_wells, "files/soakdb-export.csv", prop)
		return response

def serve_2d(request, pk):
	compound = Compounds.objects.get(pk=pk)

	mol = Chem.MolFromSmiles(compound.smiles)
	#img = Draw.MolToImage(mol)
	#response = HttpResponse(content_type="image/png")

	#img.save(response, "PNG")
	##############
	
	d2d = Draw.MolDraw2DSVG(300, 300)
	d2d.DrawMolecule(mol)
	d2d.FinishDrawing()
	svg_string = d2d.GetDrawingText()
	response2 = HttpResponse(svg_string, content_type="image/svg+xml")
	
	return response2


def serve_histogram(request, obj_type, pk, attr):
	if obj_type=="library":
		obj = Library.objects.get(pk=pk)
	if obj_type=="preset":
		obj = Preset.objects.get(pk=pk)
	if obj_type=="subset":
		obj = LibrarySubset.objects.get(pk=pk)

	g = get_histogram(obj, obj_type, attr)
	if g==204:
		print("no content")
		return HttpResponse(status=204)
	
	return HttpResponse(g)
	
def selection_histogram(request, attr):
	if request.method =="POST":
		#attr = request.POST["attr"]
		try:
			libs = [int(s) for s in request.POST["libs"].split(",")]
		except(ValueError):
			libs = 0
		try:
			subs = [int(s) for s in request.POST["subs"].split(",")]
		except(ValueError):
			subs = 0
		response = get_selection_histogram(libs, subs, attr)
		return HttpResponse(response)
	else:
		return HttpResponse("<div>Loading...</div>")
		 
def dummy(request):
	return render(request, "webSoakDB_backend/dummy.html")

def formatting(request):
	return render(request, "webSoakDB_backend/formatting-help.html")

def redirect_to_login(request):
	return HttpResponseRedirect('accounts/login/')

def dashboard(request):
	if request.user.is_staff:
		return render(request, "webSoakDB_backend/dashboard.html", {'user' : request.user})
	return HttpResponseRedirect('/selection/')

def all_histograms(request, obj_type, pk):
	if obj_type=="library":
		obj = Library.objects.get(pk=pk)
	if obj_type=="preset":
		obj = Preset.objects.get(pk=pk)
	
	g = get_all_histograms(obj, obj_type)
	return HttpResponse(g)
