from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse, HttpResponseRedirect
from .forms import LibraryPlateForm, ExternalLibraryForm, SubsetForm
from datetime import date
from .validators import data_is_valid, selection_is_valid
from .helpers import upload_plate, upload_subset, import_full_libraries, import_library_parts, export_form_is_valid
from django.core.files.storage import FileSystemStorage
from API.models import Library, LibraryPlate, LibrarySubset, Proposals
import mimetypes
from slugify import slugify

def redirect_to_login(request):
	return HttpResponseRedirect('accounts/login/')

def dashboard(request):
	if request.user.is_staff:
		return render(request, "webSoakDB_backend/dashboard.html", {'user' : request.user})
	return HttpResponseRedirect('/selection/')
		

def upload_user_library(request):
	print('fired upload user library')
	today = str(date.today())
	
	if request.method == "POST":
		form = ExternalLibraryForm(request.POST, request.FILES)
		if form.is_valid():
			print('form valid')
			log = []
			fs = FileSystemStorage()
			source = request.FILES["data_file"]
			filename = fs.save(source.name, source)
			if data_is_valid(filename, log):
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
			filename = fs.save(source.name, source)
			if selection_is_valid(filename, log, library_id):
				print('Valid data; uploading to db')
				
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
				print('Invalid data - no upload')
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
			error_str = ["Application error:  cannot generate download file. Please report to developers. You should never see this message unless you manually edit the HTML in the browser"]
			return render(request, "webSoakDB_backend/error_log.html", {'error_log': error_str})
		
		prop = request.POST.get('proposal', False)
		source_wells = import_full_libraries(prop) + import_library_parts(prop, request.POST)
		
		with open('files/soakdb-export.csv', 'r+') as f:
			f.truncate(0)
			for c in source_wells:
				line = c.library_plate.barcode + ',' + c.well + ',' + c.library_plate.library.name + ',' + c.compound.smiles + ',' + c.compound.code + "\n"
				f.write(line)
			f.close()
		
		with open('files/soakdb-export.csv', 'r+') as f:
			filename = prop + '-' + 'soakDB-source-export-' + str(date.today()) + '.csv'
			response = HttpResponse(f, content_type='text/csv')
			response['Content-Disposition'] = "attachment; filename=%s" % filename
			return response
		 
def dummy(request):
	return render(request, "webSoakDB_backend/dummy.html");

def formatting(request):
	return render(request, "webSoakDB_backend/formatting-help.html");
