from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse
from .forms import LibraryPlateForm, ExternalLibraryForm, SubsetForm
from datetime import date
from .validators import data_is_valid, selection_is_valid
from .helpers import upload_plate, upload_subset
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponseRedirect
from API.models import Proposals, Library, LibraryPlate, LibrarySubset

def upload_user_library(request):
	today = str(date.today())
	
	if request.method == "POST":
		form = ExternalLibraryForm(request.POST, request.FILES)
		if form.is_valid():
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
				user_plate = LibraryPlate.objects.create(library = user_lib, name = name, current = True, last_tested = today)
				
				#upload the compound data for the new library plate
				upload_plate(filename, user_plate)
				fs.delete(filename)
				
				#add new library to user's proposal
				proposal = Proposals.objects.get(name=proposal_name)
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
				proposal = Proposals.objects.get(name=proposal_name)
				proposal.subsets.add(subset)
				proposal.save()
			
				fs.delete(filename)
				return render(request, "webSoakDB_backend/upload_success.html")
			else:
				print('Invalid data - no upload')
				fs.delete(filename)
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})


def dummy(request):
	return render(request, "webSoakDB_backend/dummy.html");

def formatting(request):
	return render(request, "webSoakDB_backend/formatting-help.html");
