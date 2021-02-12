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
				
				#create new Library and LibraryPlate ojects
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
			print('source: ', source)
			library_id = form.cleaned_data['lib_id']
			print("library_id: ", library_id)
			filename = fs.save(source.name, source)
			if selection_is_valid(filename, log, library_id):
				
				#data to be submitted
				name = form.cleaned_data['name']
				proposal_name = form.cleaned_data['proposal']
				origin = "User selection for proposal " + proposal_name
				print("library_id: ", library_id)
				print("proposal_name: ", proposal_name)
				print("origin: ", origin)
				#create new LibrarySubset objects
				library = Library.objects.get(id=library_id)
				subset = LibrarySubset.objects.create(library = library, name = name, origin = origin)
				
				#upload the compound data for the new subset
				upload_subset(filename, library_id, name, origin)
				fs.delete(filename)
				
				#add new subset to user's proposal
				proposal = Proposals.objects.get(name=proposal_name)
				proposal.subsets.add(subset)
				proposal.save()
			
				fs.delete(filename)
				return render(request, "webSoakDB_backend/upload_success.html")
			else:
				fs.delete(filename)
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})

#copy for comparison
def upload_library_plate(request):
	from datetime import date
	today = str(date.today())
	
	if request.method == "POST":
		form = LibraryPlateForm(request.POST, request.FILES)
		if form.is_valid():
			#create a new library plate
			library = Library.objects.get(id = form.cleaned_data['library'])
			plate_name = form.cleaned_data['name']
			current = form.cleaned_data['current']
			newplate = LibraryPlate.objects.create(library = library, name = plate_name, current = current, last_tested=today)
			
			#create SourceWell objects for the new library plate
			source = request.FILES["data_file"]
			fs = FileSystemStorage()
			filename = fs.save(source.name, source)
			upload_plate(filename, newplate)
			fs.delete(filename)
			
			return HttpResponseRedirect('testing')
	else:
		form = LibraryPlateForm()
		
				
	return render(request, "playground/page.html", {
		"title": "Testing page",
		"form": form,
		})
