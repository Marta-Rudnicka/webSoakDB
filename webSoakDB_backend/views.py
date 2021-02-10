from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse
from .forms import LibraryPlateForm, ExternalLibraryForm
from datetime import date
from .library_plate_validator import data_is_valid
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponseRedirect


# Create your views here.

log = ['<p>Line 7: COMPOUND CODE ERROR: Missing compound code.</p>', 
"<p>Line 9: CONCENTRATION ERROR. Concentration value '100A' is not a number!</p>", 
'<p>Line 10: SMILES STRING ERROR. Invalid SMILES string: \'CN(C)(C)C\'<br/><code class="rdkit-err">RDKit ERROR: [21:33:50] Explicit valence for atom # 1 N, 4, is greater than permitted</code></p>',
"<p>line 13: WELL NAME ERROR. Invalid well name: 'A103'</p>",
"<p>Line 27: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>"]


def library_upload_errors(request):
	return render(request, "webSoakDB_backend/error_log.html")
	#return render(request, "webSoakDB_backend/error_log.html", {
	#"error_log": log,
	#})

def testing(request):
	today = str(date.today())
	
	if request.method == "POST":
		form = LibraryPlateForm(request.POST, request.FILES)
		if form.is_valid():
			log = []
			source = request.FILES["data_file"]
			fs = FileSystemStorage()
			filename = fs.save(source.name, source)
			if data_is_valid(filename, log):
				fs.delete(filename)
				return render(request, "webSoakDB_backend/page.html", {
					"title": "Testing page",
					"form": form,
				})	
			else:
				fs.delete(filename)
				print(log)
				return render(request, "webSoakDB_backend/error_log.html", {'error_log': log})
				#return HttpResponseRedirect('library_upload_errors/', kwargs={"error_log" : log})
				
			
			'''
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
			'''
			#return HttpResponseRedirect('library_upload_errors', "error_log" = log)
	else:
		form = LibraryPlateForm()
		
				
	return render(request, "webSoakDB_backend/page.html", {
		"title": "Testing page",
		"form": form,
		})	

@csrf_exempt
def data_test(request):
	if request.method == "POST":
		print('fired data_test')
		print('request.POST', request.POST)
		print('request.FILES', request.FILES)
	
	return HttpResponse(status=204)

#copy
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
