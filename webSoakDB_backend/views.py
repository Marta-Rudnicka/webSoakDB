from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse

# Create your views here.

@csrf_exempt
def data_test(request):
	if request.method == "POST":
		print('fired data_test')
		print('request.POST', request.POST)
	
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
