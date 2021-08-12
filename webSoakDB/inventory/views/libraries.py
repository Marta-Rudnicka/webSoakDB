from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.contrib.admin.views.decorators import staff_member_required
from webSoakDB_stack.settings import MEDIA_ROOT
import shutil
from . import inv_helpers as h
from .. import forms as forms
from API.models import Library


@staff_member_required
def libraries(request):
	libraries = Library.objects.filter(public=True).order_by("name")
	library_form = forms.LibraryForm()
	form_dict = {}
	for lib in libraries:
		form_dict[lib] = forms.LibraryForm(initial={"name": lib.name, "for_industry": lib.for_industry})
		
	return render(request, "inventory/libraries.html", {"libraries": libraries, "library_form": library_form, "form_dict" : form_dict})


@staff_member_required
def add_library(request):
	form = forms.LibraryForm(request.POST)
	
	if request.method == "POST":	
		if form.is_valid():
			name = form.cleaned_data['name']
			for_industry = form.cleaned_data['for_industry']
			Library.objects.create(name=name, for_industry=for_industry, public=True)

			return HttpResponseRedirect('../libraries')

@staff_member_required
def edit_library(request):
	form = forms.LibraryForm(request.POST)
	
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

def library_deletion_error(request):
	return render(request, "inventory/library_deletion_error.html");
