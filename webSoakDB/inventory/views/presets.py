from tools.data_storage_classes import SubsetCopyWithAvailability
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.core.files.storage import FileSystemStorage
from django.contrib.admin.views.decorators import staff_member_required
from webSoakDB_stack.settings import MEDIA_ROOT
import shutil
from . import inv_helpers as h
from .. import forms as forms
from tools.validators import selection_is_valid
from tools.histograms import update_histograms
from tools.uploads_downloads import upload_subset, current_library_selection
from API.models import Library, LibrarySubset, Preset

@staff_member_required
def presets(request):

	presets = Preset.objects.all().order_by("name")
	all_libs = current_library_selection(False)
	new_preset_form = forms.NewPresetForm(libs=([("", "Select library...")] + all_libs))

	#make copy of presets data, and add information about availability of the compounds
	presets_copy = presets 
	form_dict = {}
	
	for p in presets:
		#generate lists of valid library choices for the preset
		old_libs_set = {(s.library.id, s.library.name) for s in p.subsets.all()}
		new_libs = [("", "Select library...")] + list(set(all_libs) - old_libs_set)
		old_libs = [("", "Select library...")] + list(old_libs_set)
		
		form_dict[p] = forms.EditPresetForm(old_libs=old_libs, new_libs=new_libs, initial={
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
def add_preset(request):
	if request.method == "POST":
		form = forms.NewPresetForm(data=request.POST, files=request.FILES, libs=current_library_selection(True))
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
				return render(request, "inventory/errors.html", {'error_log': log})
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
		form = forms.EditPresetForm(data=request.POST, files=request.FILES, old_libs=old_libs, new_libs=new_libs)
		#https://django-gotchas.readthedocs.io/en/latest/forms.html <-- explanation for the weird kwargs; 
		
		if form.is_valid():
			print('VALID FORM')
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
				return render(request, "inventory/errors.html", {'error_log': log})
			
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
			print('INVALID FORM')
			return render(request, "inventory/errors.html", {"form_errors": form.errors, "non_field_errors": form.non_field_errors})			

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

def preset_availability(request, pk):
	preset = Preset.objects.get(pk=pk)
	preset = h.fake_preset_copy(preset)
	return render(request, "inventory/preset-availability.html", {'preset': preset})

def subset_availability(request, pk):
	subset = LibrarySubset.objects.get(pk=pk)
	subset_copy = SubsetCopyWithAvailability(subset)
	return render(request, "inventory/subset_availability.html", {'subset': subset_copy})
