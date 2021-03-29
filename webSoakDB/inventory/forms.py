from django import forms
from API.models import Library, Preset
#from django.forms import formset_factory

libs = [("", "Select library...")] + [(library.id, library.name) for library in Library.objects.filter(public=True)]
presets = [("", "Select preset...")] + [(preset.id, preset.name) for preset in Preset.objects.all()]


class DateInput(forms.DateInput):
	input_type= 'date'

class LibraryForm(forms.Form):
	name = forms.CharField(label='Library name', required=True)
	for_industry = forms.BooleanField(label='For industry', required=False)

class PlateUpdateForm(forms.Form):
	library = forms.ChoiceField(choices=libs, required=True)
	barcode = forms.CharField(label='Barcode')
	current = forms.BooleanField(label='Current', required=False )

class LibraryPlateForm(PlateUpdateForm): 
	plate_map = forms.FileField(label='Upload plate map:', required=True)

class AddLibToPreset(forms.Form):
	library = forms.ChoiceField(choices=libs, required=True)
	compound_list= forms.FileField(label='Upload compound list', required=True)
		
class NewPresetForm(forms.Form):
	name = forms.CharField(label='Name', required=True )
	description = forms.CharField(widget=forms.Textarea, required=False)
	new_library = forms.ChoiceField(choices=libs, label='Add new library to preset:', required=True)
	new_compound_list= forms.FileField(label='Upload compounds list for the new library', required=True)
	#field_order = ['name', 'description', 'library', 'compound_list']


class EditPresetForm(NewPresetForm):
	#new_libs = [('', '...'), ('test1', 'test1'),('test2','test2')]
	#old_libs = [('', '...'), ('test3', 'test3'),('test4','test4')]
	
	def __init__(self, old_libs=[], new_libs=[], *args, **kwargs):
		super(EditPresetForm, self).__init__(*args, **kwargs)
		self.fields['new_library'] = forms.ChoiceField(choices = new_libs, label='Add new library to preset:', required=False)
		self.fields['new_compound_list'] = forms.FileField(label='Upload compounds list for the new library', required=False)
		self.fields['edited_library'] = forms.ChoiceField(choices = old_libs, label='Update selection from:', required=False)
		self.fields['updated_compound_list'] = forms.FileField(label='Upload new compound list for selected library', required=False)
		self.fields['deleted_library'] = forms.ChoiceField(choices = old_libs, label='Delete selection from:', required=False, widget=forms.Select(attrs={'class': 'delete-button'}))

	'''new_library = forms.ChoiceField(choices = new_libs, label='Add new library to preset:', required=False)
	new_compound_list = forms.FileField(label='Upload compounds list for the new library', required=False)
	edited_library = forms.ChoiceField(choices = old_libs, label='Update selection from:', required=False)
	updated_compound_list = forms.FileField(label='Upload new compound list for selected library', required=False)
	deleted_library = forms.ChoiceField(choices = old_libs, label='Delete selection from:', required=False, widget=forms.Select(attrs={'class': 'delete-button'}))'''

#class PresetLibForm(AddLibToPreset):
#	preset = forms.ChoiceField(choices=presets, required=True)
#	field_order = ['preset', 'library', 'compound_list']
	
