from django import forms
from API.models import Library, Preset


presets = [("", "Select preset...")] + [(preset.id, preset.name) for preset in Preset.objects.all()]
#libs = [("", "Select library...")] + [(library.id, library.name) for library in Library.objects.filter(public=True)]
libs = []

'''Note on overriding __init__ : __init__ functions are overriden here to be able to dynamically generate choices
for the ChoiceFields. The reason why there are the oher fields in the __init__ functions is ordering: 
when some of the fields are declared inside __init__, and some outside of it, field_order doesn't work properly '''

class DateInput(forms.DateInput):
	input_type= 'date'

class LibraryForm(forms.Form):
	name = forms.CharField(label='Library name', required=True)
	for_industry = forms.BooleanField(label='For industry', required=False)

class PlateUpdateForm(forms.Form):
	def __init__(self, libs, *args, **kwargs):
		super(PlateUpdateForm, self).__init__(*args, **kwargs)
		self.fields['library'] = forms.ChoiceField(choices = libs, required=True)
		self.fields['barcode'] = forms.CharField(label='Barcode')
		self.fields['current'] = forms.BooleanField(label='Current', required=False )

class LibraryPlateForm(PlateUpdateForm): 
	plate_map = forms.FileField(label='Upload plate map:', required=True)
	
class AddLibToPreset(forms.Form):
	compound_list= forms.FileField(label='Upload compound list', required=True)
		
class NewPresetForm(forms.Form):
	def __init__(self, libs=[], *args, **kwargs):
		super(NewPresetForm, self).__init__(*args, **kwargs)
		self.fields['name'] = forms.CharField(label='Name', required=True )
		self.fields['description'] = forms.CharField(widget=forms.Textarea, required=False)
		self.fields['new_library'] = forms.ChoiceField(choices=libs, label='Add the first library to preset:', required=True)
		self.fields['new_compound_list'] = forms.FileField(label='Upload compounds list for the new library', required=True)


class EditPresetForm(NewPresetForm):
	def __init__(self, old_libs=[], new_libs=[], *args, **kwargs):
		super(EditPresetForm, self).__init__(*args, **kwargs)
		self.fields['new_library'] = forms.ChoiceField(choices = new_libs, label='Add new library to preset:', required=False)
		self.fields['new_compound_list'] = forms.FileField(label='Upload compounds list for the new library', required=False)
		self.fields['edited_library'] = forms.ChoiceField(choices = old_libs, label='Update selection from:', required=False)
		self.fields['updated_compound_list'] = forms.FileField(label='Upload new compound list for selected library', required=False)
		self.fields['deleted_library'] = forms.ChoiceField(choices = old_libs, label='Delete selection from:', required=False, widget=forms.Select(attrs={'class': 'delete-button'}))
