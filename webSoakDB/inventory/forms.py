from django import forms
#from API.models import Library, LibraryPlate, Preset, PlateOpening, Project
from .views.inv_helpers import make_plate_name

libs = []
#plates = [(plate.id, make_plate_name(plate) ) for plate in LibraryPlate.objects.filter(library__public=True)]
#plates = [("", "Select plate...")] + sorted(plates, key=lambda x: x[1])
plates = []
'''Note on overriding __init__ : __init__ functions are overriden here to be able to dynamically generate choices
for the ChoiceFields. The reason why there are the other fields in the __init__ functions is ordering: 
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
		self.fields['name'] = forms.CharField(label='Name', required=False)
		self.fields['current'] = forms.BooleanField(label='Current', required=False )

class LibraryPlateForm(PlateUpdateForm): 
	plate_map = forms.FileField(label='Upload plate map:', required=True)
	
class AddLibToPreset(forms.Form):
	compound_list= forms.FileField(label='Upload compound list', required=True)
		
class NewPresetForm(forms.Form):
	def __init__(self, libs=[], *args, **kwargs):
		super(NewPresetForm, self).__init__(*args, **kwargs)
		self.fields['new_preset_name'] = forms.CharField(label='Name', required=True )
		self.fields['description'] = forms.CharField(widget=forms.Textarea, required=False)
		self.fields['new_preset_library'] = forms.ChoiceField(choices=libs, label='Add the first library to preset:', required=True)
		self.fields['new_preset_compound_list'] = forms.FileField(label='Upload compounds list for the new library', required=True)


class EditPresetForm(forms.Form):
	def __init__(self, old_libs=[], new_libs=[], *args, **kwargs):
		super(EditPresetForm, self).__init__(*args, **kwargs)
		self.fields['name'] = forms.CharField(label='Name', required=True )
		self.fields['description'] = forms.CharField(widget=forms.Textarea, required=False)
		self.fields['new_library'] = forms.ChoiceField(choices = new_libs, label='Add new library to preset:', required=False)
		self.fields['new_compound_list'] = forms.FileField(label='Upload compounds list for the new library', required=False)
		self.fields['edited_library'] = forms.ChoiceField(choices = old_libs, label='Update selection from:', required=False)
		self.fields['updated_compound_list'] = forms.FileField(label='Upload new compound list for selected library', required=False)
		self.fields['deleted_library'] = forms.ChoiceField(choices = old_libs, label='Delete selection from:', required=False, widget=forms.Select(attrs={'class': 'delete-button'}))

class DTMapForm(forms.Form):
	dt_map = forms.FileField(label='Upload mapping file', required=False)

class OldProjectForm(forms.Form):
	proposal = forms.CharField(label="Proposal", required=True)

class UnsavedSubsetForm(forms.Form):
	def __init__(self, libs=[], *args, **kwargs):
		super(UnsavedSubsetForm, self).__init__(*args, **kwargs)
		self.fields['library'] = forms.ChoiceField(choices = libs, label='Library:', required=True)
		self.fields['compound_list'] = forms.FileField(label='Upload compounds list', required=True)

class NewProjectForm(forms.Form):
	proposal = forms.CharField(label="Proposal", required=True)
	title = forms.CharField(label="Name", required=False)
	industry_user = forms.BooleanField(label="Industry user", required=False)

class FindCompoundForm(forms.Form):
	string = forms.CharField(label="SMILES string or code", required=False)

class FindPlateForm(forms.Form):
	barcode = forms.CharField(label="Barcode", required=True)

class ComparePlatesForm(forms.Form):
	plate1 = forms.ChoiceField(label="Plate 1", choices = plates, required=True)
	plate2 = forms.ChoiceField(label="Plate 2", choices = plates, required=True)


class AddStaffUserFrom(forms.Form):
	username_staff = forms.CharField(label="Federal ID*")
	first_name_staff = forms.CharField(label="First Name", required=False)
	last_name_staff = forms.CharField(label="Last Name", required=False)
	email_staff = forms.CharField(label="E-mail", required=False)


class AddPowerUserFrom(forms.Form):
	username_pu = forms.CharField(label="Federal ID")
	first_name_pu = forms.CharField(label="First Name", required=False)
	last_name_pu = forms.CharField(label="Last Name", required=False)
	email_pu = forms.CharField(label="E-mail", required=False)