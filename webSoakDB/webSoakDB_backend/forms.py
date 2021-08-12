from django import forms
from API.models import Library, Project
from tools.uploads_downloads import current_library_selection
	
libs = current_library_selection(False)
libraries = []

class DateInput(forms.DateInput):
	input_type= 'date'
	
class LibraryPlateForm(forms.Form): 
	library = forms.ChoiceField(choices=libs)
	name = forms.CharField(label='Plate name')
	current = forms.BooleanField(label='Set to current plate', required=False )
	data_file = forms.FileField(label='Upload compound data:')

class ExternalLibraryForm(forms.Form):
	project = forms.IntegerField(label="project")
	name = forms.CharField(label='Library name')
	data_file = forms.FileField(label='Upload compound data:')
	

class SubsetForm(forms.Form):
	project = forms.IntegerField(label="project")
	lib_id = forms.ChoiceField(label='Select library', choices=libs)
	name = forms.CharField(label='name')
	data_file = forms.FileField(label='Upload your selection:')
