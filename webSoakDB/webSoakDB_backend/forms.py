from django import forms
from API.models import Library, Proposals
from .helpers import create_lib_selection
from django.forms import formset_factory
	
libs = create_lib_selection()
proposals = list([item.id, item.proposal] for item in Proposals.objects.all())
libraries = list([item.id, item.name] for item in Library.objects.filter(public=True))

class DateInput(forms.DateInput):
	input_type= 'date'
	
class LibraryPlateForm(forms.Form): 
	library = forms.ChoiceField(choices=libs)
	name = forms.CharField(label='Plate name')
	current = forms.BooleanField(label='Set to current plate', required=False )
	data_file = forms.FileField(label='Upload compound data:')

class ExternalLibraryForm(forms.Form):
	proposal = forms.CharField(label="proposal")
	name = forms.CharField(label='Library name')
	data_file = forms.FileField(label='Upload compound data:')
	

class SubsetForm(forms.Form):
	proposal = forms.CharField(label="proposal")
	lib_id = forms.ChoiceField(label='Select library', choices=libraries)
	name = forms.CharField(label='name')
	data_file = forms.FileField(label='Upload your selection:')
