from django.shortcuts import render, get_object_or_404
from .models import Library, LibraryPlate, Compounds, SourceWell, Protein, Proposals, Preset, CrystalPlate
from django.http import HttpResponseRedirect
from django.views.generic import ListView #RetrieveAPIView
from rest_framework import generics
from rest_framework.permissions import AllowAny
from django.views.decorators.csrf import csrf_exempt
from .serializers import (LibrarySerializer, 
							SourceWellSerializer, 
							CurrentPlateSerializer, 
							PresetSerializer, 
							CrystalPlateSerializer, 
							ProposalListSerializer, 
							LibraryPlateSerializer, 
							LibraryPlateStatSerializer, 
							ProposalUpdateSerializer)

#TODO: change persmissions everywhere; temporarily [AllowAny] for early stages of testing

#general-purpose generic endpoints

class LibraryList(generics.ListAPIView):
	queryset = Library.objects.all()
	serializer_class = LibrarySerializer
	permission_classes = [AllowAny]

class LibraryDetail(generics.RetrieveAPIView):
	queryset = Library.objects.all()
	serializer_class = LibrarySerializer
	permission_classes = [AllowAny]
	#lookup_field = 'name'

class InHouseLibraryList(generics.ListAPIView):
	queryset = Library.objects.filter(public=True)
	serializer_class = LibrarySerializer
	permission_classes = [AllowAny]

class AllPlateList(generics.ListAPIView):
	queryset = LibraryPlate.objects.all()
	serializer_class = CurrentPlateSerializer
	permission_classes = [AllowAny]

class PresetList(generics.ListAPIView):
	queryset = Preset.objects.all()
	serializer_class = PresetSerializer
	permission_classes = [AllowAny]

class CrystalsInPlates(generics.ListAPIView):
	queryset = CrystalPlate.objects.all()
	serializer_class = CrystalPlateSerializer
	permission_classes = [AllowAny]

class ProposalList(generics.ListAPIView):
	queryset = Proposals.objects.all()	
	#lookup_field = "name"
	serializer_class = ProposalListSerializer
	permission_classes = [AllowAny]

class ProposalDetail(generics.RetrieveUpdateAPIView):
	queryset = Proposals.objects.all()	
	lookup_field = "name"
	serializer_class = ProposalListSerializer
	permission_classes = [AllowAny]


class ProposalPlateList(generics.ListAPIView):
#lists current library plates for libraries selected for the proposal
#both in-house and user libraries
	def get_queryset(self):
		self.plate_list = []
		self.proposal = get_object_or_404(Proposals, name=self.kwargs['name'])
		self.libraries = self.proposal.libraries.all()
		for library in self.libraries:
			self.plates = LibraryPlate.objects.filter(library=library, current=True)
			for plate in self.plates:
				self.plate_list.append(plate)
		return self.plate_list		
	serializer_class = LibraryPlateStatSerializer	
	permission_classes = [AllowAny]
	lookup_field = 'name'




class PlateCompoundList(generics.ListAPIView):
	
	def get_queryset(self):
		self.library = get_object_or_404(Library, name=self.kwargs['library'])
		self.plate = get_object_or_404(LibraryPlate, name=self.kwargs['plate'], library__name=self.kwargs['library'])
		return self.plate.compounds.all()
	
	serializer_class = SourceWellSerializer
	permission_classes = [AllowAny]

class CurrentPlateList(generics.ListAPIView):
	def get_queryset(self):
		libs = Library.objects.filter(public=True)
		plates = []
		for lib in libs:
			c = lib.plates.filter(current=True)
			try:
				plates.append(c[0])
			except IndexError:
				print('No current plate for ', lib)
		return plates
	serializer_class = CurrentPlateSerializer
	permission_classes = [AllowAny]


class LibCurrentPlatesStatList(generics.ListAPIView):
	#list all current plates, with details about compounds (for stats)
	def get_queryset(self):
		self.library = get_object_or_404(Library, id=self.kwargs['pk'])
		print('self.library: ', self.library)
		return self.library.plates.filter(current=True);
	serializer_class = LibraryPlateStatSerializer	
	permission_classes = [AllowAny]


class UpdateProposalSelection(generics.RetrieveUpdateAPIView):
	queryset = Proposals.objects.all()	
	lookup_field = "name"
	serializer_class = ProposalUpdateSerializer
	permission_classes = [AllowAny]
	
class CurrentPlatesForLib(generics.ListAPIView):
	def get_queryset(self):
		library = get_object_or_404(Library, id=self.kwargs['pk'])
		return LibraryPlate.objects.filter(library=library, current=True)
	serializer_class = CurrentPlateSerializer
	permission_classes = [AllowAny]
