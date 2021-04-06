from django.shortcuts import render, get_object_or_404
from .models import Library, LibraryPlate, Compounds, SourceWell, Protein, Proposals, Preset, CrystalPlate, LibrarySubset
from django.http import HttpResponseRedirect
from django.views.generic import ListView #RetrieveAPIView
from rest_framework import generics
from rest_framework.permissions import AllowAny
from django.views.decorators.csrf import csrf_exempt
from .serializers import (LibrarySerializer, 
							SourceWellStatSerializer,
							CurrentPlateSerializer, 
							PresetSerializer, 
							CrystalPlateSerializer, 
							ProposalListSerializer, 
							LibraryPlateSerializer, 
							ProposalUpdateSerializer,
							LibrarySubsetStatSerializer,
							CompoundsStatSerializer)

import itertools

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

class PlateDetail(generics.RetrieveAPIView):
	queryset = LibraryPlate.objects.all()
	serializer_class = LibraryPlateSerializer
	permission_classes = [AllowAny]

class InHouseLibraryList(generics.ListAPIView):
	queryset = Library.objects.filter(public=True)
	serializer_class = LibrarySerializer
	permission_classes = [AllowAny]

class PresetList(generics.ListAPIView):
	queryset = Preset.objects.all()
	serializer_class = PresetSerializer
	permission_classes = [AllowAny]

class PresetDetail(generics.RetrieveAPIView):
	queryset = Preset.objects.all()
	serializer_class = PresetSerializer
	permission_classes = [AllowAny]


class CrystalsInPlates(generics.ListAPIView):
	queryset = CrystalPlate.objects.all()
	serializer_class = CrystalPlateSerializer
	permission_classes = [AllowAny]

class ProposalList(generics.ListAPIView):
	queryset = Proposals.objects.all()	
	serializer_class = ProposalListSerializer
	permission_classes = [AllowAny]

class ProposalDetail(generics.RetrieveUpdateAPIView):
	queryset = Proposals.objects.all()	
	lookup_field = "name"
	serializer_class = ProposalListSerializer
	permission_classes = [AllowAny]

class PlateCompoundList(generics.ListAPIView):
	
	def get_queryset(self):
		self.plate = get_object_or_404(LibraryPlate, pk=self.kwargs['pk'])
		
		return self.plate.compounds.all()
	
	serializer_class = SourceWellStatSerializer
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
		plates = self.library.plates.filter(current=True);
		wells = list(itertools.chain.from_iterable([x.compounds.all() for x in plates]))
		compounds = [well.compound for well in wells]
		return compounds

	#serializer_class = LibraryPlateStatSerializer	
	serializer_class = CompoundsStatSerializer	
	permission_classes = [AllowAny]

class LibPlatesList(generics.ListAPIView):
	#list all plates, with details about compounds (for stats)
	def get_queryset(self):
		self.library = get_object_or_404(Library, name=self.kwargs['library'])
		return self.library.plates.all();
	serializer_class = LibraryPlateSerializer	
	permission_classes = [AllowAny]


class SubsetStatList(generics.RetrieveAPIView):
	#get subset data with details about compounds (for stats)
	queryset = LibrarySubset.objects.all()
	serializer_class = LibrarySubsetStatSerializer	
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
