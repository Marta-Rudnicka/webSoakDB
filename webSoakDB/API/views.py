from django.shortcuts import get_object_or_404
from .models import Library, LibraryPlate, Proposals, Preset, LibrarySubset
from django.views.generic import ListView #RetrieveAPIView
from rest_framework import generics
from rest_framework.permissions import AllowAny
import json
from django.http import JsonResponse
from .serializers import (LibrarySerializer, 
							SourceWellStatSerializer,
							CurrentPlateSerializer, 
							PresetSerializer, 
							ProposalListSerializer, 
							LibraryPlateSerializer, 
							ProposalUpdateSerializer,
							LibrarySubsetStatSerializer,
							LibrarySubsetSerializer,
						)

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

class PresetDetail(generics.RetrieveAPIView):
	queryset = Preset.objects.all()
	serializer_class = PresetSerializer
	permission_classes = [AllowAny]

class ProposalList(generics.ListAPIView):
	queryset = Proposals.objects.all()	
	serializer_class = ProposalListSerializer
	permission_classes = [AllowAny]

class ProposalDetail(generics.RetrieveUpdateAPIView):
	permission_classes = [AllowAny]
	
	queryset = Proposals.objects.all()	
	lookup_field = "proposal"
	serializer_class = ProposalListSerializer
	

class PlateCompoundList(generics.ListAPIView):
	
	def get_queryset(self):
		self.plate = get_object_or_404(LibraryPlate, pk=self.kwargs['pk'])
		
		return self.plate.compounds.filter(active=True)
	
	serializer_class = SourceWellStatSerializer
	permission_classes = [AllowAny]

class LibPlatesList(generics.ListAPIView):
	#list all plates, with details about compounds (for stats)
	def get_queryset(self):
		self.library = get_object_or_404(Library, name=self.kwargs['library'])
		return self.library.plates.all();
	serializer_class = LibraryPlateSerializer	
	permission_classes = [AllowAny]

class SubsetDetail(generics.RetrieveAPIView):
	queryset = LibrarySubset.objects.all()
	serializer_class = LibrarySubsetSerializer
	permission_classes = [AllowAny]
	
class SubsetStatList(generics.RetrieveAPIView):
	#get subset data with details about compounds (for stats)
	queryset = LibrarySubset.objects.all()
	serializer_class = LibrarySubsetStatSerializer	
	permission_classes = [AllowAny]

class UpdateProposalSelection(generics.RetrieveUpdateAPIView):
	authentication_classes = []

	queryset = Proposals.objects.all()	
	lookup_field = "proposal"
	serializer_class = ProposalUpdateSerializer
	permission_classes = [AllowAny]
	
class CurrentPlatesForLib(generics.ListAPIView):
	def get_queryset(self):
		library = get_object_or_404(Library, id=self.kwargs['pk'])
		return LibraryPlate.objects.filter(library=library, current=True)
	serializer_class = CurrentPlateSerializer
	permission_classes = [AllowAny]



def current_library_options(request):
	'''In-house libraries with information about: 1) the number of available
	compounds in current plates; 2) single-library presets related to them'''
	class LibCopy:
		def __init__(self, library, current_plate, size, presets):
			self.id = library.id
			self.name = library.name
			self.current_plate = current_plate
			self.size = size
			self.presets = presets
	
	class PresetCopy:
		def __init__(self, preset):
			self.id = preset.id
			self.name = preset.name
			self.description = preset.description
			self.subsets = [s.id for s in preset.subsets.all()]
			self.library = preset.subsets.all()[0].library.id
			self.size = preset.subsets.all()[0].size()
			
	queryset = Library.objects.filter(public=True)
	options = []
	single_lib_presets = [p for p in Preset.objects.all() if p.subsets.count() == 1]
	presets_copy = [PresetCopy(p) for p in single_lib_presets]
	
	preset_dict = {p: p.library for p in presets_copy}
	
	for l in queryset:
		presets = [preset.__dict__ for (preset, lib) in preset_dict.items() if lib == l.id]
		current_plates = LibraryPlate.objects.filter(library=l, current = True)
		size = sum([p.compounds.filter(active=True).count() for p in current_plates])
		lib = LibCopy(l, current_plates[0].id, size, presets)
		options.append(lib.__dict__)
	
	return JsonResponse(options, safe=False)

def preset_list(request):
	class PresetCopy:
		def __init__(self, preset):
			self.id = preset.id
			self.name = preset.name
			self.description = preset.description
			self.subsets = [s.id for s in preset.subsets.all()]
			self.library = preset.subsets.all()[0].library.id
			self.size = sum([s.size() for s in p.subsets.all()])
	
	queryset = Preset.objects.all()
	presets = []
	
	for p in queryset:
		pc = PresetCopy(p)
		presets.append(pc.__dict__)
	
	return JsonResponse(presets, safe=False)
