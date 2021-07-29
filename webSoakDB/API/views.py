from django.db.models import query
from django.shortcuts import get_object_or_404, redirect
from .models import Library, LibraryPlate, Project, Preset, LibrarySubset
from rest_framework import generics, viewsets #, status
from rest_framework.permissions import AllowAny
from django.http import JsonResponse
from .serializers import (
	LibrarySerializer, 
	SourceWellStatSerializer,
	CurrentPlateSerializer, 
	PresetSerializer, 
	ProjectListSerializer, 
	LibraryPlateSerializer, 
	ProjectUpdateSerializer,
	LibrarySubsetStatSerializer,
	LibrarySubsetSerializer,
	ProjectCompoundsSerializer
)

#TODO: change persmissions everywhere; temporarily [AllowAny] for early stages of testing

#general-purpose generic endpoints

def choose_plate_view(request, pk, project_id):
	#decide whether to use a view that requires authorization

	plate = LibraryPlate.objects.get(pk=pk)
	if plate.library.public:
		return redirect('/api/public_library_plates/' + str(pk) + '/')
	else:
		return redirect('/api/project_compounds/' + str(project_id) +  '/')

class ProjectCompoundsViewSet(viewsets.ReadOnlyModelViewSet):
	#returns the project with all the libraries and their details
	#linked by foreign keys - used for getting details of a non-public library plate
	queryset = Project.objects.all()
	serializer_class = ProjectCompoundsSerializer

class PublicLibraryViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Library.objects.filter(public=True)
    serializer_class = LibrarySerializer

class LibraryViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Library.objects.all()
    serializer_class = LibrarySerializer

class PublicLibraryPlateViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = LibraryPlate.objects.filter(library__public=True)
    serializer_class = LibraryPlateSerializer

class PlateWithCompoundsViewSet(viewsets.ReadOnlyModelViewSet):
	def get_queryset(self):
		self.plate = get_object_or_404(LibraryPlate, pk=self.kwargs['pk'])
		if self.plate.library.public:	
			return self.plate.compounds.filter(active=True)
		else:
			#this should not happen in normal use, but just in case it will not show any data
			return None
	serializer_class = SourceWellStatSerializer

class ProjectViewSet(viewsets.ReadOnlyModelViewSet):
	queryset = Project.objects.all()
	serializer_class = ProjectListSerializer

class PresetDetail(generics.RetrieveAPIView):
	queryset = Preset.objects.all()
	serializer_class = PresetSerializer
	permission_classes = [AllowAny]

class LibPlatesList(generics.ListAPIView):
	#list all plates, with details about compounds (for stats)
	def get_queryset(self):
		self.library = get_object_or_404(Library, pk=self.kwargs['pk'])
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

class UpdateProjectSelection(generics.RetrieveUpdateAPIView):
	authentication_classes = []
	queryset = Project.objects.all()
	serializer_class = ProjectUpdateSerializer
	permission_classes = [AllowAny]
	
class CurrentPlatesForLib(generics.ListAPIView):
	def get_queryset(self):
		library = get_object_or_404(Library, id=self.kwargs['pk'])
		return LibraryPlate.objects.filter(library=library, current=True)
	serializer_class = CurrentPlateSerializer
	permission_classes = [AllowAny]

def current_library_options(request):
	#In-house libraries with information about: 1) the number of available
	#compounds in current plates; 2) single-library presets related to them
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
		if current_plates.count() > 0:
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
