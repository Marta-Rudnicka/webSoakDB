from django.urls import path
from . import views


urlpatterns = [
	#general purpose generic endpoints:
	path("library_list/", views.LibraryList.as_view(), name="library_list"),
	path("library_detail/<int:pk>/", views.LibraryDetail.as_view(), name="library_detail"),
	path("compounds/<int:pk>/", views.PlateCompoundList.as_view(), name="api_lib"),
	path("in_house_library_list/", views.InHouseLibraryList.as_view(), name="in_house_library_list"),
	path("plates_list/<str:library>/", views.LibPlatesList.as_view(), name="all_library_plates_list"),
	path("plate_detail/<int:pk>/", views.PlateDetail.as_view(), name="plate_detail"),
	#path("preset_list/", views.PresetList.as_view(), name="preset_list"),
	path("preset_list/", views.preset_list, name="preset_list"),
	path("preset_detail/<int:pk>/", views.PresetDetail.as_view(), name="preset_detail"),
	path("proposals/", views.ProposalList.as_view(), name="proposals"),
	path("proposals/<str:proposal>/", views.ProposalDetail.as_view(), name="proposal_detail"),
	path("current_library_options/", views.current_library_options, name="current_library_options"),
	
	#specific-purpose custom endpoints:
	
	#for Library with id=pk: list current plates, nested deep so you can access compound.compound.properties.<$attribute>
	#path("current_plates_stats/<int:pk>/", views.LibCurrentPlatesStatList.as_view(), name="current_plates_stats"),
	
	#for LibrarySubset with id=pk, list all compounds nested deep so you can access compound.compound.properties.<$attribute>
	path("subset_stats/<int:pk>/", views.SubsetStatList.as_view(), name="subset_stats"),
	
	#for Library with id=pk, list all current plates (no individual compound data)
	path("current_plates_list/<int:pk>/", views.CurrentPlatesForLib.as_view(), name="library_detail"),
	
	#list all in-house libraries with one current plate each
	#path("library_selection_list/", views.CurrentPlateList.as_view(), name="library_selection_list"),
	
	#to update library or subset selection by passing it an array of ints (representing ids)
	path("update_proposal_selection/<str:proposal>/", views.UpdateProposalSelection.as_view(), name="update_proposal_selection"),	

	#path("crystals_list", views.CrystalsInPlates.as_view(), name="crystals_list"), #needs debug


]
