from django.urls import path
from . import views


urlpatterns = [
	#general purpose views:
	path("library_list/", views.LibraryList.as_view(), name="library_list"),
	path("library_detail/<str:name>/", views.LibraryDetail.as_view(), name="library_detail"),
	path("proposal_plates/<str:name>/", views.ProposalPlateList.as_view(), name="proposal_detail"),
	path("compounds/<str:library>/<str:plate>", views.PlateCompoundList.as_view(), name="api_lib"),
	path("in_house_library_list", views.InHouseLibraryList.as_view(), name="in_house_library_list"),
	path("library_selection_list", views.CurrentPlateList.as_view(), name="library_selection_list"),
	path("library_plate_list", views.AllPlateList.as_view(), name="library_plate_list"),
	path("preset_list", views.PresetList.as_view(), name="preset_list"),
	#path("crystals_list", views.CrystalsInPlates.as_view(), name="crystals_list"), #needs debug
	path("proposals/", views.ProposalList.as_view(), name="proposals"),
]
