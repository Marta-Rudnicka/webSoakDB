from django.urls import include, path
from . import views
from rest_framework import routers

router = routers.SimpleRouter()
router.register(r'libraries', views.LibraryViewSet)
router.register(r'public_library_plates', views.PublicLibraryPlateViewSet, basename='library_plates')
router.register(r'projects', views.ProjectViewSet)
router.register(r'project_compounds', views.ProjectCompoundsViewSet)

urlpatterns = router.urls

urlpatterns = [
	path("", include(router.urls)),
	path("plate_compounds/<int:pk>/", views.PlateWithCompoundsViewSet.as_view({'get': 'list'})),
	path("library_plates/<int:pk>/<int:project_id>/", views.choose_plate_view),
	path("subset_detail/<int:pk>/", views.SubsetDetail.as_view()),
	path("plates_list/<int:pk>/", views.LibPlatesList.as_view(), name="all_library_plates_list"),
	path("preset_list/", views.preset_list, name="preset_list"),
	path("preset_detail/<int:pk>/", views.PresetDetail.as_view(), name="preset_detail"),
	path("current_library_options/", views.current_library_options, name="current_library_options"),
	path("subset_stats/<int:pk>/", views.SubsetStatList.as_view(), name="subset_stats"),
	
	#for Library with id=pk, list all current plates (no individual compound data)
	path("current_plates_list/<int:pk>/", views.CurrentPlatesForLib.as_view(), name="library_detail"),
	
	#to update library or subset selection by passing it an array of ints (representing ids)
	path("update_proposal_selection/<int:pk>/", views.UpdateProjectSelection.as_view(), name="update_proposal_selection"),

]