from django.urls import path
from . import views

urlpatterns = [
	path("", views.index, name="index"),
	
	path("libraries/", views.libraries, name="libraries"),
	path("add-library/", views.add_library, name="add-library"),
	path("edit-library/", views.edit_library, name="edit-library"),
	path("delete-library/", views.delete_library, name="delete-library"),
	path("library-deletion-error/", views.library_deletion_error, name="library-deletion-error"),
	
	path("plates/", views.plates, name="plates"),
	path("add-plate/", views.add_plate, name="add-plate"),
	path("update-plate/<int:pk>/", views.update_plate, name="update-plate"),
	path("delete-plate/", views.delete_plate, name="delete-plate"),
	path("delete-multiple-plates/", views.delete_multiple_plates),
	
	path("deactivate-compounds/", views.deactivate_compounds, name="deactivate-compounds"),
	path("edit-plate/", views.edit_plate, name="edit-plate"),
	path("track-usage/<int:pk>/<str:date>/<str:mode>/", views.track_usage, name="track-usage"),
	path("dispense-testing-map/", views.dispense_testing_map, name="dispense-testing-map"),
	path("open-plate/", views.open_plate),
	
	path("presets/", views.presets, name="presets"),
	path("add-preset/", views.add_preset, name="add-preset"),
	path("edit-preset/", views.edit_preset, name="edit-preset"),
	path("delete-preset/", views.delete_preset, name="delete-preset"),
	
	path("proposal/", views.proposal, name="proposal"),
	path("subset-map/", views.get_subset_map, name="subset-map"),
	path("locate-compounds/", views.locate_compounds, name="locate-compounds")
]
