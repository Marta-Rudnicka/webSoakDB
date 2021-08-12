from inventory.views import libraries
from django.urls import path
from inventory.views import libraries
from inventory.views import other
from inventory.views import plates
from inventory.views import presets
from inventory.views import queries
#import views

urlpatterns = [
		
	path("libraries/", libraries.libraries, name="libraries"),
	path("add-library/", libraries.add_library, name="add-library"),
	path("edit-library/", libraries.edit_library, name="edit-library"),
	path("delete-library/", libraries.delete_library, name="delete-library"),
	path("library-deletion-error/", libraries.library_deletion_error, name="library-deletion-error"),
	path("", other.index, name="index"),

	
	path("plates/", plates.plates, name="plates"),
	path("add-plate/", plates.add_plate, name="add-plate"),
	path("update-plate/<int:pk>/", plates.update_plate, name="update-plate"),
	path("delete-plate/", plates.delete_plate, name="delete-plate"),
	path("delete-multiple-plates/", plates.delete_multiple_plates),
	
	path("deactivate-compounds/", plates.deactivate_compounds, name="deactivate-compounds"),
	path("deactivate-compounds-manually/", plates.deactivate_compounds_manually),
	path("activate-single-compound/", plates.activate_single_compound),
	path("edit-plate/", plates.edit_plate, name="edit-plate"),
	path("track-usage/<int:pk>/<str:date>/<str:mode>/", plates.track_usage, name="track-usage"),
	path("dispense-testing-map/", plates.dispense_testing_map, name="dispense-testing-map"),
	path("open-plate/", plates.open_plate),
	
	path("presets/", presets.presets, name="presets"),
	path("add-preset/", presets.add_preset, name="add-preset"),
	path("edit-preset/", presets.edit_preset, name="edit-preset"),
	path("delete-preset/", presets.delete_preset, name="delete-preset"),
	path("preset-availability/<int:pk>/", presets.preset_availability),

	path("proposal/", other.proposal, name="proposal"),
	path("add-visit/", other.add_visit),
	path("projects/", other.projects, name="projects"),
	path("create-project/", other.add_project, name="add_project"),
	path("subset-map/", queries.get_subset_map, name="subset-map"),
	path("locate-compounds/", queries.locate_compounds, name="locate-compounds"),

	path("compound-lookup/<int:pk>/", queries.compound_lookup, name="compound-lookup"),

	path("find-single-compound/", queries.find_single_compound),
	path("find-library-plate/", queries.find_plate),
	path("compound-availability/<int:pk>/", queries.subset_availability),
	path("browse-data/", queries.browse_data),
	path("compare-plates/", queries.compare_plates),
	path("users/", other.manage_users),
	path("add-staff-member/", other.add_staff_member),
	path("remove-staff-member/", other.remove_from_staff)
	
	
]
