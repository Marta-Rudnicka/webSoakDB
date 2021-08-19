# API

The API app provides the views that `webSoakDB_frontend` uses to load the data from the database. It consist of views, serializers and contains the `admin.py` file, which covers models used by all the apps. It is built mostly with Django REST framework ( https://www.django-rest-framework.org/ ), and most of it is quite straightforward for someone familiar with DRF, so it will not be discussed in the documentation. Only the less obvious aspects are described.

To protect access sensitive data, the API app uses `ispyb_dja` ( https://github.com/xchem/ispyb_dja/ )

## Public vs. private library plates and the API views

Most of user-generated data is related to the `Project` model, which links to `IspybAuthorization` objects used for authentcation/authorization of the users. However the `LibraryPlate` model stores the data of both the public (XChem in-house) library plates and private ones (brought in by the user). The publicly available plates are not related to any project or user and need not be authenticated.

This is why there are completely different mechanisms used for loading library plate data.

When the user accessess a page listing all the compounds ( `SourceWell` objects) in a plate via the React app, the query is first directed to the `choose_plate_view()`. This view checks whether the plate is public or private and redirects the call to the URL that will lead to the appropriate view. The view is passed two URL parameters: the plate id, and the project id, and it passes them on depending on which view is used next

### Public library plates: `PublicLibraryPlateViewSet` and `PlateWithCompoundsViewSet`

These are standard viewsets. `choose_plate_view()` gives them the plate id (through the URL parameters), but not project id since the plates it shows are not related to any particular project. `PublicLibraryPlateViewSet` provides a standard detail view of a library plate with data about the plate itself. Just in case someone tries to use this view to access private data, it filters out all the private plates (when `queryset` variable is declared). `PlateWithCompoundsViewSet` is used to get a list of compounds (`SourceWell` objects) in the plate. To prevent somehow using this view to inspect a private plate, the `get_queryset` function checks if the plate is public first and returns `None` if someone attempts to access a private plate with it.

### Private library plates and `ProjectCompoundsViewSet`

A `LibraryPlate` instance  is not linked to any `IspybAuthorization` instance, so it is impossible to use `ispyb_dja` on this model. Instead, `ProjectCompoundsViewSet` provides an API endpoint that serves the data as an instance of `Project`. `Project` has a many-to-many field linking to `Library` instances selected to the project, and libraries are linked to `LibraryPlate` objects (private libraries have only one plate each). Thus, `ProjectCompoundsViewSet` loads the same data as `PublicLibraryPlateViewSet`, but in a more convoluted manner, following one foreign key after another until the structure is deep enough to access all the data necessary. There is no need to have a separate views for plate data and source well data, since the source well data is accessed via the plate data and the view provides both. The React component that use this endpoint filters out the unnecessary data and only display the ones related to the relevant plate (for more details see "loading API data in `<CompoundLookupPlate>`" in the `webSoakDB_frontend` documentation).

To achieve the necessary nesting, several serializers are used:
```
ProjectCompoundsSerializer (used in `ProjectCompoundsViewSet` )
	|
	|- LibraryInProjectSerializer (for `libraries` field)
		|
		|-LibraryPlateStatSerializer (for `plates` field)
			|
			|-LibrarySerializer ( for `library` field)
			|-SourceWellStatSerializer ( for `compounds` field)
				(this serialisers had `depth=1` attribute, which makes it serialise the following two models)
				* `SourceWell`
				* `Compounds` (a foreign key in `SourceWell`)

```
## Custom views (not using Django REST framework)

While DRF makes it very easy to create typical API endpoints, it is not very convenient for special cases, which is why two of the views use standard Django `JsonResponse`. While DRF views are all class-based, non-DRF views are function-based.

### `current_library_options()`

This view processes the library data for the list of libraries presented to the user for selection in order to produce additional information not directly stored in the database. The processing includes:
- counting all available compounds (`SourceWells` where `active=True`) in all the current plate(s) of the library
- finding presets that only include compounds from one library and associating them with that library 

To process that data, the view uses two classes (just regular Python classes), which store data copied from the ORM models and, plus the additional information as extra attributes. Then, the resulting `LibCopy` object is serialized into a JSON response. The Library data sent into the enpoint is the following:

- id - copied from the Library model
- name - copied from the Library model
- current_plate - one current plate (found by the view) - it is used to provide the link to this plate from the interface (from that link the user can then access links to the other plates, therefore only one plate is added to the extra data)
- presets - list of presets that include only this library (found by view) - used to list these presets directly under the library in the interface (one of the requirements from the users)

### `preset_list()`

This view works similarly, processing some of the data that are not directly available in the database, and adding the results to the preset data. The extra processing needed here is:
- adding up compounds from all the subsets in the preset
- extracting just the list of ids of the subsets belonging to the preset 

The resulting data is:
- id = preset.id
- name = preset.name
- self.description = preset.description
- subsets - list of subset ids (used in the React app to differentiate between user-uploaded subsets, and those belonging to presets) - made to avoid downloading all the preset data
- size - the number of all compounds in the preset - displayed in the front end app; while the size of the single-library presets comes from the endpoint produced by `current_library_options()`, the sizes of presets spanning multiple libraries come from this endpoint


