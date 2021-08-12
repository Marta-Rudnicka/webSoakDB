#API

The API app provides views used by `webSoakDB_frontend` to load the data from the database. It consist of views, serializers and contains the `admin.py` file, which covers models used by all the apps. It is built using Django REST framework, and most of it is quite straightforward and will not be discussed in the documentation. It has a few trickier cases, which will be discussed here.

To access sensitive data, the API app uses `ispyb_dja` [link]

##Public vs. private library plates and the API views

Most of user-generated data is related to the `Project` model, which links to `IspybAuthorization` objects used for authentcation/authorization of the users. However the `LibraryPlate` model stores the data of both the public (XChem in-house) library plates and private ones (brought in by the user). The publicly available plates are not related to any project or user and need not be authenticates.

This is why there are completely different mechanisms used for loading library plate data.

When the user accessess a page listing all the compounds ( `SourceWell` objects) in a plate via the React app, the query is first directed to the `choose_plate_view()`. This view checks whether the plate is public or private and redirects the call to the appropriate URL. The view is passed two URL parameters: the plate id, and the project id, and it passes them on depending on which view is used next

###Public library plates: `PublicLibraryPlateViewSet` and `PlateWithCompoundsViewSet`

These are standard viewsets. `choose_plate_view()` gives them the plate id (through the URL parameters), but not project id since the plates it shows are not related to any particular project. `PublicLibraryPlateViewSet` provides a standard detail view of a library plate with data about the plate itself. Just in case someone tries to use this view to access private data, it filters out all the private plates (when `queryset` variable is declared). `PlateWithCompoundsViewSet` is used to get a list of compounds (`SourceWell` objects) in the plate. To prevent somehow using this view to inspect a private plate, the `get_queryset` function checks if the plate is public first and returns `None` if someone attempts to access a private plate with it.

###Private library plates and `ProjectCompoundsViewSet` 

A `LibraryPlate` instance  is not linked to any `IspybAuthorization` instance, so it is impossible to use `ispyb_dja` on this model. Instead, `ProjectCompoundsViewSet` provides an API endpoint that serves the data as an instance of `Project`. `Project` has a many-to-many field linking to `Library` instances selected to the project, and libraries are linked to `LibraryPlate` objects (private libraries have only one plate each). Thus, `ProjectCompoundsViewSet` loads the same data as `PublicLibraryPlateViewSet`, but in a more convoluted manner, following one foreign key after another until the structure is deep enough to access all the data necessary. There is no need to have a separate views for plate data and source well data, since the source well data is accessed via the plate data and the view provides both. The React component that use this endpoint filters out the unnecessary data and only display the ones related to the relevant plate (more details in the `webSoakDB_frontend` documentation)
