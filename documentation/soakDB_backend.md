# `webSoakDB_backend` reference

This is a reference for code included in the `webSoakDB_backend` app. It is recommended to read the README.md first.

## Files:
- `helpers.py` - contains various helper functions used in `views.py` and `forms.py`
- `validators.py` - contains functions that validate the CSV files uploaded by users in the `webSoakDB_frontend` app
- `histograms.py` - contains functions responsible for producing histograms of molecular properties

## Data uploads:

Using the ReactJS frontend, users can upload the plate maps of their own fragment libraries or upload their cherry-picking lists. Users submit some data and CSV files through a form. While the page itself doesn’t use use Django `Form` objects to render the HTML forms (because these are rendered using ReactJS, not Django templates), `Form` objects are still used by the view to take advantage of their in-build functionalities of validating and sanitising the data.
After basic form validation, the CSV file goes through the validation process using functions from `validators.py`.  (SMILES strings are validated using RDKit). The validating function also create an error log: each  time an error is encountered in the file, the function adds an appropriate string to the log. If the upload is valid, new objects are created in the database and the user is redirected to a page informing about a successful upload. If not, no changes in the database are made, and the user is redirected to a page that prints out all the error messages from the error log.

The information on accepted data formats are available on the "formatting help" page (URL: `/uploads/formatting_help/`, template: `websoakDB_backend/templates/webSoakDB_backend/formatting-help.html`).

### library uploads
Uploads are handled by the `upload_user_library()` view. When a user uploads own library data:
  - new `Library` object is created
  - new `LibraryPlate` object pointing the the new `Library` is created
  - new `SourceWell` objects pointing to the new `LibraryPlate` are created
While creating a new `SourceWell` object, first the function searches the database for a `Compounds` objects with matching compound code and SMILES string. If it is found, the `SourceWell` object points to that `Compounds` object; otherwise, a new `Compounds` object is created first. Thus, there should be no duplicate `Compounds` objects in the database (the same functions are used for uploading in-house library plates, where it is more of an issue).

### cherry-picking list uploads
Uploads are handled by the `upload_user_subset()` view, creating a new `LibrarySubset` object. As mentioned before, if the subset or library plate are created by an ordinary user in the compound selection interface, references to the created `Library` and `LibraryPlate` objects are added to the relevant `Proposal` object.
Note: the function validating the cherry-picking list not only checks if the formatting is correct, but also whether the compound belongs to the selected library (to prevent errors).

__Notes on views:__
- `export_selection_for_soakdb()`
This function generates a downloadable CSV file formatted to be easily copied into SoakDB spreadsheet. The rules of copying the data from `SourceWell` and `Compounds` objects are the same as the ones described in  *Inventory data vs. experimental data* section, the only difference is that instead of copying the values into newly created `SPACompound` objects, the function prints them into a CSV file.
- `serve_2d`: serves a PNG image of a 2D molecular structure generated using RDKit; the image is generated based on SMILES string

#### Histograms
The application produces histograms of molecular properties in a selected collection (library, preset, subset, whole user selection). The properties in question are computed from the SMILES string when the Compounds object is first created in the database, and they are stored as the attributes of the Compounds object. The graphs are stored or rendered as HTML/JavScript code by `bokeh` library. There are three main functions creating historams in `histograms.py`, and two views serving them.
__non-view functions__
- `get_histogram()`: used for producing a histogram on the fly in order to serve it in an HTTP response. It takes a database object, a string determining the type of the object ("library", "preset" or "subset"), and the string representing the attribute for which the distriution will be computed. It returns either the histogram's source code, or the integer 204. 204 is returned when there is no SMILES string stored for the Compounds object, and therefore the attributes describing the molecular properties have the value of `None`. 
- `update_histograms()`: produces histograms of all stored molecular properties for a collection and saves them in the `/media/` folder (overwrites old ones if they already exist). It is used by `webSoakDB_backend` and `inventory` views that modify the compound availability data to ensure that cached histograms are always up-to-date. Note: `bokeh` has dedicated functions used for saving the generated graphs as files, but it is not used here. The files generated that way cause CORS errors when embedded in the page. Instead, `upadate_histograms()` uses the sam method of generating the source code as `get_histogram()` and saves it into the file using ordinary Python file operations.
- `get_selection_histogram()` : used for producing a histogram for the whole user selection on the fly (included collections that have not been saved in the proposal). In the process of creating a list of all the compounds involved, duplicated* compounds and compounds with no available property data are removed. The function accepts two lists of ints, which are ids of selected libraries and subsets, and a string representing the attribute to compute the graph for. It returns the source code of the graph as a string, (unless both lists it receives are empty and it produces a 'No compounds selected' message.

*duplicates could occur when user selects multi-library presets and does not notice there is some overlap with a selected library

__views__
- `serve_histogram()`: uses `get_histogram()` to produce a histogram of a collection on the fly and serves it as a HTTP response. If histogram is impossible to create because there are no SMILES strings and molecular properties unavailable, it returns HTTP status code 204. This view is used for user's own uploads only; the graphs for in-house libraries and presets are cached in the `/media/` folder and accessed via URLs produced with Django media storage features
- `selection_histogram()`: uses `get_selection_histogram()` to produce and serve the graphs for user's whole selection based on submitted form data; the form data contains the molecular property in question and the ids of selected libraries and library subsets stored in the state of the front end application