# `webSoakDB_backend` reference

This is a reference for code included in the `webSoakDB_backend` app. It is recommended to read the README.md first.

## Introduction
This app functions as a back-end to the `webSoakDB_frontend` and mainly handles user uploads and downloads. In some cases, `webSoakDB_frontend` links to `inventory` views though.

## Files:
- `helpers.py` - contains various helper functions used in `views.py` and `forms.py`
- `forms.py` - standard Django file; the only notable thing here is the fact that the forms are used only for validation
## Data uploads:

Using the ReactJS front end, users can upload the plate maps of their own fragment libraries or upload their cherry-picking lists. Users submit some data and CSV files through a form. While the page itself doesn’t use use Django `Form` objects to render the HTML forms (because these are rendered using ReactJS, not Django templates), `Form` objects are still used by the view to take advantage of their in-built functionalities of validating and sanitizing the data.
After basic form validation, the CSV file goes through the validation process using functions from `tools/validators.py`.  (SMILES strings are validated using RDKit). The validating function also creates an error log: each  time an error is encountered in the file, the function adds an appropriate string to the log. If the upload is valid, new objects are created in the database and the user is redirected to a page informing about a successful upload. If not, no changes in the database are made, and the user is redirected to a page that prints out all the error messages from the error log.

The information on accepted data formats are available on the "formatting help" page (URL: `/uploads/formatting_help/`, template: `websoakDB_backend/templates/webSoakDB_backend/formatting-help.html`).

When the user uploads a library plate or a subset, it is automatically added to the user's `Project`.

### library uploads
Uploads are handled by the `upload_user_library()` view. When a user uploads own library data:
  - new `Library` object is created (where `public=False`)
  - new `LibraryPlate` object pointing the the new `Library` is created
  - new `SourceWell` objects pointing to the new `LibraryPlate` are created
While creating a new `SourceWell` object, first the function searches the database for a `Compounds` objects with matching compound code and SMILES string. If it is found, the `SourceWell` object points to that `Compounds` object; otherwise, a new `Compounds` object is created first. Thus, there should be no duplicate `Compounds` objects in the database (the same functions are used for uploading in-house library plates, where it is more of an issue).

### cherry-picking list uploads
Uploads are handled by the `upload_user_subset()` view, creating a new `LibrarySubset` object. 
Note: the function validating the cherry-picking list not only checks if the formatting is correct, but also whether the compound belongs to the selected library (to prevent errors).

__Notes on views:__
- `export_selection_for_soakdb()`
This function generates a downloadable CSV file formatted to be easily copied into SoakDB spreadsheet. The rules of copying the data from `SourceWell` and `Compounds` objects are the same as the ones described in  [*Inventory data vs. experimental data*](https://github.com/Marta-Rudnicka/webSoakDB/blob/main/documentation/README.md#invvsexp) section, the only difference is that instead of copying the values into newly created `SPACompound` objects, the function prints them into a CSV file.
- `serve_2d()`: serves a PNG image of a 2D molecular structure generated using RDKit; the image is generated based on SMILES string

## Data downloads

There are three types of collection that are shown in detail in the selection app: library plate, library subset, and preset. For each type, the user can download two lists of compounds: one with basic information (code and smiles, and if available, well and concentration), and one with all the molecular properties produced from the smiles string. Thus, there are 6 similar views whose job is to produce and serve the CSV file with the data: `download_current_plate_map()`, `download_plate_map_with_properties`, `download_subset()`, `download_subset_with_properties()`, `download_preset()`, `download_preset_with_properties()`. Each view works the same way:
1. It retrieves the desired collection from the database
2. It produces the file header (there are some helper assisting for that)
3. It generates a descriptive filename for the download
4. It passes everything to `serve_csv_compound_list()`, which produces the final HTTP response returned by the view

### `tools.uploads_downloads.serve_csv_compound_list()`
This view produces a CSV file and returns it in a HTTP response. First, it adds the header to the file, and then, depending on the type of collection it is passed, it launches an appropriate function to produce the data rows.  There are different function for each collection type, because the available data is a bit different: only library plate produce well names, and only presets provide library names for the compounds.

The `include_details` boolean, passed first to `serve_csv_compound_list()`, and then to the lower-level helpers, determines whether the molecular properties is added to the file or not.

While all the data rows are generated, the CSV file is closed and sent in a request.

### `export_selection_for_soakdb()`
This function generates a downloadable CSV file formatted to be easily copied into SoakDB spreadsheet. The rules of copying the data from `SourceWell` and `Compounds` objects are the same as the ones described in  *Inventory data vs. experimental data* section, the only difference is that instead of copying the values into newly created `SPACompound` objects, the function prints them into a CSV file.

#### Histograms
The application produces histograms of molecular properties in a selected collection (library, preset, subset, whole user selection). The properties in question are computed from the SMILES string when the Compounds object is first created in the database, and they are stored as the attributes of the Compounds object. The graphs are stored or rendered as HTML/JavScript code by `bokeh` library, and most of the functions related to them are in `tools/histograms.py` [(docs)](https://github.com/Marta-Rudnicka/webSoakDB/blob/main/documentation/tools.md#histograms)
__views__
- `serve_histogram()`: uses `get_histogram()` to produce a histogram of a collection on the fly and serves it as a HTTP response. If the histogram is impossible to create because there are no SMILES strings and molecular properties are unavailable, it returns HTTP status code 204. This view is used for user's own uploads only; the graphs for in-house libraries and presets are cached in the `/media/` folder and accessed via URLs produced with Django media storage features
- `selection_histogram()`: uses `get_selection_histogram()` to produce and serve the graphs for user's whole selection based on submitted form data; the form data contains the molecular property in question and the ids of selected libraries and library subsets stored in the state of the front end application
## Others


 - `serve_2d()`: serves a PNG image of a 2D molecular structure generated using RDKit; the image is generated based on SMILES string
