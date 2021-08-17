# `webSoakDB_backend` reference

This is a reference for code included in the `webSoakDB_backend` app. It is recommended to read the README.md first.

## Files:
- `helpers.py` - contains various helper functions used in `views.py` and `forms.py`

## Data uploads:

Using the ReactJS frontend, users can upload the plate maps of their own fragment libraries or upload their cherry-picking lists. Users submit some data and CSV files through a form. While the page itself doesn’t use use Django `Form` objects to render the HTML forms (because these are rendered using ReactJS, not Django templates), `Form` objects are still used by the view to take advantage of their in-build functionalities of validating and sanitising the data.
After basic form validation, the CSV file goes through the validation process using functions from `tools.validators.py`.  (SMILES strings are validated using RDKit). The validating function also create an error log: each  time an error is encountered in the file, the function adds an appropriate string to the log. If the upload is valid, new objects are created in the database and the user is redirected to a page informing about a successful upload. If not, no changes in the database are made, and the user is redirected to a page that prints out all the error messages from the error log.

The information on accepted data formats are available on the "formatting help" page (URL: `/uploads/formatting_help/`, template: `websoakDB_backend/templates/webSoakDB_backend/formatting-help.html`).

### library uploads
Uploads are handled by the `upload_user_library()` view. When a user uploads own library data:
  - new `Library` object is created (where `public=False`)
  - new `LibraryPlate` object pointing the the new `Library` is created
  - new `SourceWell` objects pointing to the new `LibraryPlate` are created
While creating a new `SourceWell` object, first the function searches the database for a `Compounds` objects with matching compound code and SMILES string. If it is found, the `SourceWell` object points to that `Compounds` object; otherwise, a new `Compounds` object is created first. Thus, there should be no duplicate `Compounds` objects in the database (the same functions are used for uploading in-house library plates, where it is more of an issue).

### cherry-picking list uploads
Uploads are handled by the `upload_user_subset()` view, creating a new `LibrarySubset` object. If the subset or library plate are created by an ordinary user in the compound selection interface, references to the created `Library` and `LibraryPlate` objects are added to the relevant `Proposal` object.
Note: the function validating the cherry-picking list not only checks if the formatting is correct, but also whether the compound belongs to the selected library (to prevent errors).

__Notes on views:__
- `export_selection_for_soakdb()`
This function generates a downloadable CSV file formatted to be easily copied into SoakDB spreadsheet. The rules of copying the data from `SourceWell` and `Compounds` objects are the same as the ones described in  *Inventory data vs. experimental data* section, the only difference is that instead of copying the values into newly created `SPACompound` objects, the function prints them into a CSV file.
- `serve_2d`: serves a PNG image of a 2D molecular structure generated using RDKit; the image is generated based on SMILES string

#### Histograms

Some histograms are cached in the `/media` folder and displayed in iframes, and some are generated on the fly by one of the two views below. To learn more about how the histograms are produced, see the documentation for `histograms.py` in `tools.md`.

__views__
- `serve_histogram()`: uses `get_histogram()` to produce a histogram of a collection on the fly and serves it as a HTTP response. If the histogram is impossible to create because there are no SMILES strings and molecular properties unavailable, it returns HTTP status code 204. This view is used for user's own uploads only.

- `selection_histogram()`: uses `get_selection_histogram()` to produce and serve the graphs for user's whole selection based on submitted form data; the form data contains the molecular property in question and the ids of selected libraries and library subsets stored in the state of the front end application
