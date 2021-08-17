#`tools` directory

This directory stored helpers used in different apps, grouped by the functionalities.

## `data_storage_classes.py'

This file contains classes based on the models used in the app, which are used to manipulate the data from the model instances without affecting the database. For example, `SourceWellCopy` is used to display the data stored in `SourceWell` model instances at a certain day in the past - the current `SourceWell` data is copied into `SourceWellCopy` attributes, and then `SourceWellCopy` is changed into its historical state without affecting the `SourceWell`s.

`SubsetCopyWithAvailability` not only stores the data, but also contains methods used in generating the extra data associated to each copy of `LibrarySubset`. More in details on this class in `inventory.md` under "Generating availability information".

## `validators.py`

This file contains functions mainly used in validating CSV files uploaded by the user.

## `uploads_downloads.py`

This file contains functions used for:
* parsing various CSV files from the user
* generating downloadable files for the user

## `histogramps.py`

The React application presents histograms of molecular properties in a selected collection (library, preset, subset, whole user selection). The properties in question are computed from the SMILES string when the Compounds object is first created in the database, and they are stored as the attributes of the Compounds object. The graphs are stored or rendered as HTML/JavScript code by `bokeh` library ( [bokeh docs](https://docs.bokeh.org/en/latest/) [bokeh histogram example in docs](https://docs.bokeh.org/en/latest/docs/gallery/histogram.html))

`histogramps.py` contains the functions used for creating and updating histograms of molecular properties. While histograms themselves are only displayed in `webSoakDB_frontend`, the functions are used in other apps, because some of the inventory operations might require updating cached histograms.

### Main functions

- `get_histogram()`: used for producing a histogram on the fly in order to serve it in an HTTP response. It takes a database object, a string determining the type of the object ("library", "preset" or "subset"), and the string representing the attribute for which the distriution will be computed. It returns either the histogram's source code, or the integer 204. 204 is returned when there is no SMILES string stored for the Compounds object, and therefore the attributes describing the molecular properties have the value of `None`. 
- `update_histograms()`: produces histograms of all stored molecular properties for a collection and saves them in the `/media/` folder (overwrites old ones if they already exist). It is used by `webSoakDB_backend` and `inventory` views that modify the compound availability data to ensure that cached histograms are always up-to-date. Note: `bokeh` has dedicated functions used for saving the generated graphs as files, but it is not used here. The files generated that way cause CORS errors when embedded in the page. Instead, `upadate_histograms()` uses the same method of generating the source code as `get_histogram()` and saves it into the file using ordinary Python file operations.

- `get_selection_histogram()` : used for producing a histogram for the whole user selection on the fly (included collections that have not been saved in the proposal). In the process of creating a list of all the compounds involved, duplicated* compounds and compounds with no available property data are removed. The function accepts two lists of ints, which are ids of selected libraries and subsets, and a string representing the attribute to compute the graph for. It returns the source code of the graph as a string, (unless both lists it receives are empty and it produces a 'No compounds selected' message instead)

* duplicates could occur when user selects multi-library presets and does not notice there is some overlap with a selected library

## `set_up.py`

This file contains functions and lists of tuples used in creating a test data for unit tests (used in the `setUp` method in a `TestCase` instance).

## `compounds.py`

This file contains functions used in creating and searching `Compounds` objects.

**Background**

While testing the staging application on the full dataset used in the lab at the time, it turned out that the relationship between the compound code and the SMILES strings is more complicated thab I believed in the beginning. I assumed that one chemical compound (and consequently, one SMILES string) would have one code within one library. I did predict that there would be some instances where two different libraries happen to share a compound, and have different code for it, but that would not affect the compound searched performed by the application (that are almost always for a specific library).

It turned out that DSI-Poised library is quite inconsistent in naming and, SMILES string notation. After some investigation, it turned out that mainly there were cases where under one code, there were two different salts of the same compounds, and cases where different salts had different codes, which caused a lot of confusion. What is more, data supplied for different batches of the library, or plates using different solvent sometimes used a different convention for notating the SMILES string of the same compound, which at first made an impression that these were two different compounds.

Finally, some changes had to be introduced in uploading and locating compounds. At the time of writing this documentation, handling compounds works as following:

* when compound data is uploaded for any purpose, it is first processed with RDKit: the salt part is removed, and the rest of the string is canonicalized (converted to the same style of the notation). This is done by the `standardize_smiles()` function

* in the database model, a Compounds instance has a unique combination of SMILES string and code; it is caused by the initial idea of how compound codes are supposed to work, but should probably be changed in future. Objects that are the same chemical compound but have different codes are considered different instances od the Compounds models (as it was initially assumed they would belong to different libraries). In the lab practice, the compounds should be considered the same if they have the same (standardized) SMILES string and belong to the same library (because the information about the salt part is removed during the creation of the Compounds object, it follows that different salts will be considered the same compound for the purpose of the experiment)

* to satisfy the practical application, locating compounds in a particulat plate (for example for a cherry-picking list) has two stages:
	- 1. the application tries to find a source well linked to the instance of `Compounds` it is looking for
	- 2. if no such well is found, the application tries to look for a well with a differnt Compounds object, that has the same SMILES string
