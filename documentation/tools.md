# `tools` directory

This directory stored helpers used in different apps, grouped by the functionalities.
[data_storage_classes.py](#dsc)
[uploads_downloads.py](#ud)
[histograms.py](#histograms)
[set_up.py](#setup)
[compounds.py](#compounds)
[validators.py](#validators)

## `data_storage_classes.py`<a name="dsc"></a>

This file contains classes based on the models used in the app, which are used to manipulate the data from the model instances without affecting the database. For example, `SourceWellCopy` is used to display the data stored in `SourceWell` model instances at a certain day in the past - the current `SourceWell` data is copied into `SourceWellCopy` attributes, and then `SourceWellCopy` is changed into its historical state without affecting the `SourceWell`s.

`SubsetCopyWithAvailability` not only stores the data, but also contains methods used in generating the extra data associated to each copy of `LibrarySubset`. More in details on this class in `inventory.md` under "Generating availability information".

## `uploads_downloads.py`<a name="ud"></a>

This file contains functions used for:
* parsing various CSV files from the user
* generating downloadable files for the user

Some of the functions are discussed in relevant sections of this documentation.

## `histograms.py`<a name="histograms"></a>

The React application presents histograms of molecular properties in a selected collection (library, preset, subset, whole user selection). The properties in question are computed from the SMILES string when the Compounds object is first created in the database, and they are stored as the attributes of the Compounds object. The graphs are stored or rendered as HTML/JavScript code by `bokeh` library ( [bokeh docs](https://docs.bokeh.org/en/latest/) ; [bokeh histogram example in docs](https://docs.bokeh.org/en/latest/docs/gallery/histogram.html))

`histogramps.py` contains the functions used for creating and updating histograms of molecular properties. While histograms themselves are only displayed in `webSoakDB_frontend`, the functions are used in other apps, because some of the inventory operations might require updating cached histograms.

### Main functions

- `get_histogram()`: used for producing a histogram on the fly in order to serve it in an HTTP response. It takes a database object, a string determining the type of the object ("library", "preset" or "subset"), and the string representing the attribute for which the distriution will be computed. It returns either the histogram's source code, or the integer 204. 204 is returned when there is no SMILES string stored for the Compounds object, and therefore the attributes describing the molecular properties have the value of `None`. 
- `update_histograms()`: produces histograms of all stored molecular properties for a collection and saves them in the `/media/` folder (overwrites old ones if they already exist). It is used by `webSoakDB_backend` and `inventory` views that modify the compound availability data to ensure that cached histograms are always up-to-date. Note: `bokeh` has dedicated functions used for saving the generated graphs as files, but it is not used here. The files generated that way cause CORS errors when embedded in the page. Instead, `upadate_histograms()` uses the same method of generating the source code as `get_histogram()` and saves it into the file using ordinary Python file operations.

- `get_selection_histogram()` : used for producing a histogram for the whole user selection on the fly (included collections that have not been saved in the proposal). In the process of creating a list of all the compounds involved, duplicated* compounds and compounds with no available property data are removed. The function accepts two lists of ints, which are ids of selected libraries and subsets, and a string representing the attribute to compute the graph for. It returns the source code of the graph as a string, (unless both lists it receives are empty and it produces a 'No compounds selected' message instead)

* duplicates could occur when user selects multi-library presets and does not notice there is some overlap with a selected library

## `set_up.py`<a name="setup"></a>

This file contains functions and lists of tuples used in creating a test data for unit tests (used in the `setUp` method in a `TestCase` instance).

## `compounds.py`<a name="compounds"></a>

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

## `validators.py`<a name="validators"></a>

This description applies to both the `validators.py` file in webSoakDB and XChemSPA. This module contains functions responsible for validating upload files and, mostly in XChemSPA, validating forms.

### The structure of the code validating uploaded files
For each aspect of the data validity, there is a separate helper function checking for it, and for each value, there is a function testing all the aspects of validity.

For example, if a well name must both follow a certain regex pattern and be unique in the file, there will be one function that checks for the pattern, one function that checks for uniqueness, and one function that calls both of them to decide whether the well name is valid.
On the next level, there are functions that check if the whole input row of a csv file is correct, by checking each value with the appropriate higher-level function. Those functions are then called on each row of the file within the function that validates the whole file (which also checks if, e.g., the file is in the correct format).
For some values, there is only one thing to check, in which case it directly checks for the only relevant property and does not call any lower-level functions.
Each of the functions described above returns a boolean informing whether the passed value, file line or file fulfills the validity requirements checked by it.
The general code structure pattern can be illustrated by the following  (simplified) pseudo-code:
```
file_is_valid(file)`:
	check_file_format(file)
	open_file_and_find_formatting(file)
	for each line in the file:
		row_is_valid(line)
		
row_is_valid(row):
	value_1_id_valid(row[0])
	value_2_id_valid(row[1])
	value_3_id_valid(row[2])

value_1_id_valid(value):
	check_the_only_relevant_property

value_2_id_valid(value):
	check_property_1(value)
	check_property_2(value)

check_property_1(value):
	do_something

check_property_2(value):
	do_something_else	

etc.
```
In the files themselves, the top-level functions are written first, and the lower-level functions follow below them.
### logging and reporting errors
A view which handles uploading a CSV file launches the relevant top-level validation function imported from `validators.py` . Each of those top-level validators is passed an empty list referred to as `error_log` throughout the whole code.  This list is then passed around between different-level helper functions, and may be modified during the process. The lower-level functions are also passed the number of the line from which the value they validate comes.
The `validators.py` files also include the function `update_error_log()`, which takes a list and a string, encloses the string inside the HTML `<p>` tags and appends it to the list.
Each of the lowest-level validation functions, upon encountering an invalid value, produces as string that informs what kind of error was encountered, on which line,  which value is at fault, what values are expected (where applicable) etc. This string, together with `error_log` , is passed to `update_error_log()`. After validating the whole file, `error_log` contains the record of all the error messages produced on the way. If the file is valid, `error_log` is empty and is not used for anything, while the view proceeds to create new records in the database. If any invalid values have been found, the `error_log` is passed to a special template, which displays all its collected messages. The user can then use those messages to fix the issues with the file and try uploading it again.
This pattern, however, is not followed in XChemSPA with functions that validate CSV files generated by the lab equipment: the messages saved in `error_log` are not as detailed. This is because in this case, there can be two causes for a file being invalid: either the user uploaded the wrong file by mistake, or there is an error in the XChemSPA itself that prevents it from handling the file correctly. In both cases pointing out every single value that is wrong in the file would not help the user fix the problem.
### Validating forms
Some functions in `validators.py` (mainly in XChemSPA) check the whole form, not just the file. This is more for the sake of testing the development process than validating the user input - in most cases invalid form would be caused by an error in the application.
