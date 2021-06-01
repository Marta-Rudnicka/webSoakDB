# INVENTORY

This is a reference for code included in the `inventory` app. It is recommended to read the README.md first.

`inventory` is a conventional Django app that uses templates to render pages, with no ReactJS components and only some minimal JavaScript.  Most of the views either render pages listing database items, or process form data to create or update them. Here, Django forms are used both to render the form HTML, and to sanitise/validate the data. The inventory has two more complicated features which are discussed below detail: tracking compound availability at different dates, and dispense testing interface.

### Files:
- `inv_helpers.py` - contains various helper functions utilised by this app’s views
- `dt.py` - contains helper functions used by the dispense testing interface

### Notes on views

- `libraries()` and `presets()` (urls: `inventory/plates` and `inventory/presets`) - these views list all the in-house libraries / presets in a table. For every item, two forms are rendered: one for editing, and the other for deleting the item. These forms are by default hidden and are only shown after clicking the “edit/delete” button (achieved by a simple JavaScript script manipulating the CSS classes of the forms). Each form contains a hidden input submitting the id of the item that it modifies. Besides, both pages have one form each for creating a new item.

- `presets()` (url: inventory/presets - the preset list includes the lists of all compounds belonging to that preset. However, the view does not simply display the values stored in the database. For every Compounds object referenced in the LibrarySubset belonging to a preset, the view checks the availability of the compound in the current library plate. If it is not available, it searches the other plates of that library, and lists the plates in which the compound is still available. The relevant values copied from Compounds objects, as well the availability data, are stored as a list of temporary objects and provided as the context to the page template.
**Note on updating creating or presets:** the form creating a preset only allows for adding compounds from one library (adding one LibrarySubset). To create a preset covering multiple libraries, user first needs to create the preset with one subset, and the edit it to add more. If a preset already contains a selection from a particular library and the user wants to update this selection, the old selection is completely overwritten by the new one. The old LibrarySubset is first deleted from the database, then new one is created based on the uplodaded file and added to the Preset. E.g. if user wants to add three compounds to the list, the new uploaded list needs to contain all the old compounds + the three new ones.

- `plates()` - this view lists all library plates for XChem in-house libraries and provides links to pages processing particular plate. Above the table with the plates list, there is a `<select>` element that allows for filtering the plates by library -  this is achieved using a simple JavaScript, which is in the template file.

### `Advanced features:
	
#### Tracking plate usage across time
	
Each library plate in the list at `inventory/plates` has a link labelled “track usage”, which allows ‘travelling back in time’ to see which compounds in a plate were available at a given date. The page  with the information also contains a visual representation of the library plate.

- `track_usage()` (url: `inventory/track-usage/<int:pk>/<str:date>/<str:mode>/`)- this view copies values from the relevant SourceWell objects into new temporary objects, sets their `active` attribute to how it was on the inspected date (based on the records stored in SWStatuschange objects), and feeds the modified data into the template.

url args:
- pk: the id of the tracked plate
- date: the date to be inspected in the format yyyy-mm-dd
- mode: the desired layout of the page, the argument can be either “general-view”, which shows a various additional information about the plate, and “graphic-view”, which only provides the visual representation and buttons to navigate to other dates
	
E.g. `inventory/track-usage/8/2021-02-12/graphic-view/` shows what compounds were available on 12th February 2021 in plate with the id of 8, in the mode showing only the visual representation.

#### Notes on the visual representation:

**Appearance:**
The library plate is represented as an HTML table where green cells represent wells with available compounds, and red cells represent wells where the unavailable compounds used to be. The cells also display the name of the well. If a well is not included in the plate map, the cell representing it is empty. When users places the cursor over over a cell, the cursor disappears and the cell expands to make the well name easier to read - this is implemented in CSS.

**The HTML table vs. library plates (the physical objects)**
There are three types of library plates used in the lab; two of them (384 LDV and 384PP) have 384 wells organized in 24 columns (1-24) and 16 rows (A-P), and one (1536LDV) has 1536 wells organised into 48 columns (1-48) and 32 rows (A-AF). Therefore, the HTML table representing the plate can have two sizes with the corresponding number of rows and columns. The well names are combinations of the row and column names, e.g. ‘A13’, ‘AB01’. Note: in case of 384LDV, the actual physical layout of the wells is slightly different from what the table looks like, so it should be treated more like a schematic diagram of than an attempt to mimic the appearance of the plate (unlike dispense testing interface).

The view determines the size of the table to use based on well names: if it finds a name that is too “high” for a 384-well plate like S32, it renders the table with 1536 wells/cells. This means that if the library plate is physically stored in a 1536LDV plate, but only happens to use wells available in 384-well plates (i.e. the top left quarter of the plate), it will be represented as a 384-well table anyway.

**How the representation is rendered:**
The view renders an empty HTML table of the appropriate size, where every cell has the well name as its id. In a hidden <section>, it also renders a `<div>` element representing each `SourceWell` object belonging to the plate, with the class determined by whether it is available or not, and the well name stored in 'data-well' HTML attribute. On page load, the `arrange()` JavaScript function places each `<div>` representing a `SourceWell` in the corresponding table cell by matching the `data-well` attribute with the cell `id`. The function is in the template itself. Cells representing unused wells stay unchanged.

**The dates**
The link to the `track_usage()` page in the plate list (produced by the `plates()` view) directs to the page showing the current state of the plate in question. The date inspected is the last time the availability of any of the compounds in the plate (`SourceWell` objects) was updated (based on the last_tested attribute of the plate). From the tracking page, there are two ways of moving onto a different date:
-a drop-down list of dates - this produces a selection of all dates on which there were changes in availability recorded
-a date input allowing for selecting an arbitrary date
When a date is selected by any of those input elements, it launches the `redirect()` script (included in the template), which redirects to the appropriate URL for his date and the same display style as the previous page

**Display modes (general-view vs. graphic-view):**
In both modes, the rendered HTML is almost the same. The <main> element in the page has a different id depending on the mode, and there are different CSS rules written for <main> and its child elements depending on id (mostly, some of the elements have ```display: none``` attribute in the graphic-view). The mode also affects the URLs to which user is redirected and the link to switching modes.

## DISPENSE TESTING

Dispense testing is a process in which each compound in a library plate is dispensed onto a crystallisation plate in order to check if it still gets dispensed or if the plate has run out of it. The inventory app provides a visual interface to easily record missing compounds. The user marks empty wells on the image of a destination plate, and based on the submitted file, the application finds out what compounds were supposed to go there, and marks those compounds as unavailable.  There are two views involved in the process:

- `dispense_testing_map()` - a view that processes the file that that maps source wells (the wells in the library plate) to destination wells (the wells in the crystallisation plate into which those compounds are dispensed) and then renders the interface for marking compounds as unavailable
- `deactivate_compounds()` - a view that processes the data submitted through the dispense testing page


### Mapping source wells to destination wells
The page rendered by `update_plate()` contains a form in which the user submits the mapping file. The data from the form is submitted to the `dispense_testing_map()` view.

#### File format and well naming conventions in the crystallisation plate:
The submitted file must be a CSV file where the relevant columns have the headers “source well” and “destination well” (the script detecting the headers is not case sensitive and discards leading and trailing white space, so a name like “ Destination Well “ is also recognised). The order of the columns or the presence of other columns does not matter. Thus, **the dispense testing feature can re-use the file that is used by Echo dispenser for the same test**.
The script accepts two naming conventions of the destination wells: one that is used in the Echo file, and one that follows the marking that are physically on the plate (that I will call “human readable” here, and that are described that way in the source code too). Thus, if there is ever a need to manually create the file, or if dispensers used start recognising the human-readable naming system, there will be no need to manually convert names into the Echo system to create an input file for the dispense testing feature.

Dispense testing uses only one type of crystallisation plate: SWISSCI 3 Lens Crystallisation Plate (specs and photos available here: https://swissci.com/wp-content/uploads/2020/03/3-Lens-plate.pdf) 
In the plate, the wells are organised into groups of three with one reservoir, like in the ASCII art below (where a, c, and d are wells, and r is a resevoir).
	
```
----------
|(a) [r] | 
|        |
|(c) (d) |
----------
```
	
These groups are placed in 12 columns (marked 1-12) and 8 rows (marked A-H). The markings “a” “c” and “d” are not on the plate itself, but are commonly used. Well names are composed of the row name, column name and the letter a, c, or d specifying the position in the group, e.g. B7c, F12a or G3d.
The naming used in the Echo input files does not take groups into account and considers each column and row separately, which make 24 columns (1-24) and 16 rows (A-P), giving names such as G16, P3 etc. The conversion between the Echo naming system and the “human readable” naming system follows a pattern like this:
	
```
A1a --> A1  --------------------------
A1c --> B1  | (A1a / A1)    [r / A2]|
A1d --> B2  |                       |
A2a --> A3  | (A1c / B1) (A1d / B2) |
A2c --> B3  |-----------------------|
```
	
etc.

After the file upload, the script creates a dictionary matching source well to destination plate. If the uploaded file uses Echo-style well names, they get converted to “human-readable” names first.
	

#### The interface to mark unavailable compounds
Once the dictionary mapping the wells is created, `dispense_testing_map()` renders a page with the dispense testing interface. The main element of the page is a form rendered as an HTML table representing the crystallisation plate, which roughly mimics its layout and appearance, with individual wells represented by checkboxes, and reservoirs represented by an icon of an ‘x’ in a square (they do not “do” anything, they are there to provide a closer visual representation of the plate). The checkboxes representing inactive SourceWells are checked on load.

	
In case of errors in the upload file, error messages are displayed instead of the input form. 

During the test, users inspect the destination plate to check for wells into which no compound was dispensed. If they find such a well, they mark it on the plate with a marker. Then, in the interface, they check all the boxes representing the marked wells and submit the data. In rare cases, a compound that was not dispensed in the previous test may appear - in that case users uncheck the pre-checked box. The image below shows a photograph of a fragment of the real plate next to the data from that plate correctly entered into the form:
As the compounds are selected for deactivation, they also appear in a table listing them, including the code and source plate (JS script in the template). There is also a list of compounds that were already unavailable.

#### Rendering the form:
First, the table is rendered with <label>s but no checkboxes: each group of 3 wells and a reservoir is a <div> inside a table cell. Each label’s id is the “human-readable” style name of the well. The checkboxes are rendered separately, with HTML `id`s containing the name of the appropriate destination well, and the `data-` attributes containing information about the compound. The `arrange()` script, which launches on page load, places them in the appropriate positions in the table/form based on ids of the label and the checkbox (script in the template). The form also includes a hidden input field with a list of all the ids of the already inactive SourceWells.

#### Saving the data:
After submitting the form, the data is directed to `deactivate_compounds()`. The data includes:
* plate id
* ids of all the SourceWells from which the compounds were not dispensed ('not dispensed')
* ids of all the SourceWell from which compounds were not expected to be dispensed, i.e. which were already inactive before the test ('already missing')
For every SourceWell from 'not dispensed', if it is active, its `active` attribute is set to `False`, `deactivation_date` is set the current date, and a new SWStatuschange is created*. If it is already inactive, it is ignored.
Then the application checks if all the compounds from 'already missing' list are present in the 'not dispensed' list. If not (i.e. if any inactive SourceWell turned out to still have some compound in it), they are marked back as active, the `deactivation_date` is set to `None`, and a SwStatuschange is created*. 

Note: the application cannot simply compare the 'not dispensed' list to the list of all the inactive SourceWells in the plate taken directly from the database. In case of a larger library, it takes two or three crystallisation plates to test all the compounds. A SourceWell may be absent from the 'not dispensed' list not because it was dispensed, but because it was not tested at all in this part of the test.

**Note on status changes:**
The application will not create more than one SWStatuschange on the same day. If, on the same day, a SourceWell is set to active, and then it is set back to inactive or vice versa, the application treats it as an act of correcting a mistake, and undoes the last change. This mean that the last (today's) SWStatuschange is deleted, and `deactivation_date` as well as `active` are restored to the state from before both changes. This avoids ambiuity in tracking the changes over time. If the same change were made twice in a row (which not should happen in normal running of the application), the second one would just be ignored. 

In the end, the the `last_tested` attribute of the plate is changed to the current time stamp, and id the plate tested is a current plate, histograms of the corresponding library are updated. The user is redirected back to 'update/delete' plate page.
