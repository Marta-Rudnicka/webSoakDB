# XChemSPA webSoakDB - documentation

This directory contains an explanation of how WebSoakDB works and how it relates to the lab inventory in real life. It is meant to help navigate the codebase and give some context that will make it easier to understand and maintain the code. Basic Django features and files like `views.py` are not explained here - only things characteristic to this specific application.

This file contains the general information about the structure of the application, the apps and the ORM models used in it. It should serve as a starting point for understanding the rest of the documentation. 

Other files discuss various parts of the application by app (in the meaning of a Django app) and are meant to be reference for particular features.

## Introduction: 
webSoakDB is a part of XChemSPA. XChemSPA is a system for managing XChem experiments from the moment of selecting the compounds until the crystals are mounted and ready for data collection. On the back-end, it is composed of two Django/ReactJS applications: webSoakDB and XChemSPA. webSoakDB is used to manage the compounds, and XChemSPA is used to manage the experiments themselves.

webSoakDB is composed of two main parts: 
- the compound selection part available to XChem users with a fedID and a registered proposal. It provides the interface to browse fragment libraries provided by XChem, upload the data of users’ own libraries, and save the selection to be later used in the experiment

-inventory management part available only to XChem staff members. It provides an interface for the database of compounds, libraries, library plates etc, with tools to track the usage of library plates and easily update compound availability

### Technologies: 
webSoakDB is based on cookiecutter-xchem-stack template (see the template documentation for the list of technologies included): https://github.com/xchem/cookiecutter-xchem-stack

In addition to the tools included in the template, it also uses:
    • RDKit cheminformatics tools: https://www.rdkit.org/
    • Django REST Framework for creating API: https://www.django-rest-framework.org/
    • bokeh Python library for creating graphs: https://docs.bokeh.org/en/latest/index.html
Note on the CSS: the application uses Bootstrap, but some of the Bootstrap rules are overridden, and there is a lot of custom CSS involved

### List of Django apps in the project
- `webSoakDB_frontend` - a single-page ReactJS application for selecting compounds for an experiment and browsing the available XChem in-house libraries; contains only one view responsible for launching the ReactJS app. This is the application ordinary users interact with when they select compounds for their experiments.
- `API` - contains models, serializers and views for the API consumed by webSoakDB_frontend
- `webSoakDB_backend` - contains views and helper functions responsible for handling uploads and downloads of files in webSoakDB_frontend; it also includes some views redirecting users after login and views rendering minor static pages (e.g. file formatting guide for users); some of the functions from webSoakDB_backend are re-used in inventory
- `webSoakDB_stack` - contains configuration and setting for the whole project
- `inventory` - a backend Django application with a little plain (no ReactJS) JavaSctipt that allows the staff to manage the inventory; all the staff-only functionalities are included in this app

### ORM MODELS AND HOW THEY RELATE TO THE ITEMS IN THE LAB 

**Library**
Represents a fragment library but contains no information about its contents. A library should be thought of as a collection of library plates with some metadata about them.
If a fragment library is available in more than one solvent, each version is treated as a separate library, e.g. DSI Poised plates using DMSO and DSI Poised plates using ethylene glycol belong to two different libraries
Notes on attributes:
- `public`: if `public==True`, the library is an XChem in-house library; it is visible to all academic users and can be managed through the inventory management app; if `public==False`, it means it was created by a user and is only visible to that user and staff members; it does not show up in the inventory application along with in-house libraries; the attribute is set automatically when the library is created
- `for_industry`: if `for_industry==True`, it is suitable for industry users and shows up as an option to select in the selection application; all user-submitted libraries have this attribute set to `True`
- `plates`: a list of `LibraryPlate` objects with a foreign key pointing at the given instance of Library (implemented as a `related_name` in `LibraryPlate`)

**LibraryPlate:**
Represents a particular physical library plate. 
Notes on attributes:
- `current`: if `current=True`, the plate is the one that is currently used for experiments. Once too many compounds in the plate are used up, a new plate is opened and used as the current plate. The old plate is still stored in the lab and can be occasionally used, e.g. for cherry-picking compounds. If a library is very large, it may not physically fit into one plate, in which case there will be more than one plate where `current=True`. All user-provided libraries are automatically set to current; for in-house libraries, the attribute is manually set in the inventory management application
- `compounds`: lists `SourceWell` objects associated with a plate (implemented as a `related_name` in `SourceWell`)

**Compounds:**
Represents a chemical compound belonging to a library (it has no foreign key reference to a Library object, though); contains the data of the compound, but not its location in any library plate. A Compounds object has a unique combination of code and SMILES attribute (i.e. if two libraries use the same chemical compound under different codes, it will be treated as two separate Compounds objects; the same with two different SMILES strings under the same code)
Notes on attributes:
- `code`: a string identifying the compound used by the library supplier
- `smiles`: a SMILES string, can be left blank (industry users often keep it secret)
- molecular properties attributes: computed using RDKit from the SMILES string during data upload; if SMILES string is not provided, they stay null
- `locations`: list of `SourceWell` objects pointing to this instance of the `Compounds` object (implemented as a `related_name` in `SourceWell`) - in other words, all the places where the compound can be found

**SourceWell:**
An object that links a particular compound to a particular physical location (a well in a library plate). It includes foreign keys to `Compounds` and `LibraryPlate` objects, the name of the well where the compound is located, current availability, and (optional) concentration.
Notes on attributes:
- `active` - if `active=True`, there is still some compound left in that particular location and it can be used in an experiment; when it runs out, active should be set to `False`. Note: the location of a compound can differ between plates belonging to the same library

**LibrarySubset:**
Represents a collection of compounds (as `Compounds` objects) belonging to the same library; used to store information about user’s cherry-picked lists of compounds or preset selections.
Notes on attributes:
- `origin`: an automatically generated string informing the user what kind of subset is it; used when displaying the details of the user’s selection

**Preset:**
A selection of compounds created and saved by XChem staff for a specific purpose and available to users. It can contain one on more library subsets.

**Proposal:**
Stores the data about the proposal. Before doing anything in the compound selection interface, the user needs to select the proposal that will be managed. Then all the selected libraries, presets and uploads saved are recorded in the Proposal object related to it.
Notes on attributes:
- `libraries`: contains references to all the libraries selected for the proposal. If the user uploads the data of their own library, it is automatically added to this list
- `subsets`: contains references to all the subsets selected for the proposal. If a user uploads a cherry-picking list of compounds belonging to an in-house library, it is automatically added to this lists. If a user selects a preset, all the subsets included to this preset are added to this list. Preset objects themselves aren’t references by the Proposal object.

**SWStatuschange**
Source well status change. A `SWStatuschange` object is created each time a SourceWell is set to from active to inactive (or the other way round in rare cases). It is used for tracking compound usage over time.

**PlateOpening**
An instance of opening a library plate (removing plastic film that sits on top of the wells). A PlateOpening object is autimatically created when a plate undergoes dispense test or is used in an experiment (TODO), but it can also be manually created.

#### Collections of Compounds objects vs. collections of SourceWell objects

Some elements of the user interface are the same for libraries, presets and cherry-picking lists, but the compoiund data for these collections are processed differently. 
Any kind of statistics for "a library" are really statistics for the current plate(s) of that library, since this is the plate that would be used in the experiment. The compounds involved in calculations are Compounds objects referenced by all the active SourceWell objects belonging to the current plate(s).

Statistics for "presets" or "cherry-picking lists" are different because there is no default plate to use. When using only a selection of compounds from a library, the lab staff may decide to use one of the older plates and save the current one for full library screens. Accordingly, the collections of compounds used both in presets and cherry-picking lists are LibrarySubset objects, which are not related to any particular LibraryPlate, so they directly reference Compounds objects (instead of referencing SourceWell objects). Since without knowing the plate, it is impossible to tell which compounds are available, all the compounds in LibraryPreset objects are taken into account in the calculations.

In features that require specifying the location of the compound, the application requires the user to select the library plate first, and then it looks for the SourceWells in that plate that contain the selected compounds.


#### Inventory data vs. experimental data [not implemented yet]**

Inventory objects are not directly used while managing an experiment. When an experiment is started, the user “imports” the selection. In that process, the up-to-date data of the selected compounds is copied into `SPACompound` objects, and these are the ones that are later used. The copying rules are:
- for user’s own library data, the relevant values are just copied from what the user has uploaded
- for an in-house library, the data is taken from the current plate(s). The import process finds all active SourceWell objects belonging to the current LibraryPlate(s)The relevant data from the  `SourceWell` objects as well as `Compounds` objects related to them are copied into `SPACompound` objects
- for a cherry-picking list from an in-house library, the user needs to select the library plate that is going to be used. For a preset, the user has to choose the plate for every library included in the preset. Then, the import process finds which 'SourceWell' objects in that plate contain the desired compounds, and copies the data from those objects.

## Terms used in this documentation:
- **plate map** : a file storing the locations (wells) of particular compounds within a particular library plate
- **cherry-picking list** : sometimes users do not want to use the whole library in an experiment, but only select some of the compounds from it; a cherry-picking list is a list of desired compounds, or a file containing this list
- **collection** : any library object that references a set of compounds: a LibraryPlate, a LibrarySubset, a Preset, or Proposals object

# The `API` Django app
(There is not much to explain about this app, therefore it does not have a separate file)

This API is not a general-purpose web API to interact with the XChem inventory database; it is written to serve the specific needs of `webSoakDB_frontend` application. The views in `views.py` are mostly class-based to take advantage of the generic view classes provided by Django REST Framework(DRF) (https://www.django-rest-framework.org/api-guide/generic-views/) but there are also a couple of views that do not use DRF and provide the data using a JsonResponse (in case of these particular datasets, not using DRF was simpler). All views in `views.py` are web API enpoints. `serializers.py` contains DRF serializers, which is how code is typically organised when working with DRF.

Note on DRF endpoints: when accessed directly from the browser, they display JSON data in a very readable way, which can be useful in debugging.

[TODO - set appropriate permissions views]