# `webSoakDB_frontend`

This is a reference for code included in the `webSoakDB_frontend` app. It is recommended to read the README.md first.

This app is a single-page ReactJS application used for selecting compounds for an experiment. It provides information about XChem in-house libraries and presets, an interface to select them for an experiment, an interface to upload user's own data. The React code is located in `webSoakDB_frontend/src`; most of it in `webSoakDB_frontend/src/components`.

# Pages:
(not actual pages, but what seems like traditional web pages to the user)
- Select/change proposal: a page where the proposal is selected. Links to any other page will redirect back proposal selection if no proposal has been chosen.
- Home: a home page with links; also contains links to experiment management application (inside <Home> component)
- Compound selection: provides a list of in-house libraries and presets, upload forms for user's own data, and a table showing distribution of selected properties in a selected collection of compounds (library, preset, subset or the whole selection) (inside <Picker> component)
- Summary: a table listing all selected item with a possibility to remove an item from selection. For items uploaded by the user, this is the only place to remove them (insisde <Summary> component)
- collection details: a table listing all the compounds in the selected library plate, library subset or preset

# File names and locations

The component containing the whole application is in `webSoakDB_frontend/src/components/App.js`. For each "page" in the application, there is a subdirectory in the  `components` directory, whith the components used in that page. The file with the capitalised name is the one containing the "top-level" component. Mostly, one file covers one React component and is named after that component. In some cases, a file may contain the code for both a component and its children components - in that case, the file is named after the parent component.

The `main` subdirectory contains child components used in <App>.

``webSoakDB_frontend/src/actions/stat_functions.js` contains some helper functions used in various components.

## Navigation

The app uses actual URLs for navigation (generated and managed by react-router (https://reactrouter.com/)). Most of the code responsible for navigation is located in App.js. <Link> components replace the HTML <a> elements for navigating inside the application, and <Route> and <PrivateRoute> components determine what components are rendered under which URLs. <Route> simply renders a component, while PrivateRoute only renders its component if user has selected the proposal to manage (otherwise it redirects to the home page).

All the <Link>s in the navbar fire the `confirmLeaving()` method. It asks user for confirmation before leaving a page if there are any unsaved changes made.

The `/compounds/.../.../` link leads to <PlateLookup> component, but it does not directly load a component. It loads <UrlTranslator> that determines the appropriate props based on the URL arguments, and it is <UrlTranslator> that loads <PlateLookup>.

If there is an <a> HTML element in the app, it means the link leads outside of the ReactJS application.

**React router vs. Django views**

All the URLs used by the single-page application are also defined in the `urls.py` file for webSoakDB_frontend. Thus navigating directly to any of the URLs renders the React app (instead of e.g. throwing a 404 error when the user uses the 'go back' function in the browser). Because one of the URL patterns uses arguments, while other patterns don't, there are two views defined in `views.py`. They both only render the page containing the React application, and the only difference is that one expects to get some URL arguments and the other doesn't.

## App.js

Apart from managing the navigation within the application, <App> holds the proposal data and tracks changes to it (information about the compound selection is stored in the Porposals object in the database). It also contains the method that saves selection in the database.
The details of changes to selection made since the last save are not stored in <App>. It only tracks whether any changes have been made for the sake of the `confirmLeaving()` method.

#### State
- `proposal` - stores the data of the proposal being managed
- `unsavedChanges` - a boolean that is true when what the user selected in the interface is different than the selection saved in the database

#### Methods

- `logIn()` - downloads the proposal data from the database #TODO: rename it
- `updateSelection()` - update Proposals.libraries and Proposals.subsets in the database; the arguments for this function are either sent from <Picker> or <Summary>; in case of failure to save data, informs user about it
- `refreshAfterUpload()` - forces <App> to reload proposal data; triggered when user uploads a library or a cherry-picking list
- `trackUnsavedChanges()` - sets `unsavedChanges` to the appropriate value; triggered inside <Picker> and <Summary> methods or by `updateSelection()`
- `confirmLeaving()` - produces a confirm window each time user tries to leave a sub-page without saving changes

## src/components/main/

### Icons.js

This file stores Bootstrap icons (https://icons.getbootstrap.com/) as reusable components (where the size and function to fire on click are passed to them as props).

### proposal_selection.js

This component contains the menu to select the proposal and send it to <App>.

### url_translator.js

Translates url arguments to props for <PlateLookup> (see: Navigation)

## src/components/home/

This directory only contains the hone page component.

## src/components/picker

This directory contains components used in the compound selection page. The main component loaded on this sub-page is <Picker>, which contains three main parts:
- An interactive list of available in-house libraries and presets that allows for selecting and unselecting them, and look up additional information about them (directly in the <Picker> component)
- a side panel with upload forms where user can upload own data or download the details of the whole selection as a CSV file (inside the <Uploads> component )
- a table showing distribution of molecular properties of selected collections of compounds as histograms (in the <GraphTable> component)

## <Picker> component (in Picker.js)

This component handles the selection of items made by the user.
#### State:
- `selectedLibIds` : an array of ids of libraries that are currently selected using the interface checkboxes; if there are already some libraries selected for the proposal when the page loads, these libraries are automatically 'ticked'
- `initialLibs`: an array of ids of libraries in the proposal; when the page loads, `selectedLibIds` and `initialLibs` are the same; when user checks or unchecks a checkbox for a library, `selectedLibIds` change while initialLibs stay the same

- `selectedSubsetIds`: an array of subset ids analogous to selectedLibIds; when user checks or unchecks a preset, all the subsets within the preset are added or removed from the array
- `initialSubsets`: analogous to initialLibs

When user saves a selection, the values of `selectedLibIds` and `selectedSubsetIds` are sent in the database and saved as Proposals.libraries and Proposals.subsets. The <App> component the reloads the proposal data, and the new `initialLibs` and `initialSubsets` values are set.

- `presets`: [] TODO

- `currentLibOptions`: : this attribute stores all the data about available libraries downloaded from the database, which is used to create the list of library options

- `inHouseCompoundCount`: the sum of compounds in selected libraries and single-library presets

#### Methods

- `componentDidMount()` - downloads the current data of the in-house libraries and presets
- `componentDidUpdate()` - makes sure to reset the library and subset lists stored in the state after saving a selection, track changes and update `inHouseCompoundCount`
- `handleChangeLib()` - adds and removes library ids from `selectedLibIds` when the user checks or unchecks a library
- `handleChangePreset()` - adds and removes subset ids from `selectedSubsetIds` when the user checks or unchecks a preset
- `saveChanges()` - sends data to <App>'s `updateSelection()` in order to update the proposal in the database
- `selected()` - checks if a preset based on whether one of its subset is selected 

**helpers**

- detectUnsavedChanges() - checks if there are any differences between libraries and presets currently selected in the interface, and the selection saved in the database (by comparing `selectedLibIds` against `initialLibs` and `selectedSubsetIds` against `initialSubsets` ); its output is sent to <App>'s `trackUnsavedChanges()`
- `updateInHouseCompoundCount()` - sums up compounds in selected libraries and compounds and updates the state

#### Selection area

Each available in-house library is represented by a <LibraryOption> component, and each preset is represented by a <PresetOption>. Presets with compounds from only one library are listed under that library; presets with more subsets are listed under "OTHER PRESETS". Each item for selection has a checkbox, name of the item, number of compounds in parentheses, and a <ChevronDown> icon. The number of compounds for libraries is the number of active compounds in the current plate(s). For presets, it is just the numbers of all compounds selected for the preset, since the library plate is not pre-determined. Items already selected for the proposal are pre-checked on load. User selects or rejects items by checking and unchecking corresponding checkboxes; and each instance of user checking or unchecking a box inside <LibraryOption> or <PresetOption> triggers `handleChangeLib()` or `handleChangePreset()` in <Picker>. These functions also trigger updating `inHouseCompoundCount`. Note: inHouseCompoundCount only sums up compounds from in-house libraries and single-library presets. Including multi-library presets would carry a high risk of accidentally duplicating compounds in the selection. The number would either be unreliable, or the application would need to download all the compound data to find possible duplicates.

The background colour of the <section> containing the selection options changes when there are unsaved changes detected (controlled by `changeStatus` variable in the <Picker>'s `render()` function). When user saves a new selection or undoes the changes, the colour goes back to the original one (currently the 'neutral' colour is grey, and the colour signalling unsaved changes is light blue.)

The <ChevronDown> icon in both <...Option> components triggers creating a <Details> component, which contains additional information about the item, a link to the list of compounds, and histograms showing the distribution of various molecular properties in the library or preset. <Details> element is closed by clickind on the <XIcon> inside the component. In both <LibraryOption> and <PresetOption>, creating and destroying this component is controlled by the `details` property in the state, with `showDetails()` method triggered by clicking on <ChevronDown>, and `hideDetails()` methods triggered by clicking on <XIcon>.

Each histogram showed in <Details> is a <Graph> component, which is an <iframe> containing the html/js code displaying the histogram. The component will be discussed in detail below.

#### <Uploads>
The Uploads area contains three forms:
- a form to upload user's own library plate map in <OwnLibraryForm> component
- a form to upload user's cherry-picking list from a particular library in <CherryPickForm> component
- a from that triggers downloading compound data as a CSV file in <ExportForm> component

It also incudes lists of already uploaded libraries and cherry-picking lists under the corresponding form, and a link to the static page with the detailed information on how the uploaded files should be formatted.

Each form uses a <CSRFToken> component to generate the token necessary for Django to accept the form.

##### <OwnLibraryForm>
Submits data to the `upload_user_library()` view in webSoakDB_backend.

**`submit()`**
This method overrides the form's submit action in order to preserve the information about unsaved changes. Uploading a new plate map creates a new `Library` object in the database and adds it to `Proposal.libraries`. This forces <App> to re-upload the proposal data (by triggering `refreshAfterUpload` in <App>) so the information shown by the interface is up-to-date - for example, the newly uploaded library will now show up under the form. The side effect of it is resetting the `unsavedChanges` state attribute in <App>, which is not a desirable behaviour if there are still some unsaved changes made in the selection area. Therefore the `submit()` method ensures that the value of `unsavedChanges` is restored to whatever it was before submitting the form and reloading the proposal data.

##### <CherryPickForm>

Submits data to the `upload_user_subset()` view in webSoakDB_backend.

**`componentDidMount()`**
Uploads a list of in-house library in order to created <option>s in the <select> input for the form. Cherry-picking lists are checked against the selected library to make sure the selection is correct.

**`submit()`**
Works the same as in <OwnLibraryForm>.

##### <ExportForm>

This form submits data to `export_selection_for_soakdb()` view, which produces a CSV file and starts a download. If the proposal includes only whole libraries, the form includes only one button and one hidden input submitting the proposal string (Proposals.proposal). The compound data are automatically copied fron the current plates.
If the proposal includes subsets (either from presets, or user-uploaded cherry-picking lists), the form also includes <SubsetSelect> components, which provide an interface to select library plates for each of the subsets. Each <SubsetSelect> component covers one library: it produces a <select> inputs with all available library plates as <options>. There are as many <select>s as current plates in the library, i.e. if a library has three current plates, the application assumes that the whole library takes up three plates, and produces three <select> elements (most libraries have just one). <SubsetSelect> class is defined in the same file as <ExportForm>.

The exported file only includes saved selection.

### <GraphTable>

This component provides user with histograms of selected molecular properties in selected collections of compounds (libriaries, presets, cherry-picking lists and the whole selection). The collections represented depend on what is currently selected in the front-end interface, without the need to save it in the database.

The main parts of the component are:
- a list of properties for which graphs are shown, with checkboxes to select/unselect them
- a table showing the graphs for all the currently selected collections
- (the last row of the table provides graphs for all the selected compounds together and works slightly diffrently from the other rows)

By default, there are five properties selected. For every library, preset and cherry-picking list in the user selection (as stored in <Picker>'s `selectedLibIds` and `selectedSubsetIds` ), a histogram showing distribution of each selected property is rendered in the table. Each time user adds or removes a collection from the selection (e.g. by checking a checkbox in the library list) or adds or removes a property to show, the graphs are re-rendered.

The graphs for the whole selection are not automatically rendered and the user needs to push a button to generate them. If the selection of compounds or properties change, the last row of the table is reset and graphs need to be regenerated. This is to avoid displaying outdated graphs and other errors.

#### State
- `show` - stores the list of properties for which the histograms are displayed; it is an array of strings which are the same as the names of corresponding columns in the database (or attributes of the Compounds models); the array is modified when user checks or unchecks the checkboxes in the properties list
- `presets` - the list of preset selected by user
- `other_subsets` - the list of subsets that are selected but do not belong to a preset (i.e. user's own cherry-picking lists)

#### Methods
- `componentDidUpdate()` - updates presets and other_subsets when the selection changes
- `manageProperties()` - adds or removes properties from the `state.show` when user checks or unchecks corresponding checkboxes
- `sortSubsets()` - finds out which subsets in `selectedSubsetIds` (sent as props from <Picker>) belong to presets, and which are user-submitted cherry-picking list; produces a list od preset ids (used to set `state.presets`), and list of id of subset that do not belong to presets (used to set `state.other_subsets`)

On rendering a <GraphTable>:
- a list of properies with checkboxes is rendered based on the `properties_dict` variable.  `properties_dict` is an object that maps short strings used as attribute name in the Compounds model to longer strings with the full name of the property.
- table headers are rendered based on what properties are currently selected and stored in `state.show`
- table rows with graphs are rendered for all the selected libraries (as <GraphTr> components)
- table rows with graphs are rendered for all the selected presets (as <GraphTr> components)
- table rows with graphs are rendered for all the selected cherry-picking lists (as <GraphTr> components)
- a table row for the whole selection is rendered (as <SummaryGraphs> component)

### <GraphTr>
A <GraphTr> component is a table row in the <GraphTable>. It is defined in the same file as <GraphTable>.
It creates a HTML table row element; for every property selected by the user (stored in <GraphTable>'s `show`), it renders a table cell with a <Graph> component in it.

#### State
- `collection` - the data of the compound collection for which the graphs are rendered
- `public` - stores the information on whether the collection is "public" (an XChem in-house library or a preset) or is something the user submitted; sometimes this information is passed throught props, and sometimes need to be determined after uploading the collection data from the database; this information is there to be passed to <Graph> components, where it is used

#### Methods

- `componentDidMount()` - uploads collection data from the database based on props passed to the component. In case of presets, the data is already available in props and no API call is made.


### <Graph>

Each <Graph> component is an iframe. The histograms used by the application are rendered as HTML and JavaScript code; the `src` of the iframe is either the URL of the media file where the graph is stored (using the default Django file storage system, in the `/media/` directory), or the URL of the endpoint generating the graph on the fly (using `serve_histogram()` view). In-house libraries and presets have their graphs stored, while everything else is generated on the fly. The `props.public` property determines which method is used to display the graph.
If a user uploads a plate map without SMILES string, it is impossible to compute the properties shown in the graphs. In that case `serve_histogram()` returns HTTP status code 204 (no content) and the component displays 'Data unavailable' message instead of the graph.

(The reason for caching some of the graphs was poor perfomance while generating them, the main bottleneck being the database query necessary to gather relevant data. If in the production environment this proves not to be a propblem, <Graph> component can be simplified: `props.public` and the cached filed in the /media/ folder will no longer be needed, and all the iframe sources should lead to the `serve_histogram()` view.)

### <SummaryGraphs>

This component produces the last row of the histogram table, where the user can generate histograms for all the components in the selection. These graphs require sending data to the back end, therefore the components rendering them work differently than the ones showing single collections.

When the component mounts, no graphs are rendered. When users clicks on the "Compute" button, in each data cell the short text is replaced with a <SelectionGraph> component. <SummaryGraphs> component can also render <Explanation> component, which works similarly to <Details> component in the libraries and presets lists, and diplays information on what is and is not taken into account when graphs for the whole selection are created.

<SelectionGraph> and <Explanation> are defined in the same file and <SummaryGraphs>

#### State
- `submit` : boolean triggering creation of <SelectionGraph>s and submission of form data to the back end
- `explanation` : bollean deciding whether <Explanation> info box should be rendered or not

#### Methods
- `showExplanation()`: causes rendering <Explanation> component
- `hideExplanation()`: hides <Explanation> component
- `componentDidUpdate()`: each time user changes the selection of libraries/presets of properties to display, this method sets `state.submit` to false. This means that if `state.submit` was `true` before and graphs were displayed, the graphs will disappear.
- `submit()`: sets `state.submit` to true, which triggers creation of new <SelectionGraph> components

### <SelectionGraph>
This component contains an iframe showing the output from the `selection_histogram()` view, and an invisible form sending the current selection data to that view. Form submission is triggered when <SelectionGraph> mounts.
#### Methods
- `submit_form()` - creates HTML attributes necessary for connecting the iframe and form to the same backend view and submits the data
#### `render()` variables
- `libs` and `subs` are values POSTed to `selection_histogram()`. They are the same as <Picker>'s `selectedLibIds` and `selectedSubsetIds` with one exception: if any of those values is an empty array, it is replaced with a 0.

## Summary page

This page lists selected and uploaded collections of compounds in tables. It also provides an interface to remove them from the selection.

### <Summary>

This is the component containing the whole subpage.

#### State
- `selectedLibs`: the same function as `selectedLibIds` in <Picker>, but stores objects instead of only ids
- `selectedSubsets`: the same function as `selectedLibs`but for subsets
- `libClass`: CSS class of the table listing the complete library, responsible for changing the background colour when there are some unsaved changes in the selection
- `subsetClass` : the same as `libClass` but for the table listing incomplete libraries

#### Methods
- `resetSelection()`: sets `selectedLibs` and `selectedSubsets` to the libraries and subsets that are currently saved in the database (if any collection was removed but unsaved, it will re-appear in the table after this method is triggered)
- `removeLibrary()`: removes library from `selectedLibs` (without affecting the database), sets `unsavedChanged` in <App>'s state to `true` and changes the background of the table with libraries
- `removeSubset()` performs analogous actions from preset
- `undoChanges()` resets the state to reflect the selection currently saved in the database
- `saveChanges()` triggers `updateSelection()` in <App> 

### <LibraryInTable>
Produces a table row(s) for one library in the 'Whole libraries' table. If there is more that one current plate in the library, it produces a row for each plate.

#### State
- `plates` all the current plates in the libraries

#### Methods
- `componentDidMount()`: downloads current plate data for `plates`

#### `render()` variables / parts of a table row

If there are more than one current plate in the library, some of the table cells are shared for all the plates (and span all the plate rows) and some are individual for each plate. The plate-specific cells are stored in the `plateCells` variable. A `plateCells` element then becomes a part of a variable storing the whole <tr>: either a `firstRow`, which contains the shared table cells, or a `nextRow`, which is only a <tr> containing `plateCells`.

## Collection lookup page

This page renders a table list all the compounds in the collections with details. The table columns can be shown and hidden using buttons. Since there are differences of what kind of data is available in what kind of collection, there are some differenced in how the lookup page is rendered depending on the collection.

###  <PlateLookup>
This is the top-level component containing the rest of the page.

#### Props
The props passed to <PlateLookup> contain the information about the type of the collection listed in the component. Based on these props, <PlateLookup> "knows" what API endpoint to use to download data and which elemnts in the page should be rendered. Also, the id of the collection is passed as props.

The props are created by <UrlTranslator> (in `/src/components/main`). This component is a child of the <Route> component from `react-router`. It uses URL parameters from <Route> and baased on them, decides what props to pass to its child <PlateLookup> element.

#### State 
- `collection`: stores all the data of the collection object uploaded from the database
- `compounds`: stores the data of the compounds in the collection
- `subsets`: stored id of subsets belonging to a preset; stay null if the collection is not a preset
- `display` : decides which columns are visible in the table; each attribute in this object is a boolean, and showing and hiding columns is achieved by modifying the attribute this object; e.g. if `display.show_code` is set to false, it hides the column containing the compound code and shows a button needed to show this column

#### Methods

- `toggleDisplay()`: changes the visibility of a selected column (and related buttons)
- `componentDidUpdate()`: makes sure to re-load the data when the collection changes (useful when user follows the link from one library plate to another); makes sure to upload compound data from the API after the list of subset is downloaded (used for presets) 
- `getDataFromAPI()`: (what the name indicates); because in each kind of collection the data is organised differently, each collection needs not only a different URL pattern for an API endpoint, but also a slightly different method of extracting the necessary data; that is why each collection type has a separate method managing the download process and `getDataFromApi()` only decides which function to use based on props; each method downloads the data of the collection itself, and data of the compounds in that collection
- `getPlateData()`: sets `state.compounds` to a list of active SourceWell objects from the plate (unlike other collections, which use Compounds objects instead); sets `state.collection.name` to the plate's barcode
- `getCherryPickingListData()`: sets `state.compounds` to the list of Compounds objects in the LibrarySubset, sets `state.collection.name` to the LibrarySubset's `origin` attribute
- `getPresetdata()` : first downloads the list of LibrarySubset objects related to the Preset, and then adds the compounds (Compounds objects) from each preset to `state.compounds`; sets `collection.library` to an empty object, `collection.library.name` to the name of the preset and `collection.name` to the description of the preset; semantically it makes little sense, but it causes the preset name and description show up in the same places where for other collections, the library name and collection name are displayed
- `getSubsetCompounds()`: helper function used in downloading preset and cherry-picking list compounds; downloads data of all Compounds objects in the subset and adds them to `state.compounds` (i.e. does not overwrite the existing array, but appends new elements to it); when called from `getPresetdata()`, it is passed the library name as an argument and it adds it as an extra attribute to each compound

#### Differences in rendering the interface for different kinds of collections

<PlateLookup> and its children render slightly different elements depending on the type of type of collection, which is one of the reasons why their `render()` methods contain so many conditional statements. The differences in what is rendered are listed below:

**Library plate**
- the table has extra columns with well name and concentration of the compound (and there are buttons to show or hide it)
- there is a <PlateList> component which provides links to other library plates in the same library
- there is an <ExportBar> element, which contains an interface to download the plate map of the library plate
- the objects in `state.compounds` are SourceWell objects, not Compound objects, so to access compound-specific attributes, the `compound` attribute of the source well has to be accessed first (e.g. a SMILES string for compound `compound` in `state.compounds` is stored as `compound.compound.smiles` instead of `compound.smiles`)

**Library subset**
No extra elements

**Preset**
- the table has a column with library name (and there are buttons to show and hide them)

### <ExportBar>

Download link to a CSV plate map of a library plate.
TODO: more formats and dowload options for other types of collections

### <PlateList>

Shows only for a library plate; downloads the list of plates in the library from the API, and creates links to each plate/

### <TableHeader>

Renders the <thead> element with its children. It renders <th> elements for all possible table columns, but based on the values passed inside the `props.display` object, it gives them classes that will either make them visible or not. Each <th> element except for "Row no." also has a <button> that hides the column in the table (by triggering `toggleDisplay()`)

### <TableRow>

Produces a table row with the data for one compound; for cells showing floats, the numbers are rounded. Visibiity for individual cells works the same as in <TableHeader>. The conditional statements in `render()` ensure that the cells are appropriate for the collection type, and that the component looks for the values in the appropriate "place" within the data structre (see: "Differences in rendering the interface for different kinds of collections")
[TODO: refactor this component or maybe the whole <PlateLookup> - too many conditional statements, look too complicated; maybe add SourceWell attributes to compound objects after uploading]

### <StructurePic>

A component containing the image of the 2D structure of the compound molecule and a magnifying glass icon (<ZoomInIcon> component). The image itself is generated by the `serve_2d()` view in `webSoakDB_backend` and rendered using <LazyLoadImage> component, which delays loading of the image until the user scrolls down to it (https://www.npmjs.com/package/react-lazy-load-image-component). When the cursor is placed over the image, it expands and the magnifying glass icon disappears (all implemented using CSS, including the animations).