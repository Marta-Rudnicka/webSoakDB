# Authorization issues

The information about in-house library inventory is public, while the staff-only inventory app is Django-only. Therefore most of the DRF views do not require authorization.
The only two views that do are:
- ProjectCompoundsViewSet()
- ProjectViewSet()

which provide the data relating to what the user selected for a specific project, and provide access to private library plate data (more info in the API docs).

Since `ISpyBSafeQuerySet` model inherits from DRF's  `ReadOnlyModelViewSet`, I first wrote these views as sub-classes of `viewsets.ReadOnlyModelViewSet`, debugged them and made sure they worked. Then, I switched `viewsets.ReadOnlyModelViewSet` to `ISpyBSafeQuerySet`, modified some fields according to the instructions in `ispyb_dja` and unsuccessfully tried to access `ProjectViewSet` the same way I did before.

## Issue 1 - I don't get any data from the API
So far I only tried to debug ProjectViewSet. After an ordinary, non-staff user logs into the app, they are redirected to a page that prompts them to select a project from the drop-down list (staff members go there after choosing to manage an experiment). ProjectViewSet provides the data for that list. If user tries to access any other page without selecting the project, they will be redirected back to that page. 

[The front-end code responsible for the site is in: 
`webSoakDB_frontend/src/components/main/proposal_selection.js` ;
to modify the URL used to access the data go to the `componentDidMount()` method ]

In the central system, there are two in-house visits I am authorized to access with my fedID. There are  IspybAuthorization instances with `proposal_visit` set to the visits to which I'm authorized in CAS, and there are Project instances, for which those IspybAuthorizations are added to the `auth` field

**Expected behaviour:**
When I log in with my FedID and go to the proposal selection page, the drop down list from which I choose my proposal should contain the two projects I am authorized to in the central system

**What happens**
I get an empty list or errors

**What I tried**
I ignored the user interface and tried accessing the API endpoint directly at: http://localhost:8000/api/projects/
The site was just showing an empty list:
># Project List

>**GET**  /api/projects/

>**HTTP 200 OK**
**Allow:** GET, HEAD, OPTIONS
**Content-Type:** application/json
**Vary:** Accept 
[]

I fiddled around with a few factors to see what works
- I tried using various configurations of filter_fields and URL query parameters
- I tried using changing Project.auth from many-to-many field to a foreign key to see what happens (the example SentitiveData class in ispyb_dja docs uses a foreign key field from IspybAuthorization)
- despite the docs saying to ignore the User field, I tried creating my local User instance with the username being my FedID and adding it to IspybAuthorization.users
With all the fiddling around, I either got empty lists or  "KeyError: 'ISPYB_USER'"


## Issue 2 - authorization of PATCH requests
When decides to save a selection of libraries (or presets) for an experiment, it is processed as a PATCH request going to the  `UpdateProjectSelection()` view. Such a request should be authorized, but the ISpyBSafeQuerySet model inherits from DRF's  ReadOnlyModelViewSet, which does not handle PATCH or any other HTTP request method that could be used for here. Therefore, I expect ISpyBSafeQuerySet does not either.  (So far, I have no way of testing it because I can't get beyond proposal selection page after I start using the authorization features). The most secret data at this stage, i.e. the contents of user's own libraries, is entered to the database using a back-end script, but still there should be no technical possibility to just POST, PUT or PATCH data into someone else's project. While the API in webSoakDB  mainly for browsing data, XChemSPA is mainly for producing and saving data, so it needs to handle a lot of PUTs and POSTs.

[For the front-end bit that sends the PATCH request, see `updateSelection()` method in:
`webSoakDB_frontend/src/components/App.js`  ]
