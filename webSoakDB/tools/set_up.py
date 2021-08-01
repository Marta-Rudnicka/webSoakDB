
from API.models import (Library, 
    LibraryPlate, 
    LibrarySubset, 
    Compounds, 
    SWStatuschange, 
    SourceWell, 
    Preset, 
    Project,
    PlateOpening,
    IspybAuthorization
    )
import datetime

libraries_data = [("lib1", True, True), ("lib2", False, True), ("lib3", True, False)]
plates_data = [
    ("xyz", 0, True), ("xyz2", 0, False), ("xyz3", 0, False),
    ("abc1", 1, True), ("largest_small", 1, False), ("smallest_large", 1, False),
    ("large_row", 1, False), ("large_column", 1, False), ("empty", 1, False), ("abc", 2, True)
    ]
compounds_data = [
    ("code1", "CC"), ("code2", "CCO"), ("code3", "CCI"), ("code4", "CCF"), 
    ("code5", ""),  ("code6", ""),  ("code7", ""), 
    ("code8", "CCC"), ("code9", "CCCF"), ("code10", "CCOF"), ]
source_wells_data = [
    (0, 0, "A01", 30, True), (1, 0, "B01", 30, True), (2, 0, "C01", 30, True), (3, 0, "D01", 30, True),
    (0, 1, "A01", 30, False, datetime.date(2021, 2, 10)), (1, 1, "B01", 30, True), (2, 1, "C01", 30, False, datetime.date(2020, 10, 28)), (3, 1, "D01", 30, True), 
    (0, 2, "A01", 30, True), (1, 2, "B01", 30, True), (2, 2, "C01", 30, False), (3, 2, "D01", 30, True),
    (4, 3 ,"A01", 30, True), (5, 3, "B01", 30, True), (6, 3, "C01", 30, False),
    (7, 4, "P24", 30, True), (7, 5, "R25", 30, True), (7, 6, "R02", 30, True), (7, 7, "F25", 30, True),
    (7, 9, "X1", 30, True), (8, 9, "x2", 30, True), (9, 9, "X3", 30, True),
    ]

status_change_data = [
    (6, False, datetime.date(2020, 10, 28)), (10, False, datetime.date(2020, 10, 28)), 
    (4, False, datetime.date(2020, 10, 28)), (4, True, datetime.date(2020, 12, 2)), (4, False, datetime.date(2021, 2, 10)),
    ]

plate_opening_data = [ (1, "2020-10-28", "reason"), (2, "2020-10-28", "reason"),]
subsets_data = [("lib1-s1", 0, [1, 2]), ("lib1-s2", 0, [0, 3]), ("lib3-s1", 2, [7]), ("lib3-s2", 2, [7, 9]), ("test-missing", 0, [0, 2])  ]

auths_data = [("project-str1", "project1-1"), ("project-str2", "project2-2")]
projects_data = [(0, [0], [0, 1]), (1, [2], [])]

def set_up_libraries(data):
    libs = []
    for t in data:
        l =  Library.objects.create(name=t[0], public=t[1], for_industry=t[2])
        libs.append(l)
    return libs

def set_up_library_plates(data, library_data):
    libraries = Library.objects.all()
    if libraries.count() == 0:
        libraries = set_up_libraries(library_data)
    plates = []
    for t in data:
        p = LibraryPlate.objects.create(barcode=t[0], library=libraries[t[1]], current=t[2])
        plates.append(p)
    return plates

def set_up_compounds(data):
    compounds = []
    for t in data:
        c = Compounds.objects.create(code=t[0], smiles=t[1])
        compounds.append(c)
    return compounds

def set_up_source_wells(data, plates_data, library_data, compounds_data):
    plates = LibraryPlate.objects.all()
    if plates.count() == 0:
        plates = set_up_library_plates(plates_data, library_data)
    compounds = set_up_compounds(compounds_data)
    sws = []
    for t in data:
        sw = SourceWell.objects.create(compound=compounds[t[0]], library_plate=plates[t[1]], well=t[2], concentration=t[3], active=t[4])
        try:
            sw.deactivation_date = t[5]
            sw.save()
        except IndexError:
            pass
        sws.append(sw)
    return sws

def set_up_subsets(data, library_data, compounds_data):
    libraries = Library.objects.all()
    compounds = Compounds.objects.all()
    if libraries.count() == 0:
        libraries = set_up_libraries(library_data)
    if compounds.count() == 0:
        compounds = set_up_compounds(compounds_data)
    subs = []
    for t in data:
        l =  LibrarySubset.objects.create(name=t[0], library=libraries[t[1]])
        for index in t[2]:
            l.compounds.add(compounds[index])
            l.save()
        subs.append(l)
    return subs

def set_up_auths(data):
    for t in data:
        IspybAuthorization.objects.create(project=t[0], proposal_visit=t[1])

def set_up_projects(data, auth_data, library_data, subset_data, compounds_data):
    set_up_auths(auths_data)
    auths = IspybAuthorization.objects.all()
    subsets = LibrarySubset.objects.all()
    libraries = Library.objects.all()
    if subsets.count() == 0:
        subsets = set_up_subsets(subset_data, library_data, compounds_data)
    if libraries.count() ==0:
        libraries = Library.objects.all()
    projects = []
    for t in data:
        p =  Project.objects.create()
        p.auth.add(auths[t[0]])
        for index in t[1]:
            p.libraries.add(libraries[index])
            p.save()
        for index in t[2]:
            p.subsets.add(subsets[index])
            p.save()        
        projects.append(p)
    return projects

def set_up_status_changes(data, source_wells_data, plates_data, library_data, compounds_data):
    source_wells = SourceWell.objects.all()
    if source_wells.count() == 0:
        source_wells = set_up_source_wells(source_wells_data, plates_data, library_data, compounds_data)
    
    changes = []
    for t in data:
        sch = SWStatuschange.objects.create(source_well = source_wells[t[0]], activation=t[1], date=t[2])
        changes.append(sch)
    return changes

def set_up_openings(data, plates_data, library_data):
    plates = LibraryPlate.objects.all()
    if plates == 0:
        plates = set_up_library_plates(plates_data, library_data)
    
    openings = []
    for t in data:
        o = PlateOpening.objects.create(plate = plates[t[0]], date=t[1], reason=t[2])
        openings.append(o)
    return openings
