from API import LibrarySubset, SourceWell, Compounds

def select_plates(subset):
    plate_selection = get_initial_plate_set(subset)
    reduced = remove_redundant_current_plates(plate_selection)
    plate_selection = reduced["new_plate_data"]
    potentially_redundant_plates = reduced["potentially_redundant"]
    plate_selection = remove_remaining_redundant_plates(potentially_redundant_plates, plate_selection)
    return plate_selection



def get_initial_plate_set(subset):
    '''find a reasonable set of plates that covers as much of the seletion as possible'''
    class PlateData():
        def __init__(self, id, compounds, current):
            self.id = id
            self.compounds = compounds
            self.current = current

    desired_compounds = set([c for c in subset.compounds.all()])
    
    #sort by the number of available locations, so you start search from the "rarest" compound
    sorted_compounds = sorted(desired_compounds, key=lambda c: len(c.locations.filter(active=True)), reverse=False)

    needed_plates = set()
    plate_data = []

    for c in sorted_compounds:

        #find plates where c is available
        c_plates = [sw.plate for sw in c.locations.all()]

        #if c is not available in any plate in needed_plates, add its plates to needed_plates
        if len(c_plates.intersection(needed_plates)) == 0:
            for p in c_plates:
                needed_plates.add(p)
                #save information about plate
                plate_compounds = set([sw.compound for sw in p.compounds.filter(active=True)])
                selected_from_plate = plate_compounds.intersection(desired_compounds)

                if selected_from_plate == desired_compounds and not p.current: #an old plate that covers all selection
                    plate_data = [PlateData(p.id, selected_from_plate, p.current))]
                    return plate_data

                plate_data.append(PlateData(p.id, selected_from_plate, p.current))
    
    return plate_data

def remove_redundant_current_plates(plate_data):
    '''Removes current plates that have no unique desired compounds, and saves the list of other potentially redundant plates. 
    It is preferable to use up old plates for subsets, so in the first run, current plates should be eliminated from the
    selection if possible'''

    potentially_redundant = []
    
    for p in plate_data:
        if is_redundant(p, plate_data): 
            if p.current:
                plate_data.remove(p)
            else:
                potentially_redundant.add(p)
    
    #sort by number of desired compounds, so the plates with the fewest compounds can be kicked out first
    potentially_redundant = sorted(potentially_redundant, key=lambda p: len(p.compounds), reverse=False)

    return {"new_plate_data" : plate_data, "potentially_redundant" : potentially_redundant }

def remove_remaining_redundant_plates(potentially_redundant, plate_data):
    '''removes any other redundant plates left'''
    for plate in potentially_redundant:
        if is_redundant(plate, plate_data):
            plate_data.remove(plate)
    
    return plate_data

def is_redundant(plate, plate_data):
    other_plates = plate_data.remove(plate)
    other_compounds = set()

    #get all compounds available in all the other plates
    for p in other_plates:
        other_compounds.add(p.compounds)
        
    #if p has no unique compounds...
    if plate.compounds.issubset(other_compounds):
        return True
    else:
        return False