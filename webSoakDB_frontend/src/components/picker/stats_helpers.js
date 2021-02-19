import { deepCopyObjectArray, getAttributeArray, mean, addUniqueCompounds } from  '../../actions/stat_functions.js';

export const descriptor_names = ["mol_wt", "tpsa", "log_p", "heavy_atom_count", "heavy_atom_mol_wt", "nhoh_count" , "no_count", "num_h_acceptors", "num_h_donors", "num_het_atoms", "num_rot_bonds", "num_val_electrons", "ring_count"]

const arrayNames = descriptor_names.map(name => name + "Array");

export const dict = descriptor_names.map(function(name, index) {
  return [name, arrayNames[index]];
});

export function get_stats(collection, colType, dictionary){
	let stats = {}
	let compounds;
	compounds = collection;
	
	const properties = getAttributeArray(compounds, "properties");
	collection.forEach(compound => {
		dictionary.forEach(item =>{
			stats[item[1]] = getAttributeArray(properties, item[0]);
			stats[item[0]] = mean(stats[item[1]]).toFixed(2);
		});
	});
	return stats;
}

export function updateAllSelection(libraryId, compounds, allSelection){
	allSelection.libraries.add(parseInt(libraryId));
	allSelection.compounds = addUniqueCompounds(allSelection.compounds, compounds);
}
