import React from 'react';

export const properties_dict = {
	"mol_wt" : "Molecular weight" ,
	"tpsa" : "TPSA" ,
	"log_p" : "LogP" ,
	"num_val_electrons" : "Number of valence electrons" ,
	"num_h_acceptors" : "Number of hydrogen bond acceptors" ,
	"num_h_donors" : "Number of hydrogen bond donors" ,
	"num_het_atoms" : "Number of heteroatoms" ,
	"num_rot_bonds" : "Number of rotable bonds" ,
	"ring_count" : "Number of rings" ,
	"heavy_atom_count" : "Number of heavy atoms" ,
	"heavy_atom_mol_wt" : "Average molecular weight of heavy atoms" ,
	"nhoh_count" : "Number of NH and OH" ,
	"no_count" : "Number of nitrogens and oxygens" ,
}

export const properties_list = [
	"mol_wt", 
	"tpsa", 
	"log_p", 
	"num_val_electrons", 
	"num_h_acceptors", 
	"num_h_donors", 
	"num_het_atoms", 
	"num_rot_bonds",
	"ring_count", 
	"heavy_atom_count",
	"heavy_atom_mol_wt",
	"nhoh_count",
	"no_count"]

export default properties_dict
