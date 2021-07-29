import React from 'react';
import axios from 'axios';
import CompoundLookupCherryPick from './CompoundLookupCherryPick';
import TableRowPreset from './table_row_preset';
import ExportBar from './export_bar';

class CompoundLookupPreset extends CompoundLookupCherryPick {
	
	constructor(props) {
		super(props);
		this.toggleDisplay = this.toggleDisplay.bind(this);
		this.state = {
			collection: null,
			compounds: [],
			subsets: [],
			display : {
				show_well: false,
				show_library: true,
				show_code: true, 
				show_smiles: true, 
				show_structure: true, 
				show_concentration: false, 
				show_mol_wt: true, 
				show_tpsa: true,
				show_logp: true,
				show_logp: false,
				show_heavy_atom_count: false,
				show_heavy_atom_mol_wt: false,
				show_nhoh_count: false,
				show_no_count: false,
				show_num_h_acceptors: false,
				show_num_h_donors: false,
				show_num_het_atoms: false,
				show_num_rot_bonds: false,
				show_num_val_electrons: false,
				show_ring_count: false,
			}
		};
	}
	
	getDataFromAPI() {
			
		const apiUrl = '/api/preset_detail/' + this.props.id + '/';
			
		let subsets = [];
		axios.get(apiUrl)
			.then(res => {
			const subsets = res.data.subsets;
			this.setState({subsets});
			let collection = res.data;
			collection.library = {};
			collection.library.name = collection.name;
			collection.name = collection.description;
			
			this.setState({collection});
		});
		
		subsets.forEach(s => {
			compounds = this.getSubsetCompounds(s.id, s.library.name)
		})
	}

	getRows(){
		let rows = <tr><td colSpan="6" className="large-text">Loading compounds...</td></tr>;
		if (this.state.compounds.length > 0){
			rows = this.state.compounds.map((compound, index) => {
				return <TableRowPreset
					key = {index +1}
					counter = {index + 1}
					compound={compound}
					display = {this.state.display}
				/>
			});
		}
		return rows;
	}

	getExtraButton(){
	    return (
	    	<button 
		    	key="0" 
				className={this.state.display.show_library ? "hidden" : "small-button"}
				onClick={event => this.toggleDisplay(show_library)}
			>
				Show Library
			</button>
		);
	}

	getExportBar(){
		return <ExportBar url="preset" id={this.props.id} label= "compound list"/>;
	}

}

export default CompoundLookupPreset;

