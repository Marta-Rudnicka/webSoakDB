import React from 'react';
import CompoundLookupPlate from './CompoundLookupPlate.js';
import TableRowCherryPick from './table_row_cherry_pick.js';
import ExportBar from './export_bar.js';
import axios from 'axios';
import {display_options} from './display_options.js';
import { deepCopyObjectArray } from  '../../actions/stat_functions.js';

class CompoundLookupCherryPick extends CompoundLookupPlate {
	
	constructor(props) {
		super(props);
		this.toggleDisplay = this.toggleDisplay.bind(this);
		this.state = {
			collection: null,
			compounds: [],
			subsets: [],
			display : {
				show_well: false,
				show_library: false,
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
	
	componentDidUpdate(prevProps, prevState) {
		if (prevProps !== this.props) {
			this.getDataFromAPI()
	}
		
		if (prevState.subsets !== this.state.subsets){
			this.state.subsets.forEach(subset =>{
				this.getSubsetCompounds(subset.id, subset.library.name);
			});
		}
	}
	
	getDataFromAPI() {
			
		this.getSubsetCompounds(this.props.id, null)
			
		const apiUrl = '/api/subset_detail/' + this.props.id + '/';
		
		axios.get(apiUrl)
			.then(res => {
			let collection = res.data;
			collection.name = collection.origin
			this.setState({ collection });
		});
	}

	getSubsetCompounds(id, libName){
		const apiUrl = '/api/subset_stats/' + id + '/';
		
		axios.get(apiUrl)
			.then(res => {
				const oldCompounds = deepCopyObjectArray(this.state.compounds);
				let compounds = res.data.compounds;
				if (libName){
					compounds.forEach(compound => compound.library = libName);
				}
				compounds.push(...oldCompounds);
				this.setState({ compounds: compounds });
		});
	}

	getDisplayButtons() {
		const buttons = display_options.map((option, index) => {
			if (option[0] !== 'show_well' && option[0] !=='show_concentration'){
			return (
				<button 
					key={index} 
					className={this.state.display[option[0]] ? "hidden" : "small-button"} 
					onClick={() => this.toggleDisplay(option[0])}>Show {option[1]}
				</button>
			);}
		});
		return buttons;
	}

	getPlateList(){
		return null;
	}

	getExportBar(){
		return <ExportBar url="subset" id={this.props.id} label= "compound list"/>;
	}

	getRows(){
		let rows = null;
		if (this.state.compounds.length > 0){
			rows = this.state.compounds.map((compound, index) => {
				return <TableRowCherryPick
					key = {index +1}
					counter = {index + 1}
					compound={compound}
					display = {this.state.display}
				/>
			});
		}
		return rows;
	}
	
}

export default CompoundLookupCherryPick;