import React from 'react';
import ExportBar from './export_bar.js';
//import DataTable from './data_table.js';
import PlateList from './plate_list.js';
import TableHeader from './table_header.js';
import TableRow from './table_row.js';
//import {descriptor_names} from '../picker/stats_helper';
import axios from 'axios';
import {display_options} from './display_options.js';
import { deepCopyObjectArray, getAttributeArray, mean } from  '../../actions/stat_functions.js';


class PlateLookup extends React.Component {
	
	constructor(props) {
		super(props);
		this.toggleDisplay = this.toggleDisplay.bind(this);
		this.state = {
			compounds: [],
			subsets: [],
			display : {
				show_well: this.props.lookup_args.is_a_plate,
				show_library: this.props.lookup_args.is_a_preset,
				show_code: true, 
				show_smiles: true, 
				show_structure: false, 
				show_concentration: this.props.lookup_args.is_a_plate, 
				show_mol_wt: !this.props.lookup_args.is_a_plate, 
				show_tpsa: !this.props.lookup_args.is_a_plate, 
				show_logp: !this.props.lookup_args.is_a_plate,
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
	
	toggleDisplay(display_option){
		let displayCopy = Object.assign({}, this.state.display);
		displayCopy[display_option] = !this.state.display[display_option];
		this.setState({display : displayCopy});
	}
	
	componentDidMount() {
		this.uploadDataFromAPI()
		   		
	}
	
	componentDidUpdate(prevProps, prevState) {
		if (prevProps.lookup_args !== this.props.lookup_args) {
			uploadDataFromAPI()
		}
		
		if (prevState.subsets !== this.state.subsets){
			this.state.subsets.forEach(subset =>{
				this.uploadSubset(subset.id, subset.library.name);
			});
		}
	}
	
	uploadSubset(id, libName){
		const apiUrl = 'api/subset_stats/' + id + '/';
		
		const compounds = axios.get(apiUrl)
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
	
	uploadDataFromAPI() {
				
		//for library plate
		if (this.props.lookup_args.is_a_plate) {
			const library = this.props.lookup_args.library.name;
			const plate = this.props.lookup_args.collection.name;
			const apiUrl = 'api/compounds/' + library + '/' + plate;
			
			axios.get(apiUrl)
				.then(res => {
				const compounds = res.data;
				this.setState({ compounds });
		  });
		}
		//for a single cherrypicking list
		else if (!this.props.lookup_args.is_a_preset) {
			const compounds = this.uploadSubset(this.props.lookup_args.collection.id, null)
		}
		//for a preset
		else {
			const apiUrl = 'api/preset_detail/' + this.props.lookup_args.collection.id + '/';
			
			let subsets = [];
			
			axios.get(apiUrl)
				.then(res => {
				subsets = res.data.subsets;
				this.setState({subsets})
			});
		}
	}
	
	
	render() {
		let library = this.props.lookup_args.collection.library;
		if (this.props.lookup_args.is_a_preset){
			library = 'Preset';
		}
		
		const collection = this.props.lookup_args.collection;
		let current =""
		if (this.props.lookup_args.is_a_plate && this.props.lookup_args.is_a_preset){
			console.log('ERROR in this.props : both is_a_plate and is_a_preset are set to true!!! May cause unexpected behaviour')
		
		}
		
		if (this.props.lookup_args.collection.current){
			current = "(current)"; 
			}
		
		//Create show buttons for hidden columns. For subsets and presets, don't make show buttons for well and concentration.
		let buttons = display_options.map((option, index) => {
			if ((this.props.lookup_args.is_a_plate) || (option[0] !== 'show_well' && option[0] !=='show_concentration')){
				return <button key={index} className={this.state.display[option[0]] ? "hidden" : "small-button"} onClick={event => this.toggleDisplay(option[0])}>Show {option[1]}</button>
			}
		});
		
		let extra_button = null;
		if (this.props.lookup_args.is_a_preset){
			extra_button = <button key="0" className={this.state.display.show_library ? "hidden" : "small-button"} onClick={event => this.toggleDisplay(show_library)}>Show Library</button>
		}
		
		let plateList = null;
		if (this.props.lookup_args.is_a_plate){
			plateList = <PlateList library={library} showPlate={this.props.showPlate} />;
		}
		
		let rows = null;
		
		if (this.state.compounds.length > 0){
			rows = this.state.compounds.map((compound, index) => {
				return <TableRow
					key = {index +1}
					counter = {index + 1}
					compound={compound}
					display = {this.state.display}
					lookup_args = {this.props.lookup_args}
				/>
			});
		}
		
		return (
		<div id="plate-lookup">
				<h1>{library.name} </h1>
				<h2>{collection.name} {current}</h2>	
			<main>
				<div className="sidebar-div">
					<div>
					<h3>Show more table columns:</h3>
						{extra_button}
						{buttons}
					</div>
					
					<ExportBar />
					{plateList}
					
					
				</div>
				<table className="datatable">
					<caption>
						Compound list ({this.state.compounds.length} items)
					</caption>
					<TableHeader
						display = {this.state.display}
						onButtonClick = {this.toggleDisplay}
						lookup_args={this.props.lookup_args}
					/>
					<tbody id="datatable-body">
						{rows}
					</tbody> 
				</table>		
			</main>
		</div>
		);
	}
}

export default PlateLookup;

