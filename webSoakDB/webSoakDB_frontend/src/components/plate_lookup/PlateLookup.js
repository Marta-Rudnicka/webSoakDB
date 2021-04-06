import React from 'react';
import ExportBar from './export_bar.js';
import PlateList from './plate_list.js';
import TableHeader from './table_header.js';
import TableRow from './table_row.js';
import axios from 'axios';
import {display_options} from './display_options.js';
import { deepCopyObjectArray } from  '../../actions/stat_functions.js';
import {Link } from "react-router-dom";


class PlateLookup extends React.Component {
	
	constructor(props) {
		super(props);
		this.toggleDisplay = this.toggleDisplay.bind(this);
		this.state = {
			collection: null,
			compounds: [],
			subsets: [],
			display : {
				show_well: this.props.is_a_plate,
				show_library: this.props.is_a_preset,
				show_code: true, 
				show_smiles: true, 
				show_structure: false, 
				show_concentration: this.props.is_a_plate, 
				show_mol_wt: !this.props.is_a_plate, 
				show_tpsa: !this.props.is_a_plate,
				show_logp: !this.props.is_a_plate,
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
	
	toggleDisplay(display_option){
		let displayCopy = Object.assign({}, this.state.display);
		displayCopy[display_option] = !this.state.display[display_option];
		this.setState({display : displayCopy});
	}
	
	componentDidMount() {
		this.uploadDataFromAPI()
		   		
	}
	
	componentDidUpdate(prevProps, prevState) {
		if (prevProps !== this.props) {
			this.uploadDataFromAPI()
	}
		
		if (prevState.subsets !== this.state.subsets){
			this.state.subsets.forEach(subset =>{
				this.uploadSubset(subset.id, subset.library.name);
			});
		}
	}
	
	uploadSubset(id, libName){
		const apiUrl = '/api/subset_stats/' + id + '/';
		
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
		if (this.props.is_a_plate) {
			let apiUrl = '/api/compounds/' + this.props.id + '/';
			
			axios.get(apiUrl)
				.then(res => {
				const compounds = res.data;
				this.setState({ compounds });
			});
			
			apiUrl = '/api/plate_detail/' + this.props.id + '/';
			
			axios.get(apiUrl)
				.then(res => {
				const collection = res.data;
				this.setState({ collection });
			});
			
			
		}
		//for a single cherrypicking list
		else if (!this.props.is_a_preset) {
			const compounds = this.uploadSubset(this.props.id, null)
			
			const apiUrl = '/api/subset_stats/' + this.props.id + '/';
			
			axios.get(apiUrl)
				.then(res => {
				const collection = res.data;
				console.log('Setting collection')
				this.setState({ collection });
			});
		}
		//for a preset
		else {
			console.log('Found preset');
			const apiUrl = '/api/preset_detail/' + this.props.id + '/';
			
			let subsets = [];
			
			axios.get(apiUrl)
				.then(res => {
				const subsets = res.data.subsets;
				this.setState({subsets});
				const collection = res.data;
				this.setState({collection});
			});
		}
	}
	
	
	render() {
		let collection = null;
		let name = null;
		let current ="";
		if (this.state.collection){
			collection = this.state.collection;
			if (this.state.collection.library){
				name = collection.library.name; 
			}
			if (this.state.collection.current){
				current = "(current)"; 
			}
		}
		
		if (this.props.is_a_preset){
			name = 'Preset';
		}
		
		
		if (this.props.is_a_plate && this.props.is_a_preset){
			console.log('ERROR in this.props : both is_a_plate and is_a_preset are set to true!!! May cause unexpected behaviour')
		
		}
		
		
		
		//Create show buttons for hidden columns. For subsets and presets, don't make show buttons for well and concentration.
		let buttons = display_options.map((option, index) => {
			if ((this.props.is_a_plate) || (option[0] !== 'show_well' && option[0] !=='show_concentration')){
				return <button key={index} className={this.state.display[option[0]] ? "hidden" : "small-button"} onClick={event => this.toggleDisplay(option[0])}>Show {option[1]}</button>
			}
		});
		
		let extra_button = null;
		if (this.props.is_a_preset){
			extra_button = <button key="0" className={this.state.display.show_library ? "hidden" : "small-button"} onClick={event => this.toggleDisplay(show_library)}>Show Library</button>
		}
		
		let plateList = null;
		if (this.props.is_a_plate && collection && collection.library){
			plateList = <PlateList library={collection.library} />;
		}
		
		let rows = null;
		
		if (this.state.compounds.length > 0){
			rows = this.state.compounds.map((compound, index) => {
				return <TableRow
					key = {index +1}
					counter = {index + 1}
					compound={compound}
					display = {this.state.display}
					is_a_plate = {this.props.is_a_plate}
					is_a_preset = {this.props.is_a_preset}
				/>
			});
		}
		
		return (
		<div id="plate-lookup">
			
				<h1>{name ? name : ""} </h1>
				<h2>{collection ? collection.name : ""} {current}</h2>	
	
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
				<table data-toggle="table" data-pagination="true" data-search="true" className="table table-bordered table-hover" id="table">
					<caption>
						Compound list ({this.state.compounds.length} items)
					</caption>
					<TableHeader
						display = {this.state.display}
						onButtonClick = {this.toggleDisplay}
						is_a_preset={this.props.is_a_preset}
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

