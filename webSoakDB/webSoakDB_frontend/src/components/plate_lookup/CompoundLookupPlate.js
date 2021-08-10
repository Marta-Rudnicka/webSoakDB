import React from 'react';
import ExportBar from './export_bar.js';
import PlateList from './plate_list.js';
import TableHeader from './table_header.js';
import TableRowPlate from './table_row_plate.js';
import axios from 'axios';
import {display_options} from './display_options.js';

class CompoundLookupPlate extends React.Component {
	
	constructor(props) {
		super(props);
		this.toggleDisplay = this.toggleDisplay.bind(this);
		this.state = {
			collection: null,
			compounds: [],
			subsets: [],
			display : {
				show_well: true,
				show_library: false,
				show_code: true, 
				show_smiles: true, 
				show_structure: true, 
				show_concentration: true, 
				show_mol_wt: false, 
				show_tpsa: false,
				show_logp: false,
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
		this.getDataFromAPI();	
	}
	
	componentDidUpdate(prevProps, prevState) {
		if (prevProps !== this.props) {
			this.getDataFromAPI()
		}
	}
	
	getDataFromAPI() {
		let apiUrl = '/api/library_plates/' + this.props.id + '/' + this.props.project_id + '/';
		
		axios.get(apiUrl)
			.then(res => {
			//public and private plate data are taken from different API endpoints with 
			//different data structures
			if (res.data.barcode){
				this.managePublicPlateData(res.data)
			}
			else{
				this.managePrivatePlateData(res.data)
			}
			//const collection = res.data;
			//collection.name = collection.barcode;
			//this.setState({ collection });
		});
		
		apiUrl = '/api/plate_compounds/' + this.props.id + '/';
			
		axios.get(apiUrl)
			.then(res => {
			const compounds = res.data;

			//private plates don't receive any new data at this point
			if(res.data){ 
				this.setState({ compounds });
			}
		});
		
		
	}

	managePublicPlateData(data){
		const collection = data;
		if (! collection.name){
			collection.name = collection.barcode;
		}
		else {
			collection.name = collection.name + ' (' + collection.barcode + ')';
		}
		this.setState({ collection });
	}

	managePrivatePlateData(data){
		console.log('fired managePrivatePlateData with: ', data)
		//find the desired library in the project data
		//(private libraries have only one plate, so only the first plate 
		//in each libraries.plates array needs checking)
		const lib = data.libraries.find(library => library.plates[0].id === parseInt(this.props.id))
		let collection = lib.plates[0]
		collection.name = collection.barcode;
		let compounds = collection.compounds
		this.setState({ collection, compounds });
	}

	getDisplayButtons() {
		const buttons = display_options.map((option, index) => {
			return (
				<button 
					key={index} 
					className={this.state.display[option[0]] ? "hidden" : "small-button"} 
					onClick={() => this.toggleDisplay(option[0])}>Show {option[1]}
				</button>
			);
		});
		return buttons;
	}

	getExtraButton(){
		return null; //to be overwritten in sub-classes
	}

	getPlateList() {
		if (this.state.collection){
			return <PlateList library={this.state.collection.library} />;
		}
		else {
			return null;
		}
	}

	getExportBar() {
		return <ExportBar url="plate-map" id={this.props.id} label= "plate map"/>;
	}

	getPageHeaders(){
		let collection = null;
		let name = "Loading data...";
		let current = null;
		if (this.state.collection){
			collection = this.state.collection;
			name = collection.library.name;
			if (this.state.collection.current){
				current = "(current)"; 
			}
		}
		return {collection : collection, name : name, current : current}
	}

	getRows(){
		let rows = <tr><td colSpan="6" className="large-text">Loading compounds...</td></tr>;
		if (this.state.compounds.length > 0){
			rows = this.state.compounds.map((compound, index) => {
				return <TableRowPlate
					key = {index +1}
					counter = {index + 1}
					compound={compound}
					display = {this.state.display}
				/>
			});
		}
		return rows;
	}

	render() {
		const headerStrings = this.getPageHeaders();
		const collection = headerStrings.collection;
		const name = headerStrings.name;
		const current = headerStrings.current;
		const buttons = this.getDisplayButtons();
		const extra_button = this.getExtraButton();
		const plateList = this.getPlateList();
		const export_bar = this.getExportBar(); 
		const rows = this.getRows();
		
		
		
		return (
		<div id="plate-lookup">
			
				<h1>{name} </h1>
				<h2>{collection ? collection.name : ""} {current}</h2>	
	
			<main>
				<div className="sidebar-div">
					<div>
					<h3>Show more table columns:</h3>
						{extra_button}
						{buttons}
					</div>
					{export_bar}
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

export default CompoundLookupPlate;
