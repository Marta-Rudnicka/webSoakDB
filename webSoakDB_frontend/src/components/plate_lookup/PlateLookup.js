import React from 'react';
import ExportBar from './export_bar.js';
import DataTable from './data_table.js';
import PlateList from './plate_list.js';
import TableHeader from './table_header.js';
//import {descriptor_names} from '../picker/stats_helper';
import axios from 'axios';
import {display_options} from './display_options.js';



class PlateLookup extends React.Component {
	
	constructor(props) {
		super(props);
		this.toggleDisplay = this.toggleDisplay.bind(this);
		this.state = {
			compounds: [],
			show_well: true, 
			show_code: true, 
			show_smiles: true, 
			show_structure: false, 
			show_concentration: true, 
			show_mw: true, 
			show_tpsa: false, 
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
		};
	}
	
	toggleDisplay(display_option){
		this.setState({[display_option]: !this.state[display_option]});
	}
	
	componentDidMount() {
		const library = this.props.library.name;
		const plate = this.props.plate.name;
		const apiUrl = 'api/compounds/' + library + '/' + plate;
		
		axios.get(apiUrl)
			.then(res => {
			const compounds = res.data;
			this.setState({ compounds });
      });     		
	}
	
	componentDidUpdate(prevProps, prevState) {
		if (prevProps.lookup_args !== this.props.lookup_args) {
			this.componentDidMount()
		}
	}
	
	
	render() {
	//	console.log('this.props: ', this.props)
		const library = this.props.plate.library;
		const plate = this.props.plate;
	//	console.log('library: ', library, 'plate: ', plate)
		let current =""
		
		if (this.props.current){
			current = "(current)"; 
			}
		
		const buttons = display_options.map((option, index) => {
			return <button key={index} className={this.state[option[0]] ? "hidden" : "small-button"} onClick={event => this.toggleDisplay(option[0])}>Show {option[1]}</button>
			});
		
		const cols =  display_options.map((option, index) => {
			return <col key={index}  className={this.state[option[0]] ? "" : "hidden"} />
			});
		
		let display = {}
		
		display_options.forEach(item => {
			display[item[0]] = this.state[item[0]]
		});
		
		console.log('display: ', display) 
		
		return (
		<div id="plate-lookup">
				<h1>Library: {library.name} </h1>
				<h2>Plate: {plate.name} {plate.current}</h2>	
			<main>
				<div className="sidebar-div">
					<div>
					<h3>Show more table columns:</h3>
						{buttons}
					</div>
					
					<ExportBar />
					<PlateList library={library} showPlate={this.props.showPlate} />
					
				</div>
				<table className="datatable">
					<caption>
						Compound list ({this.state.compounds.length} items)
					</caption>
					<colgroup>
						<col/>
						{cols}
					</colgroup>
					<TableHeader
						display = {display}
						onButtonClick = {this.toggleDisplay}
					/>
					<DataTable 
						display = {display}
						compounds = {this.state.compounds}
						show_well = {this.state.show_well}
						show_code = {this.state.show_code}
						show_smiles = {this.state.show_smiles}
						show_structure = {this.state.show_structure}
						show_concentration = {this.state.show_concentration} 
						show_mw = {this.state.show_mw}
						show_tpsa = {this.state.show_tpsa}
						show_logp = {this.state.show_logp}
						show_heavy_atom_count = {this.state.show_heavy_atom_count}
						show_heavy_atom_mol_wt = {this.state.show_heavy_atom_mol_wt}
						show_nhoh_count = {this.state.show_nhoh_count}
						show_no_count = {this.state.show_no_count}
						show_num_h_acceptors = {this.state.show_num_h_acceptors}
						show_num_h_donors = {this.state.show_num_h_donors}
						show_num_het_atoms = {this.state.show_num_het_atoms}
						show_num_rot_bonds = {this.state.show_num_rot_bonds}
						show_num_val_electrons = {this.state.show_num_val_electrons}
						show_ring_count = {this.state.show_ring_count}
					/> 
				</table>		
			</main>
		</div>
		);
	}
}

export default PlateLookup;

