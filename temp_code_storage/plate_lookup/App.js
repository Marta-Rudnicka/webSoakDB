import React from 'react';
import './plate_details.css';
//import ButtonBar from './button_bar.js';
import ExportBar from './export_bar.js';
import DataTable from './data_table.js';
import PlateList from './plate_list.js';
import TableHeader from './table_header.js';

//const proposal = "Test proposal"

class App extends React.Component {
	
	constructor(props) {
		super(props);
		this.toggleDisplay = this.toggleDisplay.bind(this);
	//	this.setToFalse = this.setToFalse.bind(this);
		this.state = {
			show_well: true, 
			show_code: true, 
			show_smiles: true, 
			show_structure: false, 
			show_concentration: true, 
			show_mw: true, 
			show_p3: false, 
			show_p4: false,
		};
	}
	
	toggleDisplay(display_option){
		this.setState({[display_option]: !this.state[display_option]});
	}
	
	render() {
		const library = this.props.plate.library.name;
		const name = this.props.plate.name;
		
		let current = "";
		if (this.props.plate.current){
			current = "current"; 
			}
		
		const display_options = [
			["show_well", "Show Well"],
			["show_code", "Show Compound Code"],
			["show_smiles", "Show SMILES"],
			["show_structure", "Show 2D Structure"],
			["show_concentration", "Show Concentration"],
			["show_mw", "Show Molecular Weight"],
			["show_p3", "Show [Property3]"],
			["show_p4", "Show [Property4]"]
		];
		
		const buttons = display_options.map( option => {
			return <button className={this.state[option[0]] ? "hidden" : "small-button"} onClick={event => this.toggleDisplay(option[0])}>{option[1]}</button>
			});
		
		const cols =  display_options.map( option => {
			return <col className={this.state[option[0]] ? "" : "hidden"} />
			});
		
		return (
		<div id="all">
				<h1>Library: {library} </h1>
				<h2>Plate: {name} ({current})</h2>	
			<main>
				<div className="sidebar-div">
					<div>
					Show more table columns:
						{buttons}
					</div>
					
					<ExportBar />
					<PlateList current = {this.props.plate.current} library ={library} />
					
				</div>
				<table className="datatable">
					<caption>
						<h2>Compound list ({this.props.compounds.length} items)</h2>
					</caption>
					<colgroup>
						<col/>
						{cols}
					</colgroup>
					<TableHeader
						show_well = {this.state.show_well}
						show_code = {this.state.show_code}
						show_smiles = {this.state.show_smiles}
						show_structure = {this.state.show_structure}
						show_concentration = {this.state.show_concentration} 
						show_mw = {this.state.show_mw}
						show_p3 = {this.state.show_p3}
						show_p4 = {this.state.show_p4}
						onButtonClick = {this.toggleDisplay}
					/>
					<DataTable 
						compounds = {this.props.compounds}
						show_well = {this.state.show_well}
						show_code = {this.state.show_code}
						show_smiles = {this.state.show_smiles}
						show_structure = {this.state.show_structure}
						show_concentration = {this.state.show_concentration} 
						show_mw = {this.state.show_mw}
						show_p3 = {this.state.show_p3}
						show_p4 = {this.state.show_p4}
					/> 
				</table>		
			</main>
		</div>
		); 
	}
}

export default App;

