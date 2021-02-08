import React from 'react';
import ExportBar from './export_bar.js';
import DataTable from './data_table.js';
import PlateList from './plate_list.js';
import TableHeader from './table_header.js';
import axios from 'axios';


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
			show_p3: false, 
			show_p4: false,
		};
	}
	
	toggleDisplay(display_option){
		this.setState({[display_option]: !this.state[display_option]});
	}
	
	componentDidMount() {
		const library = this.props.library;
		const plate = this.props.plate;
		const apiUrl = 'api/compounds/' + library + '/' + plate;
		
		axios.get(apiUrl)
			.then(res => {
			const compounds = res.data;
			this.setState({ compounds });
      });     		
	}
	
	render() {
		const library = this.props.library;
		const name = this.props.plate;
		let current =""
		
		if (this.props.current){
			current = "(current)"; 
			}
		
		const display_options = [
			["show_well", "Well"],
			["show_code", "Compound Code"],
			["show_smiles", "SMILES"],
			["show_structure", "2D Structure"],
			["show_concentration", "Concentration"],
			["show_mw", "Molecular Weight"],
			["show_p3", "[Property3]"],
			["show_p4", "[Property4]"]
		];
		
		const buttons = display_options.map((option, index) => {
			return <button key={index} className={this.state[option[0]] ? "hidden" : "small-button"} onClick={event => this.toggleDisplay(option[0])}>Show {option[1]}</button>
			});
		
		const cols =  display_options.map((option, index) => {
			return <col key={index}  className={this.state[option[0]] ? "" : "hidden"} />
			});
		
		return (
		<div id="plate-lookup">
				<h1>Library: {library} </h1>
				<h2>Plate: {name} {current}</h2>	
			<main>
				<div className="sidebar-div">
					<div>
					Show more table columns:
						{buttons}
					</div>
					
					<ExportBar />
					<PlateList current = {true} library ={library} />
					
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
						compounds = {this.state.compounds}
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

export default PlateLookup;

