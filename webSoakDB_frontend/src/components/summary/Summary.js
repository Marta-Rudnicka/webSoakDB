import React from 'react';
//import './home.css';
import CollectionRow from './collection_row.js';
import TableHeader from './th.js';
import LibraryInTable from './library_rows.js';
import { deepCopyObjectArray, getAttributeArray, mean } from  '../../actions/stat_functions.js';
import axios from 'axios';

class Summary extends React.Component {
	
	constructor(props) {
		super(props);
		this.removeLibrary = this.removeLibrary.bind(this);
		this.state = { 	selectedLibs: [],
						selectedLibIds: [],
					};
	}
	
	
	componentDidMount() {
		const apiUrl = 'api/proposals/' + this.props.proposalName;;
		
		axios.get(apiUrl)
			.then(res => {
			const selectedLibs = res.data.libraries;
			this.setState({ selectedLibs, selectedLibIds : getAttributeArray(selectedLibs, "id") });
      });
	}
	
	removeLibrary(id){
		const libs = deepCopyObjectArray(this.state.selectedLibs)
		const found = libs.find(object => object.id === parseInt(id));
		libs.splice(libs.indexOf(found), 1);
		this.setState({selectedLibs : libs});
	}
	
	undoChages(){
		this.componentDidMount();
	}
	
	render() {
		
		let libraries = []
		
		if(this.state.selectedLibs){
			libraries = this.state.selectedLibs.map((library, index) => {
				return <LibraryInTable
						key={index} 
						library={library} 
						handleClick={this.removeLibrary}
						showPlate={this.props.showPlate} />		
			});
		}
		
		return (
			<div id="all">
				<h1>Selected Compounds for {this.props.proposalName} </h1>
				<section>
				<h2>Whole libraries</h2>
					<table className="summary-table">
						<TableHeader compoundsDescription="Available"/>
						<tbody>
							{libraries}
						</tbody>
					</table>
				</section>
				<section>
					<h2>Selections from libraries (cherrypicked compounds)</h2>
					<table className="summary-table">
						<TableHeader compoundsDescription="Selected"/>
						<tbody>
						
						</tbody>
					</table>
				</section>
				<div>
				
					<button onClick={event => this.props.updateLibrarySelection(getAttributeArray(this.state.selectedLibs, "id"))}>Save changes </button>
					<button onClick={event => this.undoChages()}>Undo all changes </button>
				</div>
			</div>
		); 
	}
}

export default Summary;
