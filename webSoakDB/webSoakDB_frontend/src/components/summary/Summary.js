import React from 'react';
//import './home.css';
import CollectionRow from './collection_row.js';
import TableHeader from './th.js';
import LibraryInTable from './library_rows.js';
import SubsetTable from './subset_table.js';
import { deepCopyObjectArray, getAttributeArray, mean } from  '../../actions/stat_functions.js';
import axios from 'axios';

class Summary extends React.Component {
	
	constructor(props) {
		super(props);
		this.removeLibrary = this.removeLibrary.bind(this);
		this.removeSubset = this.removeSubset.bind(this);
		this.state = { 	selectedLibs: [],
						selectedSubsets: [],
						libClass : null,
						subsetClass : null,
					};
	}
	
	
	componentDidMount() {
		const apiUrl = 'api/proposals/' + this.props.proposal.name;
		
		axios.get(apiUrl)
			.then(res => {
			const selectedLibs = res.data.libraries;
			this.setState({ selectedLibs, selectedLibIds : getAttributeArray(selectedLibs, "id") });
			const selectedSubsets = res.data.subsets;
			this.setState({ selectedSubsets, selectedSubsetIds : getAttributeArray(selectedSubsets, "id") });
      });
	}
	
	removeLibrary(id){
		const libs = deepCopyObjectArray(this.state.selectedLibs)
		const found = libs.find(object => object.id === parseInt(id));
		libs.splice(libs.indexOf(found), 1);
		this.setState({selectedLibs : libs});
		this.props.trackUnsavedChanges(true);
		this.setState({libClass: "changed"});
	}
	
	removeSubset(id){
		const subsets = deepCopyObjectArray(this.state.selectedSubsets)
		const found = subsets.find(object => object.id === parseInt(id));
		subsets.splice(subsets.indexOf(found), 1);
		this.setState({selectedSubsets : subsets});	
		this.props.trackUnsavedChanges(true);
		this.setState({subsetClass: "changed"})
	}
	
	undoChages(){
		this.componentDidMount();
		this.props.trackUnsavedChanges(false);
		this.setState({libClass : null, subsetClass : null});
	}
	
	saveChanes(){
		this.props.updateLibrarySelection(getAttributeArray(this.state.selectedLibs, "id"));
		this.props.updateSubsetSelection(getAttributeArray(this.state.selectedSubsets, "id"));
		this.props.trackUnsavedChanges(false);
		this.setState({libClass : null, subsetClass : null});
	
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
				<h1>Selected Compounds for {this.props.proposal.name} </h1>
				<main id="summary-main">
				<section className={this.state.libClass}>
					<h2>Whole libraries</h2>
					<table className="table table-bordered" id="table">
						<thead>
							<tr>
								<th>Library</th>
								<th>Origin</th>
								<th>Currently used plate</th>
								<th>Available <br/>compounds</th>
								<th>Compounds</th>
								<th>Remove</th>
							</tr>
						</thead>	
						<tbody>
							{libraries}
						</tbody>
					</table>
				</section>
				<section className={this.state.subsetClass}>
					<SubsetTable 
						selectedSubsets={this.state.selectedSubsets} 
						lookup_args={this.props.lookup_args} 
						handleClick={this.removeSubset} 
						showPlate={this.props.showPlate}
					/>
				</section>
				<div>
				
					<button onClick={event => this.saveChanes()}>Save changes </button>
					<button onClick={event => this.undoChages()}>Undo unsaved changes </button>
				</div>
			</main>
			</div>
			
		); 
	}
}

export default Summary;