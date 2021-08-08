import React from 'react';
import LibraryInTable from './library_rows.js';
import SubsetTable from './subset_table.js';
import { removeFromArray, getAttributeArray } from  '../../actions/stat_functions.js';

class Summary extends React.Component {
	
	constructor(props) {
		super(props);
		this.removeLibrary = this.removeLibrary.bind(this);
		this.removeSubset = this.removeSubset.bind(this);
		this.state = { 	selectedLibs: this.props.proposal.libraries,
						selectedSubsets: this.props.proposal.subsets,
						libClass : null,
						subsetClass : null,
		};
	}

	resetSelection(){
		this.setState({
			selectedLibs: this.props.proposal.libraries, 
			selectedSubsets: this.props.proposal.subsets
		  })
	}

	removeLibrary(id){
		const lib = this.state.selectedLibs.find(object => object.id === parseInt(id));
		const newSelectedLibs = removeFromArray(this.state.selectedLibs, [lib]);
		this.setState({selectedLibs : newSelectedLibs});
		this.props.trackUnsavedChanges(true);
		this.setState({libClass: "changed"});
	}
	
	removeSubset(id){
		const sub = this.state.selectedSubsets.find(object => object.id === parseInt(id));
		const newSelectedSubsets = removeFromArray(this.state.selectedSubsets, [sub]);
		this.setState({selectedSubsets: newSelectedSubsets})
		this.props.trackUnsavedChanges(true);
		this.setState({subsetClass: "changed"})
	}
	
	undoChanges(){
		this.resetSelection();
		this.props.trackUnsavedChanges(false);
		this.setState({libClass : null, subsetClass : null});
	}
	
	saveChanges(){
		this.props.updateSelection(getAttributeArray(this.state.selectedLibs, "id"), getAttributeArray(this.state.selectedSubsets, "id"));
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
						project={this.props.proposal}
						/>		
			});
		}
		
		return (
		  <div id="all">
			<h1>Selected Compounds for {this.props.proposal.auth[0].proposal_visit} </h1>
			<main id="summary-main">
			<section id="whole" className={this.state.libClass}>
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
			<section id="subsets" className={this.state.subsetClass}>
			  <SubsetTable 
				selectedSubsets={this.state.selectedSubsets} 
				handleClick={this.removeSubset}
			  />
			  <div>
				<button onClick={event => this.saveChanges()}>Save changes </button>
				<button onClick={event => this.undoChanges()}>Undo unsaved changes </button>
			  </div>
			</section>
		  </main>
		  </div>
		); 
	}
}

export default Summary;
