import React from 'react';
import LibraryOption from './library_option.js';
import axios from 'axios';

import { deepCopyObjectArray, getAttributeArray, mean, shareAllElements } from  '../../actions/stat_functions.js';

class Libraries extends React.Component {

	constructor(props) {
		super(props);
		this.state = {
			currentLibPlates: [],
		};
	}
	
	componentDidMount() {
		const apiUrl = 'api/library_selection_list/';
		
		axios.get(apiUrl)
			.then(res => {
			const currentLibPlates = res.data;
			this.setState({ currentLibPlates });
      });
	}
	
	isInProposal(plate, libArray){
		if (libArray){
			const plateLibId = plate.library.id;
			const arrayIds = [];
			libArray.map(lib => arrayIds.push(lib.id));
			return arrayIds.includes(plateLibId);
		}
		return false;
	}	
	
	render(){
		
		const proposalLibs =  this.props.proposal.libraries;
		
		const libraries = this.state.currentLibPlates.map((plate, index) => { 
			return <LibraryOption 
				key={index} 
				plate={plate} 
				showPlate={this.props.showPlate}
				handleCheckboxChange = {this.props.handleChange}
				defaultChecked={this.isInProposal(plate, proposalLibs)}
				/>;
		});
		
		const selectedLibIds = getAttributeArray(this.props.selectedLibs, "id");
		const selectionHasNotChanged = shareAllElements(getAttributeArray(this.props.proposal.libraries, "id"), selectedLibIds);
			
		return (
		<section id="libraries">
			<h2>XChem in-house fragment libraries</h2>
			
			<form id="libform" >
				<div id="libs">
					{libraries}
				</div>
				<button type="submit" onClick={event => this.props.handleSubmit(selectedLibIds)} disabled={selectionHasNotChanged}>Save changes in your library selection</button>
			</form>
		</section>
		)
	}
}

export default Libraries;
