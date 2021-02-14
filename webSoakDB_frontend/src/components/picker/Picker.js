import React from 'react';
import Libraries from './libraries.js';
import LibraryOption from './library_option.js';
import Presets from './presets.js';
import Uploads from './uploads.js';
import Stats from './stats.js';
import Graphs from './graphs.js';
import axios from 'axios';

import { deepCopyObjectArray, getAttributeArray, mean, shareAllElements } from  '../../actions/stat_functions.js';

const presets = [
	{
		"id": 1,
		"name": "Preset 1",
		"description": "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed nec tempus libero, quis porta mi. ",
	},
	{
		"id": 2,
		"name": "Preset 2",
		"description": "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed nec tempus libero, quis porta mi. ",
	}
]

const proposal = 'Test string';

class Picker extends React.Component {
	
	constructor(props) {
		super(props);
		this.state = {
			presets: presets,
			//libSelection: =
			//selectedLibIds: this.props.libSelection,
			selectedLibIds: getAttributeArray(this.props.proposal.libraries, "id"),
			initialLibs: getAttributeArray(this.props.proposal.libraries, "id"),

			selectedSubsetIds: [],
			currentLibPlates: [],
		}
		this.handleChange = this.handleChange.bind(this);
		this.updateSelection = this.updateSelection.bind(this);
	}
	
	componentDidMount() {
		const apiUrl = 'api/library_selection_list/';
		
		axios.get(apiUrl)
			.then(res => {
			const currentLibPlates = res.data;
			this.setState({ currentLibPlates });
      });
      
	}
	componentDidUpdate(prevProps, prevState) {
		if (prevProps.proposal !== this.props.proposal) {
			this.setState({selectedLibIds : getAttributeArray(this.props.proposal.libraries, "id"), initialLibs: getAttributeArray(this.props.proposal.libraries, "id")});
		}
	}
	
	removeLibraryFromSelected(id){
		const selectedLibIdsCopy = this.state.selectedLibIds.slice(0, this.state.selectedLibIds.length);
		const found = selectedLibIdsCopy.find(item => item === parseInt(id));
		selectedLibIdsCopy.splice(selectedLibIdsCopy.indexOf(found), 1);
		this.setState({selectedLibIds : selectedLibIdsCopy});
	}
	
	addLibraryToSelected(id){
		const selectedLibIdsCopy = this.state.selectedLibIds.slice(0, this.state.selectedLibIds.length);
		selectedLibIdsCopy.push(parseInt(id));
		this.setState({selectedLibIds : selectedLibIdsCopy});
	}
	
	handleChange(event){
		if(event.target.checked === true){
			this.addLibraryToSelected(event.target.value);
		}
		else{
			this.removeLibraryFromSelected(event.target.value);
		}
	}
	
	updateSelection(){
		this.props.updateLibrarySelection(this.state.selectedLibIds, 'Picker');
	}
	
	render() {
		
		const libraries = this.state.currentLibPlates.map((plate, index) => { 
			return <LibraryOption 
				key={index} 
				plate={plate} 
				showPlate={this.props.showPlate}
				handleCheckboxChange = {this.handleChange}
				defaultChecked={this.state.selectedLibIds.includes(plate.library.id)}
				/>;
		});
		
		const proposalLibs =  this.props.proposal.libraries;
		//const selectionHasNotChanged = shareAllElements(getAttributeArray(this.props.proposal.libraries, "id"), this.state.selectedLibIds);			
		const selectionHasNotChanged = shareAllElements(this.state.initialLibs, this.state.selectedLibIds);			
		
		return (
		<div id="picker">
			<h1>Select compounds for {this.props.proposal.name}</h1>
			<main id="main-picker">
				<section id="libraries">
					<h2>XChem in-house fragment libraries</h2>
					
					<form id="libform" >
						<div id="libs">
							{libraries}
						</div>
						<button type="submit" onClick={event => this.updateSelection()} disabled={selectionHasNotChanged}>Save changes in your library selection</button>
					</form>
				</section>
				
				<Presets presets={this.state.presets}  proposal={this.props.proposal}/>
				<Uploads proposal={this.props.proposal} changeMainPage={this.props.changeMainPage}/>
				<Stats proposal={this.props.proposal} selectedLibIds={this.state.selectedLibIds}/>
			
			</main>
		</div>
		); 
	}

}

export default Picker;
