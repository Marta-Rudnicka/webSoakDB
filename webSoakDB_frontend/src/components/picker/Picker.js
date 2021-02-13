import React from 'react';
import Libraries from './libraries.js';
import Presets from './presets.js';
import Uploads from './uploads.js';
import Stats from './stats.js';
import Graphs from './graphs.js';
import axios from 'axios';

import { deepCopyObjectArray, getAttributeArray, mean } from  '../../actions/stat_functions.js';

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
			selectedLibs: this.props.proposal.libraries,
		}
		this.handleChange = this.handleChange.bind(this);
	}
	
	removeLibraryFromSelected(id){
		const libs = deepCopyObjectArray(this.state.selectedLibs)
		const found = libs.find(object => object.id === parseInt(id));
		libs.splice(libs.indexOf(found), 1);
		this.setState({selectedLibs : libs});
	}
	
	addLibraryToSelected(id){
		const libs = deepCopyObjectArray(this.state.selectedLibs);
		const apiUrl = 'api/library_detail/' + id;

		axios.get(apiUrl)
			.then(res => {
			const lib = res.data;
			libs.push(lib);
			this.setState({selectedLibs: libs});			
      });
		
	}
	
	handleChange(event){
		if(event.target.checked === true){
			this.addLibraryToSelected(event.target.value);
		}
		else{
			this.removeLibraryFromSelected(event.target.value);
		}
	}
	
	render() {
		
		return (
		<div id="picker">
			<h1>Select compounds for {this.props.proposal.name}</h1>
			<main id="main-picker">
				<Libraries showPlate={this.props.showPlate}  proposal={this.props.proposal} handleChange={this.handleChange} selectedLibs={this.state.selectedLibs} handleSubmit={this.props.updateLibrarySelection}/>
				<Presets presets={this.state.presets}  proposal={this.props.proposal}/>
				<Uploads proposal={this.props.proposal} changeMainPage={this.props.changeMainPage}/>
				<Stats proposal={this.props.proposal} selectedLibs={this.state.selectedLibs}/>
			</main>
		</div>
	); 
}

}

export default Picker;
