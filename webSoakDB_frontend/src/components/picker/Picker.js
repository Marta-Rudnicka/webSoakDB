import React from 'react';
import Libraries from './libraries.js';
import LibraryOption from './library_option.js';
import Presets from './presets.js';
import Uploads from './uploads.js';
import Stats from './stats.js';
import Graphs from './graphs.js';
import PresetOption from './preset_option.js';
import axios from 'axios';

import { deepCopyObjectArray, getAttributeArray, mean, shareAllElements } from  '../../actions/stat_functions.js';

const proposal = 'Test string';

class Picker extends React.Component {
	
	constructor(props) {
		super(props);
		this.state = {
			//presets: presets,
			selectedLibIds: getAttributeArray(this.props.proposal.libraries, "id"),
			initialLibs: getAttributeArray(this.props.proposal.libraries, "id"),
			selectedSubsetIds: getAttributeArray(this.props.proposal.subsets, "id"),
			initialSubsets: getAttributeArray(this.props.proposal.subsets, "id"),
			currentLibPlates: [],
			presets: [],
		}
		this.handleChangeLib = this.handleChangeLib.bind(this);
		this.handleChangePreset = this.handleChangePreset.bind(this);
		this.updateSelectionLibs = this.updateSelectionLibs.bind(this);
	}
	
	componentDidMount() {
		const libApiUrl = 'api/library_selection_list/';
		
		axios.get(libApiUrl)
			.then(res => {
			const currentLibPlates = res.data;
			this.setState({ currentLibPlates });
      });
      
      const presetApiUrl = 'api/preset_list/';
		
		axios.get(presetApiUrl)
			.then(res => {
			const presets = res.data;
			this.setState({ presets });
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
	
	removePresetFromSelected(id){
		const removedPreset = this.state.presets.find(preset => preset.id === parseInt(id))
		const selectedSubsetIdsCopy = this.state.selectedSubsetIds.slice(0, this.state.selectedSubsetIds.length);	
		removedPreset.subsets.forEach(subset => {
			const found = selectedSubsetIdsCopy.find(item => item === parseInt(id));
			selectedSubsetIdsCopy.splice(selectedSubsetIdsCopy.indexOf(found), 1);
		});
		this.setState({selectedSubsetIds : selectedSubsetIdsCopy});
	}
	
	addPresetToSelected(id){
		const addedPreset = this.state.presets.find(preset => preset.id === parseInt(id))
		const selectedSubsetIdsCopy = this.state.selectedSubsetIds.slice(0, this.state.selectedSubsetIds.length);
		selectedSubsetIdsCopy.push(...addedPreset.subsets);
		this.setState({selectedSubsetIds : selectedSubsetIdsCopy});
	}

	handleChangeLib(event){
		if(event.target.checked === true){
			this.addLibraryToSelected(event.target.value);
		}
		else{
			this.removeLibraryFromSelected(event.target.value);
		}
	}
	
	handleChangePreset(event){
		if(event.target.checked === true){
			this.addPresetToSelected(event.target.value);
		}
		else{
			this.removePresetFromSelected(event.target.value);
		}
	}
	
	updateSelectionLibs(){
		this.props.updateLibrarySelection(this.state.selectedLibIds, 'Picker');
	}
	
	updateSelectionSubsets(){
		this.props.updateSubsetSelection(this.state.selectedSubsetIds, 'Picker');
	}
	
	render() {
		
		const libraries = this.state.currentLibPlates.map((plate, index) => { 
			return <LibraryOption 
				key={index} 
				plate={plate} 
				showPlate={this.props.showPlate}
				handleCheckboxChange = {this.handleChangeLib}
				defaultChecked={this.state.selectedLibIds.includes(plate.library.id)}
				/>;
		});
		
		const presets = this.state.presets.map((preset, index) => {
		return <PresetOption
			key={index}
			id={preset.id}
			name = {preset.name}
			handleCheckboxChange = {this.handleChangePreset}
			description = {preset.description}
			/>}
		)
		
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
						<button type="submit" onClick={event => this.updateSelectionLibs()} disabled={selectionHasNotChanged}>Save changes in your library selection</button>
					</form>
				</section>
				
				<section id="presets">
					<h2>Presets</h2>
					<p>Specific-purpose compounds selections from in-house libraries</p>
					<form id="properties-form">
						{presets}
					</form>
					<button type="submit" onClick={event => this.updateSelectionSubsets()}>Submit</button>
				</section>
				
				<Uploads proposal={this.props.proposal} changeMainPage={this.props.changeMainPage}/>
				<Stats proposal={this.props.proposal} selectedLibIds={this.state.selectedLibIds}/>
			
			</main>
		</div>
		); 
	}

}

export default Picker;
