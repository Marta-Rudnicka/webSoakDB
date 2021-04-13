import React from 'react';
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
			selectedLibIds: getAttributeArray(this.props.proposal.libraries, "id"),
			initialLibs: getAttributeArray(this.props.proposal.libraries, "id"),
			selectedSubsetIds: getAttributeArray(this.props.proposal.subsets, "id"),
			initialSubsets: getAttributeArray(this.props.proposal.subsets, "id"),
			currentLibPlates: [],
			presets: [],
			//unsavedChanges: false,
			
		}
		this.detectUnsavedChanges = this.detectUnsavedChanges.bind(this);
		this.handleChangeLib = this.handleChangeLib.bind(this);
		this.handleChangePreset = this.handleChangePreset.bind(this);
		this.updateSelectionLibs = this.updateSelectionLibs.bind(this);
	}
	
	componentDidMount() {
		const libApiUrl = '/api/library_selection_list/';
		
		axios.get(libApiUrl)
			.then(res => {
			const currentLibPlates = res.data;
			this.setState({ currentLibPlates });
      });
      
      const presetApiUrl = '/api/preset_list/';
		
		axios.get(presetApiUrl)
			.then(res => {
			const presets = res.data;
			this.setState({ presets });
      });
      
	}

	componentDidUpdate(prevProps, prevState) {
		if (prevProps.proposal !== this.props.proposal) {
			console.log('Proposal changed')
			this.setState({
				selectedLibIds : getAttributeArray(this.props.proposal.libraries, "id"), 
				initialLibs: getAttributeArray(this.props.proposal.libraries, "id"),
				selectedSubsetIds: getAttributeArray(this.props.proposal.subsets, "id"),
				initialSubsets: getAttributeArray(this.props.proposal.subsets, "id"),
				});
		}
	}
	
	removeLibraryFromSelected(id){
		const selectedLibIdsCopy = this.state.selectedLibIds.slice(0, this.state.selectedLibIds.length);
		const found = selectedLibIdsCopy.find(item => item === parseInt(id));
		selectedLibIdsCopy.splice(selectedLibIdsCopy.indexOf(found), 1);
		this.setState({selectedLibIds : selectedLibIdsCopy});
		this.handleUnsavedChanges(selectedLibIdsCopy, null)
	}
	
	addLibraryToSelected(id){
		const selectedLibIdsCopy = this.state.selectedLibIds.slice(0, this.state.selectedLibIds.length);
		selectedLibIdsCopy.push(parseInt(id));
		this.setState({selectedLibIds : selectedLibIdsCopy});
		this.handleUnsavedChanges(selectedLibIdsCopy, null);
	}
	
	removePresetFromSelected(id){
		const removedPreset = this.state.presets.find(preset => preset.id === parseInt(id))
		const selectedSubsetIdsCopy = this.state.selectedSubsetIds.slice(0, this.state.selectedSubsetIds.length);	
		removedPreset.subsets.forEach(subset => {
			const found = selectedSubsetIdsCopy.find(item => item === parseInt(id));
			selectedSubsetIdsCopy.splice(selectedSubsetIdsCopy.indexOf(found), 1);
		});
		this.setState({selectedSubsetIds : selectedSubsetIdsCopy});
		this.handleUnsavedChanges(null, selectedSubsetIdsCopy);
	}
	
	addPresetToSelected(id){
		const addedPreset = this.state.presets.find(preset => preset.id === parseInt(id))
		const selectedSubsetIdsCopy = this.state.selectedSubsetIds.slice(0, this.state.selectedSubsetIds.length);
		addedPreset.subsets.forEach(subset => selectedSubsetIdsCopy.push(subset.id));
		this.setState({selectedSubsetIds : selectedSubsetIdsCopy});
		this.handleUnsavedChanges(null, selectedSubsetIdsCopy);
	}

	handleChangeLib(event, unsaved){
		if(event.target.checked === true){
			this.addLibraryToSelected(event.target.value);
		}
		else{
			this.removeLibraryFromSelected(event.target.value);
		}
	}
	
	handleChangePreset(event, unsaved){
		if(event.target.checked === true){
			this.addPresetToSelected(event.target.value);
		}
		else{
			this.removePresetFromSelected(event.target.value);
		}
	}
	
	updateSelectionLibs(){
		console.log('fired updateSelectionLibs');
		event.preventDefault();
		this.handleUnsavedChanges(this.state.initialLibs, this.state.selectedSubsetIds);
		this.props.updateLibrarySelection(this.state.selectedLibIds);//, 'Picker');
		
	}
	
	updateSelectionSubsets(){
		this.handleUnsavedChanges(this.state.selectedLibIds, this.state.initialSubsets);
		this.props.updateSubsetSelection(this.state.selectedSubsetIds);
		
	}
	
	presetAreadySelected(preset){
		if (this.state.selectedSubsetIds.includes(preset.subsets[0].id)){
			return true;
		}
		else{
			return false;
		}
	}
	
	handleUnsavedChanges(newLibs, newSubsets){
	/*detects unsaved changes and passes them to trackUnsavedChanges
	 *if caller function does not change library selection, pass null for newLibs 
	 *if caller function is supposed to ignore changes in library selection (e.g.
	 * while saving changes in the db), pass this.state.initialLibIds for newLibs
	 * (and the same for subset selection and newSubsets arg)*/
		console.log('fired handleUnsavedChanges')
	 	this.props.trackUnsavedChanges(this.detectUnsavedChanges(newLibs, newSubsets));
	}
	
	detectUnsavedChanges(newLibs, newSubsets){
	/*determine if newLibs and newSubsets are different from the saved 
	 * selection of compounds; if both args are null, compares current selection
	 * to the saved one */
	  	
		let currentLibs = this.state.selectedLibIds;
		let currentSubsets = this.state.selectedSubsetIds;
		if (newLibs){
			currentLibs = newLibs;
		}
		if (newSubsets){
			currentSubsets = newSubsets;
		}		
		if (shareAllElements(currentLibs, this.state.initialLibs) && shareAllElements(currentSubsets, this.state.initialSubsets)){
			return false;
		}
		else {
			return true;
		}
	}
	
	
	render() {
		let sameLibs
		if (!shareAllElements(this.state.initialLibs, this.state.selectedLibIds)){
			sameLibs = 'changed';
		}
		let sameSubsets
		if (!shareAllElements(this.state.initialSubsets, this.state.selectedSubsetIds)){
			sameSubsets = 'changed';
		}
		const unsaved = !(sameLibs && sameSubsets)
		
		const libraries = this.state.currentLibPlates.map((plate, index) => { 
			return <LibraryOption 
				key={index} 
				plate={plate}
				handleCheckboxChange = {this.handleChangeLib}
				defaultChecked={this.state.selectedLibIds.includes(plate.library.id)}
				unsaved={unsaved}
				/>;
		});
		
		const presets = this.state.presets.map((preset, index) => {
			return <PresetOption
				key={index}
				preset = {preset}
				handleCheckboxChange = {this.handleChangePreset}
				defaultChecked={this.presetAreadySelected(preset)}
				unsaved={unsaved}
				/>}
		)
		
		const proposalLibs =  this.props.proposal.libraries;
		
		
		
		let publicSubsets = [];
		this.state.presets.forEach(preset => publicSubsets.push(...preset.subsets));
		
		return (
		<div id="picker">
			<h1>Select compounds for {this.props.proposal.proposal}</h1>
			<main id="main-picker">
				<section id="libraries" className={sameLibs}>
					<h2>XChem in-house fragment libraries</h2>
					
					<form id="libform" >
						<div id="libs">
							{libraries}
						</div>
						<button type="submit" onClick={event => this.updateSelectionLibs(event)}>Save changes in your selection</button>
					</form>
				</section>
				
				<section id="presets" className={sameSubsets}>
					<h2>Presets</h2>
					<p>Specific-purpose compounds selections from in-house libraries</p>
					<form id="preset-form">
						<div id="pres">
							{presets}
						</div>
						<button type="submit" onClick={event => this.updateSelectionSubsets()}>Save changes in your selection</button>
					</form>
				</section>
				
				<Uploads proposal={this.props.proposal}
					publicSubsets={publicSubsets} 
					refreshAfterUpload={this.props.refreshAfterUpload}
					detectUnsavedChanges={this.detectUnsavedChanges}
					proposal={this.props.proposal}
					trackUnsavedChanges={this.props.trackUnsavedChanges}
					presets={this.state.presets}
				/>
				<Stats 
					proposal={this.props.proposal} 
					selectedLibIds={this.state.selectedLibIds} 
					selectedSubsetIds={this.state.selectedSubsetIds}
				/>
			
			</main>
		</div>
		); 
	}

}

export default Picker;
