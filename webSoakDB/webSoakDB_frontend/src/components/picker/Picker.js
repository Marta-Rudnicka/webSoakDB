import React from 'react';
import LibraryOption from './library_option.js';
import Uploads from './uploads.js';
import GraphTable from './graph_table.js';
import PresetOption from './preset_option.js';
import axios from 'axios';

import { 
	removeFromArray, 
	getAttributeArray, 
	shareAllElements, 
	getSubsetIds,
	} 
	from  '../../actions/stat_functions.js';

const proposal = 'Test string';

class Picker extends React.Component {
  
  constructor(props) {
    super(props);
    this.state = {
      selectedLibIds: getAttributeArray(this.props.proposal.libraries, "id"),
      initialLibs: getAttributeArray(this.props.proposal.libraries, "id"),
      selectedSubsetIds: getAttributeArray(this.props.proposal.subsets, "id"),
      initialSubsets: getAttributeArray(this.props.proposal.subsets, "id"),
      presets: [],
      currentLibOptions: [],
      inHouseCompoundCount: 0,      
    }
    
    this.detectUnsavedChanges = this.detectUnsavedChanges.bind(this);
    this.handleChangeLib = this.handleChangeLib.bind(this);
    this.handleChangePreset = this.handleChangePreset.bind(this);
    this.saveChanges = this.saveChanges.bind(this);
    this.selected = this.selected.bind(this);
  }
  
  componentDidMount() {
        
    const presetApiUrl = '/api/preset_list/';
    
    axios.get(presetApiUrl)
      .then(res => {
      const presets = res.data;
      this.setState({ presets });
      });
      
    const libApiUrl = '/api/current_library_options/';
    
    axios.get(libApiUrl)
      .then(res => {
      const currentLibOptions = res.data;
      this.setState({ currentLibOptions })
      });
  }

  componentDidUpdate(prevProps, prevState) {
    if (prevProps.proposal !== this.props.proposal) {
      this.setState({
        selectedLibIds : getAttributeArray(this.props.proposal.libraries, "id"), 
        initialLibs: getAttributeArray(this.props.proposal.libraries, "id"),
        selectedSubsetIds: getAttributeArray(this.props.proposal.subsets, "id"),
        initialSubsets: getAttributeArray(this.props.proposal.subsets, "id"),
        });
    }
    if (prevState.selectedLibIds !== this.state.selectedLibIds
      || prevState.presets !== this.state.presets
      || prevState.selectedSubsetIds != this.state.selectedSubsetIds
      || prevState.currentLibOptions !== this.state.currentLibOptions
      )
      {
        this.updateInHouseCompoundCount();
        this.props.trackUnsavedChanges(this.detectUnsavedChanges());
      }
  }
  
  handleChangeLib(event){
  /*add or remove library from temporary (unsaved) selection on checkbox change*/
    const id = parseInt(event.target.value);
    
    if(event.target.checked === true){
      const newSelectedLibIds = this.state.selectedLibIds.slice(0, this.state.selectedLibIds.length);
      newSelectedLibIds.push(parseInt(id));
      this.setState({selectedLibIds : newSelectedLibIds});
    }
    else{
      const newSelectedLibIds = removeFromArray(this.state.selectedLibIds, [id]);
      this.setState({selectedLibIds : newSelectedLibIds});
    }
  }
  
  handleChangePreset(event){
  /*add or remove subsets from temporary (unsaved) selection on preset checkbox change*/
    const id = parseInt(event.target.value);
    
    if(event.target.checked === true){
      const newSelectedSubsetIds = this.state.selectedSubsetIds.slice(0, this.state.selectedSubsetIds.length);
      newSelectedSubsetIds.push(...getSubsetIds(this.state, id))
      this.setState({selectedSubsetIds : newSelectedSubsetIds});
    }
    else{
      const newSelectedSubsetIds = removeFromArray(this.state.selectedSubsetIds, getSubsetIds(this.state, id));
      this.setState({selectedSubsetIds : newSelectedSubsetIds});
    }
  }

  saveChanges(){
    //event.preventDefault();
    this.props.updateSelection(this.state.selectedLibIds, this.state.selectedSubsetIds);
    location.reload();
  }


  resetSelection(){
    //event.preventDefault();
    this.props.updateSelection([], []);
    location.reload();
  }

  selected(preset){
    return this.state.selectedSubsetIds.includes(preset.subsets[0])
  }

  detectUnsavedChanges(){
    if (!shareAllElements(this.state.initialLibs, this.state.selectedLibIds)){
      return true;
    }
    if (!shareAllElements(this.state.initialSubsets, this.state.selectedSubsetIds)){
      return true;
    }
    return false;
  }
  
  updateInHouseCompoundCount(){
    let count = 0;
    this.state.currentLibOptions.forEach(lib =>{
      if (this.state.selectedLibIds.includes(lib.id)){
        count = count + lib.size;
      }
    });
    this.state.presets.forEach(preset =>{
      const subs = getSubsetIds(this.state, preset.id)
      if (this.state.selectedSubsetIds.includes(subs[0])){
        count = count + preset.size;
      }
    });
    this.setState({inHouseCompoundCount: count});
  }
  
  render() {
    let changeStatus = "";
    if (this.props.unsavedChanges){
      changeStatus = 'changed';
    }
    
    const libraries = this.state.currentLibOptions.map((lib, index) => {
      const checked =  this.state.selectedLibIds.includes(lib.id)
      return <LibraryOption 
        key={index} 
        lib={lib}
        handleCheckboxChange = {this.handleChangeLib}
        handleChangePreset = {this.handleChangePreset}
        defaultChecked={this.state.selectedLibIds.includes(lib.id)}
        selected = {this.selected}
        />;
    });
    
    const otherPresets = this.state.presets.map((preset, index) => {
      if (preset.subsets.length > 1){
        return <li key={index}><PresetOption
          preset = {preset}
          handleCheckboxChange = {this.handleChangePreset}
          defaultChecked={this.selected(preset)}
          /></li>}
      }
    )
  
    
    let publicSubsets = [];
    this.state.presets.forEach(preset => publicSubsets.push(...preset.subsets));
    
    return (
    <div id="picker">
      <h1>Select compounds for {this.props.proposal.proposal}</h1>
      <main id="main-picker">
        <section id="libraries" className={changeStatus}>
		  <h2>Libraries available:</h2>
		  <div id="libform" >
			  <div id="libs">
			    {libraries}
			  </div>
			<p className="bottom-right"><strong>Total number of compounds selected: {this.state.inHouseCompoundCount}</strong></p>
          </div>
          <div id="pres">
            <p>OTHER PRESETS</p>
            <ul>
              {otherPresets}
            </ul>
          </div>
          <button type="submit" onClick={event => this.saveChanges(event)}>Save changes in your selection</button>
          <button type="submit" onClick={event => this.resetSelection(event)}>Reset selection (demo only)</button>         
        </section>
                
        <Uploads proposal={this.props.proposal}
          publicSubsets={publicSubsets} 
          refreshAfterUpload={this.props.refreshAfterUpload}
          detectUnsavedChanges={this.detectUnsavedChanges}
          proposal={this.props.proposal}
          trackUnsavedChanges={this.props.trackUnsavedChanges}
          presets={this.state.presets}
        />
        <section id="stats">
          <GraphTable
            parentState = {this.state}
          />
        </section>
      </main>
    </div>
    ); 
  }

}

export default Picker;
