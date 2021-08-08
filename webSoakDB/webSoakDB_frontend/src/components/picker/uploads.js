import React from 'react';
import OwnLibraryForm from './own_lib_form.js';
import CherryPickForm from './cherrypick_form.js';
import ExportForm from './export_form.js';

class Uploads extends React.Component {
  
  /*
  constructor(props){
    super(props)
    
    this.state = {
      proposal : this.props.proposal,
    }
  }
  
  componentDidUpdate(prevProps, prevState) {
    if (prevProps.proposal !== this.props.proposal) {
      this.setState({proposal : this.props.proposal});
    }
  }
*/
  getPublicSubsetIds(){
    let publicSubsetIds = []
    this.props.presets.forEach(preset => {
      publicSubsetIds.push(...preset.subsets);
      //preset.subsets.forEach(subs => publicSubsetIds.push(...subs));
    });
    return publicSubsetIds;
  }

  getUserSubsets(){
    //list selected subsets uploaded by the user (not from a preset)
    const publicSubs = this.getPublicSubsetIds();
    const allSelected = this.props.proposal.subsets;
    const userSubsets = allSelected.filter(subset => !publicSubs.includes(subset.id));
    return userSubsets;
  }

  render(){
    const uploadedLibs = this.props.proposal.libraries.filter(lib => !lib.public)
    const uploadedLibsHeader = uploadedLibs.length > 0 ? <h3>Already uploaded:</h3> : null;
    const liblist = uploadedLibs.map((lib, index) => <li key={index}>{lib.name}</li>);
    let uploadedSubsets = this.getUserSubsets();
    const uploadedSubsetsHeader = uploadedSubsets.length > 0 ? <h3>Already uploaded:</h3> : null;
    const subsetlist = uploadedSubsets.map((subset, index) => <li key={index}>{subset.name}</li>);
    
    return (
    <section id="upload">
      <div className="form-container">
        <h2>Upload your own library</h2>
          <OwnLibraryForm 
            proposal={this.props.proposal} 
            refreshAfterUpload={this.props.refreshAfterUpload}
            detectUnsavedChanges={this.props.detectUnsavedChanges}
            trackUnsavedChanges={this.props.trackUnsavedChanges}
            showOverlay={this.props.showOverlay}
            waiting = {this.props.waiting}
            updateWaitingStatus = {this.props.updateWaitingStatus}
          />
            {uploadedLibsHeader}
            <ul className="old-uploads">
              {liblist}
            </ul>
      </div>
      <div className="form-container">
        <hr />
        <h2>Upload cherrypicking list</h2>
        <CherryPickForm 
          libs={this.props.libs}
          proposal={this.props.proposal} 
          refreshAfterUpload={this.props.refreshAfterUpload} 
          detectUnsavedChanges={this.props.detectUnsavedChanges}
          trackUnsavedChanges={this.props.trackUnsavedChanges}
          showOverlay={this.props.showOverlay}
          waiting = {this.props.waiting}
          updateWaitingStatus = {this.props.updateWaitingStatus}
        />
            {uploadedSubsetsHeader}
            <ul className="old-uploads">
            {subsetlist}
            </ul>
            
      </div>
      <div id="formatting-help">
        <hr />
        <a href="/uploads/formatting_help/" target="_blank">CSV formatting help</a>
      </div>
      <div id="exports">
        <hr />
        <ExportForm proposal={this.props.proposal}/>
      </div>
    </section>
    )
  }
}

export default Uploads
