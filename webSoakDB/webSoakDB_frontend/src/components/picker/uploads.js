import React from 'react';
import OwnLibraryForm from './own_lib_form.js';
import CherryPickForm from './cherrypick_form.js';
import ExportForm from './export_form.js';

class Uploads extends React.Component {
  
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

  render(){
    let uploadedLibs = []
    this.state.proposal.libraries.forEach(lib => {
      if (!lib.public){
        uploadedLibs.push(lib);
        }
      });
    
    const liblist = () => { 
      if (uploadedLibs.length > 0){
        return(
        <React.Fragment> 
          <h3>Already uploaded:</h3>
          <ul className="old-uploads">
             {uploadedLibs.map((lib, index) => <li key={index}>{lib.name}</li>)}
          </ul>
        </React.Fragment>
        ); 
        }
      else {
        return "";
      }
    };
    
    let publicSubsets = []
    this.props.presets.forEach(preset => {
      preset.subsets.forEach(subset=> publicSubsets.push(subset.id));
    });
    
    let uploadedSubsets = []
    
    this.state.proposal.subsets.forEach(subset => {
      if (!publicSubsets.includes(subset.id)){
        uploadedSubsets.push(subset);
        }
      });
      
    const subsetlist = () => { 
      if (uploadedSubsets.length > 0){
        return(
        <React.Fragment> 
          <h3>Already uploaded:</h3>
          <ul className="old-uploads">
             {uploadedSubsets.map((subset, index) => <li key={index}>{subset.name}</li>)}
          </ul>
        </React.Fragment>
        ); 
        }
      else {
        return "";
      }
    };
    
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
          />
            {liblist()}      
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
        />
            {subsetlist()}
      </div>
      <div id="formatting-help">
        <hr />
        <a href="/uploads/formatting_help/" target="_blank">CSV formatting help</a>
      </div>
      <div id="exports">
        <hr />
        <ExportForm proposal={this.state.proposal}/>
      </div>
    </section>
    )
  }
}

export default Uploads
