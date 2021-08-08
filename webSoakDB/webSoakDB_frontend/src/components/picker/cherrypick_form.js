import React from 'react';
import CSRFToken from './csrf.js';
import OwnLibraryForm from './own_lib_form.js';

class CherryPickForm extends OwnLibraryForm {
  
  triggerFormSubmission(changesBeforeSubmission){
    document.forms["own_subset"].requestSubmit();
    const initialSubNumber = this.props.proposal.subsets.length;
    this.refreshUntilAdded(initialSubNumber, changesBeforeSubmission);
  }
  
  checkForChanges(initialNumber){
    if (this.props.proposal.subsets.length !== initialNumber){
      this.props.updateWaitingStatus(false); //this will trigger stopReloading()
    }
  }
  
  render(){
    let options = null;
    if (this.props.libs){
      options = this.props.libs.map((lib, index) => 
        <option key={index} value={lib.id}>{lib.name}</option>
      );
    }
    
    const name = this.state.name;
    const file = this.state.file;
    const library = this.state.library;
    return (
    <form className="compound-upload" method="post" action="/uploads/subset_upload_form/" id="own_subset" encType="multipart/form-data" target="_blank">
      <CSRFToken />
      <input type="hidden" name="project" value={this.props.proposal.id} />
      <label htmlFor="id_lib_id">Select library:</label> 
      <select name="lib_id" id="id_lib_id" value={library}>
        <option>...</option>
      {options}
      </select>
  
      <label htmlFor="id_name">Name your selection: </label> 
      <input type="text" name="name" required id="id_name_chpck" value={name}/>

      <label htmlFor="id_data_file">Upload your selection:</label> 
      <input type="file" name="data_file" required id="id_data_file" value={file} />
      <button type="submit" onClick={ ()=> this.submit()}>Upload and add to your selection</button>
    </form>
    );
  }
}

export default CherryPickForm;
