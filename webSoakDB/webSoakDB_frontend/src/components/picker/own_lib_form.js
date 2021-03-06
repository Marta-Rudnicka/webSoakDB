import React from 'react';
import CSRFToken from './csrf.js';

class OwnLibraryForm extends React.Component {
  
  submit() {
    const unsavedChanges = this.props.detectUnsavedChanges();
    document.forms["own_lib"].requestSubmit();
    setTimeout(() => {this.props.refreshAfterUpload();}, 2000);
    this.props.trackUnsavedChanges(unsavedChanges);
  }
  
  render(){
    return (
      <form className="compound-upload" method="post" action="/uploads/library_upload_form/" id="own_lib" encType="multipart/form-data" >
      <CSRFToken />
      <input type="hidden" name="proposal" value={this.props.proposal.proposal} />
        <label htmlFor="id_name">Enter library name: </label> 
        <input type="text" name="name" required id="id_name" />
        
        <label htmlFor="id_data_file_lib">Upload plate map:</label> 
        <input type="file" name="data_file" required id="id_data_file_lib" />
        <button type="submit" onClick={() => this.submit()}>Upload and add to your selection</button>
      </form>
    )
  }
}

export default OwnLibraryForm;
