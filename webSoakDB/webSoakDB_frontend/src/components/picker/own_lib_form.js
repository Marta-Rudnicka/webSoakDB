import React from 'react';
import CSRFToken from './csrf.js';

class OwnLibraryForm extends React.Component {

  constructor(props){
    super(props)
    this.state = {
      name: "",
      file: "",
    }
  }
  
  changeName(e){
    this.setState({name: e.target.value})
  }

  changeFile(e){
    this.setState({file: e.target.value})
  }

  submit() {
    if (!window.confirm("Are you sure you want to submit the data despite missing security features in the application demo?")){
      event.preventDefault();
    }
    if (!(this.state.name && this.state.file)){
      return;
    }
    else{
      this.props.showOverlay();
      const unsavedChanges = this.props.detectUnsavedChanges();
      document.forms["own_lib"].requestSubmit();
      setTimeout(() => {this.props.refreshAfterUpload();}, 2000);
      this.props.trackUnsavedChanges(unsavedChanges);
    }
  }
  
  render(){
    const name = this.state.name;
    const file = this.state.file;
    return (
      <form className="compound-upload" method="post" action="/uploads/library_upload_form/" id="own_lib" encType="multipart/form-data" >
      <CSRFToken />
      <input type="hidden" name="proposal" value={this.props.proposal.proposal} />
        <label htmlFor="id_name">Enter library name: </label> 
        <input type="text" name="name" required id="id_name" value={name} onChange={(e) => this.changeName(e)}/>
        
        <label htmlFor="id_data_file_lib">Upload plate map:</label> 
        <input type="file" name="data_file" required id="id_data_file_lib"  value={file} onChange={(e) => this.changeFile(e)} />
        <button type="submit" onClick={() => this.submit()}>Upload and add to your selection</button>
      </form>
    )
  }
}

export default OwnLibraryForm;
