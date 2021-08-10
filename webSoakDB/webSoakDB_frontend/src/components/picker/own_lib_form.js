import React from 'react';
import CSRFToken from './csrf.js';

class OwnLibraryForm extends React.Component {

  constructor(props){
    super(props)
    this.state = {
      //waiting: false, //true while data is reloaded from API after submitting the form
      intervalIds: [],
      reloadCount: 0,
    }
  }
  
  componentDidUpdate(prevProps, prevState){
    if (!this.props.waiting && this.state.intervalIds.length > 0 
      || this.state.reloadCount > 20 ){
      this.stopReloading();
    }

    /*Stop reloading data when:
    1) API data suggests a new library (or subset in case of CherryPickForm) was added to the project
    2) the data was reloaded 30 times
    If the uploaded file is invalid, no new item is added to the project, therefore case 2) is added to 
    eventually stop reloads */
  }

  componentWillUnmount(){
    this.stopReloading();
  }
  
  triggerFormSubmission(unsavedChanges){
    document.forms["own_lib"].requestSubmit();
    const initialLibNumber = this.props.proposal.libraries.length;
    this.refreshUntilAdded(initialLibNumber, unsavedChanges);
  }

  submit() {
      event.preventDefault(); //prevents submitting the form twice
      const changesBeforeSubmission = this.props.detectUnsavedChanges();
      //remember if there were any other unsaved changes before submission

      this.triggerFormSubmission(changesBeforeSubmission);
  }

  refreshUntilAdded(initialLibNumber, changesBeforeSubmission){
    /*reload proposal data from the API every 4s until the newly uploaded 
    library appears in project.libraries. Purpose: the new library will 
    eventually appear in the "already uploaded" list in the UI */
    
    this.props.updateWaitingStatus(true);
    const intervalId = setInterval(()=> {this.reloadData(initialLibNumber, changesBeforeSubmission)}, 4000);
    const intervalIdsCopy = [...this.state.intervalIds, intervalId];
    this.setState({intervalIds: intervalIdsCopy});
  }

  reloadData(initialNumber, changesBeforeSubmission){
    const newReloadCount = this.state.reloadCount + 1;
    this.setState({reloadCount: newReloadCount})
    this.props.refreshAfterUpload();
    //set back to the state before form submission
    this.props.trackUnsavedChanges(changesBeforeSubmission);
    this.checkForChanges(initialNumber);
    
  }

  checkForChanges(initialNumber){
    if (this.props.proposal.libraries.length !== initialNumber){
      this.props.updateWaitingStatus(false);
    }
  }
  
  stopReloading(){
    this.state.intervalIds.forEach( id => clearInterval(id));
    this.setState({intervalIds : [], reloadCount : 0});
    
  }

  render(){
    const name = this.state.name;
    const file = this.state.file;
    return (
      <form className="compound-upload" method="post" action="/uploads/library_upload_form/" id="own_lib" encType="multipart/form-data" target="_blank" >
      <CSRFToken />
      <input type="hidden" name="project" value={this.props.proposal.id} />
        <label htmlFor="id_name">Enter library name: </label> 
        <input type="text" name="name" required id="id_name" value={name}/>
        
        <label htmlFor="id_data_file_lib">Upload plate map:</label> 
        <input type="file" name="data_file" required id="id_data_file_lib"  value={file} />
        <button type="submit" onClick={() => this.submit()}>Upload and add to your selection</button>
      </form>
    )
  }
}

export default OwnLibraryForm;