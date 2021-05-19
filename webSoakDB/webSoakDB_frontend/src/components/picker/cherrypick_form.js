import React from 'react';
import axios from 'axios';
import CSRFToken from './csrf.js';

class CherryPickForm extends React.Component {
  
  constructor(props){
    super(props);
    this.state = {
      libs: [],
      }
  }
  
  componentDidMount() {
    const apiUrl = '/api/in_house_library_list/';
    
    axios.get(apiUrl)
      .then(res => {
      const libs = res.data;
      this.setState({ libs });
      });
  }
  
  submit() {
    const unsavedChanges = this.props.detectUnsavedChanges(null, null);
    document.forms["own_subset"].requestSubmit();
    setTimeout(() => {this.props.refreshAfterUpload()}, 2000);
    this.props.trackUnsavedChanges(unsavedChanges);
  }
  
  render(){
    const options = this.state.libs.map((lib, index) => 
       <option key={index} value={lib.id}>{lib.name}</option>
    );
    
    return (
    <form className="compound-upload" method="post" action="/uploads/subset_upload_form/" id="own_subset" encType="multipart/form-data">
      <CSRFToken />
      <input type="hidden" name="proposal" value={this.props.proposal.proposal} />
      <label htmlFor="id_lib_id">Select library:</label> 
      <select name="lib_id" id="id_lib_id">
        <option>...</option>
      {options}
      </select>
  
      <label htmlFor="id_name">Name your selection: </label> 
      <input type="text" name="name" required id="id_name_chpck" />

      <label htmlFor="id_data_file">Upload your selection:</label> 
      <input type="file" name="data_file" required id="id_data_file" />
      <button type="submit" onClick={ ()=> this.submit()}>Upload and add to your selection</button>
    </form>
    );
  }
}

export default CherryPickForm;
