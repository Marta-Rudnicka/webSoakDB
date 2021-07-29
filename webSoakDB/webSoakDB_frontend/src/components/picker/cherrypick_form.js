import React from 'react';
import axios from 'axios';
import CSRFToken from './csrf.js';

class CherryPickForm extends React.Component {
  
  constructor(props){
    super(props);
    this.state = {
      //libs: [],
      name: "",
      file: "",
      library: "",
      }
  }
  
  /*
  componentDidMount() {
    const apiUrl = '/api/in_house_library_list/';
    
    axios.get(apiUrl)
      .then(res => {
      const libs = res.data;
      this.setState({ libs });
      });
  }*/
  

  changeName(e){
    this.setState({name: e.target.value})
  }

  changeFile(e){
    this.setState({file: e.target.value})
  }

  changeLibrary(e){
    this.setState({library: e.target.value})
  }
  submit() {
    if (!(this.state.name && this.state.file && this.state.library)){
      return;
    }
    this.props.showOverlay();
    const unsavedChanges = this.props.detectUnsavedChanges(null, null);
    document.forms["own_subset"].requestSubmit();
    setTimeout(() => {this.props.refreshAfterUpload()}, 2000);
    this.props.trackUnsavedChanges(unsavedChanges);
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
      <select name="lib_id" id="id_lib_id" value={library} onChange={(e) => this.changeLibrary(e)}>
        <option>...</option>
      {options}
      </select>
  
      <label htmlFor="id_name">Name your selection: </label> 
      <input type="text" name="name" required id="id_name_chpck" value={name} onChange={(e) => this.changeName(e)}/>

      <label htmlFor="id_data_file">Upload your selection:</label> 
      <input type="file" name="data_file" required id="id_data_file" value={file} onChange={(e) => this.changeFile(e)} />
      <button type="submit" onClick={ ()=> this.submit()}>Upload and add to your selection</button>
    </form>
    );
  }
}

export default CherryPickForm;
