import React from 'react';

import CSRFToken from './csrf.js';

class OwnLibraryForm extends React.Component {
	render(){
		return (
			<form className="compound-upload" method="post" action="uploads/library_upload_form/" id="own_lib" encType="multipart/form-data" target="_blank">
			<CSRFToken />
			<input type="hidden" name="proposal" value={this.props.proposal.name} />
				<label htmlFor="id_name">Enter library name: </label> 
				<input type="text" name="name" required id="id_name" />
				
				<label htmlFor="id_data_file">Upload compound data:</label> 
				<input type="file" name="data_file" required id="id_data_file" />
				<button type="submit">Upload and add to your selection</button>
			</form>
		)
	}
}

export default OwnLibraryForm;
