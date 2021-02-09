import React from 'react';

class OwnLibraryForm extends React.Component {
	render(){
		return (
			<form className="compound-upload" method="post" action="data_test/" id="own_lib" encType="multipart/form-data">
				<p>
					<label htmlFor="id_name">Enter library name: </label> 
					<input type="text" name="name" required id="id_name" />
				</p>
				<p>				
					<label htmlFor="id_data_file">Upload compound data:</label> 
					<input type="file" name="data_file" required id="id_data_file" />
				</p>
				<button type="submit">Upload and add to your selection</button>
			</form>
		)
	}
}

export default OwnLibraryForm;
