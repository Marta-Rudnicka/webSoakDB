import React from 'react';

class OwnLibraryForm extends React.Component {
	render(){
		return (
			<form method="post" action="/playground/upload_user_library" id="own_lib" encType="multipart/form-data">
				<p>
					<label htmlFor="id_name">Library name:</label> 
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
