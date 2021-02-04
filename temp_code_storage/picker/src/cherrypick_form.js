import React from 'react';

class CherryPickForm extends React.Component {
	render(){
		const options = this.props.libs.map((lib) => 
			 <option value="{lib.library.id}">{lib.library.name}</option>
		)
		
		return (
			<form method="post" action="/playground/upload_subset" id="own_subset" enctype="multipart/form-data">
				<p>
					<label for="id_lib_id">Select library:</label> 
					<select name="lib_id" id="id_lib_id">
					{options}
					</select>
				</p>
				<p>
					<label for="id_data_file">Upload your selection:</label> 
					<input type="file" name="data_file" required id="id_data_file" />
				</p>
				<button type="submit">Upload and add to your selection</button>
			</form>
		)
	}
}

export default CherryPickForm;
