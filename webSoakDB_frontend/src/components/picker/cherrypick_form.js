import React from 'react';

class CherryPickForm extends React.Component {
	render(){
		const options = this.props.libs.map((lib, index) => 
			 <option key={index} value="{lib.library.id}">{lib.library.name}</option>
		)
		
		return (
			<form method="post" action="/playground/upload_subset" id="own_subset" encType="multipart/form-data">
				<p>
					<label htmlFor="id_lib_id">Select library:</label> 
					<select name="lib_id" id="id_lib_id">
					{options}
					</select>
				</p>
				<p>
					<label htmlFor="id_data_file">Upload your selection:</label> 
					<input type="file" name="data_file" required id="id_data_file" />
				</p>
				<button type="submit">Upload and add to your selection</button>
			</form>
		)
	}
}

export default CherryPickForm;
