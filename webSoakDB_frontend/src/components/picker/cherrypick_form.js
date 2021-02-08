import React from 'react';
import axios from 'axios';

class CherryPickForm extends React.Component {
	
	constructor(props){
		super(props);
		this.state = {
			libs: [],
			}
	}
	
	componentDidMount() {
		const apiUrl = 'api/in_house_library_list/';
		
		axios.get(apiUrl)
			.then(res => {
			const libs = res.data;
			this.setState({ libs });
      });
	}
	
	render(){
		const options = this.state.libs.map((lib, index) => 
			 <option key={index} value={lib.id}>{lib.name}</option>
		);
		
		return (
			<form method="post" action="/playground/upload_subset" id="own_subset" encType="multipart/form-data">
				<p>
					<label htmlFor="id_lib_id">Select library:</label> 
					<select name="lib_id" id="id_lib_id">
						<option>...</option>
					{options}
					</select>
				</p>
				<p>
					<label htmlFor="id_data_file">Upload your selection:</label> 
					<input type="file" name="data_file" required id="id_data_file" />
				</p>
				<button type="submit">Upload and add to your selection</button>
			</form>
		);
	}
}

export default CherryPickForm;
