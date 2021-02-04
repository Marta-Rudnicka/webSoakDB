import React from 'react';
import OwnLibraryForm from './own_lib_form.js';
import CherryPickForm from './cherrypick_form.js';

class Uploads extends React.Component {


	render(){
		return (
		<section id="upload">
			<h2>Upload your own library</h2>
				<OwnLibraryForm />
			<hr />
			<h2>Upload cherrypicking list</h2>
			<CherryPickForm libs={this.props.libs}/>
		</section>
		)
	}
}

export default Uploads
