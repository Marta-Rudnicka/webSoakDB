import React from 'react';
import OwnLibraryForm from './own_lib_form.js';
import CherryPickForm from './cherrypick_form.js';

class Uploads extends React.Component {

	render(){
		let uploadedLibs = []
		this.props.proposal.libraries.map(lib => {
			if (!lib.public){
				uploadedLibs.push(lib);
				}
			});
		
		const liblist = () => { 
			if (uploadedLibs.length > 0){
				return(
				<React.Fragment> 
					<h3>Already uploaded:</h3>
					<ul className="old-uploads">
						 {uploadedLibs.map((lib, index) => <li key={index}>{lib.name}</li>)}
					</ul>
				</React.Fragment>
				); 
				}
			else {
				return "";
			}
		};
		
		
		return (
		<section id="upload">
			<h2>Upload your own library</h2>
				<OwnLibraryForm proposal={this.props.proposal}/>
				{liblist()}
			<hr />
			<h2>Upload cherrypicking list</h2>
			<CherryPickForm libs={this.props.libs}/>
		</section>
		)
	}
}

export default Uploads
