import React from 'react';
import LibraryOption from './library_option.js';
import axios from 'axios';

class Libraries extends React.Component {

	constructor(props) {
		super(props);
		this.state = {
			libs: [],
		};
	}
	
	componentDidMount() {
		const apiUrl = 'api/library_selection_list/';
		
		axios.get(apiUrl)
			.then(res => {
			const libs = res.data;
			this.setState({ libs });
      });
	}
	
	isInProposal(plate, libArray){
		if (libArray){
			const plateLibId = plate.library.id;
			const arrayIds = [];
			libArray.map(lib => arrayIds.push(lib.id));
			return arrayIds.includes(plateLibId);
		}
		return false;
	}
	render(){
		
		const proposalLibs =  this.props.proposal.libraries;
		
		const libraries = this.state.libs.map((lib, index) => { 
			return <LibraryOption 
				key={index} 
				plate={lib} 
				showPlate={this.props.showPlate}
				defaultChecked={this.isInProposal(lib, proposalLibs)}
				/>;
		});
		
		return (
		<section id="libraries">
			<h2>XChem in-house fragment libraries</h2>
			
			<form id="libform" >
				<div id="libs">
					{libraries}
				</div>
				<button type="submit">Add selected to your collection</button>
			</form>
		</section>
		)
	}
}

export default Libraries;
