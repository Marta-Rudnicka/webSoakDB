import React from 'react';
import LibraryOption from './library_option.js';

class Libraries extends React.Component {

	

	render(){
		
		const libraries = this.props.libs.map((lib, index) => <LibraryOption key={index} plate={lib}/>)
		
		return (
		<section id="libraries">
			<h2>XChem in-house fragment libraries</h2>
			
			<form id="libform" >
				<label>SELECT A LIBRARY TO USE IN YOUR EXPERIMENT:</label>	
				<br />
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
