import React from 'react';

class LibraryOption extends React.Component {

	render(){
		const lib = this.props.plate.library;
		const plate = this.props.plate;
		
		return (
			<div className="library">
				<input type="checkbox" 
					value={this.props.plate.library.id} 
					name="lib_ids" defaultChecked={this.props.defaultChecked} 
					onChange={event => this.props.handleCheckboxChange(event)} 
				/>
				<label htmlFor={this.props.plate.library.id}>{lib.name} </label>
				<br />
				<span className="pseudo-link" onClick={event => this.props.showPlate(lib, plate)} >See compounds</span>
			</div>
		)
	}
}

export default LibraryOption;
