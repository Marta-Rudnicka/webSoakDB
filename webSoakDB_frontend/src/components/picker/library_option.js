import React from 'react';

class LibraryOption extends React.Component {

	render(){
		const libname = this.props.plate.library.name;
		const platename = this.props.plate.name;
		
		return (
			<div className="library">
				<input type="checkbox" value={this.props.plate.library.id} name="lib_ids" defaultChecked={this.props.defaultChecked} onChange={event => this.props.handleCheckboxChange(event)} />
				<label htmlFor={this.props.plate.library.id}>{libname} </label>
				<br />
				<span className="pseudo-link" onClick={event => this.props.showPlate(libname, platename, true)} >See compounds</span>
			</div>
		)
	}
}

export default LibraryOption;
