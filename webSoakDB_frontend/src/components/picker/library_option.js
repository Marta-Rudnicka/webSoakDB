import React from 'react';

class LibraryOption extends React.Component {

	render(){
		const libname = this.props.plate.library.name;
		const platename = this.props.plate.name;
		
		return (
			<p className="library">
				<input type="checkbox" value="this.props.plate.id" name="lib_ids" defaultChecked={this.props.defaultChecked}/>
				<label htmlFor="{this.props.plate.id}">{libname} </label>
				<br />
				<span className="pseudo-link" onClick={event => this.props.showPlate(libname, platename, true)} >See compounds</span>
			</p>
		)
	}
}

export default LibraryOption;
