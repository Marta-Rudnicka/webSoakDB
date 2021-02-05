import React from 'react';

class LibraryOption extends React.Component {

	render(){
		const libname = this.props.plate.library.name
		const platename = this.props.plate.name
		const url = 'localhost:8000/playground/' + libname + '/' + platename
		return (
			<p className="library">
				<input type="checkbox" value="{this.props.plate.id}" name="lib_ids" />
				<label htmlFor="{this.props.plate.id}">{libname} </label>
				<br />
				<a href={url}>See compounds</a>
			</p>
		)
	}
}

export default LibraryOption;
