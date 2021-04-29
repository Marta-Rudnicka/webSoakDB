import React from 'react';
import {Link} from "react-router-dom";


class LibraryOption extends React.Component {

	render(){
		const lib = this.props.plate.library;
		const plate = this.props.plate;
		
		return (
			<div className="library">
				<input type="checkbox" 
					value={this.props.plate.library.id} 
					name="lib_ids" defaultChecked={this.props.defaultChecked} 
					onChange={event => this.props.handleCheckboxChange(event, this.props.unsaved)} 
				/>
				<label htmlFor={this.props.plate.library.id}>{lib.name} </label>
				<br />
				<Link to={"../../compounds/plate/" + this.props.plate.id + "/"} >See compounds</Link>
			</div>
		)
	}
}

//<span className="pseudo-link" onClick={event => this.props.showPlate(lib, plate, true, false)} >See compounds</span>
export default LibraryOption;
