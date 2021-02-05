import React from 'react';

class PresetOption extends React.Component {


	render(){
		const id = this.props.id
		const name = this.props.name
		const desc = this.props.description
		return (
			<div className="preset">
				<input type="checkbox" value="{id}" name="preset_ids" />
				<label htmlFor="{id}">{name} </label>
					<div className="show-on-hover">Description...
						
					</div>
					<p>{desc}</p>
			</div>
		)
	}
}

export default PresetOption;
