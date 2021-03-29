import React from 'react';

class PresetOption extends React.Component {


	render(){
		const id = this.props.preset.id
		const name = this.props.preset.name
		const desc = this.props.preset.description
		return (
			<div className="preset">
				<input 
					type="checkbox" 
					value={id} name="preset_ids" 
					onChange={event => this.props.handleCheckboxChange(event, this.props.unsaved)} 
					defaultChecked={this.props.defaultChecked}/>
														   			 
				<label htmlFor="{id}">{name} </label>
				<br/>
				<span className="pseudo-link" onClick={() => this.props.showPlate(null, this.props.preset, false, true)}>See compounds</span>
					<div className="show-on-hover">Description...
						
					</div>
					<p>{desc}</p>
			</div>
		)
	}
}

export default PresetOption;