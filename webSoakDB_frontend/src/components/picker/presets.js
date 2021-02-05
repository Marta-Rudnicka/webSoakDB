import React from 'react';
import PresetOption from './preset_option.js';

class Presets extends React.Component {
	
	render(){
		const presets = this.props.presets.map((preset, index) => 
		<PresetOption
			key={index}
			id={preset.id}
			name = {preset.name}
			description = {preset.description}
		/>
	)
		return (
		<section id="presets">
			<h2>Select a preset</h2>
			<form id="properties-form">
				{presets}
			</form>
			<button>Submit</button>
		</section>
		)
	}
}

export default Presets;
