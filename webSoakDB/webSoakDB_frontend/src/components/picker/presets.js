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
			<h2>Presets</h2>
			<p>Specific-purpose compounds selections from in-house libraries</p>
			<form id="properties-form">
				{presets}
			</form>
			<button>Submit</button>
		</section>
		)
	}
}

export default Presets;
