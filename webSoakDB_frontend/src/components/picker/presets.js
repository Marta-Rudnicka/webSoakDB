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
			<p>Collections of compounds that do not cover whole libraries. 
			They can be parts of larger libraries intended for experiments 
			with a smaller number of crystals, or sets of compounds 
			cherry-picked from multiple libraries based on properties</p>
			<form id="properties-form">
				{presets}
			</form>
			<button>Submit</button>
		</section>
		)
	}
}

export default Presets;
