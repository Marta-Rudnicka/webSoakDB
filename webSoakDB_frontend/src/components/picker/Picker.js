import React from 'react';
import Libraries from './libraries.js';
import Presets from './presets.js';
import Uploads from './uploads.js';
import Stats from './stats.js';
import Graphs from './graphs.js';
import axios from 'axios';

const presets = [
	{
		"id": 1,
		"name": "Preset 1",
		"description": "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed nec tempus libero, quis porta mi. ",
	},
	{
		"id": 2,
		"name": "Preset 2",
		"description": "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed nec tempus libero, quis porta mi. ",
	}
]

const proposal = 'Test string';

class Picker extends React.Component {
	
	constructor(props) {
		super(props);
		this.state = {
					presets: presets,
		}
	}
	render() {
		return (
		<div id="picker">
			<h1>Select compounds for {this.props.proposal.name}</h1>
			<main id="main-picker">
				<Libraries showPlate={this.props.showPlate}  proposal={this.props.proposal}/>
				<Presets presets={this.state.presets}  proposal={this.props.proposal}/>
				<Uploads proposal={this.props.proposal}/>
				<Stats />
				<Graphs />
			</main>
		</div>
	); 
}

}

export default Picker;
