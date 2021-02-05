import React from 'react';
//import './style.css';
import Libraries from './libraries.js';
import Presets from './presets.js';
import Uploads from './uploads.js';
import Stats from './stats.js';
import Graphs from './graphs.js';


//const libs = fetch('http://127.0.0.1:8000/playground/library_list')
//  .then(response => response.json())
//  .then(data => console.log(data));

const libs = [
    {
        "library": {
            "id": 1,
            "name": "DSI_Poised_DMSO",
            "for_industry": true,
            "public": true
        },
        "name": "test_plate"
    },
    {
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate1"
    },
    {
        "library": {
            "id": 3,
            "name": "FragLite",
            "for_industry": false,
            "public": true
        },
        "name": "test plate"
    },
    {
        "library": {
            "id": 7,
            "name": "EU_openscreen",
            "for_industry": false,
            "public": true
        },
        "name": "plate 1 - 9510050007054"
    }
]


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
//const libs = "string"
//function App() {
class Picker extends React.Component {
	
	constructor(props) {
		super(props);
		this.state = {libs: libs, 
					presets: presets,}
	}

	render() {
		return (
		<div id="all">
			<h1>Select compounds - proposal: {proposal}</h1>
			<main id="main-picker">
			
				<Libraries libs={this.state.libs}/>
				<Presets presets={this.state.presets}/>
				<Uploads libs={this.state.libs} />
				<Stats />
				<Graphs />
			</main>
		</div>
	); 
}

}

export default Picker;
