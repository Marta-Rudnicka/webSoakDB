import React from 'react';
import axios from 'axios';
import StatHeaders from './stat_headers.js';

import { deepCopyObjectArray, getAttributeArray, mean } from  '../../actions/stat_functions.js';
import { descriptor_names, get_stats, updateAllSelection, dict } from './stats_helpers.js';

class Stats extends React.Component {
	constructor(props) {
		super(props);
		this.state = { 	
			selectedPlates: [],
			};
	}
	
	componentDidMount(){
		let plates = [];
		this.props.selectedLibs.map(library => {
			const apiUrl = 'api/current_plates_stats/' + library.id + '/';		
			axios.get(apiUrl)
				.then(res => {
					const addedPlate = res.data;
					plates = deepCopyObjectArray(this.state.selectedPlates);
					plates.push(...addedPlate);
					this.setState({selectedPlates: plates});
			});
		});
	}
	
	componentDidUpdate(prevProps) {
 
		if (this.props.selectedLibs !== prevProps.selectedLibs) {
			this.setState({selectedPlates: []});
			this.componentDidMount();
		}
	}
	

	render(){
		const allSelection = {libraries : new Set(), plates : new Set (), compounds : [] };

		const selectedPlates = this.state.selectedPlates;
		const rows = selectedPlates.map((plate, index) => {
			
			//get sums and means for each property			
			const stats = get_stats(plate.compounds, dict);
			const NumberOfCompounds = plate.compounds.length;
			
			//add the plate data to the stats of the whole selection
			updateAllSelection(plate.library.id, plate.id, plate.compounds, allSelection)
			
			//create cells with mean values
			const stat_cells = descriptor_names.map((string, index) => {
					return <td key={index}>{stats[string]}</td>
				});
			
			//create table row for the specific library plates
			return <tr key={index}>
						<td>{plate.library.name}</td>
						<td>{plate.name}</td>
						<td>{NumberOfCompounds}</td>
						<td>{NumberOfCompounds}</td>
						{stat_cells}
					</tr>
			});
			
			allSelection.stats = get_stats(allSelection.compounds, dict)
			
			const all_stat_cells = descriptor_names.map((string, index) => {
				return <td key={index}>{allSelection.stats[string]}</td>
			});
			
		return (
		<section id="stats">
			<h2>Statistics for the current selection</h2>
			<table id="library-stats">
				<thead>
					<tr>
						<th>Library</th>
						<th>Plate</th>
						<th>Compounds</th>
						<th>Selected <br />compounds</th>
						<StatHeaders />						
					</tr>
				</thead>
				<tbody>
					{rows}
					<tr>
						<td colSpan="17"><h3>ALL SELECTION</h3></td>
					</tr>
					<tr>
						<td><strong>Libraries</strong></td>
						<td><strong>Plates</strong></td>
						<td colSpan="2"><strong>Selected compounds</strong></td>
						<StatHeaders />
					</tr>
					<tr>
						<td>{allSelection.libraries.size}</td>
						<td>{allSelection.plates.size}</td>
						<td colSpan="2">{allSelection.compounds.length}</td>
						{all_stat_cells}
					</tr>
				</tbody>
			</table>
		</section>
		)
	}
}

export default Stats
