import React from 'react';
import axios from 'axios';

import { deepCopyObjectArray, getAttributeArray, mean } from  '../../actions/stat_functions.js';

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
		const allSelectionStats = {libraries : new Set(), plates : 0, compounds : 0, mws : []};
		
		function updateAllSelectionStats(library, NumberOfCompounds, mwArray){
			allSelectionStats.libraries.add(library.name);
			allSelectionStats.plates ++;
			allSelectionStats.compounds = allSelectionStats.compounds + NumberOfCompounds;
			allSelectionStats.mws.push(...mwArray);
		}
		
		const selectedPlates = this.state.selectedPlates;
		const rows = selectedPlates.map((plate, index) => {
			
			const sourceWells = getAttributeArray(plate.compounds, "compound");
			const mwArray = getAttributeArray(sourceWells, "molecular_weight");
			const NumberOfCompounds = plate.compounds.length;
			const mw = mean(mwArray).toFixed(4);
		
			updateAllSelectionStats(plate.library, NumberOfCompounds, mwArray);
			
			//create table rows for specific library plates
			return <tr key={index}>
						<td>{plate.library.name}</td>
						<td>{plate.name}</td>
						<td>{NumberOfCompounds}</td>
						<td>{NumberOfCompounds}</td>
						<td>{mw}</td>
						<td>TODO</td>
						<td>TODO</td>
					</tr>
			});
		
			allSelectionStats.mw = (allSelectionStats.mws.length > 0) ? mean(allSelectionStats.mws).toFixed(4) : "N/A";
		
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
						<th>Mean molecular<br/>weight</th>
						<th>Stat 2</th>
						<th>Stat 3</th>
					</tr>
				</thead>
				<tbody>
					{rows}
					<tr>
						<td colSpan="7"><h3>ALL SELECTION</h3></td>
					</tr>
					<tr>
						<td><strong>Libraries</strong></td>
						<td><strong>Plates</strong></td>
						<td colSpan="2"><strong>Selected compounds</strong></td>
						<td><strong>Mean molecular<br/>weight</strong></td>
						<td>Stat 2</td>
						<td>Stat 3</td>
					</tr>
					<tr>
						<td>{allSelectionStats.libraries.size}</td>
						<td>{allSelectionStats.plates}</td>
						<td colSpan="2">{allSelectionStats.compounds}</td>
						<td>{allSelectionStats.mw}</td>
						<td>TODO</td>
						<td>TODO</td>
						
					</tr>
				</tbody>
			</table>
		</section>
		)
	}
}

export default Stats
