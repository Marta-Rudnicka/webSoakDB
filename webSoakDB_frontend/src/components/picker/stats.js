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
			selectedSubsets: [],
			content: [],
			};
	}
	
	calculate(){
		/* downloads compound data and generates the contents of the stats table */
		console.log('fired calculate()')
		let plates = [];
		let subsets = [];
		this.props.selectedLibIds.forEach(id => {
			const apiUrl = 'api/current_plates_stats/' + id + '/';		
			axios.get(apiUrl)
				.then(res => {
					const addedPlate = res.data;
					plates = deepCopyObjectArray(this.state.selectedPlates);
					plates.push(...addedPlate);
					this.setState({selectedPlates: plates});
					this.setState({content: this.generateContent()});
			});
		});
		
		this.props.selectedSubsetIds.forEach(id => {
			console.log('id: ', id)
			const apiUrl = 'api/subset_stats/' + id + '/';		
			axios.get(apiUrl)
				.then(res => {
					const addedSubset = res.data;
					console.log('addedSubset: ', addedSubset )
					console.log('addedSubset: ', addedSubset )
					subsets = deepCopyObjectArray(this.state.selectedSubsets);
					subsets.push(addedSubset);
					this.setState({selectedSubsets: subsets});
					this.setState({content: this.generateContent()});
			});
		});
		
		/* to improve performance with larger data sets */
		this.setState({selectedPlates: []});
		this.setState({selectedSubsets: []});
		
	}
	
	componentDidUpdate(prevProps) {
 
		if (this.props.selectedLibs !== prevProps.selectedLibs) {
			this.setState({selectedPlates: []});
			this.componentDidMount();
		}
	}
	
	generateContent(){
		console.log('fired generateContent()')
		const allSelection = {libraries : new Set(), plates : new Set (), compounds : [] };
		
		//FULL LIBRARIES
		const selectedPlates = this.state.selectedPlates;
	
		const libRows = selectedPlates.map((plate, index) => {
			
			//get sums and means for each property			
			const stats = get_stats(plate.compounds, "plate", dict);
			const NumberOfCompounds = plate.compounds.length;
			
			//add the plate data to the stats of the whole selection
			updateAllSelection(plate.library.id, plate.id, getAttributeArray(plate.compounds, "compound"), allSelection)
			
			//create cells with mean values
			const stat_cells = descriptor_names.map((string, index) => {
					return <td key={index}>{stats[string]}</td>
				});
			
			//create table row for the specific library plate
			return <tr key={index}>
						<td>{plate.library.name}</td>
						<td>{plate.name}</td>
						<td>{NumberOfCompounds} (all)</td>
						{stat_cells}
					</tr>
		});
		
		//SUBSETS
		const selectedSubsets = this.state.selectedSubsets;
		const subsetRows = selectedSubsets.map((subset, index) => {
			
			//get sums and means for each property			
			const stats = get_stats(subset.compounds, "subset", dict);
			//const selectedCompounds = subset.compounds.length;
			
			//add the plate data to the stats of the whole selection
			updateAllSelection(subset.library.id, null, subset.compounds, allSelection)
			
			//create cells with mean values
			const stat_cells = descriptor_names.map((string, index) => {
					return <td key={index}>{stats[string]}</td>
				});
			
			//create table row for the specific library plate
			return <tr key={index} className="subset-row">
						<td>{subset.library.name}</td>
						<td>N/A</td>
						<td>{subset.compounds.length}</td>
						{stat_cells}
					</tr>
		});
			
		//generate page content for overall selection
		allSelection.stats = get_stats(allSelection.compounds, "all", dict)
			
		const all_stat_cells = descriptor_names.map((string, index) => {
			return <td key={index}>{allSelection.stats[string]}</td>
		});
			
		const sums = <React.Fragment><td>{allSelection.libraries.size}</td><td>{allSelection.plates.size}</td><td>{allSelection.compounds.length}</td></React.Fragment>
		
		return [libRows, subsetRows, sums, all_stat_cells]
	
	}
	

	render(){
			
			const libRows = this.state.content[0];
			const subsetRows = this.state.content[1];
			const sums = this.state.content[2];
			const all_stats = this.state.content[3];
						
		return (
		<section id="stats">
			<h2>Statistics for the current selection</h2>
			<button id="stat-calculator" onClick={event => this.calculate()}>(Re)calculate</button>
			<table id="library-stats">
				<thead>
					<tr>
						<th>Library</th>
						<th>Plate</th>
						<th>Selected <br />compounds</th>
						<StatHeaders />						
					</tr>
				</thead>
				<tbody>
					{libRows}
					{subsetRows}
					<tr>
						<td colSpan="16"><h3>ALL SELECTION</h3></td>
					</tr>
					<tr>
						<td><strong>Libraries</strong></td>
						<td><strong>Plates</strong></td>
						<td><strong>Selected compounds</strong></td>
						<StatHeaders />
					</tr>
					<tr>
						{sums}
						{all_stats}
					</tr>
				</tbody>
			</table>
		</section>
		)
	}
}

export default Stats
