import React from 'react';
import axios from 'axios';
import StatHeaders from './stat_headers.js';

import { deepCopyObjectArray, getAttributeArray, mean } from  '../../actions/stat_functions.js';
import { descriptor_names, get_stats, updateAllSelection, dict } from './stats_helpers.js';

class Stats extends React.Component {
	constructor(props) {
		super(props);
		this.state = {
			libraries: [], 	
			libraryStats: [],
			subsetStats: [],
			content: [],
			disabled : false,
			};
	}
	
	calculate(){
		/* downloads compound data and generates the contents of the stats table */
		this.setState({disabled : true});
		let libraries = [];
		let subsets = [];
		this.state.libraries.forEach(library => {	
			const apiUrl = '/api/current_plates_stats/' + library.id + '/';
			axios.get(apiUrl)
				.then(res => {
					const compounds = res.data;
					let stats =  JSON.parse(JSON.stringify(this.state.libraryStats));
					const addedLib = {id : library.id, name : library.name, compounds : compounds}
					libraries.push(addedLib);
					this.setState({libraryStats: libraries});
					this.setState({content: this.generateContent()});
				})
		});		
		
		this.props.selectedSubsetIds.forEach(id => {
			const apiUrl = '/api/subset_stats/' + id + '/';		
			axios.get(apiUrl)
				.then(res => {
					const addedSubset = res.data;
					let subsets = deepCopyObjectArray(this.state.subsetStats);
					subsets.push(addedSubset);
					this.setState({subsetStats: subsets});
					this.setState({content: this.generateContent()});
			});
		});
	}
	
	uploadLibraryData(){
		let libs = []
		this.props.selectedLibIds.forEach(id => {
			const apiUrl = '/api/library_detail/' + id + '/';
			axios.get(apiUrl)
				.then(res => {
					const newLib = res.data;
					libs.push(newLib)
				});
			});
		this.setState({libraries : libs});
	}
	
	componentDidMount(){
		this.uploadLibraryData();
	}
	
	componentDidUpdate(prevProps) {
		if (this.props.selectedLibIds !== prevProps.selectedLibIds) {
			this.uploadLibraryData();
		}
	}
	
	generateContent(){
		const allSelection = {libraries : new Set(), compounds : [] };
		let libRows = null
		
		if(this.state.libraryStats){
			libRows = this.state.libraryStats.map((lib, index) => {
				
				//get sums and means for each property		
				const stats = get_stats(lib.compounds, "subset", dict);
				const NumberOfCompounds = lib.compounds.length;
				
				//add the plate data to the stats of the whole selection
				updateAllSelection(lib.id, lib.compounds, allSelection)
				
				//create cells with mean values
				const stat_cells = descriptor_names.map((string, index) => {
						return <td key={index}>{stats[string]}</td>
					});
				
				//create table row for the specific library plate
				return <tr key={index}>
							<td>{lib.name}</td>
							<td>{NumberOfCompounds} (all)</td>
							{stat_cells}
						</tr>
			});
		}
		
		
		//SUBSETS
		const subsetStats = this.state.subsetStats;
		let subsetRows = null;
		if(subsetStats){
			subsetRows = subsetStats.map((subset, index) => {
				
				//get sums and means for each property			
				const stats = get_stats(subset.compounds, "subset", dict);
				
				//add the plate data to the stats of the whole selection
				updateAllSelection(subset.library.id, subset.compounds, allSelection)
				
				//create cells with mean values
				const stat_cells = descriptor_names.map((string, index) => {
						return <td key={index}>{stats[string]}</td>
					});
				
				//create table row for the specific library plate
				return <tr key={index} className="subset-row">
							<td>{subset.library.name}</td>
							<td>{subset.compounds.length}</td>
							{stat_cells}
						</tr>
			});
		}
		//generate page content for overall selection
		allSelection.stats = get_stats(allSelection.compounds, "all", dict)
			
		const all_stat_cells = descriptor_names.map((string, index) => {
			return <td key={index}>{allSelection.stats[string]}</td>
		});
			
		const sums = <React.Fragment><td>{allSelection.libraries.size}</td><td>{allSelection.compounds.length}</td></React.Fragment>
		
		//if done, flush the data (improves performance)
		if (this.state.subsetStats.length === this.props.selectedSubsetIds.length && 
				this.state.libraryStats.length === this.state.libraries.length){
				this.setState({subsetStats : [], libraryStats : [], disabled : false});
		}
		
		return [libRows, subsetRows, sums, all_stat_cells]
	}
	

	render(){
		
		let className = null;
		const items = this.props.selectedLibIds.length + this.props.selectedSubsetIds.length
		const stats = this.state.libraryStats.length + this.state.subsetStats.length
		
		if(stats > 0 && items > stats){
			className = 'changed';
		}
			
		const libRows = this.state.content[0];
		const subsetRows = this.state.content[1];
		const sums = this.state.content[2];
		const all_stats = this.state.content[3];
						
		return (
		<section id="stats" className={className}>
			<h2>Statistics for the current selection (including unsaved items)</h2>
			<button id="stat-calculator" onClick={event => this.calculate()} disabled={this.state.disabled}>(Re)calculate</button>
			<table id="library-stats" className="table table-bordered" >
				<thead>
					<tr>
						<th>Library</th>
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
						<td><strong>Unique selected compounds</strong></td>
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
