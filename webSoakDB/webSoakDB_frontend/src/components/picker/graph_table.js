import React from 'react';
//import {Link} from "react-router-dom";
//import { ChevronDown } from '../icons.js';
//import Details from './details.js';
import properties_dict from './properties_dict.js'
import Graph from './graph.js';
import { removeFromArray, getSubsetIds } from  '../../actions/stat_functions.js';

class GraphTr extends React.Component {
	render(){
		
		const cells = this.props.properties.map((p, index) => {
			return <td key={index}><Graph id={this.props.col.id} type={this.props.type} property={p} caption={properties_dict[p]} /></td>
		});
		
		return(
		<tr>
			<td><strong>{this.props.col.name}</strong></td>
			{cells}
		</tr>
		);
	}
}

class GraphTable extends React.Component {
	
	constructor(props){
		super(props);
		this.state = {
			show: ["mol_wt", "log_p", "num_h_donors", "num_h_acceptors", "num_rot_bonds"],
			//libs: [],
			libs: this.selectedInHouseLibs(),
			presets: this.selectedPresets(),
		}
	}
	
	componentDidUpdate(prevProps){
		if(this.props.parentState !== prevProps.parentState ){
			this.setState({libs : this.selectedInHouseLibs(), presets : this.selectedPresets()});
		}
	}
	
	manageProperties(event){
		const p = event.target.value;
		
		if(event.target.checked === true){
			const newShow = this.state.show.slice(0, this.state.show.length);
			newShow.push(p);
			this.setState({show : newShow});
		}
		else{
			const newShow = removeFromArray(this.state.show, [p]);
			this.setState({show : newShow});
		}
	}
	
	selectedInHouseLibs(){
		const output = [];
		this.props.parentState.currentLibOptions.forEach(l =>{
			if (this.props.parentState.selectedLibIds.includes(l.id)){
				output.push(l);
			}
		})
		return output;
	}
	
	selectedPresets(){
		console.log('firing selectedPresets', this.props.parentState.presets);
		console.log('this.props.parentState', this.props.parentState);
		const output = [];
		this.props.parentState.presets.forEach(p =>{
			console.log('checking: ', p);
			const subs = getSubsetIds(this.props.parentState, p.id)
			console.log('subs: ', subs);
			if (this.props.parentState.selectedSubsetIds.includes(subs[0])){
				output.push(p);
			}
		})
		console.log('Selected presets output: ', output)
		return output;
	}
	
	render(){
		const checkboxes = Object.keys(properties_dict).map((key, index) => {
				let checked = false;
				//if (["mol_wt", "tpsa", "log_p"].includes(key)){
				if (this.state.show.includes(key)){
					checked = true;
				}
				return (
				<div key={index}>
					<input type="checkbox" value={key} id={key} onChange={e => this.manageProperties(e)} checked={checked}/>
					<label htmlFor={key}>{properties_dict[key]}</label>
				</div>
				)
			});
		
		const headers = this.state.show.map((prop, index) =>{
			return <th key={index}>{properties_dict[prop]}</th>
		});
		
		const l_rows = this.state.libs.map((lib, index) => {
			return <GraphTr key={index} type="library" col={lib} properties={this.state.show}/>
		
		});
		const s_rows = this.state.presets.map((p, index) => {
			return <GraphTr key={index} type="preset" col={p} properties={this.state.show}/>
		
		});
		
		return(
		<div>
			<h3>Select properties to include: </h3>
			<div id="properties">
				{checkboxes}
			</div>
			<table className="table">
				<thead>
					<tr>
						<th>Collection</th>
						{headers}
					</tr>
				</thead>
				<tbody>
					{l_rows}
					{s_rows}
				</tbody>
			</table>
		</div>
		
		
		)
	}

}

export default GraphTable;
