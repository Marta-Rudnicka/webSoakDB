import React from 'react';
import {display_options} from './display_options.js';
import {descriptor_names} from '../picker/stats_helpers.js';

function showHide(value) {
	if (value){
		return "";
	}
	else {
		return "hidden";
	}	
}

const stat_options = display_options.slice(5)

class TableRow extends React.Component {
	
	
	render() {
		
		let classes = {};
		
		display_options.forEach(option => {
			classes[option[0]] = showHide(this.props.display[option[0]]);
		});
		
	//	console.log('classes', classes);
		
		const stat_cells = descriptor_names.map((option, index) => {
			let value = this.props.compound.compound.properties[option];
			if (value % 1 !== 0){
				value = value.toFixed(2);
			}
			
			return <td key={index} className={classes[stat_options[index][0]]}>{value}</td>
		});
		
		//console.log(stat_cells)
		
		return (
			<tr>
				<td>{this.props.counter}</td>
				<td className={classes.show_well} >{this.props.compound.well}</td>
				<td className={classes.show_code} > {this.props.compound.compound.code}</td>
				<td className={classes.show_smiles} >{this.props.compound.compound.smiles}</td>
				<td className={classes.show_structure}>
					[PICTURE GOES HERE]
				</td>
				<td className={classes.show_concentration}>{this.props.compound.concentration}</td>
				{stat_cells}
			
			</tr>
		);
	}
}

export default TableRow;
