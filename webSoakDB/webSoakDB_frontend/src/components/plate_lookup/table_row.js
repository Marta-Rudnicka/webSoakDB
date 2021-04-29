import React, { Suspense } from 'react';
import {display_options} from './display_options.js';
import {descriptor_names} from '../picker/stats_helpers.js';
import StructurePic from './structure_pic.js';
//import {ExpandIcon} from '../icons.js';


function showHide(value) {
	if (value){
		return "";
	}
	else {
		return "hidden";
	}	
}

const stat_options = display_options.slice(5)
/*
const 

function StructurePicLoader(props) {
	return(	
		<Suspense fallback={<div>Loading image...</div>}>
			<StructurePic smiles={props.smiles} />
		</Suspense>
	);
}
*/

class TableRow extends React.Component {
	
	
	render() {
		let url_smiles = "";
		const img_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/' + url_smiles + '/PNG'
		//let img_url = 'https://cactus.nci.nih.gov/chemical/structure/' + url_smiles + '/image'
		let classes = {};
		
		display_options.forEach(option => {
			classes[option[0]] = showHide(this.props.display[option[0]]);
		});
		
		let general_cells;
		
		if (this.props.is_a_plate){
			url_smiles = this.props.compound.compound.smiles.replace('#', '%23');
			general_cells = (
				<React.Fragment>
					<td className={classes.show_well} >{this.props.compound.well}</td>
					<td className={classes.show_code} > {this.props.compound.compound.code}</td>
					<td className={classes.show_smiles} >{this.props.compound.compound.smiles}</td>
					<td className={classes.show_structure}>
						<StructurePic smiles={this.props.compound.compound.smiles} />
					</td>
					<td className={classes.show_concentration}>{this.props.compound.concentration}</td>
				</React.Fragment>
			)
		}
		else {
			url_smiles = this.props.compound.smiles.replace('#', '%23');
			general_cells = (
				<React.Fragment>
					
					<td className={classes.show_code} > {this.props.compound.code}</td>
					<td className={classes.show_smiles} >{this.props.compound.smiles}</td>
					<td className={classes.show_structure}>
						<StructurePic smiles={this.props.compound.smiles} />
					</td>
				</React.Fragment>
			)
		}
		
		const stat_cells = descriptor_names.map((option, index) => {
			let value 
			if (this.props.is_a_plate) {
		
				value = this.props.compound.compound[option];
			}
			else {
				value = this.props.compound[option];
			}
			if (value % 1 !== 0){
				value = value.toFixed(2);
			}
			
			return <td key={index} className={classes[stat_options[index][0]]}>{value}</td>
		});
		
		let library_cell = null;
		
		if (this.props.is_a_preset){
			library_cell = <td className={classes.show_library}>{(this.props.compound.library !== undefined) ? this.props.compound.library : '...'}</td>
		}
		
		return (
			<tr>
				<td>{this.props.counter}</td>
				{library_cell}
				{general_cells}
				{stat_cells}
			</tr>
		);
	}
}

export default TableRow;

//<img className="img2d" src={img_url} />
