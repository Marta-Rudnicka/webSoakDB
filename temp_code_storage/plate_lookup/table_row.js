import React from 'react';

function showHide(value) {
	if (value){
		return "";
	}
	else {
		return "hidden";
	}	
}

class TableRow extends React.Component {
	
	
	render() {
		const well = showHide(this.props.show_well);
		const code = showHide(this.props.show_code);
		const smiles = showHide(this.props.show_smiles);
		const structure = showHide(this.props.show_structure);
		const concentration = showHide(this.props.show_concentration);
		const mw = showHide(this.props.show_mw);
		const p3 = showHide(this.props.show_p3);
		const p4 = showHide(this.props.show_p4);
	
		return (
			<tr>
				<td>{this.props.counter}</td>
				<td className={well} >{this.props.compound.well}</td>
				<td className={code} > {this.props.compound.compound.code}</td>
				<td className={smiles} >{this.props.compound.compound.smiles}</td>
				<td className={structure}>
					<img alt="alttext" src="../../static/ethanol.png" />
				</td>
				<td className={concentration}> {this.props.compound.concentration}</td>
				<td className={mw}>{this.props.compound.compound.molecular_weight}</td>
				<td className={p3}>TODO</td>
				<td className={p4}>TODO</td>
			
			</tr>
		);
	}
}

export default TableRow;
