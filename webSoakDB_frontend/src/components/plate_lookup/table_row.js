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
		const tpsa = showHide(this.props.show_tpsa);
		const logp = showHide(this.props.show_logp);
		const heavy_atom_count = showHide(this.props.show_heavy_atom_count);
		const heavy_atom_mol_wt = showHide(this.props.show_heavy_atom_mol_wt);
		const nhoh_count = showHide(this.props.show_nhoh_count);
		const no_count = showHide(this.props.show_no_count);
		const num_h_acceptors = showHide(this.props.show_num_h_acceptors);
		const num_h_donors = showHide(this.props.show_num_h_donors);
		const num_het_atoms = showHide(this.props.show_num_het_atoms);
		const num_rot_bonds = showHide(this.props.show_num_rot_bonds);
		const num_val_electrons = showHide(this.props.show_num_val_electrons);
		const ring_count = showHide(this.props.show_ring_count);
	
		return (
			<tr>
				<td>{this.props.counter}</td>
				<td className={well} >{this.props.compound.well}</td>
				<td className={code} > {this.props.compound.compound.code}</td>
				<td className={smiles} >{this.props.compound.compound.smiles}</td>
				<td className={structure}>
					[PICTURE GOES HERE]
				</td>
				<td className={concentration}> {this.props.compound.concentration}</td>
				<td className={mw}>{this.props.compound.compound.properties.mol_wt.toFixed(4)}</td>
				<td className={tpsa}>{this.props.compound.compound.properties.tpsa.toFixed(2)}</td>
				<td className={logp}>{this.props.compound.compound.properties.mol_log_p.toFixed(4)}</td>
				<td className={heavy_atom_count}>{this.props.compound.compound.properties.heavy_atom_count}</td>
				<td className={heavy_atom_mol_wt}>{this.props.compound.compound.properties.heavy_atom_mol_wt.toFixed(4)}</td>
				<td className={nhoh_count}>{this.props.compound.compound.properties.nhoh_count}</td>
				<td className={no_count}>{this.props.compound.compound.properties.no_count}</td>
				<td className={num_h_acceptors}>{this.props.compound.compound.properties.num_h_acceptors}</td>
				<td className={num_h_donors}>{this.props.compound.compound.properties.num_h_donors}</td>
				<td className={num_het_atoms}>{this.props.compound.compound.properties.num_het_atoms}</td>
				<td className={num_rot_bonds}>{this.props.compound.compound.properties.num_rot_bonds}</td>
				<td className={num_val_electrons}>{this.props.compound.compound.properties.num_val_electrons}</td>
				<td className={ring_count}>{this.props.compound.compound.properties.ring_count}</td>
			
			</tr>
		);
	}
}

export default TableRow;
