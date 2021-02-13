import React from 'react';
import TableRow from './table_row.js';

class DataTable extends React.Component {
	
	render() {
		const rows = this.props.compounds.map((compound, index) => 
		<TableRow
			key = {index +1}
			counter = {index + 1}
			compound={compound}
			show_well = {this.props.show_well}
			show_code = {this.props.show_code}
			show_smiles = {this.props.show_smiles}
			show_structure = {this.props.show_structure}
			show_concentration = {this.props.show_concentration} 
			show_mw = {this.props.show_mw}
			show_tpsa = {this.props.show_tpsa}
			show_logp = {this.props.show_logp}
			show_heavy_atom_count = {this.props.show_heavy_atom_count}
			show_heavy_atom_mol_wt = {this.props.show_heavy_atom_mol_wt}
			show_nhoh_count = {this.props.show_nhoh_count}
			show_no_count = {this.props.show_no_count}
			show_num_h_acceptors = {this.props.show_num_h_acceptors}
			show_num_h_donors = {this.props.show_num_h_donors}
			show_num_het_atoms = {this.props.show_num_het_atoms}
			show_num_rot_bonds = {this.props.show_num_rot_bonds}
			show_num_val_electrons = {this.props.show_num_val_electrons}
			show_ring_count = {this.props.show_ring_count}
		/>
		);
		return (
			<tbody id="datatable-body">
				{rows}
			</tbody>
		);
	}
}

export default DataTable;
