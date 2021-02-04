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
			show_p3 = {this.props.show_p3}
			show_p4 = {this.props.show_p4}
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
