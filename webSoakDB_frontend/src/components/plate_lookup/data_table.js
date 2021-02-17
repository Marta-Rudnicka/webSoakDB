import React from 'react';
import TableRow from './table_row.js';

class DataTable extends React.Component {
	
	render() {
		const rows = this.props.compounds.map((compound, index) => {
			//console.log('compound: ', compound,)
		return <TableRow
			key = {index +1}
			counter = {index + 1}
			compound={compound}
			display = {this.props.display}
		/>}
		);
		return (
			<tbody id="datatable-body">
				{rows}
			</tbody>
		);
	}
}

export default DataTable;
