import React from 'react';

class TableHeader extends React.Component {
	render() {
		return(
		<thead>
			<tr>
				<th>Library</th>
				<th>Origin</th>
				<th>Currently used plate</th>
				<th>{this.props.compoundsDescription} <br/>compounds</th>
				<th>Compounds</th>
				<th></th>
			</tr>
		</thead>	
		); 
	}

}

export default TableHeader;
