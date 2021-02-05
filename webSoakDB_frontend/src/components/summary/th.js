import React from 'react';

//function App() {
class TableHeader extends React.Component {
	render() {
		return(
		<thead>
			<tr>
				<th>Library</th>
				<th>Currently used plate</th>
				<th>{this.props.compoundsDescription} <br/>compounds</th>
				<th>Origin</th>
				<th>Compounds</th>
				<th></th>
			</tr>
		</thead>	
		); 
	}

}

export default TableHeader;
