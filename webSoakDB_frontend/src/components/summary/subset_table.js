import React from 'react';
import axios from 'axios';

class SubsetTable extends React.Component {
	
	render() {
		const rows = this.props.selectedSubsets.map((subset, index)=>
			<tr key={index}>
				<td>{subset.name}</td>
				<td>{subset.origin}</td>
				<td>{subset.size}</td>
				<td><span className="pseudo-link" onClick={event => this.props.showPlate(subset.library, subset, false, false)}>See compounds</span></td>
				<td><button onClick={event => this.props.handleClick(subset.id)}>Remove</button></td>
			</tr>
		);
		
		if (this.props.selectedSubsets.length>0){			
			return(
				<section>
					<h2>Selections from libraries (cherrypicked compounds)</h2>
					<table className="summary-table">
						<thead>
							<tr>
								<th>Library</th>
								<th>Origin</th>
								<th>Selected <br/>compounds</th>
								<th>Compounds</th>
								<th>Remove</th>
							</tr>
						</thead>
						<tbody>
							{rows}
						</tbody>
					</table>
				</section>
			)
		}
		else {
			return null;
		}
	}
}

export default SubsetTable;



