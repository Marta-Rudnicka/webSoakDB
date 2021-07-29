import React from 'react';
import {Link } from "react-router-dom";


class SubsetTable extends React.Component {
	
	render() {
		const rows = this.props.selectedSubsets.map((subset, index)=>
			<tr key={index}>
				<td>{subset.name}</td>
				<td>{subset.origin}</td>
				<td>{subset.size}</td>
				<td><a href={"/compounds/subset/" + subset.id + "/0/"} target="_blank">See compounds</a></td>
				<td><button onClick={event => this.props.handleClick(subset.id)}>Remove</button></td>
			</tr>
		);
		
		if (this.props.selectedSubsets.length>0){			
			return(
				<React.Fragment>
					<h2>Selections from libraries (cherrypicked compounds)</h2>
					<table className="table table-bordered" id="table">
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
				</React.Fragment>
			)
		}
		else {
			return null;
		}
	}
}

export default SubsetTable;



