import React from 'react';

class ProposalTable extends React.Component {
  render() {
	  return (
		<table id="experiment-summary">
		  <tbody>
			  <tr>
			    <th>Proposal</th>
				  <th>Lab visit</th>
				  <th>Date:</th>
				  <th>Protein</th>
				  <th>Plate type:</th>
				</tr>
				<tr>
				  <td>{this.props.proposal}</td>
				  <td>{this.props.visit}</td>
				  <td>{this.props.date}</td>
				  <td>{this.props.protein}</td>
				  <td>autodetected</td>
				</tr>
			</tbody>
		</table>
		);	
	}
}

export default ProposalTable;
