import React from 'react';

class StatHeaders extends React.Component {
	render(){
		return(
			<React.Fragment>
				<th>Mean molecular<br/>weight</th>
				<th>Mean TPSA</th>
				<th>Mean LogP</th>
				<th>Mean Heavy Atom Count</th>
				<th>Mean Heavy Atom Mol. Weight</th>
				<th>Mean NH or OH Count</th>
				<th>Mean Nitrogen and Oxygen Count</th>
				<th>Mean Hydrogen Bond Acceptor Count</th>
				<th>Mean Hydrogen Bond Donor Count</th>
				<th>Mean Heteroatoms Count</th>
				<th>Mean Rotatable Bonds Count</th>
				<th>Mean Valence Electrons Count</th>
				<th>Mean Ring Count</th>						
			</React.Fragment>
		
		)}
}

export default StatHeaders;
