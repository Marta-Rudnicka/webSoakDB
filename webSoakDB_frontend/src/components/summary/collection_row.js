import React from 'react';

//function App() {
class CollectionRow extends React.Component {
	render() {
		const col = this.props.collection;
		return(
			<tr>
				<td>{col.library.name}</td>
				<td>{col.name}</td>
				<td>{col.size}</td>
				<td>{col.origin}</td>
				<td><span className="pseudo-link" onClick={event => this.props.showPlate(col.library.name, col.name, true)}>See compounds</span></td>
				<td>
					<button onClick={event => this.props.handleClick(col)}>Remove</button>
				</td>
			</tr>
		); 
	}

}

export default CollectionRow;
