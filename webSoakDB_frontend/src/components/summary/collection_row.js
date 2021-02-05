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
				<td><a href={"/playground/" + col.library.name + "/" + col.name}>See compounds</a></td>
				<td>
					<button onClick={event => this.props.handleClick(col)}>Remove</button>
				</td>
			</tr>
		); 
	}

}

export default CollectionRow;
