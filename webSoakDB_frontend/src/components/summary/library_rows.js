import React from 'react';
import axios from 'axios';

class LibraryInTable extends React.Component {
	
	constructor(props) {
		super(props);
		this.state = { plates: []};
	}
	
	componentDidMount() {
		const apiUrl = 'api/current_plates_list/' + this.props.library.id + '/';
		
		axios.get(apiUrl)
			.then(res => {
			const plates = res.data;
			this.setState({ plates });
      });
	}
	
	render() {
		let first = true;
		const library = this.props.library;
		const plates = this.state.plates;
		const showPlate = this.props.showPlate
		const handleClick = this.props.handleClick
		
		function plateCells(plate){
			return(
			<React.Fragment>
				<td>{plate.name}</td>
				<td>{plate.size}</td>
				<td><span className="pseudo-link" onClick={event => showPlate(library.name, plate.name, true)}>See compounds</span></td>
			</React.Fragment>
			)
		}
		
		function firstRow(plate, index){
			return (<tr key={index}>
				<td rowSpan={plates.length}>{library.name}</td>
				<td rowSpan={plates.length}>{library.public ? "XChem in-house library" : "User-submitted library"}</td>
				{plateCells(plate)}
				<td rowSpan={plates.length}><button onClick={event => handleClick(library.id)}>Remove</button></td>
			</tr>)
		}
		
		function nextRow(plate, index){
			return <tr key={index}>{plateCells(plate)}</tr>
		}
		
			
		const plateRows = plates.map((plate, index) => {
			if (first) {
				first = false;
				return firstRow(plate, index);				
			}
			else {
				return nextRow(plate, index);
			}
		});
		
		function empty(){
			return <tr><td colSpan="6">No libraries added to the proposal</td></tr>;
		}
		
		const content = plates.length > 0 ? plateRows : empty();

		return(
		<React.Fragment>
			{content}
		</React.Fragment>
		); 
	}

}

export default LibraryInTable;
