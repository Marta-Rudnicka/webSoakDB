import React from 'react';
import axios from 'axios'; 

class PlateList extends React.Component {
	
	constructor(props){
		super(props)
		
		this.state={
			plates: [],
		}
	}
	
	componentDidMount() {
		const apiUrl = 'api/plates_list/' + this.props.library.name;		
		axios.get(apiUrl)
			.then(res => {
			const plates = res.data;
			this.setState({ plates });
      });     		
	}
	
	render() {
		
		if (this.props.library.public && this.state.plates.length > 1){
			
			const list = this.state.plates.map((plate, index)=> {
				let current;
				
				if (plate.current){
					current = '(current)';
				}
				return <li key={index}><span className="pseudo-link" onClick={event => this.props.showPlate(plate.library, plate, true, false)}>{plate.name}</span>{current}</li>;
			});
			
			return (
				<div id="other-plates">
					<h3>Other plates for {this.props.library.name}</h3>
					{list}
					<p>(If you wish to use an old plate in your experiment, contact XChem staff)</p>
				</div>
			);
		}
		else{
			return null;
		}
	}
}

export default PlateList;
