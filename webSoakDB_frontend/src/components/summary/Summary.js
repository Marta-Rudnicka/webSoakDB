import React from 'react';
//import './home.css';
import CollectionRow from './collection_row.js';
import TableHeader from './th.js';
import axios from 'axios';


function deepCopyObjectArray(array){
	//this sucks; needs to be changed
	let output = [];
	array.forEach(object => {
		const objectDeepCopy = JSON.parse(JSON.stringify(object));
		output.push(objectDeepCopy);
	});
	return output;
}

class Summary extends React.Component {
	
	constructor(props) {
		super(props);
		this.removeCollection = this.removeCollection.bind(this);
		this.state = { 	proposalPlates: [],
						//proposal: {},
					};
	}
	
	componentDidMount() {
		const platesApiUrl = 'api/proposal_plates/' + this.props.proposal.name;;
		//const proposalApiUrl = 'api/proposals/' + this.props.proposal.name;
		
		axios.get(platesApiUrl)
			.then(res => {
			const proposalPlates = res.data;
			this.setState({ proposalPlates });
      });
      /*
		axios.get(proposalApiUrl)
			.then(res => {
			const proposal = res.data;
			this.setState({ proposal });
      });
      	*/	
	}
	
	get_plates(){
		let plates = []
			
		this.state.proposalPlates.forEach(plate =>{
			if (plate.library.public){
				plate.origin = 'XChem in-house library';
			}
			else{
				plate.origin = "User-submitted library"
			}
			plates.push(plate);
		});
		return plates;
	}
	
	get_subsets(){
		let subsets = []
			
		this.props.proposal.subsets.forEach(subset =>{
			subset.size = subset.compounds.length;
			subsets.push(subset);
		});
		return subsets;
	}
	
	removeCollection(plate){
		const plates = deepCopyObjectArray(this.state.proposalPlates)
		const found = plates.find(object => object.id === plate.id);
		plates.splice(plates.indexOf(found), 1);
		this.setState({proposalPlates : plates});
	}
	
	
	undoChages(){
		this.componentDidMount();
	}
	
	render() {
		
		const plates = this.get_plates();		
		
		const library_rows = plates.map(plate =>{
			return <CollectionRow 
				key={plate.name + '-' + plate.library} 
				collection={plate} 
				handleClick={this.removeCollection}
				showPlate={this.props.showPlate}
			/>
		});
		
		let subset_rows;
		
		if (this.props.proposal.subsets){
			const subsets = this.get_subsets();
			subset_rows = subsets.map(subset =>{
				return <CollectionRow 
					key={subset.name + '-' + subset.library} 
					collection={subset} 
					handleClick={this.removeCollection}
					showPlate={this.props.showPlate}
				/>
			});
		}
		else {
			subset_rows = [];
		}
		
		return (
			<div id="all">
				<h1>Selected Compounds for {this.props.proposal.name} </h1>
				<section>
				<h2>Whole libraries</h2>
					<table className="summary-table">
						<TableHeader compoundsDescription="Available"/>
						<tbody>
							{library_rows}
						</tbody>
					</table>
				</section>
				<section>
					<h2>Selections from libraries (cherrypicked compounds)</h2>
					<table className="summary-table">
						<TableHeader compoundsDescription="Selected"/>
						<tbody>
							{subset_rows}
						</tbody>
					</table>
				</section>
				<div>
				
					<button>Save changes </button>
					<button onClick={event => this.undoChages()}>Undo all changes </button>
				</div>
			</div>
		); 
	}
}

export default Summary;
