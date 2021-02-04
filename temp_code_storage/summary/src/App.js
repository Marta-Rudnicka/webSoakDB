import React from 'react';
import './home.css';
import CollectionRow from './collection_row.js';
import TableHeader from './th.js';

function deepCopyObjectArray(array){
	//this sucks; needs to be changed
	let output = [];
	array.forEach(object => {
		const objectDeepCopy = JSON.parse(JSON.stringify(object));
		output.push(objectDeepCopy);
	});
	return output;
}

//function App() {
class App extends React.Component {
	
	constructor(props) {
		super(props);
		this.removeCollection = this.removeCollection.bind(this);
		this.state = { 	libraries: this.props.proposal.libraries,
						ownplates: this.props.ownplates,
						proposalPlates: this.props.proposalPlates};
	}
	
	get_compound_collections(){
		let collections = []
			
		this.state.proposalPlates.forEach(plate =>{
			if (plate.library.public){
				plate.origin = 'XChem in-house library';
			}
			else{
				plate.origin = "User-submitted library"
			}
			collections.push(plate);
		});
		//TODO: get subsets from selected presets and cherrypicking lists
		return collections;
	}
	
	removeCollection(plate){
		const plates = deepCopyObjectArray(this.state.proposalPlates)
		const found = plates.find(object => object.id === plate.id);
		plates.splice(plates.indexOf(found), 1);
		this.setState({proposalPlates : plates});
	}
	
	
	undoChages(){
		console.log('fired undoChanges');
		console.log('libraries: ', this.state.libraries)
		console.log('libraries in props: ', this.props.proposal.libraries)
		
		this.setState({proposalPlates: this.props.proposalPlates});
	}
	
	render() {
		
		const collections = this.get_compound_collections()
		
		const library_rows = collections.map(collection =>{
			return <CollectionRow 
				key={collection.name + '-' + collection.library} 
				collection={collection} 
				handleClick={this.removeCollection}
			/>
		});
		
		const subset_rows = [];
		return (
			<div id="all">
				<nav>
			<a href="/playground/">Home</a> | <a href="/playground/proposal">Change proposal</a>
			</nav>
			<h1>Selected Compounds for {this.props.proposal.name} </h1>
			<section>
			<h2>Whole libraries</h2>
				<table>
					<TableHeader compoundsDescription="Available"/>
					<tbody>
						{library_rows}
					</tbody>
				</table>
			</section>
			<section>
				<h2>Selections from libraries (cherrypicked compounds)</h2>
				<table>
					<TableHeader compoundsDescription="Available"/>
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

export default App;
