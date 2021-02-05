import React from 'react';
//import './home.css';
import CollectionRow from './collection_row.js';
import TableHeader from './th.js';


const proposal = {
        "name": "test_proposal_3",
        "libraries": [
            {
                "id": 2,
                "name": "DSI_Poised_EG",
                "for_industry": true,
                "public": true
            },
            {
                "id": 3,
                "name": "FragLite",
                "for_industry": false,
                "public": true
            },
            {
                "id": 6,
                "name": "mylib_for_test_proposal_3",
                "for_industry": true,
                "public": false
            }
        ],
        "subsets": []
    };

const proposal_plates = [
    {
        "id": 3,
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate1",
        "current": true,
        "size": 308
    },
    {
        "id": 7,
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate2",
        "current": true,
        "size": 308
    },
    {
        "id": 9,
        "library": {
            "id": 3,
            "name": "FragLite",
            "for_industry": false,
            "public": true
        },
        "name": "test plate",
        "current": true,
        "size": 31
    },
    {
        "id": 6,
        "library": {
            "id": 6,
            "name": "mylib_for_test_proposal_3",
            "for_industry": true,
            "public": false
        },
        "name": "mylib_for_test_proposal_3",
        "current": true,
        "size": 19
    }
]
const current = [
    {
        "id": 3,
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate1",
        "current": true,
        "size": 308
    },
    {
        "id": 7,
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate2",
        "current": true,
        "size": 308
    },
    {
        "id": 9,
        "library": {
            "id": 3,
            "name": "FragLite",
            "for_industry": false,
            "public": true
        },
        "name": "test plate",
        "current": true,
        "size": 31
    },
    {
        "id": 6,
        "library": {
            "id": 6,
            "name": "mylib_for_test_proposal_3",
            "for_industry": true,
            "public": false
        },
        "name": "mylib_for_test_proposal_3",
        "current": true,
        "size": 19
    }
]

const ownplates = [{
        "id": 6,
        "library": {
            "id": 6,
            "name": "mylib_for_test_proposal_3",
            "for_industry": true,
            "public": false
        },
        "name": "mylib_for_test_proposal_3",
        "size": 19
    }]


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
		this.state = { 	libraries: proposal.libraries,
						ownplates: ownplates,
						proposalPlates: proposal_plates};
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
		console.log('libraries in props: ', proposal.libraries)
		
		this.setState({proposalPlates: proposal_plates});
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
				<h1>Selected Compounds for {proposal.name} </h1>
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

export default Summary;
