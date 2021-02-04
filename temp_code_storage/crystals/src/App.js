import React from 'react';
import './Crystal gallery_files/crystals.css';

import ProposalTable from './proposal_table.js';
import Sidebar from './crystal_sidebar.js';
import Navigation from './navigation.js';
import Main from './crystal_main.js';

import cryst from './Crystal gallery_files/crystal.png';


class Crystal {
	constructor(well, x, y, score, picture, status){
		this.well = well;
		this.x = x;
		this.y = y;
		this.score = score;
		this.picture = picture;
		this.status = status;
	}
}

class PlateData {
	constructor(proposal, name, crystalArray, dropVolume){
			this.proposal = proposal;
			this.name = name;
			this.crystalArray = crystalArray;
			this.dropVolume = dropVolume;
	}
}

/* TEST DATASET */

const crystal1 = new Crystal("A42", 67, 78, 5, cryst, "accepted");
const crystal2 = new Crystal("B69", 15, -68, 2, cryst, "used");
const crystal3 = new Crystal("B54", 18, -68, 7, cryst, "accepted");

const crystal4= new Crystal("A42", 67, 78, 5, cryst, "used");
const crystal5 = new Crystal("B69", 15, -68, 2, cryst, "rejected");
const crystal6 = new Crystal("B54", 18, -68, 7, cryst, "used");

const crystal7 = new Crystal("A42", 67, 78, 5, cryst, "accepted");
const crystal8 = new Crystal("B69", 15, -68, 2, cryst, "accepted");
const crystal9 = new Crystal("B54", 18, -68, 7, cryst, "accepted");

const crystalArrayPart = new PlateData("123", 'Partially used', [crystal1, crystal2, crystal3], 45);
const crystalArrayUsed = new PlateData("123", 'All used', [crystal4, crystal5, crystal6], 30);
const crystalArrayFresh = new PlateData("123", 'Fresh', [crystal7, crystal8, crystal9], 30);

const testData = [crystalArrayUsed, crystalArrayPart, crystalArrayFresh];
//Will normally be fetched from an API
/* --------------------------------- */

function deepCopyPlatesArray(array){
	//this sucks; needs to be changed
	let output = [];
	array.forEach(plate => {
		const plateDeepCopy = JSON.parse(JSON.stringify(plate));
		output.push(plateDeepCopy);
	});
	return output;
}

function changeStatusBelowScore(plateArray, score, status){
	console.log('fired changeStatusBelowScore');
	plateArray.forEach(plate =>{
		plate.crystalArray.forEach(crystal => {
			if (crystal.score < score && crystal.status !== 'used'){
				crystal.status = status;
			}
		});
	});
}

function changeStatusAboveScore(plateArray, score, status){
	plateArray.forEach(plate =>{
		plate.crystalArray.forEach(crystal => {
			if (crystal.score > score && crystal.status !== 'used'){
				crystal.status = status;
			}
		});
	});
}

//function App() {
class App extends React.Component {
	
	constructor(props) {
		super(props);
		this.state = {plates: testData};
		this.changeCrystalStatus = this.changeCrystalStatus.bind(this);
		this.handleReject = this.handleReject.bind(this);
		this.handleAccept = this.handleAccept.bind(this);
		this.filterByScore = this.filterByScore.bind(this);
	}
	
	changeCrystalStatus(plateName, well, newStatus){
		//TODO - is there a better way that overwriting all of the data?
		//might cause problems when hundreds of crystals are involved
		let platesCopy = deepCopyPlatesArray(this.state.plates);
		const modifiedPlate = platesCopy.find(plate => plate.name === plateName);
		const modifiedCrystal = modifiedPlate.crystalArray.find(crystal => crystal.well === well);
		modifiedCrystal.status = newStatus;
		this.setState({plates: platesCopy});
		}
		
	handleReject(plateName, well){
		this.changeCrystalStatus(plateName, well, 'rejected');
	}
	
	handleAccept(plateName, well){
		this.changeCrystalStatus(plateName, well, 'accepted');
	}
	
	filterByScore(rejectLimit, acceptLimit){
		let platesCopy = deepCopyPlatesArray(this.state.plates);
		if(rejectLimit){
			changeStatusBelowScore(platesCopy, rejectLimit, "rejected");
		}
		if(acceptLimit){
			changeStatusAboveScore(platesCopy, acceptLimit, "accepted");
		}
		this.setState({plates: platesCopy});
	}
		
	render() {
		return (
		<div id="all">
			<ProposalTable proposal="proposal number" visit="visit" date="date" protein="protein" />
			<Navigation />
			<h1>Crystal gallery</h1>
			<Main plates={this.state.plates} handleAccept={this.handleAccept} handleReject={this.handleReject} test={this.state.test} />
			<Sidebar plates={this.state.plates} filterByScore={this.filterByScore}/>
		</div>
	); 
}

}

export default App;
