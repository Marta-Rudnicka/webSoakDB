import React from 'react';
import PlateButtons from './plate_buttons.js';
//import hide from './Crystal gallery_files/hide.png';
import CrystalGroup from './crystal_group.js';

function getDisplayStatus(crystals){
	const accepted = crystals.filter(crystal => crystal.status === 'accepted');
	
	if (accepted.length > 0){
		return 'show-plate';
	}
	else {
		return 'hide';
	}
}

function setClassOrHide(array, className){
	if (array.length === 0){
		return 'hide';
	}
	else {
		return className;
	}
}


class Plate extends React.Component {
	
	constructor(props) {
		super(props);
		
		this.state = {
			displayStatus: getDisplayStatus(this.props.crystals),
		};
		
		this.hidePlate = this.hidePlate.bind(this);
		this.showPlate = this.showPlate.bind(this);
	}
	
	hidePlate(){
		this.setState({displayStatus: 'hide'});
	}
	
	showPlate(){
		this.setState({displayStatus: 'show-plate'});
	}
	
	
	render() {
		const usedList = this.props.crystals.filter(crystal => crystal.status === 'used');
		const acceptedList = this.props.crystals.filter(crystal => crystal.status === 'accepted');
		const rejectedList = this.props.crystals.filter(crystal => crystal.status === 'rejected');
		const plateName = this.props.name;
		
		let sectionClass;
		if (acceptedList.length===0) {
			sectionClass = 'all-used';
		}
		let buttonStatus;
		
		if (usedList.length===0) {
			buttonStatus = 'unused';
			}
		
		let usedDivClass = setClassOrHide(usedList, 'used-div');
		let acceptedDivClass = setClassOrHide(acceptedList, 'accepted-div');
		let rejectedDivClass = setClassOrHide(rejectedList, 'rejected-div');
		
				
		return (
		<section className={sectionClass} id={plateName}>
			<h2>Crystallisation plate: {plateName}</h2>
			<PlateButtons displayStatus={this.state.displayStatus} plateName={plateName} show={this.showPlate} hide={this.hidePlate} status={buttonStatus}/>
			<div>Drop volume: {this.props.dropvol}<br />
				Used: {usedList.length}<br />
				Accepted: {acceptedList.length}<br />
				Rejected: {rejectedList.length}
			</div>
			<div className={this.state.displayStatus}>
				<CrystalGroup divClass={usedDivClass} heading="Used crystals:" crystalClass="gallery" crystals={usedList} iconAction={this.props.handleReject} plateName={plateName} />
				<CrystalGroup divClass={acceptedDivClass} heading="Accepted crystals:" crystalClass="gallery" crystals={acceptedList} iconAction={this.props.handleReject} plateName={plateName}  />
				<CrystalGroup divClass={rejectedDivClass} heading="Rejected crystals:" crystalClass="gallery" crystals={rejectedList} iconAction={this.props.handleAccept} plateName={plateName} />
			</div>
		</section>
		);
	}
}

export default Plate;
