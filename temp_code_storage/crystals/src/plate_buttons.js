import React from 'react';

class PlateButtons extends React.Component {
	
	constructor(props) {
		super(props);
		this.handleShowAllDetails = this.handleShowAllDetails.bind(this);
		this.handleHideAllDetails = this.handleHideAllDetails.bind(this);
	}
	
	makeDisplayButtons(display){
		let displayButton;
		if (display === 'show-plate'){
			displayButton = <button className="hide-plate-button" onClick={() => this.hidePlate()}>Hide plate</button>;
		}
		else {
			displayButton = <button className="show-plate-button" onClick={() => this.showPlate()} >Show plate</button>;
		}
		return displayButton;
	}
	
	/* Showing and hiding of the details is implemented using CSS only;
	 * display attribute of the div containing the details is decided by 
	 * radio buttons. Therefore the two methods below only check radio buttons.*/
	handleShowAllDetails(){
		const section = document.getElementById(this.props.plateName);
		section.querySelectorAll('.show-icon img').forEach(icon => {
			if (icon.checked !== true) {
				icon.click();
			}
		});
	}
	
	handleHideAllDetails(){
		const section = document.getElementById(this.props.plateName);
		section.querySelectorAll('.hide-icon img').forEach(icon => {
			if (icon.checked !== true) {
				icon.click();
			}
		});
	}
	
	makeShowDetailsButton(display){
		if (display === 'show-plate'){
			return <button className="show-all" onClick={this.handleShowAllDetails}>Display all crystal details</button>;
		}
	}
	makeHideDetailsButton(display){
		if (display === 'show-plate'){
			return <button className="hide-all" onClick={this.handleHideAllDetails}>Hide all crystal details</button>;
		}
	}
	makeRemoveButton(status){
		if (status === 'unused'){
			return <button>Remove plate from experiment </button>;
		}
	}
	
	showPlate(){
		this.props.show();
	}
	
	hidePlate(){
		this.props.hide();
	}
	
	render() {
		const display = this.props.displayStatus;
		return(
		  <div className="plate-buttons">
			{this.makeDisplayButtons(display)}
			{this.makeShowDetailsButton(display)}
			{this.makeHideDetailsButton(display)}
			{this.makeRemoveButton(this.props.status)}
		  </div>
        );
	}
}

export default PlateButtons;
