import React from 'react';
import CrystalTile from './crystal_tile.js';
import cryst from './Crystal gallery_files/crystal.png';
import hide from './Crystal gallery_files/hide.png';
import show from './Crystal gallery_files/show.png';

//const cryst = '';

class CrystalGroup extends React.Component {
	
	constructor(props) {
    super(props);
    this.handleShowCrystals = this.handleShowCrystals.bind(this);
    this.handleHideCrystals = this.handleHideCrystals.bind(this);
    this.handleClick = this.handleClick.bind(this);
    this.state = {dummy: 'value',
				displayGroup: true,}
  }
	
	handleClick(well) {
		this.props.iconAction(well);
	}
	
	handleShowCrystals() {
		this.setState({displayGroup: true});
	}
	
	handleHideCrystals(){
		this.setState({displayGroup: false});
	}
	
	render() {
		const crystalTiles = this.props.crystals.map((crystal) =>
			<CrystalTile
			crystal = {crystal}
			key={crystal.well}
			crystalpic={cryst} 
			crystalClass={this.props.crystalClass}
			iconAction={this.props.iconAction}
			plateName={this.props.plateName}
			/>
		);
		
		let crystalGroupClass;
		let hideIconClass;
		let showIconClass;
		
		if (this.state.displayGroup) {
			crystalGroupClass = 'crystal-group';
			hideIconClass = "hide-used-crystals";
			showIconClass = 'hide';
		}
		else {
			crystalGroupClass = 'hide';
			hideIconClass = "hide";
			showIconClass = 'hide-used-crystals';
			}
			
		return(
		<div className={this.props.divClass}>
			<h3>{this.props.heading}</h3>
			<img className={showIconClass} src={show} alt="show crystals" onClick={this.handleShowCrystals} />
			<img className={hideIconClass} src={hide} alt="hide crystals" onClick={this.handleHideCrystals} />
			<div className={crystalGroupClass}>
				{crystalTiles}
			</div>
		</div>
		);
	}
}

export default CrystalGroup;
