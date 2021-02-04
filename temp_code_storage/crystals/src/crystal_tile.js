/*Component representing one crystal; contains photo of the crystal and 
 * crystal data in <div class="infobox">. The infobox is shown and hidden
 * by clicking on icons using CSS only. */

import React from 'react';
import binpic from './Crystal gallery_files/bin.png';
import recyclepic from './Crystal gallery_files/recycle.jpg';

import info from './Crystal gallery_files/show-info.png';
import hide from './Crystal gallery_files/hide.png';

class CrystalTile extends React.Component {
	
  render() {
	const well = this.props.crystal.well;
	const plate = this.props.plateName;
	let pic;
	if (this.props.crystal.status === 'accepted') {
		pic = binpic;
	}
	else {
		pic = recyclepic;
	}
	
	return (
	  <div className='gallery' id={well}>
		  <div>
			  <p className="well-name">{well}</p>	
				<img className="bin-pic" src={pic} alt="remove crystal from experiment" onClick={() => this.props.iconAction(plate, well)}/>
				
				<input id={plate + well + '-checkbox'} className="show" type="radio" name={plate + well + '-radiogroup'} hidden />
				<label className="show-icon" htmlFor={plate + well + '-checkbox'}>
					<img className="show-pic" src={info} alt="show details" />
				</label>
				<label className="hide-icon" htmlFor={plate + well + '-hide'}>	
					<img className="hide-pic" src={hide} alt="alttext" />
				</label>

				<div className="infobox">
				  <span className="right"><strong>X:</strong></span><span>{this.props.crystal.x}</span>
				  <span className="right"><strong>Y:</strong></span><span>{this.props.crystal.y}</span>
				  <span className="right"><strong>Score:</strong></span> <span>{this.props.crystal.score}</span>
				</div>
				<input id={plate + well + '-hide'} className="show" type="radio" name={plate + well + '-radiogroup'} hidden />
				<img className="main-pic" src={this.props.crystalpic} alt="crystal" /> 
			</div>
		</div>	
	);
  }

}

export default CrystalTile;
