import React from 'react';
import {Link} from "react-router-dom";
import Details from './details.js';
import PresetOption from './preset_option.js';
import { ChevronDown } from '../icons.js';

class LibraryOption extends React.Component {
	
	constructor(props){
		super(props);
		this.state = {details : false};
		this.showDetails = this.showDetails.bind(this);
		this.hideDetails = this.hideDetails.bind(this);
	}
	
	showDetails() {
		this.setState({details: true });
	}
	
	hideDetails() {
		this.setState({details: false });
	}
	
	
	render(){
		const details = this.state.details ? <Details collection={this.props.lib} type="library" hideDetails={this.hideDetails}/> : null;
		const lib = this.props.lib;	
		const presets = lib.presets.map((p, index) => {
				return <li key={index}>
							<PresetOption 
							preset={p}
							handleCheckboxChange = {this.props.handleChangePreset}
							defaultChecked={this.props.selected(p)}
							/>
						</li>
			})
		
		return (
			<div className="library">
				<input type="checkbox" 
					id={"l_" + this.props.lib.id}
					value={this.props.lib.id} 
					name="lib_ids" defaultChecked={this.props.defaultChecked} 
					onChange={event => this.props.handleCheckboxChange(event)} 
				/>
				<label htmlFor={"l_" + this.props.lib.id}>{lib.name} ({lib.size}) <ChevronDown size="16" handleClick={this.showDetails}/></label>
				{details}
				<ul>
				{presets}
				</ul>
			</div>
		)
	}
}

export default LibraryOption;
