import React from 'react';
import {Link} from "react-router-dom";
import { ChevronDown } from '../main/icons.js';
import Details from './details.js';

class PresetOption extends React.Component {
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
    const id = this.props.preset.id
    const name = this.props.preset.name
    const size = this.props.preset.size
    const details = this.state.details ? <Details collection={this.props.preset} type="preset" hideDetails={this.hideDetails}/> : null;
    return (
      <React.Fragment>
      <div>
        <input 
          id={"p_" + id}
          type="checkbox" 
          value={id} name="preset_ids" 
          onChange={event => this.props.handleCheckboxChange(event)} 
          defaultChecked={this.props.defaultChecked}
        />
                                      
		  &nbsp;{name} ({size}) <ChevronDown size="14" handleClick={this.showDetails}/>
		</div>
		{details}
      </React.Fragment>
    )
  }
}

export default PresetOption;
