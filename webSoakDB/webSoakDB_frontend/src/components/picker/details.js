import React from 'react';
import {Link} from "react-router-dom";
import { XIcon } from '../main/icons.js';
import { properties_list } from './properties_dict.js'

class Details extends React.Component {

  render(){
    const col = this.props.collection;
    const id = col.id;
    const type = this.props.type;
    let url_type = "plate"
    let url_id = null; 
    let description = null;
    if (type === "preset"){
      url_type = type;
      url_id = id;
      description = col.description;
    }
    if (type === "library"){
      url_id = col.current_plate;
    }
    const link_url = "../../compounds/" + url_type + '/' + url_id + '/';
    const graphs_url = "/media/html_graphs/" + type + '/' + id + '/';

    const graphs = properties_list.map((key, index) => {
		return <iframe key={index} src={graphs_url + key + '.html'}></iframe>
	});
    
    return (
      <div className="library-details">
        <div className="x-icon">
          <XIcon size="20" handleClick={this.props.hideDetails}/>
        </div>
        <h2>{col.name}</h2>
        <h3>({col.size} compounds {type === "library" ? "currently available" : null}) <br/>
          <Link className="plate-link" to={link_url} target="_blank">See the full list of compounds</Link>
        </h3>
        <p>{description}</p>
        <div>
          <h3>Compound properties distribution:</h3>
          <div className="graphs">
            {graphs}
          </div>
        </div>
      </div>
    )
  }
}

export default Details;

