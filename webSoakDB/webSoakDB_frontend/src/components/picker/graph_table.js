import React from 'react';
import properties_dict from './properties_dict.js'
import Graph from './graph.js';
import { removeFromArray, getSubsetIds } from  '../../actions/stat_functions.js';
import axios from 'axios';
import SummaryGraphs from './summary_graphs.js';
import { LazyLoadComponent } from 'react-lazy-load-image-component';
import { shareAllElements } from '../../actions/stat_functions.js';

class GraphTr extends React.Component {
  constructor(props){
    super(props);
    this.state = {
      collection: null,
      public: this.props.public,
    }
  }
  
  componentDidMount(){
	if (this.props.col){
		this.setState({collection: this.props.col});
	}
	else{
	  let apiUrl = ""
	  if(this.props.type==="library"){	
		  apiUrl = '/api/libraries/' + this.props.col_id + '/';
	  }
	  if(this.props.type==="subset"){	
		  apiUrl = '/api/subset_detail/' + this.props.col_id + '/';
	  }
    axios.get(apiUrl)
      .then(res => {
      const collection = res.data;
      this.setState({ collection });
      if (collection.public){
		this.setState({public: collection.public});
		}
        });
	}
  }

  render(){
    let cells = null;
    let name = null;
	
	if (this.state.collection){
      cells = this.props.properties.map((p, index) => {
		return <td key={index}><Graph id={this.state.collection.id} type={this.props.type} property={p} public={this.state.public} /></td>
	  });
	  name = this.state.collection.name;
	  if(this.props.type==="subset"){
		  name = name + "(selected from " + this.state.collection.library.name + ")"
	  }
  }

   return(
      
		<tr>
		  <td><strong>{name}</strong></td>
		  {cells}
		</tr>
	  
    );
  }
}

class GraphTable extends React.Component {
  
  constructor(props){
    super(props);
    this.state = {
      show: ["mol_wt", "log_p", "num_h_donors", "num_h_acceptors", "num_rot_bonds"],
      presets: this.sortSubsets()[0],
      other_subsets: this.sortSubsets()[1]
    }
  }
  
  
  componentDidUpdate(prevProps, prevState){
    if(
      (!shareAllElements(prevProps.parentState.selectedLibIds, this.props.parentState.selectedLibIds) 
      || !shareAllElements(prevProps.parentState.selectedSubsetIds,this.props.parentState.selectedSubsetIds))
      && !this.props.parentState.waitingForSave
      ){
        console.log('fired componentDidUpdate')
        this.setState({presets : this.sortSubsets()[0], other_subsets: this.sortSubsets()[1]});
    }
  }

  shouldComponentUpdate(nextProps, nextState){
    //prevent continuous updates while a user file is being processed
    if (nextProps.parentState.waitingForSave){
      return false;
    }
    else {
      return true;
    }
  }
  
  manageProperties(event){
    const p = event.target.value;
    
    if(event.target.checked === true){
      const newShow = this.state.show.slice(0, this.state.show.length);
      newShow.push(p);
      this.setState({show : newShow});
    }
    else{
      const newShow = removeFromArray(this.state.show, [p]);
      this.setState({show : newShow});
    }
  }
  
  sortSubsets(){
    const presets = [];
    let other = this.props.parentState.selectedSubsetIds;
    this.props.parentState.presets.forEach(p =>{
      const subs = getSubsetIds(this.props.parentState, p.id)
      if (this.props.parentState.selectedSubsetIds.includes(subs[0])){
        presets.push(p);
        other = removeFromArray(other, subs);
      }
    })
    return [presets, other];
  }
  
  render(){
    console.log('rendering GraphTable')
    const checkboxes = Object.keys(properties_dict).map((key, index) => {
        let checked = false;
        if (this.state.show.includes(key)){
          checked = true;
        }
        return (
        <div key={index}>
          <input type="checkbox" value={key} id={key} onChange={e => this.manageProperties(e)} checked={checked}/>
          <label htmlFor={key}>&nbsp;{properties_dict[key]} &nbsp;&nbsp;</label>
        </div>
        )
      });
    
    const headers = this.state.show.map((prop, index) =>{
      return <th key={index}>{properties_dict[prop]}</th>
    });
    
    const l_rows = this.props.parentState.selectedLibIds.map((lib, index) => {
      return <GraphTr key={index} type="library" col_id={lib} public={false} properties={this.state.show}/>
    
    });
    const p_rows = this.state.presets.map((p, index) => {
      return <GraphTr key={index} type="preset" col={p} public={true} properties={this.state.show}/>
    
    });
    const s_rows = this.state.other_subsets.map((sub, index) => {
      return <GraphTr key={index} type="subset" col_id={sub} public={false} properties={this.state.show}/>
    
    });
    
    const last_row = (
    <SummaryGraphs 
      parentState={this.props.parentState} 
      properties={this.state.show}
      updateWaitingStatus = {this.props.updateWaitingStatus}
      
      />);
    
    return(
    <div>
      <h2>Properties to compute: </h2>
      <div id="properties">
        {checkboxes}
      </div>
      <LazyLoadComponent>
        <table className="table">
		  <caption>Molecular parameters distribution</caption>
          <thead>
            <tr>
              <th>Collection</th>
              {headers}
            </tr>
          </thead>
          <tbody>
            {l_rows}
            {p_rows}
            {s_rows}
            {last_row}
          </tbody>
        </table>
      </LazyLoadComponent>
    </div>
    )
  }
}

export default GraphTable;
