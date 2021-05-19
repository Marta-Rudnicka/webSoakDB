import React from 'react';
import properties_dict from './properties_dict.js'
import CSRFToken from './csrf.js';
import { XIcon, ChevronDown } from '../main/icons.js';

class SelectionGraph extends React.Component {
  
  componentDidMount(){
    this.submit_form();
  }

  submit_form(){
    const form_id = this.props.attr + '_form';
    const frame_id = this.props.attr + '_frame';
    if (document.getElementById(frame_id)){
      document.forms[form_id].requestSubmit();
    }
  }
  
  render(){
    const libs = this.props.parentState.selectedLibIds
    if (!libs){
      libs = 0;
    }
    const subs = this.props.parentState.selectedSubsetIds
    if (!subs){
      subs = 0;
    }
    const frame_id = this.props.attr + '_frame';
    const form_id = this.props.attr + '_form';
    const url = '/selection-histogram/' + this.props.attr + '/';

    return(
      <div>
        <iframe name={frame_id} id={frame_id} src={url}></iframe>
        
        <form method="post" action={url} target={frame_id} id={form_id}>
          <CSRFToken />
          <input type="hidden" name="attr" value={this.props.attr} />
          <input type="hidden" name="libs" value={libs} />
          <input type="hidden" name="subs" value={subs} />
        </form>
      </div>
    )
  }

}


class Explanation extends React.Component {
  render(){
    return(
      <div className="explanation">
        <div className="x-icon">
          <XIcon size="20" handleClick={this.props.hideExplanation}/>
        </div>
        <p>Includes:</p>
        <ul>
          <li><small>XChem libraries and presets selected (including unsaved)</small></li>
          <li><small>User-submitted libraries and cherry-picking lists (unless deleted)</small></li>
        </ul>
        <p>Excludes:</p>
        <ul>
          <li><small>Duplicated compounds (e.g. if there is an overlap between a library and a preset)</small></li>
          <li><small>Libraries uploaded without SMILES strings</small></li>
        </ul>
      </div>
    )
  }
}

class SummaryGraphs extends React.Component {
  
  constructor(props){
    super(props);
    this.state = {
		submit: false, 
		explanation : false
		};
    this.showExplanation = this.showExplanation.bind(this);
    this.hideExplanation = this.hideExplanation.bind(this);
  }
  
  showExplanation() {
    this.setState({explanation: true });
  }
  
  hideExplanation() {
    this.setState({explanation: false });
  }
  
  componentDidUpdate(prevProps, prevState){
    if(prevProps !== this.props){
      this.setState({submit: false});
    }
  }
  submit(){
    this.setState({submit: true})
  }
  
  render(){
  let cells = null;
  if(this.state.submit){
    cells = this.props.properties.map((p, index) => {
      return <td key={index}><SelectionGraph parentState={this.props.parentState} attr={p} caption={properties_dict[p]} submit={this.state.submit}/></td>
    });
  }
  else{
    cells = this.props.properties.map((p, index) => {
      return <td key={index}>No computed yet</td>
    });
  }
  const explanation = this.state.explanation ? <Explanation hideExplanation={this.hideExplanation}/> : null;
    return(
    <tr>
      <td>
        <p><strong>All selected <ChevronDown size="16" handleClick={this.showExplanation}/></strong></p>
        {explanation}
        <button onClick={() => this.submit()}>Compute</button>
      </td>
      {cells}
    </tr>
    );
  }
}

export default SummaryGraphs;
