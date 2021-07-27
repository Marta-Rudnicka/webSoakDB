import React, { Component, Fragment } from 'react';

import ReactDOM from 'react-dom';
import Home from './home/Home.js';
import Picker from './picker/Picker.js';
import Summary from './summary/Summary.js';
import ProposalSelection from './main/proposal_selection.js';
import UrlTranslator from './main/url_translator.js';

import {
  BrowserRouter as Router,
  Switch,
  Route,
  Link,
  Redirect,
} from "react-router-dom";

import axios from 'axios';

class PrivateRoute extends Component {

  render(){
    let content = null;
    if (this.props.proposal){
      content = <Route path={this.props.path}> {this.props.children} </Route>;
    }
    else {
      content = <Redirect to='/selection/'/>;
    }
    
    return (
    <React.Fragment>
    {content}
    </React.Fragment>
    );
  }
}

class App extends Component {
  
  constructor(props){
    super(props);
    this.logIn = this.logIn.bind(this);
    this.updateSelection = this.updateSelection.bind(this);
    this.trackUnsavedChanges = this.trackUnsavedChanges.bind(this);
    this.refreshAfterUpload = this.refreshAfterUpload.bind(this);
    this.state = {
      proposal: null,
      unsavedChanges: false,
      }
    };
  
  
  logIn(id){
    if(event){
      event.preventDefault();  
    }
    if (id){
      const proposalApiUrl = '/api/proposals/' + id + '/';
      axios.get(proposalApiUrl)
        .then(res => {
        const proposal = res.data;
        this.setState({ proposal })
      });
    }
    else{
      this.setState({ proposal : null });
    }
  }
  
  updateSelection(libraries, subsets){
    const propId = this.state.proposal.id;
    console.log(this.state.proposal);
    
    const apiUrl='/api/update_proposal_selection/' + propId + '/';
    axios.patch(apiUrl, {libraries: libraries, subsets: subsets})
    .then(res =>{
      event.preventDefault();
      if (res.status===200){
        console.log(this.state.proposal)
        this.logIn(propId);
      }
    })
    .catch(error => {
      console.log('updateSelection: axios error:');
      console.log(error)
      alert('Failed to save changes in your selection due to server error. Please try again.');
      this.trackUnsavedChanges(true);
    })
  }
  

  refreshAfterUpload(){
    this.logIn(this.state.proposal.proposal);  
  }
   
  trackUnsavedChanges(bool){
    this.setState({unsavedChanges: bool})
  }
  
  confirmLeaving(event){
    if(this.state.unsavedChanges){
      if (!window.confirm("You have unsaved changes in your selection. Do you want to discard them and leave this page?")){
        event.preventDefault();
      }
      else{
        this.trackUnsavedChanges(false);
      }
    }
  }

    render() {
        return (
         <div id="content"> 
      <Router>
        <nav className="navbar navbar-expand-sm navbar-light bg-light">
          <ul className="nav nav-tabs">
          
          <li className="nav-item">
            <Link onClick={e => this.confirmLeaving(e)} className="navbar-brand" to="/selection/home/">Home</Link>
          </li>
          <li className="nav-item">
            <Link onClick={e => this.confirmLeaving(e)} className="navbar-brand" to="/selection/selection/">Select compounds</Link>
          </li>
          <li className="nav-item">
            <Link onClick={e => this.confirmLeaving(e)} className="navbar-brand" to="/selection/summary/">Selection summary</Link>
          </li>
          <li className="nav-item">
            <Link onClick={e => this.confirmLeaving(e)} className="navbar-brand" to="/selection/">Change proposal</Link>
          </li>
          <li className="nav-item">
            <a onClick={e => this.confirmLeaving(e)} className="navbar-brand" href="/accounts/logout/">Log out</a>
          </li>
          </ul>
        </nav>
        
        <Switch> 
          <PrivateRoute path="/selection/home/" proposal={this.state.proposal} >
            <Home proposal={this.state.proposal}/>
          </PrivateRoute>
          <PrivateRoute path="/selection/selection/" proposal={this.state.proposal}>
            <Picker 
              proposal={this.state.proposal}
              trackUnsavedChanges = {this.trackUnsavedChanges}
              updateSelection = {this.updateSelection}
              refreshAfterUpload = {this.refreshAfterUpload}
              unsavedChanges = {this.state.unsavedChanges}
              />
            </PrivateRoute>
          <PrivateRoute path="/selection/summary/" proposal={this.state.proposal}>
            <Summary 
              proposal={this.state.proposal}
              trackUnsavedChanges = {this.trackUnsavedChanges}
              updateSelection = {this.updateSelection}
            />
          </PrivateRoute>
          <Route path="/compounds/:type/:id/" children={<UrlTranslator />} proposal={this.state.proposal}/>
          
          <Route path="/selection/">
          { this.props.proposal ? <Redirect to="/selection/home/" /> : <ProposalSelection logIn={this.logIn} proposal={this.state.proposal} /> }
          </Route>  
        </Switch>
      </Router>
        </div>
        )
    }
}

export default App;

ReactDOM.render(<App />, document.getElementById('app'));
