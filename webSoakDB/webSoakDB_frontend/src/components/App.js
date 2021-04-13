import React, { Component, Fragment } from 'react';

import ReactDOM from 'react-dom';
import Home from './home/Home.js';
import Main from './main/Main.js';
import Picker from './picker/Picker.js';
import Summary from './summary/Summary.js';
import ProposalSelection from './main/proposal_selection.js';
import PlateLookup from './plate_lookup/PlateLookup.js';
import UrlTranslator from './url_translator.js';

import {
  BrowserRouter as Router,
  Switch,
  Route,
  Link,
  useParams,
  Redirect,
  useHistory,
  useLocation
} from "react-router-dom";

import { deepCopyObjectArray, getAttributeArray, mean, shareAllElements } from  '../actions/stat_functions.js';

import axios from 'axios';

class PrivateRoute extends Component {

	render(){
		let content = null;
		if (this.props.proposal){
		  content = <Route path={this.props.path}> {this.props.children} </Route>;
		}
		else {
			content = <Redirect to='/'/>;
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
		this.updateLibrarySelection = this.updateLibrarySelection.bind(this);
		this.updateSubsetSelection = this.updateSubsetSelection.bind(this);
		this.trackUnsavedChanges = this.trackUnsavedChanges.bind(this);
		this.refreshAfterUpload = this.refreshAfterUpload.bind(this);
		this.state = {
			proposal: null,
			}
		};
	
	
	logIn(name){
		if(event){
			event.preventDefault();	
		}
		if (name){
			const proposalApiUrl = '/api/proposals/' + name;
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
	
	updateLibrarySelection(idArray){
		console.log('fired updateLibrarySelection with', idArray);
		event.preventDefault()
		const apiUrl='/api/update_proposal_selection/' + this.state.proposal.proposal + '/';
		axios.patch(apiUrl, {libraries: idArray}) 
		.catch(error => {
			console.log('updateLibrarySelection: axios error:');
			console.log(error)
		})
		
		this.logIn(this.state.proposal.proposal);
	}
	
	updateSubsetSelection(idArray){
		event.preventDefault()
		const apiUrl='/api/update_proposal_selection/' + this.state.proposal.proposal + '/';
		axios.patch(apiUrl, {subsets: idArray}) 
		.catch(error => {
			console.log('updateLibrarySelection: axios error:');
			console.log(error)
		})
		
		setTimeout(() => {this.logIn(this.state.proposal.proposal)}, 1000);
	}
	
	refreshAfterUpload(){
		this.props.logIn(this.state.proposal.proposal);	
	}
	
	componentDidUpdate(prevProps, prevState) {
		if (prevState.proposal !== this.state.proposal){
			console.log('fired componentDidUpdate from App');
		}
	}
	
	
	//to give warning when user risks discarding changes
	trackUnsavedChanges(bool){
		console.log('fired trackUnsavedChanges')
		this.setState({unsavedChanges: bool})
	}

	


    render() {
		let app;
		if(this.state.proposal){
			app = <Main proposal={this.state.proposal} logIn={this.logIn} />;
		}
		else{
			app = <ProposalSelection logIn={this.logIn} />;
		}
		
        return (
         <div id="content"> 
			<Router>
				<nav className="navbar navbar-expand-sm navbar-light bg-light">
				  <ul className="nav nav-tabs">
					
					<li className="nav-item">
					  <Link className="navbar-brand" to="/home/">Home</Link>
					</li>
					<li className="nav-item">
					  <Link className="navbar-brand" to="/selection/">Select compounds</Link>
					</li>
					<li className="nav-item">
					  <Link className="navbar-brand" to="/summary/">Selection summary</Link>
					</li>
					<li className="nav-item">
					  <Link className="navbar-brand" to="/">Log in/out</Link>
					</li>
				  </ul>
				</nav>
				
				<Switch> 
				  <PrivateRoute path="/home/"proposal={this.state.proposal} >
					<Home proposal={this.state.proposal}/>
				  </PrivateRoute>
				  <PrivateRoute path="/selection/" proposal={this.state.proposal}>
					<Picker 
						proposal={this.state.proposal}
						trackUnsavedChanges = {this.trackUnsavedChanges}
						updateLibrarySelection = {this.updateLibrarySelection}
						updateSubsetSelection = {this.updateSubsetSelection}
						refreshAfterUpload = {this.refreshAfterUpload}
						/>
				  </PrivateRoute>
				  <PrivateRoute path="/summary/" proposal={this.state.proposal}>
					<Summary 
						proposal={this.state.proposal}
						trackUnsavedChanges = {this.trackUnsavedChanges}
						updateLibrarySelection = {this.updateLibrarySelection}
						updateSubsetSelection = {this.updateSubsetSelection}
					/>
				  </PrivateRoute>
				  <PrivateRoute path="/compounds/:type/:id/" children={<UrlTranslator />} proposal={this.state.proposal}/>
				  <Route path="/">
					{ this.props.proposal ? <Redirect to="/home/" /> : <ProposalSelection logIn={this.logIn} proposal={this.state.proposal} /> }
				  </Route>	
				</Switch>
			</Router>
        </div>
        )
    }
}

export default App;

ReactDOM.render(<App />, document.getElementById('app'));
