import React, { Component, Fragment } from 'react';
import ReactDOM from 'react-dom';

import Main from './main/Main.js';
import ProposalSelection from './main/proposal_selection.js';
import { deepCopyObjectArray, getAttributeArray, mean, shareAllElements } from  '../actions/stat_functions.js';

import axios from 'axios';

class App extends Component {
	
	constructor(props){
		super(props);
		this.logIn = this.logIn.bind(this);
		this.state = {
			proposal: null,
			}
		};
	
	
	logIn(name){
		console.log('Logging in');
		if(event){
			event.preventDefault();	
		}
		if (name){
			const proposalApiUrl = 'api/proposals/' + name;
			axios.get(proposalApiUrl)
				.then(res => {
				const proposal = res.data;
				this.setState({ proposal });
			});
		}
		else{
			this.setState({ proposal : null });
		}
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
         <div>  
			{app}
        </div>
        )
    }
}

export default App;

ReactDOM.render(<App />, document.getElementById('app'));
