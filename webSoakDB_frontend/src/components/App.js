import React, { Component, Fragment } from 'react';
import ReactDOM from 'react-dom';

import Main from './main/Main.js';
import ProposalSelection from './main/proposal_selection.js';

import axios from 'axios';
//import { Provider } from 'react-redux';
//import store from '../store';

class App extends Component {
	
	constructor(props){
		super(props);
		this.logIn = this.logIn.bind(this);
		this.state = {proposalName: "",
		};
	}
	
	logIn(name){
		this.setState({proposalName: name});
	}
	
    render() {
		let app;
		if(this.state.proposalName){
			app = <Main proposalName={this.state.proposalName} logIn={this.logIn}/>;
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
