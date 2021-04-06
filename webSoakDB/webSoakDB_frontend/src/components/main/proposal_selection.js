import React from 'react';
import axios from 'axios';

import { Redirect } from "react-router-dom";

class ProposalSelection extends React.Component {
	
	constructor(props){
		super(props);
		this.state = {
			proposals: [],
			chosenName: "",
			}
		this.getSelection = this.getSelection.bind(this);
	}
	
	componentDidMount() {
		const apiUrl = '/api/proposals/';
		
		axios.get(apiUrl)
			.then(res => {
			const proposals = res.data;
			this.setState({ proposals });
      });
      
      this.props.logIn(null);
	}
	
	getSelection(e){
		this.setState({chosenName: e.target.value});
	}
	
	
	render(){
		const options = this.state.proposals.map((proposal, index) => 
			 <option key={index} value={proposal.name} htmlFor="proposal_name">{proposal.name}</option>
		);
		
		let loggedIn = null;
		loggedIn = this.props.proposal ?  <Redirect to="/home/" /> : null;
		
		return (
		<div>
			<h1>Temporary logging page</h1>
			<div>This is a temporary page made for testing purposes. Normally, a user will have to log in first using their fed id, and then select the proposal.</div>
			<form>
				<legend>Choose proposal before browsing other pages</legend>
				<p>
					<label htmlFor="proposal_name">Select proposal:</label> 
					<select name="proposal_name" onChange={event => {this.getSelection(event)}}>
						<option value="">...</option>
					{options}
					</select>
				</p>
				<button onClick={event => this.props.logIn(this.state.chosenName)} >Confirm selection</button>
			</form>
			{loggedIn}
		</div>
		);
	}
}

export default ProposalSelection;
