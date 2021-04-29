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
			 <option key={index} value={proposal.proposal} htmlFor="proposal_name">{proposal.proposal}</option>
		);
		
		let loggedIn = null;
		loggedIn = this.props.proposal ?  <Redirect to="/selection/home/" /> : null;
		
		return (
		<main id="choose-proposal">
			<section>
			<h1>Proposals</h1>
			<form>
				<legend>Select proposal to manage</legend>
				<p>
					<select name="proposal_name" onChange={event => {this.getSelection(event)}}>
						<option value="">...</option>
					{options}
					</select>
				</p>
				<button onClick={event => this.props.logIn(this.state.chosenName)} >Confirm selection</button>
			</form>
			</section>
			{loggedIn}
		</main>
		);
	}
}

export default ProposalSelection;
