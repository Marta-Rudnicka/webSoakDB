import React from 'react';
import axios from 'axios';

import { Redirect } from "react-router-dom";

class ProposalSelection extends React.Component {
	
	constructor(props){
		super(props);
		this.state = {
			proposals: [],
			chosenProject: "",
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
      
     // this.props.logIn(null);
	}
	
	getSelection(e){
		this.setState({chosenProject: e.target.value});
	}
	
	getProposalString(project){
		if (! project.auth[0]) {
			return null;
		}
		const proposal_visit = project.auth[0].proposal_visit;
		const regex =  /([A-Za-z0-9_]+)(\-[0-9]+)/;
		const m = proposal_visit.match(regex);
	
		if (m===null){
			return null;
		}
		else{
			return m[1];
		}
	}
	
	render(){
		
		const options = this.state.proposals.map((proposal, index) => {
			 return (
			<option
				key={index} 
				value={proposal.id} 
				htmlFor="proposal_name"
			>
				{this.getProposalString(proposal)}: {proposal.auth[0] ? proposal.auth[0].project : null}
			</option>);
		});
		
		//let loggedIn = null;
		let loggedIn = this.props.proposal ?  <Redirect to="/selection/home/" /> : null;
		
		return ( 
			
			<main id="choose-proposal">
				<div>{loggedIn} </div>
				<section>
					<h1>Proposals</h1>
					<form>
						<legend>Select project to manage</legend>
						<p>
							<select name="proposal_name" onChange={event => {this.getSelection(event)}}>
								<option value="">...</option>
							{options}
							</select>
						</p>
						<button onClick={event => this.props.logIn(this.state.chosenProject)} >Confirm selection</button>
					</form>
			</section>
		</main>
		);
	}
}

export default ProposalSelection;