import React from 'react';
import axios from 'axios';
import { getProposalString } from '../../actions/stat_functions';
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
      
      this.props.logIn(null);
	}
	
	getSelection(e){
		this.setState({chosenProject: e.target.value});
	}
	
	
	
	render(){
		
		const options = this.state.proposals.map((proposal, index) => {
			 return (
			<option
				key={index} 
				value={proposal.id} 
				htmlFor="proposal_name"
			>
				{getProposalString(proposal)}: {proposal.auth[0] ? proposal.auth[0].project : null}
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