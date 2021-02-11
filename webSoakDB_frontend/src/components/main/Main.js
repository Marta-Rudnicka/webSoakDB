import React, { Component, Fragment } from 'react';
import ReactDOM from 'react-dom';

import axios from 'axios';

import Home from '../home/Home.js';
import Picker from '../picker/Picker.js';
import Summary from '../summary/Summary.js';
import PlateLookup from '../plate_lookup/PlateLookup.js';

class Main extends Component {
	
	constructor(props){
		super(props);
		
		this.changeMainPage = this.changeMainPage.bind(this);
		this.showPlate = this.showPlate.bind(this);

		this.updateLibrarySelection = this.updateLibrarySelection.bind(this);
		this.state = {page: <Home />};
		}
	
	componentDidMount() {
		const proposalApiUrl = 'api/proposals/' + this.props.proposalName;
		axios.get(proposalApiUrl)
			.then(res => {
			const proposal = res.data;
			this.setState({ proposal });
      });
      		
	}
	
	changeMainPage(page){
		switch(page){
			case 'Home':
				this.setState({page: <Home handleClick={this.changeMainPage}/>});
				break;
			case 'Picker':
				this.setState({page: <Picker showPlate={this.showPlate}  proposal={this.state.proposal} updateLibrarySelection={this.updateLibrarySelection} changeMainPage={this.changeMainPage}/>});
				break;
			case 'Summary':
				this.setState({page: <Summary showPlate={this.showPlate} proposalName={this.props.proposalName} updateLibrarySelection={this.updateLibrarySelection}/>});
				break;
		}
	}
	
	showPlate(library, plate, current){
		this.setState({page: <PlateLookup library={library} plate={plate} current={current} />});
	}
	
	updateLibrarySelection(idArray){
		const apiUrl='api/update_proposal_selection/' + this.props.proposalName + '/';
		
		axios.patch(apiUrl, {libraries: idArray}) 
		.catch(error => {
			console.log(error)
		})
		
		this.componentDidMount();
		this.changeMainPage('Summary');
	}
	
    render() {
		const app = this.state.page;
		
        return (
         <div>
			<nav>
				<span className="pseudo-link" onClick={event => this.changeMainPage('Home')}>Home | </span>
				<span className="pseudo-link" onClick={event => this.changeMainPage('Picker')}> Select compounds |</span>
				<span className="pseudo-link" onClick={event => this.changeMainPage('Summary')}>Selection summary | </span>
				<span className="pseudo-link" onClick={event => this.props.logIn('')}>Log out | </span>
			</nav>
			{app}
        </div>
        )
    }
}

export default Main;

