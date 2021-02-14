import React, { Component, Fragment } from 'react';
import ReactDOM from 'react-dom';

import axios from 'axios';

import Home from '../home/Home.js';
import Picker from '../picker/Picker.js';
import Summary from '../summary/Summary.js';
import PlateLookup from '../plate_lookup/PlateLookup.js';

import { deepCopyObjectArray, getAttributeArray, mean, shareAllElements } from  '../../actions/stat_functions.js';

class Main extends Component {
	
	constructor(props){
		super(props);
		
		
		this.changeMainPage = this.changeMainPage.bind(this);
		this.showPlate = this.showPlate.bind(this);

		this.updateLibrarySelection = this.updateLibrarySelection.bind(this);
		this.updateSubsetSelection = this.updateSubsetSelection.bind(this);
		this.state = {
			page: 'Home',
		}
	}
	
	componentDidUpdate(prevProps, prevState) {
		if (prevProps.libSelection !== this.props.libSelection) {
			this.setState({libSelection : this.props.libSelection});
		}
	}

	
	changeMainPage(page, safe){
		if (['Home', 'Picker', 'Summary'].includes(page)){
			if (safe){
				this.setState({page: page});
			}
			else {
				if (window.confirm("You have not saved changes to your selection. Do you want to leave page without saving?")){
					this.setState({page: page});
				}
			}
		}
		else {
			console.log('changeMainPage: invalid argument: ', page);
		}
	}
	
	mainPage(){
		switch(this.state.page){
			case 'Home':
				return <Home key="home" handleClick={this.changeMainPage}/>
				break;
			case 'Picker':
				return <Picker 
						key="picker" 
						showPlate={this.showPlate} 
						proposal={this.props.proposal} 
						updateLibrarySelection={this.updateLibrarySelection} 
						updateSubsetSelection={this.updateSubsetSelection} 
						changeMainPage={this.changeMainPage} 
						/>
				break;
			case 'Summary':
				return <Summary key="summary" showPlate={this.showPlate} proposal={this.props.proposal} libSelection={this.props.libSelection} updateLibrarySelection={this.updateLibrarySelection}/>
				break;
		}
	}
	
	showPlate(library, plate, current){
		this.setState({page: <PlateLookup library={library} plate={plate} current={current} />});
	}
	
	updateLibrarySelection(idArray, page){
		event.preventDefault()
		const apiUrl='api/update_proposal_selection/' + this.props.proposal.name + '/';
		axios.patch(apiUrl, {libraries: idArray}) 
		.catch(error => {
			console.log('updateLibrarySelection: axios error:');
			console.log(error)
		})
		
		this.props.logIn(this.props.proposal.name);
	}
	
	updateSubsetSelection(idArray, page){
		event.preventDefault()
		const apiUrl='api/update_proposal_selection/' + this.props.proposal.name + '/';
		axios.patch(apiUrl, {subsets: idArray}) 
		.catch(error => {
			console.log('updateLibrarySelection: axios error:');
			console.log(error)
		})
		
		this.props.logIn(this.props.proposal.name);
	}
	
	
    render() {
		const content = this.mainPage();
		
        return (
         <div>
			<nav>
				<span className="pseudo-link" onClick={event => this.changeMainPage('Home', true)}>Home | </span>
				<span className="pseudo-link" onClick={() => this.changeMainPage('Picker', true)}> Select compounds |</span>
				<span className="pseudo-link" onClick={event => this.changeMainPage('Summary', true)}>Selection summary | </span>
				<span className="pseudo-link" onClick={event => this.props.logIn(null)}>Log out | </span>
			</nav>
			{content}
        </div>
        )
    }
}

export default Main;
