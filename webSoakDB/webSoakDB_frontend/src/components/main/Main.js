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
		this.refreshAfterUpload = this.refreshAfterUpload.bind(this);
		this.updateLibrarySelection = this.updateLibrarySelection.bind(this);
		this.updateSubsetSelection = this.updateSubsetSelection.bind(this);
		this.trackUnsavedChanges = this.trackUnsavedChanges.bind(this);
		
		this.state = {
			page: 'Home',
			lookup_args: {},
			unsavedChanges: false,
		}
	}
	
	changeMainPage(page, safe){
		if (['Home', 'Picker', 'Summary', 'PlateLookup'].includes(page)){
			if (safe){
				this.setState({page: page});
			}
			else {
				if (window.confirm("You have some unsaved changes in your selection. Do you want to discard them and leave this page?")){
					this.setState({page: page});
					this.trackUnsavedChanges(false);
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
				return <Home key="home" handleClick={this.changeMainPage} proposal={this.props.proposal}/>
				break;
			case 'Picker':
				return <Picker 
						key="picker" 
						showPlate={this.showPlate} 
						proposal={this.props.proposal} 
						updateLibrarySelection={this.updateLibrarySelection} 
						updateSubsetSelection={this.updateSubsetSelection} 
						changeMainPage={this.changeMainPage}
						trackUnsavedChanges={this.trackUnsavedChanges}
						initialLibs={this.state.initialLibs}
						initialSubsets={this.state.initialSubsets}
						refreshAfterUpload={this.refreshAfterUpload}
						/>
				break;
			case 'Summary':
				return <Summary key="summary" showPlate={this.showPlate} 
						proposal={this.props.proposal} 
						libSelection={this.props.libSelection}
						updateLibrarySelection={this.updateLibrarySelection}
						updateSubsetSelection={this.updateSubsetSelection}
						changeMainPage={this.changeMainPage}
						unsavedChanges={this.state.unsavedChanges}
						trackUnsavedChanges={this.trackUnsavedChanges}
						lookup_args={this.state.lookup_args}
						/>
				break;
			case 'PlateLookup':
				return <PlateLookup 
							library={this.state.lookup_args[0]} 
							plate={this.state.lookup_args[1]} 
							is_a_plate={this.state.lookup_args[2]} 
							showPlate={this.showPlate}
							lookup_args={this.state.lookup_args}/>
				break;
		}
	}
	
	showPlate(library, collection, is_a_plate, is_a_preset){
		this.setState({lookup_args: {library : library, collection : collection, is_a_plate: is_a_plate, is_a_preset : is_a_preset}});
		this.setState({page: 'PlateLookup'});
		this.trackUnsavedChanges(false);
	}
	
	updateLibrarySelection(idArray){
		event.preventDefault()
		const apiUrl='api/update_proposal_selection/' + this.props.proposal.name + '/';
		axios.patch(apiUrl, {libraries: idArray}) 
		.catch(error => {
			console.log('updateLibrarySelection: axios error:');
			console.log(error)
		})
		
		this.props.logIn(this.props.proposal.name);
	}
	
	updateSubsetSelection(idArray){
		event.preventDefault()
		const apiUrl='api/update_proposal_selection/' + this.props.proposal.name + '/';
		axios.patch(apiUrl, {subsets: idArray}) 
		.catch(error => {
			console.log('updateLibrarySelection: axios error:');
			console.log(error)
		})
		
		setTimeout(() => {this.props.logIn(this.props.proposal.name)}, 1000);
	}
	
	refreshAfterUpload(){
		this.props.logIn(this.props.proposal.name);	
	}
	
	//to give warning when user risks discarding changes
	trackUnsavedChanges(bool){
		this.setState({unsavedChanges: bool})
	}
	
    render() {
		const content = this.mainPage();
		
        return (
		<div id="all-container">
         <div>
			<nav>
				<span className="pseudo-link" onClick={event => this.changeMainPage('Home', !this.state.unsavedChanges)}>Home | </span>
				<span className="pseudo-link" onClick={() => this.changeMainPage('Picker', !this.state.unsavedChanges)}> Select compounds |</span>
				<span className="pseudo-link" onClick={event => this.changeMainPage('Summary', !this.state.unsavedChanges)}>Selection summary | </span>
				<a href="/inventory/">Inventory | </a>
				<span className="pseudo-link" onClick={event => this.props.logIn(null)}>Log out | </span>
			</nav>
			{content}
        </div>
       </div>
        )
    }
}

export default Main;
