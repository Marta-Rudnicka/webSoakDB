import React, { Component, Fragment } from 'react';
import ReactDOM from 'react-dom';

import Home from './home/Home.js';
import Picker from './picker/Picker.js';
import Summary from './summary/Summary.js';

import { Provider } from 'react-redux';
import store from '../store';

class App extends Component {
	
	constructor(props){
		super(props);
		this.state = {page: <Home />}
		}
	
	changePage(page){
		console.log('Fired changePage, input: ', page);
		switch(page){
			case 'Home':
				this.setState({page: <Home />});
				break;
			case 'Picker':
				this.setState({page: <Picker />});
				break;
			case 'Summary':
				this.setState({page: <Summary />});
				break;
		}
	}
	
    render() {
		const app = this.state.page
        return (
         <div>       
			<nav>
				<span onClick={event => this.changePage('Home')}>Home | </span>
				<span onClick={event => this.changePage('Picker')}> Select compounds |</span>
				<span onClick={event => this.changePage('Summary')}>Selection summary | </span>
				<span onClick={event => console.log('Logout')}>Log out | </span>
			</nav>
			{app}
        </div>
        )
    }
}

ReactDOM.render(<App />, document.getElementById('app'));
