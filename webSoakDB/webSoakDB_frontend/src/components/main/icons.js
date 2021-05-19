import React, { Component } from 'react';

export class ZoomInIcon extends Component {
	
	render() {
		 return (
			<svg xmlns="http://www.w3.org/2000/svg" width={this.props.size} height={this.props.size} onClick={() => this.props.handleClick()} fill="currentColor" className="bi bi-zoom-in" viewBox="0 0 16 16">
			  <path fillRule="evenodd" d="M6.5 12a5.5 5.5 0 1 0 0-11 5.5 5.5 0 0 0 0 11zM13 6.5a6.5 6.5 0 1 1-13 0 6.5 6.5 0 0 1 13 0z"/>
			  <path d="M10.344 11.742c.03.04.062.078.098.115l3.85 3.85a1 1 0 0 0 1.415-1.414l-3.85-3.85a1.007 1.007 0 0 0-.115-.1 6.538 6.538 0 0 1-1.398 1.4z"/>
			  <path fillRule="evenodd" d="M6.5 3a.5.5 0 0 1 .5.5V6h2.5a.5.5 0 0 1 0 1H7v2.5a.5.5 0 0 1-1 0V7H3.5a.5.5 0 0 1 0-1H6V3.5a.5.5 0 0 1 .5-.5z"/>
			</svg>
		)
	}
}

export class ChevronDown extends Component {
	render() {
		 return (
			<svg xmlns="http://www.w3.org/2000/svg" width={this.props.size} height={this.props.size} onClick={() => this.props.handleClick()} fill="currentColor" className="bi bi-chevron-down" viewBox="0 0 16 16">
			  <path fillRule="evenodd" d="M1.646 4.646a.5.5 0 0 1 .708 0L8 10.293l5.646-5.647a.5.5 0 0 1 .708.708l-6 6a.5.5 0 0 1-.708 0l-6-6a.5.5 0 0 1 0-.708z"/>
			</svg>
		)
	}
}

export class XIcon extends Component {
	render() {
		 return (
			<svg xmlns="http://www.w3.org/2000/svg" width={this.props.size} height={this.props.size} onClick={() => this.props.handleClick()} fill="currentColor" className="bi bi-x-square" viewBox="0 0 16 16">
			  <path d="M14 1a1 1 0 0 1 1 1v12a1 1 0 0 1-1 1H2a1 1 0 0 1-1-1V2a1 1 0 0 1 1-1h12zM2 0a2 2 0 0 0-2 2v12a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V2a2 2 0 0 0-2-2H2z"/>
			  <path d="M4.646 4.646a.5.5 0 0 1 .708 0L8 7.293l2.646-2.647a.5.5 0 0 1 .708.708L8.707 8l2.647 2.646a.5.5 0 0 1-.708.708L8 8.707l-2.646 2.647a.5.5 0 0 1-.708-.708L7.293 8 4.646 5.354a.5.5 0 0 1 0-.708z"/>
			</svg>
		)
	}
}
