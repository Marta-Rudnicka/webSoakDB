import React from 'react';

class ExportBar extends React.Component {
	
	render() {
		return (
			<div>
				<h3>Export data as:</h3>
				<a href={'/downloads/' + this.props.url + '/' + this.props.id + '/'} download><button className="small-button" id="csv-export">CSV {this.props.label}</button></a>
				<a href={'/downloads/' + this.props.url + '-properties/' + this.props.id + '/'} download><button className="small-button" id="csv-export">CSV {this.props.label} with mol. properties</button></a>
			</div>
		);
	}
}

export default ExportBar;
