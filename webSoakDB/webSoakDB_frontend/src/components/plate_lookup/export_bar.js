import React from 'react';

class ExportBar extends React.Component {
	
	render() {
		return (
			<div>
				<h3>Export data as:</h3>
				<a href={'/downloads/plate-map/' + this.props.id + '/'} download><button className="small-button" id="csv-export">CSV plate map</button></a>
			</div>
		);
	}
}

export default ExportBar;
