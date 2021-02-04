import React from 'react';

class ExportBar extends React.Component {
	
	render() {
		return (
			<div>
				Export data as:
				<button className="small-button" id="csv-export">csv file</button>
				<button className="small-button" id="other-export">some other format</button>
			</div>
		);
	}
}

export default ExportBar;
