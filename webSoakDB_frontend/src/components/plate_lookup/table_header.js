import React from 'react';
import {display_options} from './display_options';

function showHide(value) {
	if (value){
		return "";
	}
	else {
		return "hidden";
	}	
}

class TableHeader extends React.Component {
	

	
	render() {
	//console.log('display in TableHeader: ', this.props.display)
		const ths = display_options.map((option, index) => {
			return	<th key={index} className={this.props.display[option[0]] ? "" : "hidden"}>{option[1]} 
						<br />
						<button className="in-table" onClick={event => this.props.onButtonClick([option[0]])}>Hide</button>
					</th>
			});
		
		return (
		<thead>
			<tr>
				<th className="row-no">Row no.</th>
				{ths}
			</tr>
		</thead>
		);
	}
}

export default TableHeader;
