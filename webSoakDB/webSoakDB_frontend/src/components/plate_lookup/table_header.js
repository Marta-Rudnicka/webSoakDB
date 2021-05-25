import React from 'react';
import {display_options} from './display_options';

class TableHeader extends React.Component {
	
	render() {
		let extra = null;
		if (this.props.is_a_preset){
			extra = <th className={this.props.display.show_library ? "library" : "hidden"}>
						Library	<br />
						<button className="in-table" onClick={event => this.props.onButtonClick(show_library)}>Hide</button>
					</th>;
		}
		
		let ths = display_options.map((option, index) => {
			return	<th key={index} className={this.props.display[option[0]] ? option[0] : "hidden"}>
						{option[1]}
						<br />
						<button className="in-table" onClick={event => this.props.onButtonClick([option[0]])}>Hide</button>
					</th>
			});
		
		return (
		<thead className="undefined sticky-header">
			<tr>
				<th className="row-no">Row no.</th>
				{extra}
				{ths}
			</tr>
		</thead>
		);
	}
}

export default TableHeader;
