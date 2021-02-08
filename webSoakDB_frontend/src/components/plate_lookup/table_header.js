import React from 'react';

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
		const display_options = [
			["show_well", "Well"],
			["show_code", "Compound Code"],
			["show_smiles", "SMILES"],
			["show_structure", "2D Structure"],
			["show_concentration", "Concentration"],
			["show_mw", "Molecular Weight"],
			["show_p3", "[Property3]"],
			["show_p4", "[Property4]"]
		];
		
		const ths = display_options.map((option, index) => {
			return	<th key={index} className={this.props[option[0]] ? "" : "hidden"}>{option[1]} 
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
