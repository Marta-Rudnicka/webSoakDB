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
		const well = showHide(this.props.show_well);
		const code = showHide(this.props.show_code);
		const smiles = showHide(this.props.show_smiles);
		const structure = showHide(this.props.show_structure);
		const concentration = showHide(this.props.show_concentration);
		const mw = showHide(this.props.show_mw);
		const p3 = showHide(this.props.show_p3);
		const p4 = showHide(this.props.show_p4);
		
		const display_options = [
			["show_well", "Show Well"],
			["show_code", "Show Compound Code"],
			["show_smiles", "Show SMILES"],
			["show_structure", "Show 2D Structure"],
			["show_concentration", "Show Concentration"],
			["show_mw", "Show Molecular Weight"],
			["show_p3", "Show [Property3]"],
			["show_p4", "Show [Property4]"]
		];
		
		const ths = display_options.map( option => {
			return	<th className={this.props[option[0]] ? "" : "hidden"}>{option[1]} 
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
	
	/*	return (
		
			<thead>
				<tr>
					<th className="row-no">Row no.</th>
					<th className={well}>Well 
						<br />
						<button className="in-table" id="hide-well" >Hide</button>
					</th>
					<th className={code}>Code
						<br />
						<button className="in-table" id="hide-code">Hide</button>
					</th>
					<th className={smiles}>SMILES
						<br />
						<button className="in-table" id="hide-smiles">Hide</button>
					</th>
					<th className={structure}>2D structure
						<br />
						<button className="in-table" id="hide-c2d">Hide</button>
					</th>
					<th className={concentration}>Concentration
						
						<br />
						<button className="in-table" id="hide-concentration">Hide</button>		
					</th>
					<th className={mw}>Molecular weight
						
						<br />
						<button className="in-table" id="hide-p2">Hide</button>
					</th>
					<th className={p3}>[Property3]
					
						<br />
						<button className="in-table" id="hide-p3">Hide</button>
					</th>
					<th className={p4}>[Property4]
						<br />
						<button className="in-table" id="hide-p4">Hide</button>
					</th>
				</tr>
			</thead> 
		
		);
		*/
	}
}

export default TableHeader;
