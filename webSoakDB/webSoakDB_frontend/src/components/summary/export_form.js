import React from 'react';
import axios from 'axios';
import CSRFToken from '../picker/csrf.js';

class SubsetSelect extends React.Component {
	constructor(props) {
		super(props);
		this.state = { library: null};
	}
	
	
	componentDidMount() {
		const apiUrl = '/api/library_detail/' + this.props.subset.library.id + '/';
		
		axios.get(apiUrl)
			.then(res => {
			const library = res.data;
			this.setState({library});
		});
	}
	
	countCurrent(){
		if (this.state.library){
			let c = 0;
			this.state.library.plates.forEach(p => {
				if(p.current){
					c++;
				}
			});
			return c;
		}
	}
	
	render(){
		
		let options = null;
		let selects = null;
		if (this.state.library){
			options = this.state.library.plates.map((plate, index) =>
				<option key={index} value={plate.id}>{plate.name} {plate.current ? '(current)' : ''}</option>
			);
			if (this.countCurrent() === 1){
				selects = (<p>
								<select key="1" name={this.props.subset.library.id} required>
									<option value="">Select plate...</option>
									{options}
								</select>
							</p>
						);
			}
			else {
				selects = [];
				for (let i = 1; i <= this.countCurrent(); i++){
					const sel = (<p>
						<select key ={i} name={this.props.subset.library.id + '-' + i} required>
							<option value="">Select plate...</option>
							{options}
						</select>
						</p>);
					selects.push(sel);				
				}
			}
		}
		
		return(
		<tr>
			<td>
				<label>Select library plate(s) for {this.props.subset.library.name} ({this.props.subset.origin})</label>
			</td>
			<td>
				{selects}
			</td>
		</tr>
		
		);
	}
}

class ExportForm extends React.Component {
	
	render() {
		
		const subsets = this.props.proposal.subsets.map((subset, index) =>{
			return <SubsetSelect key={index} subset={subset} />;
		});
		
		
		return (
			<form method="post" action="/downloads/export-for-soakdb/">
				<CSRFToken />
				<input type="hidden" name="proposal" value={this.props.proposal.name} />
				<legend>Export selection for SoakDB </legend>
				{ this.props.proposal.subsets.length > 0 ? <p><strong>Select library plates for cherry-picking lists and subsets</strong></p> : null }
				<table className="table">
					<tbody>
						{subsets}
					</tbody>
				</table>
				<button type="submit">Generate file</button>
			</form>
		); 
	}
}

export default ExportForm;


