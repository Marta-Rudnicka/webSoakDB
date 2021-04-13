import React from 'react';
import { Link } from "react-router-dom";

class Home extends React.Component {
	
	constructor(props) {
		super(props);
	}
	render() {
		return (
		<div id="home">
			<h1>XChemSPA home - {this.props.proposal.proposal}</h1>
			<main>
			<section id="prep">
				<h2>Prepararation</h2>
				<ul>
					<li>
						<p><Link to="/selection/">Selection of compounds</Link></p>
						<p>the place to select libraries, cherry-pick compounds, and upload the data of your own library</p>
					</li>
					<li>
						<p><Link to="/summary/">Summary</Link></p>
						<p> the place to see and manage your protein and compound data</p>
					</li>
					<li>
						<p><a href="/inventory/">Inventory management</a></p>
						<p>The place to manage XChem in-house libraries, plates etc.</p>
					</li>
				</ul>
			</section>
			
			<section id="exp">
				<h2>Experiment</h2>
				<ul>
					<li>
						<p><a href="/dummy" target="_blank" >Solvent characterization application</a> </p><p> software used for managing the solvent characterisation process in the lab</p>
					</li>
					<li>
						<p><a href="/dummy" target="_blank">Compound screen application</a> </p><p> software used for managing the experiment in the lab</p>
					</li>
					<li>
						<p><a href="/dummy" target="_blank">Data summary</a> </p><p> the place to see and download all the data produced in your experiment (at any stage)</p>
					</li>
				</ul>
			</section>
			</main>
		</div>

		); 
	}
}

export default Home;
