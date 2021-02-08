import React from 'react';

const proposal = "Test proposal x"

class Home extends React.Component {
	
	constructor(props) {
		super(props);
		this.state = {proposal: proposal};
	}
	render() {
		return (
		<div id="home">
			<h1>SoakDB home - {this.state.proposal}</h1>
			<main>
			<section id="prep">
				<h2>Prepararation</h2>
				<ul>
					<li>
						<p><span className="pseudo-link" onClick={() => this.props.handleClick("Picker")}>Selection of compounds</span> </p>
						<p>the place to select libraries, cherry-pick compounds, and upload the data of your own library</p>
					</li>
					<li>
						<p><span className="pseudo-link" onClick={() => this.props.handleClick("Summary")}>Summary</span> </p>
						<p> the place to see and manage your protein and compound data</p>
					</li>
				</ul>
			</section>
			
			<section id="exp">
				<h2>Experiment</h2>
				<ul>
					<li>
						<p><a href="/playground/dummy">Solvent characterization application</a> </p><p> software used for managing the solvent characterisation process in the lab</p>
					</li>
					<li>
						<p><a href="/playground/dummy">Compound screen application</a> </p><p> software used for managing the experiment in the lab</p>
					</li>
					<li>
						<p><a href="/playground/dummy">Data summary</a> </p><p> the place to see and download all the data produced in your experiment (at any stage)</p>
					</li>
				</ul>
			</section>
			</main>
		</div>

		); 
	}
}

export default Home;
