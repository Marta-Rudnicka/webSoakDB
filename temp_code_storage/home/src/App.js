import React from 'react';
import './home.css';

const proposal = "Test proposal x"
//function App() {
class App extends React.Component {
	
	constructor(props) {
		super(props);
		this.state = {proposal: proposal};
	}
	render() {
		return (
		<div id="all">
			<nav>
		<a href="/playground/">Home</a> | <a href="/playground/proposal">Change proposal</a>
		</nav>

			<h1>SoakDB home - {this.state.proposal}</h1>
			<main>
			<section id="prep">
				<h2>Prepararation</h2>
				<ul>
					<li>
						<p><a href="/playground/picker">Selection of compounds</a> </p>
						<p>the place to select libraries, cherry-pick compounds, and upload the data of your own library</p>
					</li>
					<li>
						<p><a href="/playground/summary">Summary</a> </p>
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
						<p><a href="/playground/all">Data summary</a> </p><p> the place to see and download all the data produced in your experiment (at any stage)</p>
					</li>
				</ul>
			</section>
			</main>
		</div>

		); 
	}

}

export default App;
