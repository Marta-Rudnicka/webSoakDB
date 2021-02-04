import React from 'react';


class Navigation extends React.Component {
	render() {
		return (
			<nav>
				<a href="solvents/source.html">Go to solvent characterization | </a><a href="all.html"> See all in one table (read only)</a>
				<br/><hr/>
				<a href="source"> Source compounds | </a>
				<a href="crystals">Crystals | </a>
				<a href="batches">Batches | </a>
				<a href="soak">Soak | </a>
				<a href="cryo">Cryo | </a>
				<a href="harvesting">Harvesting | </a>
				<a href="collection">Data collection</a>
			</nav>
		);
	}
}

export default Navigation;
