import React from 'react';
import pic1 from './dsi-ClogP.png';
import pic2 from './dsi-Fsp3.png';
import pic3 from './dsi-mol-weight.png';
import pic4 from './dsi-TPSA.png';

class Graphs extends React.Component {

	render(){
		return (
			<section id="graphs">
				<h2>Summary</h2>
				<div>
					<img alt="graph" src={pic1} />
					<img alt="graph" src={pic2} />
					<img alt="graph" src={pic3} />
					<img alt="graph" src={pic4} />
				</div>
			</section>
		)
	}
}

export default Graphs;
