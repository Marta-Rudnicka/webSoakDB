import React from 'react';

class Graph extends React.Component{
	render(){
		const path = './public/' + this.props.type + '/' + this.props.id + '/' + this.props.property +'.svg';
		return(
		<div>
			<figcaption>{this.props.caption}</figcaption>
			<img src={require(`${path}`).default} className="hist-small" alt="failed to load image"/>
		</div>
		);
	}
}

export default Graph;
