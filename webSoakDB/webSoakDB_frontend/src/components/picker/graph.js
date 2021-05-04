import React from 'react';

class Graph extends React.Component{
	render(){
		const path = '/media/images/graphs/' + this.props.property + '/' + this.props.id + '.svg';
		return(
		<div>
			<figcaption>{this.props.caption}</figcaption>
			
			<img src={path} className="hist-small" alt="failed to load image"/>
		</div>
		);
	}
}

export default Graph;
