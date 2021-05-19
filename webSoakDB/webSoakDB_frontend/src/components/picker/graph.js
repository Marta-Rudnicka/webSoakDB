import React from 'react';
import axios from 'axios';

class Graph extends React.Component{
  constructor(props){
    super(props)
    this.state = {
		status: null,
		path: null,
	} 
  }
  
  componentDidMount(){
	this.load();
  }
  
  componentDidUpdate(prevProps, prevState){
	if(prevProps !== this.props){
		this.load();
    }
  }
  
  load(){
	let path;
	if(this.props.public){
		path ='/media/html_graphs/' + this.props.type + '/' + this.props.id + '/' + this.props.property + '.html';
	}
	else{
		path = '/histogram/' + this.props.type + '/' + this.props.id + '/' + this.props.property + '/';
	}
	this.setState({path: path});
	
	axios.get(path)
    .then(res =>{
	this.setState({status: res.status});
	});
  }
    
  render(){
    
    let content = <iframe src={this.state.path} loading="lazy"/>;
    
    if (this.state.status===204){
		content = <p>Data unavailable</p>;
	}
    
    return(
		<div>{content}</div>
    );
  }
}

export default Graph;
