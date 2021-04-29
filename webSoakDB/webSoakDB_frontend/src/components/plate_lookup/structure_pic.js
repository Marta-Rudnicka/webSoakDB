import React from 'react';
import { LazyLoadImage } from 'react-lazy-load-image-component';
import { ZoomInIcon } from '../icons.js';


class StructurePic extends React.Component {
	
	constructor(props){
		super(props);
		this.state = {
			url_suffix: "",
			attempts : 0,
			url_start : 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/', 
			url_end: '/PNG?'}
	}
	
	
	//adds new query string to the end of the URL to force re-upload
	reload() {
		let att = this.state.attempts;
		
		if(att < 11){
			att ++;
			setTimeout(this.setState({url_suffix: Date.now(), attempts: att}), 2000);
		}
		else if(att < 13){ //try different service
			att ++;
			this.setState({url_start : 'https://cactus.nci.nih.gov/chemical/structure/',
							url_end : '/image?',
							attempts: att})
		}
	}
	
		
	render() {
		let url_smiles = this.props.smiles.replace(/#/g, '%23');
		url_smiles = url_smiles.replace(/\//g, '%2F');
		const img_url = this.state.url_start + url_smiles + this.state.url_end + this.state.url_suffix;
	  
		return ( 
		<div className="img2d">
			<LazyLoadImage
			  alt="Loading image..."
			  src={img_url}
			  onError={() => this.reload()}
			/>
			&nbsp;
			<ZoomInIcon size="20"/>
		</div>);
	 }
}

export default StructurePic;

