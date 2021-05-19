import React from 'react';
import { LazyLoadImage } from 'react-lazy-load-image-component';
import { ZoomInIcon } from '../main/icons.js';


class StructurePic extends React.Component {
		
	render() {
		const img_url = '/molecule/' + this.props.id + '/'
		return ( 
		<div className="img2d">
			<LazyLoadImage
			  alt="Loading image..."
			  src={img_url}
			/>
			&nbsp;
			<ZoomInIcon size="20"/>
		</div>);
	 }
}

export default StructurePic;

