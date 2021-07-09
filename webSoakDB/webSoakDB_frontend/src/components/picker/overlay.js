import React from 'react';

class Overlay extends React.Component {

    constructor(props){
        super(props);
        this.state = {dots: '.'};
    }
	
    componentDidMount(){
        this.animateDots();
    }

	animateDots(){
		let text = '. ';
		for (let i = 0; i < 30; i++){
			setTimeout(()=>{
				text = text + '. ';
				this.setState({dots : text});
			},i * 1000);
		}
	}

    render(){
        const dots = this.state.dots;
        return(
            <div id="overlay">
                <div id="overlay-message">
                    <p>Processing data. Please wait.</p>
                    <p id="dots">{dots}</p>
                </div>
            </div>
        );
    }
}

export default Overlay;