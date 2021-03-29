//document.addEventListener('DOMContentLoaded', arrange());
/*		
function arrange(){
	const spans = document.querySelectorAll('.c');
	spans.forEach(span => {
		const re = new RegExp('[A-Z]{1,2}[1-9]$');
		let well = span.dataset.well
		
		//insert leading zero where needed
		if (re.test(well)){
			well = well.slice(0, well.length-1) + '0' + well[well.length-1];
		}
		
		const cell = document.getElementById(well)
		if (cell) {
			cell.appendChild(span);
		}
	});
}		

console.log('loaded');

function test(){
	console.log('fired test');
	setTimeout(() => arrange(), 100);
	
	console.log('exiting test');
}

document.addEventListener('DOMContentLoaded', test());
*/
