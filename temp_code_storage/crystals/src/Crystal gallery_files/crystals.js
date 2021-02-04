/* Defines actions of buttons and icons in the /compound_screen/crystals
 * path. Uses functions defined in crystals_functions.js	*/

document.addEventListener('DOMContentLoaded', () => {
	updateTotals();
	
	console.log('loaded crystals.js');
	//rejecting or recycling a crystal	
	document.querySelectorAll('.bin-pic').forEach(image => {
		image.onclick = () => {
			const parentTile = image.parentElement.parentElement;
			
			hideTile(image);
			
			//show new element in the other div
			if (parentTile.parentElement.className === 'accepted') {
				showTile(parentTile, '.rejected');
				updateCounter(parentTile, '.rejected');			
			}
			else {
				showTile(parentTile, '.accepted');
				updateCounter(parentTile, '.accepted');
			}
			
		}
	})
			
	//show new plate
	document.getElementById('texrank').onclick = () => {
		document.querySelectorAll('section').forEach(section => {
			section.style.display = 'block';
		})
		document.getElementById('summary').hidden = false;
		document.getElementById('crystal-form').style.display = 'grid';
	}
		
	//show all crystal details in the plate
	document.querySelectorAll('section').forEach(section => {
		let checker = new MassInputChecker(section, '.show-icon', '.show-all');
		checker.addListeners();
	})
	
	//hide all crystal details in the plate
	document.querySelectorAll('section').forEach(section => {
		let checker = new MassInputChecker(section, '.hide-icon', '.hide-all');
		checker.addListeners();
	})
	
	//hide and show crystals in a used plate
	document.querySelectorAll('.all-used').forEach(plate => {
		let group = new Hideable(plate, '.plate-body', 'block', '.hide-plate', '.show-plate', 'inline-block' );
		group.addListeners();		
	})	
	
	//show and hide used crystals in a partially used plate 
	document.querySelectorAll('.used-div').forEach(div => {
		let group = new Hideable(div, '.used', 'flex', '.hide-used-crystals', '.show-used-crystals', 'inline-block');
		group.addListeners();
	})	
	
})
