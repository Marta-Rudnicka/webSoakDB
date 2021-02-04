//small functions to make setTimeout easier	
function hideParentTile(element){
	element.parentElement.parentElement.className = 'hidden';
}
	
function show(element){
	element.style.display = 'block';
}
	
function showDiv(element){
	element.className = 'gallery';
}
	
function showHidden(element){
	element.hidden = false;
}
	
function clearStyle(element) {
	element.style="";
}

//update info about how many crystals are accepted, and how many rejected
function updateTotals() {
	var totalAccepted = 0;
	var totalRejected = 0;
	
	document.querySelectorAll('.counter-accepted').forEach(span => {
		totalAccepted = totalAccepted + parseInt(span.innerHTML);
		})
	document.querySelectorAll('.counter-rejected').forEach(span => {
		totalRejected = totalRejected + parseInt(span.innerHTML);
		})
		
	document.getElementById('total-accepted').innerHTML = totalAccepted;
	document.getElementById('total-rejected').innerHTML = totalRejected;
}


function updateCounter(tile, containerClass){
	const counterAccepted = tile.parentElement.parentElement.querySelector(containerClass).parentElement.querySelector('.counter-accepted');
	const counterRejected = tile.parentElement.parentElement.querySelector(containerClass).parentElement.querySelector('.counter-rejected');
	var accepted = parseInt(counterAccepted.innerHTML);
	var rejected = parseInt(counterRejected.innerHTML);		
	
	//update counters of accepted and rejected crystals in each plate
	if (containerClass === '.rejected') {
		counterAccepted.innerHTML = accepted - 1;
		counterRejected.innerHTML = rejected + 1;
	}
	else {
		counterAccepted.innerHTML = accepted + 1;
		counterRejected.innerHTML = rejected - 1;
	}
	updateTotals();
}


//hide the tile that contains image
function hideTile(image) {	
	const parentTile = image.parentElement.parentElement;
	
	//animate hiding a tile
	parentTile.querySelector('.bin-pic').hidden = true;
	
	//tile.querySelector('.info-pic').hidden = true;
	parentTile.style = 'animation-play-state: running; animation-duration: 0.4s';
	parentTile.querySelector('.main-pic').style = 'animation-play-state: running; animation-duration: 0.4s';
	setTimeout(hideParentTile, 400, image);
	setTimeout(clearStyle, 400, parentTile);
	setTimeout(clearStyle, 400, parentTile.querySelector('.main-pic'));
}

//show corresponding tile in the accepted/rejected container
function showTile(tile, containerClass) {
	wellName = tile.id.replace('rejected', '');	
			
	//find corresponding tile to show	
	const galleryTiles = tile.parentElement.parentElement.querySelector(containerClass).querySelectorAll('div');
	var newTile = null;
	galleryTiles.forEach( div => {
		if (div.id.includes(wellName)) {
			newTile = div;
			return;
		}
	})
		
	//animate showing a new tile
	showDiv(newTile);
	newTile.style = 'animation-play-state: running; animation-direction: reverse; animation-duration: 0.4s';
	newTile.querySelector('.main-pic').style = 'animation-play-state: running; animation-direction: reverse; animation-duration: 0.4s';
	setTimeout(clearStyle, 400, newTile);
	setTimeout(clearStyle, 400, newTile.querySelector('.main-pic'));
	setTimeout(showHidden, 420, newTile.querySelector('.bin-pic'));
	setTimeout(showHidden, 420, newTile.querySelector('.show-pic'));
	newTile.className = "gallery";
}
