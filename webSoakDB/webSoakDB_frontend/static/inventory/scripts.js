
console.log('LOADED');

function submitDelete(id, object_name){
	str = "Are you sure you want to delete " + object_name + "?"
	if(confirm(str)){
		document.forms['delete-' + id].requestSubmit();
	}
	else{
		event.preventDefault();
	}
}

function findAndCheck(str){
	const feedback = document.getElementById("finder-feedback");
	const row =  document.getElementById('w-' + str.toUpperCase());
	if (!row){
		feedback.innerHTML = "No such well in the plate map: <strong>" + str.toUpperCase() + "</strong>";
	}
			
	if (row.className === "inactive-row"){
		feedback.innerHTML = "<strong>" + str.toUpperCase() + "</strong> is already inactive";
	}
	else {
		const input = row.querySelector('input');
		input.checked = true;
		moveToTop(input.id);
		feedback.innerHTML = "Selected <strong>" + str.toUpperCase() + "</strong> for deactivation";
	}
}
		
function selectToDeactivate(){
	const value = document.getElementById("well").value;
	findAndCheck(value);
	document.getElementById("well").value="";
}
		
function moveToTop(id){
	const box = document.getElementById(id);
	if (box.checked===true) {
		const firstLine = document.getElementById('first-line');
		const row = box.closest('tr');
		const parentNode = firstLine.parentNode
		setTimeout(function() {parentNode.insertBefore(row, firstLine.nextElementSibling)}, 300);
	}
}

function showRows(className){
	const hid = className + ' hidden';
	console.log("show: ", className);
	document.querySelectorAll('#table tr').forEach(tr => {
		//console.log('inspecting: ', tr)
		
		if (tr.className===hid){
			//console.log('found hidden: ', tr);
			tr.className=className
		}
	});
}

function hideRows(className){
	const hid = className + ' hidden';
	console.log("hide: ", className);
	document.querySelectorAll('#table tr').forEach(tr => {
		//console.log('inspecting: ', tr)
		if (tr.className===className){
			//console.log('found visible: ', tr);
			tr.className=hid;
		}
	});
}

function filterRows(value){
	console.log('fired filterRows');
	if (value==="all"){
		showRows("active-row");
		showRows("inactive-row");
		}
	else if (value=="active-only"){
		showRows("active-row");
		hideRows("inactive-row");
	}
	else{
		showRows("inactive-row");
		hideRows("active-row");
	}

}


function show(id){
	console.log('show');
	//subTable = document.getElementById(id);
	//subTable.className = "";
	const details = String.raw`{% include "availability-details.html" with initial_class="" %}`;
	const tId = 'container-table-' + id;
	document.getElementById(tId).innerHTML = `${details}`;
	hider =  document.getElementById(id + "-hider");
	hider.className = "";
	shower = document.getElementById(id + "-shower");
	shower.className = "hidden";
  }

  function hide(id){
	console.log('hide');
	subTable = document.getElementById(id);
	subTable.className = "hidden";
	hider =  document.getElementById(id + "-hider");
	hider.className = "hidden";
	shower = document.getElementById(id + "-shower");
	shower.className = "";
  }