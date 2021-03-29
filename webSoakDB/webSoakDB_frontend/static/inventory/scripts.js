
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
