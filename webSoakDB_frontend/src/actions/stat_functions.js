//helper functions for statistics

export function test(){
	console.log('TEST');

}

export function mean(array){
	let sum = 0;
	let items = 0;
	array.map(number => {
		if (isNaN(parseFloat(number))){
			throw "Error in mean(): Not a number!"
		}
		sum = sum + number;
		items ++;
	});
	return sum/items;
}



export function getAttributeArray(array, attr){
	let outputArray = [];
	array.map(item => {
		if (item[attr] !== undefined){
			outputArray.push(item[attr]);
		}
		else {
			console.log('getAttributeArray: Undefined value of ', attr, ' for ', item, '!');
		}
	});
	return outputArray;
}

export function deepCopyObjectArray(array){
	//this sucks; needs to be changed
	let output = [];
	array.forEach(object => {
		const objectDeepCopy = JSON.parse(JSON.stringify(object));
		output.push(objectDeepCopy);
	});
	return output;
}
