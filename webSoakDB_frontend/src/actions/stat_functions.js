//helper functions for statistics
import React from 'react';

export function mean(array){
	/*returns arithmetic mean of all the numbers in an array; throws 
	 * error if fed something else than numbers */
	let sum = 0;
	let items = 0;
	array.forEach(number => {
		if (isNaN(parseFloat(number))){
			throw "Error in mean(): Not a number!"
		}
		sum = sum + number;
		items ++;
	});
	return sum/items;
}



export function getAttributeArray(array, attr){
	/* returns an array of values of the attribute $attr of objects
	 * in $array*/
	let outputArray = [];
	array.forEach(item => {
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

export function shareAllElements(array1, array2){
	/*returns true if array1 and array2 contain exactly the same elements;
	 * otherwise returns false; USE ONLY FOR ARRAYS WITH UNIQUE ELEMENTS*/
	
	if(array1.length !== array2.length){
		return false;
	}
	
	let sameElements = true;
	
	array1.forEach(element => {
		if (!array2.includes(element)){
			sameElements = false;
		}
	});
	return sameElements;
}

export function addUniqueCompounds(oldArray, newArray){
	/* Takes in arrays of SourceWell objects. Returns an array that 
	 * contains all unique Compounds objects that are attributes in 
	 * oldArray and newArray (no duplicates). Warning: that does 
	 * not necessary mean  unique chemical compounds! In rare cases 
	 * Compound objects with different ids share SMILES string
	 * */
	
	//const compounds = getAttributeArray(oldArray, "compound")
	const ids = getAttributeArray(oldArray, "id");
	let idSet = new Set();
	ids.forEach(id => idSet.add(parseInt(id)));
	newArray.forEach(compound => {
		if (!idSet.has(compound.id)){
			oldArray.push(compound);
			idSet.add(parseInt(compound.id));
		}
	});
	return oldArray;
}
