import React from 'react';
import { useParams } from "react-router-dom";
import PlateLookup from './plate_lookup/PlateLookup.js';

/* translates url parameters to appropriate props for PlateLookup*/
export default function UrlTranslator() {
	let { type, id }  = useParams();
	let is_a_plate;
	let is_a_preset;
	
	if (type==="plate"){
		is_a_plate = true;
		is_a_preset = false;
	}
	else if (type=="preset"){
		is_a_plate = false;
		is_a_preset = true;
	}
	else {
		is_a_plate = false;
		is_a_preset = false;
	}
	
	return (
		<
			PlateLookup 
			id={id} 
			is_a_plate={is_a_plate} 
			is_a_preset={is_a_preset} 
			/>
	)
	
}
