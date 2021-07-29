import React from 'react';
import { useParams } from "react-router-dom";
import CompoundLookupPlate from '../plate_lookup/CompoundLookupPlate.js';
import CompoundLookupPreset from '../plate_lookup/CompoundLookupPreset.js';
import CompoundLookupCherryPick from '../plate_lookup/CompoundLookupCherryPick.js';

/* translates url parameters and chooses the right component to display compound data*/
export default function UrlTranslator() {
	let { type, id, project_id }  = useParams();

	if (type==="plate"){
		return (
			<CompoundLookupPlate id={id} project_id={project_id}/>
		);
	}
	else if (type=="preset"){
		return (
			<CompoundLookupPreset id={id} />
		);
	}
	else {
		return (
			<CompoundLookupCherryPick id={id} />
		);
	}
}
