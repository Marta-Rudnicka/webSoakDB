//import React from 'react';
import ReactDOM from 'react-dom';
//import App from './App';

console.log('Loaded index.js')

const platedata = {
	"library": {
		"name": "Test library",
		"public": true
	},
	"name": "Test plate",
	"current": true
};

const compounds = 
[
    {
        "well": "A1",
        "compound": {
            "code": "NCL-00023819",
            "smiles": "BrC1=CNN=C1",
            "molecular_weight": 145.947960196
        },
        "concentration": 100
    },
    {
        "well": "A2",
        "compound": {
            "code": "NCL-00023818",
            "smiles": "IC1=CNN=C1",
            "molecular_weight": 193.934096096
        },
        "concentration": 100
    },
    {
        "well": "A3",
        "compound": {
            "code": "NCL-00023820",
            "smiles": "BrC1=CON=C1",
            "molecular_weight": 146.931975784
        },
        "concentration": 100
    },
    {
        "well": "A4",
        "compound": {
            "code": "NCL-00023823",
            "smiles": "NC1=NC=CC(Br)=C1",
            "molecular_weight": 171.96361026
        },
        "concentration": 100
    },
    {
        "well": "A5",
        "compound": {
            "code": "NCL-00023822",
            "smiles": "NC1=NC=CC(I)=C1",
            "molecular_weight": 219.94974616
        },
        "concentration": 100
    },
    {
        "well": "A6",
        "compound": {
            "code": "NCL-00023825",
            "smiles": "O=C1C=C(Br)C=CN1",
            "molecular_weight": 172.947625848
        },
        "concentration": 100
    },
    {
        "well": "A7",
        "compound": {
            "code": "NCL-00023824",
            "smiles": "O=C1C=C(I)C=CN1",
            "molecular_weight": 220.933761748
        },
        "concentration": 100
    },
    {
        "well": "A8",
        "compound": {
            "code": "NCL-00023830",
            "smiles": "BrC1=CC=C(S(N)(=O)=O)C=C1",
            "molecular_weight": 234.930261532
        },
        "concentration": 100
    },
    {
        "well": "A9",
        "compound": {
            "code": "NCL-00023829",
            "smiles": "IC1=CC=C(S(N)(=O)=O)C=C1",
            "molecular_weight": 282.916397432
        },
        "concentration": 100
    },
    {
        "well": "A10",
        "compound": {
            "code": "NCL-00023827",
            "smiles": "O=C1CC2=C(C=C(Br)C=C2)N1",
            "molecular_weight": 210.963275912
        },
        "concentration": 100
    },
    {
        "well": "A11",
        "compound": {
            "code": "NCL-00023828",
            "smiles": "BrC1=CC=C(NC(C)=O)C=C1",
            "molecular_weight": 212.978925976
        },
        "concentration": 100
    },
    {
        "well": "A12",
        "compound": {
            "code": "NCL-00023826",
            "smiles": "IC1=CC=C(C(N)=O)C=C1",
            "molecular_weight": 246.949411812
        },
        "concentration": 100
    },
    {
        "well": "A13",
        "compound": {
            "code": "NCL-00023832",
            "smiles": "BrC1=CN=CN=C1",
            "molecular_weight": 157.947960196
        },
        "concentration": 100
    },
    {
        "well": "A14",
        "compound": {
            "code": "NCL-00023831",
            "smiles": "IC1=CN=CN=C1",
            "molecular_weight": 205.934096096
        },
        "concentration": 100
    },
    {
        "well": "A15",
        "compound": {
            "code": "NCL-00023836",
            "smiles": "BrC1=CC(OC)=NC=C1",
            "molecular_weight": 186.963275912
        },
        "concentration": 100
    },
    {
        "well": "A16",
        "compound": {
            "code": "NCL-00023833",
            "smiles": "BrC1=CC=NC2=NC=CC=C21",
            "molecular_weight": 207.96361026
        },
        "concentration": 100
    },
    {
        "well": "A17",
        "compound": {
            "code": "NCL-00023835",
            "smiles": "BrC1=CC=C(S(C)(=O)=O)C=C1",
            "molecular_weight": 233.935012564
        },
        "concentration": 100
    },
    {
        "well": "A18",
        "compound": {
            "code": "NCL-00024670",
            "smiles": "BrC1=CC(CO)=NC=C1",
            "molecular_weight": 186.963275912
        },
        "concentration": 100
    },
    {
        "well": "A19",
        "compound": {
            "code": "NCL-00024774",
            "smiles": "OC1=C(OC)C=C(Br)C=C1",
            "molecular_weight": 201.962941564
        },
        "concentration": 100
    },
    {
        "well": "A20",
        "compound": {
            "code": "NCL-00024674",
            "smiles": "BrC1=CC(COC)=NC=C1",
            "molecular_weight": 200.978925976
        },
        "concentration": 100
    },
    {
        "well": "A21",
        "compound": {
            "code": "NCL-00024661",
            "smiles": "OC1=C(C#N)C=C(Br)C=C1",
            "molecular_weight": 196.947625848
        },
        "concentration": 100
    },
    {
        "well": "A22",
        "compound": {
            "code": "NCL-00024662",
            "smiles": "BrC1=CC(OC)=C(CO)C=C1",
            "molecular_weight": 215.978591628
        },
        "concentration": 100
    },
    {
        "well": "A23",
        "compound": {
            "code": "NCL-00024671",
            "smiles": "BrC1=CN(CCO)N=C1",
            "molecular_weight": 189.974174944
        },
        "concentration": 100
    },
    {
        "well": "A24",
        "compound": {
            "code": "NCL-00024667",
            "smiles": "BrC1=CC(O)=C(C(O)=O)C=C1",
            "molecular_weight": 215.94220612
        },
        "concentration": 100
    },
    {
        "well": "B1",
        "compound": {
            "code": "NCL-00024663",
            "smiles": "BrC1=CC(C#N)=C(OC)C=C1",
            "molecular_weight": 210.963275912
        },
        "concentration": 100
    },
    {
        "well": "B2",
        "compound": {
            "code": "NCL-00024673",
            "smiles": "BrC1=CN(CCOC)N=C1",
            "molecular_weight": 203.989825008
        },
        "concentration": 100
    },
    {
        "well": "B3",
        "compound": {
            "code": "NCL-00024890",
            "smiles": "BrC1=CN(CC(O)=O)N=C1",
            "molecular_weight": 203.9534395
        },
        "concentration": 100
    },
    {
        "well": "B4",
        "compound": {
            "code": "NCL-00024672",
            "smiles": "O=C1N(CCO)C=CC(Br)=C1",
            "molecular_weight": 216.973840596
        },
        "concentration": 100
    },
    {
        "well": "B5",
        "compound": {
            "code": "NCL-00024387",
            "smiles": "O=C1N(CCOC)C=CC(Br)=C1",
            "molecular_weight": 230.98949066
        },
        "concentration": 100
    },
    {
        "well": "B6",
        "compound": {
            "code": "NCL-00024773",
            "smiles": "O=C1N(CC(O)=O)C=CC(Br)=C1",
            "molecular_weight": 230.953105152
        },
        "concentration": 100
    },
    {
        "well": "B7",
        "compound": {
            "code": "NCL-00024665",
            "smiles": "BrC1=CC(OC)=C(CC(O)=O)C=C1",
            "molecular_weight": 243.973506248
        },
        "concentration": 100
    }
]
/*
ReactDOM.render(
  <React.StrictMode>
    <App plate = {platedata} compounds = {compounds} />
  </React.StrictMode>,
  document.getElementById('root')
);
*/
//const domContainer = document.querySelector('#root');
//ReactDOM.render(e(App), domContainer);
