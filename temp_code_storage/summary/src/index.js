import React from 'react';
import ReactDOM from 'react-dom';
import './index.css';
import App from './App';
import reportWebVitals from './reportWebVitals';

const proposal = {
        "name": "test_proposal_3",
        "libraries": [
            {
                "id": 2,
                "name": "DSI_Poised_EG",
                "for_industry": true,
                "public": true
            },
            {
                "id": 3,
                "name": "FragLite",
                "for_industry": false,
                "public": true
            },
            {
                "id": 6,
                "name": "mylib_for_test_proposal_3",
                "for_industry": true,
                "public": false
            }
        ],
        "subsets": []
    };

const proposal_plates = [
    {
        "id": 3,
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate1",
        "current": true,
        "size": 308
    },
    {
        "id": 7,
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate2",
        "current": true,
        "size": 308
    },
    {
        "id": 9,
        "library": {
            "id": 3,
            "name": "FragLite",
            "for_industry": false,
            "public": true
        },
        "name": "test plate",
        "current": true,
        "size": 31
    },
    {
        "id": 6,
        "library": {
            "id": 6,
            "name": "mylib_for_test_proposal_3",
            "for_industry": true,
            "public": false
        },
        "name": "mylib_for_test_proposal_3",
        "current": true,
        "size": 19
    }
]
const current = [
    {
        "id": 3,
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate1",
        "current": true,
        "size": 308
    },
    {
        "id": 7,
        "library": {
            "id": 2,
            "name": "DSI_Poised_EG",
            "for_industry": true,
            "public": true
        },
        "name": "plate2",
        "current": true,
        "size": 308
    },
    {
        "id": 9,
        "library": {
            "id": 3,
            "name": "FragLite",
            "for_industry": false,
            "public": true
        },
        "name": "test plate",
        "current": true,
        "size": 31
    },
    {
        "id": 6,
        "library": {
            "id": 6,
            "name": "mylib_for_test_proposal_3",
            "for_industry": true,
            "public": false
        },
        "name": "mylib_for_test_proposal_3",
        "current": true,
        "size": 19
    }
]

const ownplates = [{
        "id": 6,
        "library": {
            "id": 6,
            "name": "mylib_for_test_proposal_3",
            "for_industry": true,
            "public": false
        },
        "name": "mylib_for_test_proposal_3",
        "size": 19
    }]

ReactDOM.render(
  <React.StrictMode>
    <App proposal={proposal} current={current} ownplates={ownplates} proposalPlates = {proposal_plates}/>
  </React.StrictMode>,
  document.getElementById('root')
);

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();
