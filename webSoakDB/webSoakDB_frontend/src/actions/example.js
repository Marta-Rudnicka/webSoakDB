import axios from 'axios';
import { GET_EXAMPLE } from './types';

// Get Example
export const getExample = () => dispatch => {
    axios
        .get('/api/example/')
        .then(res =>
        {
            dispatch({
                type: GET_EXAMPLE,
                payload: res.data
            });
        })
        .catch(err => console.log(err));
};