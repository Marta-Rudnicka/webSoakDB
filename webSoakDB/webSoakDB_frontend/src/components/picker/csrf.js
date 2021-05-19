//Based on: https://www.techiediaries.com/django-react-forms-csrf-axios/

import React from 'react';

function getCookie(name) {
  let cookieValue = null;
  if (document.cookie && document.cookie !== '') {
    let cookies = document.cookie.split(';');
    for (let i = 0; i < cookies.length; i++) {
      let cookie = cookies[i].trim();
      if (cookie.substring(0, name.length + 1) === (name + '=')) {
          cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
          break;
      }
    }
  }
  return cookieValue;
}



const CSRFToken = () => {
  let csrftoken = getCookie('csrftoken');
    
  return (
    <input type="hidden" name="csrfmiddlewaretoken" value={csrftoken} />
  );
};

export default CSRFToken;
