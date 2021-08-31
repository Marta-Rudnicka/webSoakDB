# webSoakDB

This branch is meant for integrating the app with ispyb_dja authorization package.

webSoakDB is going to be a part of the larger XChemSPA application connected to XChem DB. Currently, webSoakDB is using a dummy SQLite database with some test data (based on real inventory data from the past, so can be used to simulate realistic use cases).

Detailed documentation aimed at developers, which explains what the application does and why, is in the `documentation` folder.

## Building the Stack for Devlopment Purposes (for local deployment)

The application uses RDKit, and the recommended way to use RDKit is to create a conda environment that includes this library.

[Installing Anaconda](https://docs.anaconda.com/anaconda/install/)
[Installing RDKit](https://www.rdkit.org/docs/Install.html)

After creating an environment with Python3, Node and RDKIt
```bash
conda activate <env_name>

pip install django djangorestframework django-rest-knox bokeh python-slugify ispyb_dja django_filters

npm init -y
npm i -D webpack webpack-cli @babel/core babel-loader @babel/preset-env @babel/preset-react babel-plugin-transform-class-properties
npm i react react-dom prop-types axios react-dom react-redux redux redux-devtools-extension redux-thunk remote-redux-devtools react-lazy-load-image-component react-router-dom 
npm run dev

cd webSoakDB
```

Then, open webSoakDB_stack/settings.py and add secret ISPyB environment variables to the file.
Then:

```
python ./manage.py migrate
python ./manage.py runserver

```

IMPORTANT! To access the locally running application through the browser, use the address: `localhost:8000`, not `127.0.0.1:8000` - otherwise redirection to CAS logging page will not work.

AUTHOR: Marta Rudnicka

LICENSE: MIT
