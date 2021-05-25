# webSoakDB

webSoakDB is going to be a part of the larger XChemSPA application connected to XChem DB. Currently, webSoakDB is using a dummy SQLite database with some test data (based on real inventory data from the past, so can be used to simulate realistic use cases).

Detailed documentation explaining what the application does and how it works is available in the `documentation` folder.

## Building the Stack for Devlopment Purposes (for local deployment)

The application uses RDKit, and the recommended way to use RDKit is to create a conda environment that includes this library.

Installing Anaconda: https://docs.anaconda.com/anaconda/install/
Installing RDKit: https://www.rdkit.org/docs/Install.html

After creating an environment with Python3, Node and RDKIt
```bash
conda activate <env_name>

pip install django djangorestframework django-rest-knox bokeh python-slugify 

npm init -y
npm i -D webpack webpack-cli @babel/core babel-loader @babel/preset-env @babel/preset-react babel-plugin-transform-class-properties
npm i react react-dom prop-types axios react-dom react-redux redux redux-devtools-extension redux-thunk remote-redux-devtools react-lazy-load-image-component react-router-dom 
npm run dev

cd webSoakDB
python ./manage.py migrate
python ./manage.py runserver
```
TODO: verify the build, update requirements.txt

## Building Container and Deploying to Kubernetes

[TODO] Current Dockerfile needs updating

AUTHOR: Marta Rudnicka

LICENSE: MIT