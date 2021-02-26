# webSoakDB

AUTHOR: Marta Rudnicka

LICENSE: MIT

## Building the Stack for Devlopment Purposes (for local deployment)

Install python3 and node, both should be available on DLS (I hope).

We recommend creating conda environment to handle the pythony stuff. After cookiecutting the template, cd into this directory (where `ls` yields webSoakDB, basically where the Dockerfile is aswell.

```bash
# Can use previous environment from cookiecutter
conda create --name <env_name> python=3.8
conda activate <env_name>

pip install django djangorestframework django-rest-knox

npm init -y
npm i -D webpack webpack-cli @babel/core babel-loader @babel/preset-env @babel/preset-react babel-plugin-transform-class-properties
npm i react react-dom prop-types axios react-dom react-redux redux redux-devtools-extension redux-thunk remote-redux-devtools
npm run dev

cd webSoakDB
python ./manage.py migrate
python ./manage.py runserver
```

## Building Container and Deploying to Kubernetes

We have included a Framework Dockerfile and Kubernetes yaml file that should work out of box with minor configuration.
If you require additional packages and tools please edit the Docker file as you go to save problems later on.

The kubernetes.yaml file may require some additional tweaking and you will require a namespace etc to ensure stuff can run.

For diamond users you will be allowed to run deployments from your own namespace (your FEDID UID) but you will need to register a DNS domain from Chris to be able to see it!

Additionally, we provide two scripts that will rebuild the Dockerfile and commit to the Diamond GCloud docker registry and another one which will automatically disperse and redeploy kubenetes pods.

Anyway:
To rebuild the Docker image use:

```bash
./rebuild_docker_image.sh
```

To apply the kubenetes config:

```bash
module load pollux
kubectl apply -f kubernetes.yaml
```

To redeploy pods

```bash
./redeploy_pods.sh
```
