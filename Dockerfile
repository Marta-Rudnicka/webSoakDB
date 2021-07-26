#Start with RDkit
FROM informaticsmatters/rdkit-python3-debian:Release_2021_03_2

# Install node prereqs, nodejs and yarn
# Ref: https://deb.nodesource.com/setup_14.x
# Ref: https://yarnpkg.com/en/docs/install

USER root

RUN apt-get update && apt-get install -y gnupg2

RUN \
  echo "deb https://deb.nodesource.com/node_14.x buster main" > /etc/apt/sources.list.d/nodesource.list && \
  wget -qO- https://deb.nodesource.com/gpgkey/nodesource.gpg.key | apt-key add - && \
  echo "deb https://dl.yarnpkg.com/debian/ stable main" > /etc/apt/sources.list.d/yarn.list && \
  wget -qO- https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - && \
  apt-get update && \
  apt-get install -yqq nodejs yarn && \
  pip3 install -U pip && pip3 install pipenv \
#  npm i npm@6 && \
  curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python && ln -s /root/.poetry/bin/poetry /usr/local/bin && \
  rm -rf /var/lib/apt/lists/*

# Install Python Stuff
RUN pip3 install django djangorestframework django-rest-knox psycopg2-binary seaborn python-slugify

# Move Stuff to a new place
WORKDIR /usr/app
COPY ./ /usr/app
RUN cd /usr/app/

# Configure npm and build site
RUN npm init -y
RUN npm i -D webpack webpack-cli @babel/core babel-loader @babel/preset-env @babel/preset-react babel-plugin-transform-class-properties
RUN npm i react react-dom prop-types axios react-dom react-redux redux redux-devtools-extension redux-thunk remote-redux-devtools
RUN npm run dev

# Expose port
EXPOSE 8000

# Run on Entry
CMD cd /usr/app/webSoakDB ; python3 ./manage.py migrate ; python3 ./manage.py runserver 0.0.0.0:8000
