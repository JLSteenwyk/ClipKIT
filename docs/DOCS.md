## Documentation
The documentation is built using Sphinx. Assuming you have Python 3 and [pipenv](https://pipenv.readthedocs.io/en/latest/install/#installing-pipenv)
installed, run the following:

```shell
cd docs
pipenv install --dev
```

#### Running the Development Server
You can run a local development server to preview changes using the following:

```shell
cd docs
pipenv run serve
```

#### Running Manually
You can build the docs manually by running:
```shell
cd docs
pipenv run make html
```