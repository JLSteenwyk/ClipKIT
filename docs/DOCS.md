## Documentation
The documentation is built using Sphinx.

#### Option 1: pipenv (recommended)
Assuming you have Python 3 and [pipenv](https://pipenv.readthedocs.io/en/latest/install/#installing-pipenv)
installed, run:

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

#### Option 2: pip/venv
If you prefer a plain virtual environment instead of pipenv:

```shell
python -m venv .venv
source .venv/bin/activate
python -m pip install -r docs/requirements.txt
make -C docs html
```
