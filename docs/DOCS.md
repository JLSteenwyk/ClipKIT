## Documentation
The documentation is built using Sphinx.

#### Option 1: pip/venv (recommended)
Use a plain virtual environment with the pinned docs requirements:

```shell
python -m venv .venv
source .venv/bin/activate
python -m pip install -r docs/requirements.txt
make -C docs html
```

#### Option 2: pipenv (compatibility path)
If you prefer [pipenv](https://pipenv.pypa.io/en/latest/), you can still use:

```shell
cd docs
pipenv install --dev
pipenv run make html
```

For live preview with pipenv:

```shell
cd docs
pipenv run serve
```
