# Release Automation

ClipKIT releases are automated with GitHub Actions via `.github/workflows/release.yml`.

## Supported release triggers

1. Push a version tag (recommended): `vX.Y.Z`
2. Manually run the `Release` workflow from GitHub Actions (optional `version` input)

## Behavior

The workflow will:

1. Validate that `clipkit/version.py` matches the requested/tagged version.
2. Build source and wheel distributions with `python -m build`.
3. Run `twine check` on built artifacts.
4. Create a GitHub Release and attach `dist/*` artifacts.

The workflow does not publish to PyPI. Publish the locally validated artifacts
with the manual procedure below before pushing the release tag.

## Manual PyPI release

PyPI versions and filenames are immutable. Before uploading, bump
`clipkit/version.py` and the changelog files to a version that does not already
exist on PyPI.

```shell
source venv/bin/activate
python -m pip install --upgrade build twine
rm -rf dist build *.egg-info
python -m build
twine check dist/*
twine upload dist/* -r pypi
```

Do not pass `--universal`; ClipKIT supports Python 3.10+ and should publish a
`py3-none-any` wheel.

## PyPI authentication

Use a project-scoped PyPI API token through Twine, either in the local
`pypi` entry in `~/.pypirc` or through Twine's environment variables. Never
commit a token or print it in release logs.
