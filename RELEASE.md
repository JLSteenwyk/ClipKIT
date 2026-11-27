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
5. For version tags, publish those same validated artifacts to PyPI in a separate job.

The PyPI job uses the `PYPI_API_TOKEN` Actions secret. To retry a PyPI
failure, rerun the failed job or dispatch the workflow against the existing
release tag. Dispatches from a branch build the GitHub release only.

## Manual PyPI release

Use this fallback if the automated job cannot publish. PyPI versions and filenames are immutable. Before uploading, bump
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

Do not pass `--universal`; ClipKIT supports Python 3.12+ and should publish a
`py3-none-any` wheel.

## PyPI authentication

Use a project-scoped PyPI API token. Store it as the `PYPI_API_TOKEN`
Actions secret for automated releases, or use the local `pypi` entry in
`~/.pypirc` or Twine environment variables for manual recovery. Never commit
a token or print it in release logs.

## Bioconda

After a version tag publishes to PyPI, the release workflow calls **Sync
Bioconda release**. It checks both PyPI artifact hashes, generates the recipe
from wheel and source metadata, and builds and tests the package. A separate
job uses `BIOCONDA_GITHUB_TOKEN` to open or update a PR in
`bioconda/bioconda-recipes`. Bioconda reviews the PR and publishes the package.

Add the repository Actions secret `BIOCONDA_GITHUB_TOKEN` with a dedicated
credential that can push to the Bioconda recipes fork and open and edit
upstream PRs. The default fork is `JLSteenwyk/bioconda-recipes`. Set the
repository variable `BIOCONDA_FORK` to override it. Use GitHub settings or
`gh secret set BIOCONDA_GITHUB_TOKEN --repo JLSteenwyk/ClipKIT` and enter the
credential interactively.

To validate a published release without credentials or upstream writes:

```sh
gh workflow run bioconda-release.yml --repo JLSteenwyk/ClipKIT \
  --ref main -f version=2.16.2 -f dry_run=true
```

Manual runs default to dry run. If a maintainer or autobump PR already changes
the ClipKIT recipe, the publish job reports its URL and stops before pushing.
The generated recipe remains downloadable from the run, including on Conda build
failure. A recipe that changes upstream during the build requires a rerun.

The generated noarch recipe starts with Python 3.12 and evaluates the wheel's
Python-version markers for that minor. This preserves ClipKIT's NumPy bounds in
a recipe that has one runtime dependency set. PyPI support is independent of
that Conda pin. Review the pin and mappings in `scripts/bioconda/sync.py` when
changing ClipKIT's Python or dependency policy.
