import hashlib
import tarfile
import io
import json
from pathlib import Path
import subprocess
import sys
from urllib.error import HTTPError, URLError
from zipfile import ZipFile

import pytest

pytest.importorskip("ruamel.yaml")

from scripts.bioconda import publish, sync


RECIPE = '''{% set name = "clipkit" %}
{% set version = "2.7.0" %}

package:
  name: "{{ name|lower }}"
  version: "{{ version }}"
source:
  url: "https://pypi.io/{{ name }}/{{ version }}.tar.gz"
  sha256: old
build:
  noarch: python
  number: 0
  entry_points:
    - clipkit = clipkit.cli:main
  script: "{{ PYTHON }} -m pip install . --no-deps"
  run_exports:
    - {{ pin_subpackage('clipkit', max_pin="x") }}
requirements:
  host:
    - python >=3.10
    - pip
    - setuptools
  run:
    - python >=3.10
test:
  imports:
    - clipkit
  commands:
    - clipkit --help
    - echo keep-this-test
about:
  license: MIT
extra:
  recipe-maintainers:
    - pauldg
'''
ENTRIES = {"clipkit": "clipkit.cli:main"}
SOURCE = {"url": "https://files.pythonhosted.org/clipkit.tar.gz",
          "filename": "clipkit.tar.gz", "digests": {"sha256": "a" * 64}}


def make_wheel(version="2.7.0", requirements="Requires-Dist: numpy>=1.26\nRequires-Dist: biopython>=1.82\n", python=">=3.10",
               entries="clipkit = clipkit.cli:main\n"):
    buffer = io.BytesIO()
    with ZipFile(buffer, "w") as archive:
        archive.writestr("clipkit.dist-info/METADATA", (
            f"Name: clipkit\nVersion: {version}\nRequires-Python: {python}\n{requirements}\n"
        ))
        archive.writestr("clipkit.dist-info/entry_points.txt", "[console_scripts]\n" + entries)
    return buffer.getvalue()


def update(text=RECIPE, version="2.7.0", test_changed=False):
    return sync.update_recipe(text, version, SOURCE, ENTRIES, "python >=3.11",
                              ["biopython >=1.82", "numpy >=1.26"], ["setuptools >=83", "wheel >=0.30"], test_changed)


@pytest.mark.parametrize("value", ["v2.7.0", "2.7.0"])
def test_stable_version(value):
    assert sync.release_version(value) == "2.7.0"


@pytest.mark.parametrize("value", ["2.7", "2.7.0rc1", "2.7.0.dev1", "2.7.0+local", "$(id)", "../2.7.0"])
def test_reject_invalid_version(value):
    with pytest.raises(ValueError, match="stable"):
        sync.release_version(value)


def test_metadata_uses_python_specific_release_pins():
    entries, python, deps = sync.wheel_metadata(make_wheel(requirements=(
        "Requires-Dist: biopython>=1.82; python_version < '3.13'\n"
        "Requires-Dist: biopython>=1.84; python_version >= '3.13'\n"
        "Requires-Dist: numpy>=1.26,<2.1; python_version == '3.12'\n"
        "Requires-Dist: numpy>=2.3.3; python_version >= '3.14'\n"
    )), "2.7.0")
    assert entries == ENTRIES
    assert python == "python 3.12.*"
    assert deps == ["biopython >=1.82", "numpy <2.1,>=1.26"]


def make_source():
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        pyproject = (b'[build-system]\nrequires = ["setuptools>=83", "wheel>=0.30"]\n'
                     b'build-backend = "setuptools.build_meta"\n')
        member = tarfile.TarInfo("clipkit-2.7.0/pyproject.toml")
        member.size = len(pyproject)
        archive.addfile(member, io.BytesIO(pyproject))
    return buffer.getvalue()


def test_source_supplies_build_requirements():
    assert sync.source_build_requirements(make_source(), "2.7.0") == [
        "setuptools >=83", "wheel >=0.30"
    ]


@pytest.mark.parametrize("requirement,match", [
    ("new-dependency>=1", "manual mapping"),
    ('numpy>=1; sys_platform == "linux"', "Unsupported dependency marker"),
    ("numpy[extra]>=1", "manual mapping"),
    ("numpy @ https://example.com/numpy.whl", "manual mapping"),
    ("numpy~=1.24", "Unsupported Conda"),
])
def test_reject_unmapped_requirements(requirement, match):
    with pytest.raises(ValueError, match=match):
        sync.wheel_metadata(make_wheel(requirements=f"Requires-Dist: {requirement}\n"), "2.7.0")


@pytest.mark.parametrize("kwargs,match", [
    ({"version": "2.6.0"}, "name/version"),
    ({"python": ""}, "must support recipe Python"),
    ({"requirements": ""}, "runtime dependencies"),
    ({"entries": "clipkit = os:system\n"}, "Unexpected ClipKIT console"),
    ({"entries": "clipkit = clipkit.cli:main\npk_x = os:system\n"}, "Unexpected ClipKIT console"),
])
def test_reject_invalid_wheel(kwargs, match):
    with pytest.raises(ValueError, match=match):
        sync.wheel_metadata(make_wheel(**kwargs), "2.7.0")


def test_verify_artifact_hash(monkeypatch):
    artifact = dict(SOURCE, digests={"sha256": hashlib.sha256(b"wheel").hexdigest()})
    monkeypatch.setattr(sync, "read_url", lambda url: b"wheel")
    assert sync.verified_download(artifact) == b"wheel"
    with pytest.raises(ValueError, match="Checksum mismatch"):
        sync.verified_download(SOURCE)


@pytest.mark.parametrize("url", ["http://files.pythonhosted.org/a", "https://example.org/a", "file:///tmp/a"])
def test_reject_artifact_origin(url):
    with pytest.raises(ValueError, match="HTTPS artifact"):
        sync.verified_download(dict(SOURCE, url=url))


def test_waits_for_pypi_and_both_artifacts(monkeypatch):
    wheel = {"filename": "clipkit-2.7.0-py3-none-any.whl", "packagetype": "bdist_wheel"}
    source = dict(SOURCE, packagetype="sdist")
    responses = iter([
        HTTPError("url", 404, "not found", {}, None),
        {"info": {"version": "2.7.0"}, "urls": [wheel]},
        {"info": {"version": "2.7.0"}, "urls": [wheel, source]},
    ])
    def read(url):
        response = next(responses)
        if isinstance(response, Exception):
            raise response
        return json.dumps(response).encode()
    monkeypatch.setattr(sync, "read_url", read)
    monkeypatch.setattr(sync.time, "sleep", lambda delay: None)
    assert sync.wait_for_release("2.7.0", attempts=3) == (wheel, source)


def test_wait_has_bounded_retries(monkeypatch):
    def read(url):
        raise URLError("offline")
    monkeypatch.setattr(sync, "read_url", read)
    with pytest.raises(ValueError, match="unavailable"):
        sync.wait_for_release("2.7.0", attempts=1)


def test_recipe_build_bump_preserves_jinja_and_unrelated_fields():
    result = update()
    _, parsed, _, _, version = sync.parse_recipe(result)
    assert version == "2.7.0"
    assert parsed["build"]["number"] == 1
    assert parsed["extra"]["recipe-maintainers"] == ["pauldg"]
    assert "{{ PYTHON }}" in result
    assert "{{ pin_subpackage('clipkit', max_pin=\"x\") }}" in result
    assert parsed["requirements"]["host"] == ["python >=3.11", "pip", "setuptools >=83", "wheel >=0.30"]
    assert "pk_old --help" not in parsed["test"]["commands"]
    assert "echo keep-this-test" in parsed["test"]["commands"]
    assert sync.TEST_NAME in parsed["test"]["files"]
    assert update(result) == result


def test_new_version_resets_build_number():
    result = update(update(), version="2.8.0")
    _, parsed, _, _, version = sync.parse_recipe(result)
    assert version == "2.8.0"
    assert parsed["build"]["number"] == 0


def test_test_only_change_bumps_build_once():
    result = update(update(), test_changed=True)
    assert sync.parse_recipe(result)[1]["build"]["number"] == 2
    assert update(result) == result


def test_refuse_downgrade_and_unsupported_template():
    with pytest.raises(ValueError, match="downgrade"):
        update(version="2.6.0")
    with pytest.raises(ValueError, match="Unsupported Jinja"):
        update(RECIPE + "{% if linux %}\n")


def test_detect_conflicting_pull_request(monkeypatch):
    pull = {"number": 123, "url": "https://github.com/example/123", "headRefName": "autobump"}
    def gh(*args):
        if args[0] == "pr":
            return [pull]
        return [[{"filename": "recipes/clipkit/meta.yaml"}]]
    monkeypatch.setattr(publish, "gh_json", gh)
    with pytest.raises(ValueError, match="Another ClipKIT recipe PR"):
        publish.check_pull_requests("JLSteenwyk/bioconda-recipes", "clipkit-release/2.7.0")


def test_reuse_own_pr_and_ignore_unrelated_files(monkeypatch):
    own = {"number": 1, "url": "own", "headRefName": "clipkit-release/2.7.0",
           "headRepositoryOwner": {"login": "JLSteenwyk"},
           "headRepository": {"name": "bioconda-recipes"}}
    other = {"number": 2, "url": "other", "headRefName": "other"}
    monkeypatch.setattr(publish, "gh_json", lambda *args: (
        [own, other] if args[0] == "pr" else [[{"filename": "recipes/not-clipkit/meta.yaml"}]]
    ))
    assert publish.check_pull_requests("JLSteenwyk/bioconda-recipes", "clipkit-release/2.7.0") == own


def test_publish_requires_explicit_credential(monkeypatch, tmp_path):
    monkeypatch.delenv("GH_TOKEN", raising=False)
    monkeypatch.setattr(sys, "argv", ["publish", "--artifact", str(tmp_path), "--checkout", str(tmp_path)])
    with pytest.raises(ValueError, match="BIOCONDA_GITHUB_TOKEN"):
        publish.main()


def test_cli_sync_and_second_run_are_idempotent(monkeypatch, tmp_path):
    recipe = tmp_path / "recipe"
    recipe.mkdir()
    (recipe / "meta.yaml").write_text(RECIPE)
    result = tmp_path / "result.json"
    monkeypatch.setattr(sync, "wait_for_release", lambda *args: ({"wheel": True}, SOURCE))
    monkeypatch.setattr(sync, "verified_download", lambda artifact: make_wheel() if "wheel" in artifact else make_source())
    monkeypatch.setattr(sys, "argv", ["sync", "--version", "v2.7.0", "--recipe", str(recipe), "--output", str(result)])
    sync.main()
    assert json.loads(result.read_text())["changed"] is True
    first = (recipe / "meta.yaml").read_text()
    sync.main()
    assert json.loads(result.read_text())["changed"] is False
    assert (recipe / "meta.yaml").read_text() == first
    assert (recipe / sync.TEST_NAME).exists()


def test_release_workflow_security_contract():
    from ruamel.yaml import YAML
    root = Path(__file__).resolve().parents[2]
    workflow = YAML(typ="safe").load((root / ".github/workflows/bioconda-release.yml").read_text())
    release = YAML(typ="safe").load((root / ".github/workflows/release.yml").read_text())
    assert workflow["on"]["workflow_call"]["inputs"]["version"]["required"]
    assert "publish-pypi" in release["jobs"]["sync-bioconda"]["needs"]
    assert workflow["on"]["workflow_dispatch"]["inputs"]["dry_run"]["default"] is True
    assert workflow["permissions"] == {"contents": "read"}
    assert workflow["concurrency"]["cancel-in-progress"] is False
    build = workflow["jobs"]["build"]
    assert "secrets." not in json.dumps(build)
    assert workflow["jobs"]["publish"]["needs"] == "build"
    for job in workflow["jobs"].values():
        for step in job["steps"]:
            if "uses" in step:
                assert len(step["uses"].split("@")[1]) == 40


def test_publisher_pushes_once_and_reuses_branch_without_force(monkeypatch, tmp_path):
    """Exercise real Git operations with local remotes; mock only GitHub APIs."""
    upstream = tmp_path / "upstream"
    fork = tmp_path / "fork.git"
    artifact = tmp_path / "artifact"
    artifact.mkdir()
    (artifact / "recipe").mkdir()
    subprocess.run(["git", "init", "-b", "master", str(upstream)], check=True, capture_output=True)
    def git(*args):
        return subprocess.check_output(["git", "-C", str(upstream), *args], text=True).strip()
    git("config", "user.name", "Test")
    git("config", "user.email", "test@example.com")
    recipe = upstream / "recipes" / "clipkit"
    recipe.mkdir(parents=True)
    (recipe / "meta.yaml").write_text(RECIPE)
    git("add", ".")
    git("commit", "-m", "Initial upstream")
    subprocess.run(["git", "clone", "--bare", str(upstream), str(fork)], check=True, capture_output=True)
    manifest = {"version": "2.7.0", "changed": True, "entry_points": 2,
                "source_sha256": "a" * 64, "recipe_tree": git("rev-parse", "HEAD:recipes/clipkit")}
    (artifact / "recipe" / "meta.yaml").write_text(update())
    (artifact / "recipe" / sync.TEST_NAME).write_text("print('package test')\n")
    manifest["recipe_sha256"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (artifact / "recipe").iterdir()
    }
    (artifact / "result.json").write_text(json.dumps(manifest))
    real_run = publish.run
    calls = []
    pulls = []
    def run(*args, cwd=None):
        calls.append(args)
        if args[0] == "gh":
            if args[1:3] == ("auth", "setup-git"):
                return ""
            if args[1:3] == ("pr", "create"):
                pulls.append({"number": 1, "url": "https://example.com/pr/1"})
                return pulls[0]["url"]
            if args[1:3] == ("pr", "edit"):
                return ""
            raise AssertionError(args)
        args = tuple(str(arg).replace(f"https://github.com/{publish.UPSTREAM}.git", str(upstream))
                     .replace("https://github.com/JLSteenwyk/bioconda-recipes.git", str(fork)) for arg in args)
        return real_run(*args, cwd=cwd)
    monkeypatch.setattr(publish, "run", run)
    monkeypatch.setattr(publish, "gh_json", lambda *args: {
        "fork": True, "parent": {"full_name": publish.UPSTREAM},
    })
    monkeypatch.setattr(publish, "check_pull_requests", lambda *args: pulls[0] if pulls else None)
    monkeypatch.setenv("GH_TOKEN", "test-token-not-used-for-network")
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    for index in range(2):
        monkeypatch.setattr(sys, "argv", ["publish", "--artifact", str(artifact),
                                         "--checkout", str(tmp_path / f"checkout-{index}")])
        publish.main()
    assert len(pulls) == 1
    assert sum(args[:2] == ("git", "commit") for args in calls) == 1
    assert not any("--force" in args or "--force-with-lease" in args for args in calls)
    assert any(args[:3] == ("gh", "pr", "edit") for args in calls)
    # A changed upstream recipe invalidates the previously built artifact.
    (recipe / "meta.yaml").write_text(RECIPE + "# upstream changed\n")
    git("add", ".")
    git("commit", "-m", "Change upstream recipe")
    monkeypatch.setattr(sys, "argv", ["publish", "--artifact", str(artifact),
                                     "--checkout", str(tmp_path / "stale-checkout")])
    with pytest.raises(ValueError, match="changed during the build"):
        publish.main()
