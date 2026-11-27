#!/usr/bin/env python3
"""Synchronize a Bioconda recipe from verified, published PyPI artifacts.

Does not import ClipKIT or execute code from its release artifacts.
"""

from __future__ import annotations

import argparse
import configparser
from email.parser import BytesParser
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile
import time
import tomllib
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import urlopen
from zipfile import ZipFile

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name
from packaging.version import Version
from ruamel.yaml import YAML


DEPENDENCIES = {"biopython", "numpy"}
BUILD_DEPENDENCIES = {"setuptools", "wheel"}
# A noarch recipe has one dependency set. ClipKIT's NumPy ranges vary by Python minor.
CONDA_PYTHON = "3.12"
TEST_NAME = "clipkit_release_test.py"
VERSION_LINE = re.compile(r'(?m)^{% set version = ["\']([^"\']+)["\'] %}$')


def release_version(value):
    value = value.removeprefix("v")
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", value):
        raise ValueError("Expected a stable X.Y.Z version or vX.Y.Z tag")
    return str(Version(value))


def read_url(url):
    with urlopen(url, timeout=60) as response:
        data = response.read(50 * 1024 * 1024 + 1)
    if len(data) > 50 * 1024 * 1024:
        raise ValueError(f"Download exceeds 50 MiB: {url}")
    return data


def wait_for_release(version, attempts=30, delay=30):
    if attempts < 1 or delay < 0:
        raise ValueError("Attempts must be positive and delay nonnegative")
    for attempt in range(attempts):
        try:
            release = json.loads(read_url(f"https://pypi.org/pypi/clipkit/{version}/json"))
            if release["info"]["version"] != version:
                raise ValueError("PyPI returned a different version")
            files = [f for f in release["urls"] if not f.get("yanked")]
            wheels = [f for f in files if f["filename"].endswith("-py3-none-any.whl")]
            sources = [f for f in files if f["packagetype"] == "sdist"]
            if len(wheels) == len(sources) == 1:
                return wheels[0], sources[0]
        except HTTPError as error:
            if error.code not in (404, 429, 500, 502, 503, 504):
                raise
        except URLError:
            pass
        if attempt + 1 < attempts:
            print(f"Waiting for both PyPI artifacts for {version} ({attempt + 1}/{attempts})", flush=True)
            time.sleep(delay)
    raise ValueError(f"PyPI release {version} is unavailable, incomplete, or yanked")


def verified_download(artifact):
    parsed = urlparse(artifact["url"])
    if parsed.scheme != "https" or parsed.netloc != "files.pythonhosted.org":
        raise ValueError("Expected a files.pythonhosted.org HTTPS artifact")
    digest = artifact["digests"]["sha256"]
    if not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise ValueError("Invalid SHA256 digest")
    data = read_url(artifact["url"])
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"Checksum mismatch: {artifact['filename']}")
    return data


def conda_specifiers(value):
    result = []
    for specifier in sorted(SpecifierSet(value), key=str):
        if specifier.operator not in ("==", "!=", ">=", "<=", ">", "<"):
            raise ValueError(f"Unsupported Conda version constraint: {specifier}")
        if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)*(?:\.\*)?", specifier.version):
            raise ValueError(f"Unsupported Conda version: {specifier.version}")
        result.append(str(specifier))
    return ",".join(result)


def wheel_metadata(data, version):
    with ZipFile(io.BytesIO(data)) as archive:
        paths = [p for p in archive.namelist() if p.endswith(".dist-info/METADATA")]
        if len(paths) != 1:
            raise ValueError("Expected one wheel METADATA file")
        metadata = BytesParser().parsebytes(archive.read(paths[0]))
        if canonicalize_name(metadata["Name"]) != "clipkit" or metadata["Version"] != version:
            raise ValueError("Wheel name/version does not match the requested release")
        config = configparser.ConfigParser(interpolation=None)
        config.optionxform = str
        config.read_string(archive.read(paths[0].replace("METADATA", "entry_points.txt")).decode())
    entries = dict(config["console_scripts"])
    if entries != {"clipkit": "clipkit.cli:main"}:
        raise ValueError("Unexpected ClipKIT console scripts")
    python = metadata.get("Requires-Python")
    if not python or CONDA_PYTHON not in SpecifierSet(python):
        raise ValueError(f"Wheel must support recipe Python {CONDA_PYTHON}")
    dependencies = mapped_requirements(metadata.get_all("Requires-Dist", []), DEPENDENCIES)
    if {requirement.split()[0] for requirement in dependencies} != DEPENDENCIES:
        raise ValueError("Missing runtime dependencies for recipe Python")
    return entries, f"python {CONDA_PYTHON}.*", dependencies


def mapped_requirements(values, names):
    dependencies = []
    for value in values:
        requirement = Requirement(value)
        name = canonicalize_name(requirement.name)
        if name not in names or requirement.extras or requirement.url:
            raise ValueError(f"Dependency needs manual mapping: {value}")
        if requirement.marker:
            marker = str(requirement.marker)
            if not re.fullmatch(r'python_version (?:==|!=|<=|>=|<|>) "[0-9]+\.[0-9]+"', marker):
                raise ValueError(f"Unsupported dependency marker: {value}")
            if not requirement.marker.evaluate({"python_version": CONDA_PYTHON}):
                continue
        constraint = conda_specifiers(str(requirement.specifier))
        dependencies.append(f"{name} {constraint}".strip())
    return sorted(set(dependencies))


def source_build_requirements(data, version):
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        member = archive.getmember(f"clipkit-{version}/pyproject.toml")
        if not member.isfile() or member.size > 65536:
            raise ValueError("Invalid pyproject.toml in the source archive")
        with archive.extractfile(member) as handle:
            build = tomllib.load(handle)["build-system"]
    if build.get("build-backend") != "setuptools.build_meta" or build.get("backend-path"):
        raise ValueError("Unexpected source build backend")
    return mapped_requirements(build["requires"], BUILD_DEPENDENCIES)


def parse_recipe(text):
    """Mask Jinja tokens before round-trip YAML parsing; never execute templates."""
    versions = VERSION_LINE.findall(text)
    if len(versions) != 1:
        raise ValueError("Expected one literal Jinja version assignment")
    tokens = {}

    def mask(match):
        token = f"__CLIPKIT_JINJA_{len(tokens)}__"
        tokens[token] = match.group()
        return token

    headers = re.findall(r"(?m)^{%.*%}\n?", text)
    if any(not re.fullmatch(r'{% set (?:name|version) = ["\'][^"\']+["\'] %}\n?', h) for h in headers):
        raise ValueError("Unsupported Jinja statement in recipe")
    body = re.sub(r"(?m)^{%.*%}\n?", "", text).lstrip("\n")
    body = re.sub(r"{{.*?}}", mask, body)
    if "{%" in body or "{{" in body or "{#" in body:
        raise ValueError("Unsupported Jinja template in recipe")
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.width = 120
    yaml.indent(mapping=2, sequence=4, offset=2)
    return yaml, yaml.load(body), headers, tokens, versions[0]


def update_recipe(text, version, source, entries, python, dependencies, build_dependencies, test_changed):
    yaml, recipe, headers, tokens, previous = parse_recipe(text)
    if Version(version) < Version(previous):
        raise ValueError(f"Refusing to downgrade recipe {previous} to {version}")
    if recipe["build"].get("noarch") != "python":
        raise ValueError("Expected a noarch: python recipe")
    if not isinstance(recipe["source"], dict):
        raise ValueError("Multiple recipe sources require manual review")
    # Preserve the upstream order for existing commands; append new ones deterministically.
    old_entries = recipe["build"].get("entry_points", [])
    names = [value.split(" = ")[0] for value in old_entries]
    order = list(dict.fromkeys([n for n in names if n in entries] + sorted(entries)))
    desired_entries = [f"{name} = {entries[name]}" for name in order]
    run = sorted(dependencies + [python])
    # The upstream build installs with --no-build-isolation.
    host = [python, "pip", *build_dependencies]
    test = recipe.setdefault("test", {})
    files = list(test.get("files", []))
    if TEST_NAME not in files:
        files.append(TEST_NAME)
    commands = [c for c in test.get("commands", []) if not (
        c.endswith(" --help") and c[:-7] in names and c[:-7] not in entries
    )]
    if f"python {TEST_NAME}" not in commands:
        commands.append(f"python {TEST_NAME}")
    changed = any((
        previous != version,
        recipe["source"].get("url") != source["url"],
        recipe["source"].get("sha256") != source["digests"]["sha256"],
        old_entries != desired_entries,
        sorted(recipe["requirements"]["run"]) != run,
        list(recipe["requirements"]["host"]) != host,
        list(test.get("files", [])) != files,
        list(test.get("commands", [])) != commands,
        test_changed,
    ))
    if not changed:
        return text
    # PyPI's version-templated shortcut can return 404 for a published sdist.
    recipe["source"]["url"] = source["url"]
    recipe["source"]["sha256"] = source["digests"]["sha256"]
    build = recipe["build"]
    build["number"] = 0 if previous != version else int(build["number"]) + 1
    build["entry_points"] = desired_entries
    recipe["requirements"]["host"] = host
    recipe["requirements"]["run"] = run
    test["files"] = files
    test["commands"] = commands
    output = io.StringIO()
    yaml.dump(recipe, output)
    body = output.getvalue()
    for token, value in tokens.items():
        body = body.replace(token, value)
    header = "".join(headers)
    header = VERSION_LINE.sub(lambda m: m.group().replace(previous, version), header)
    return header.rstrip() + "\n\n" + body


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--attempts", type=int, default=30)
    parser.add_argument("--delay", type=float, default=30)
    parser.add_argument("--output", type=Path, required=True, help="JSON result manifest")
    args = parser.parse_args()
    version = release_version(args.version)
    wheel, source = wait_for_release(version, args.attempts, args.delay)
    entries, python, dependencies = wheel_metadata(verified_download(wheel), version)
    build_dependencies = source_build_requirements(verified_download(source), version)
    path = args.recipe / "meta.yaml"
    before = path.read_text()
    test_path = args.recipe / TEST_NAME
    test = Path(__file__).with_name("package_test.py").read_text()
    test_changed = not test_path.exists() or test_path.read_text() != test
    after = update_recipe(before, version, source, entries, python, dependencies, build_dependencies, test_changed)
    path.write_text(after)
    test_path.write_text(test)
    result = {"version": version, "changed": before != after or test_changed,
              "entry_points": len(entries), "source_sha256": source["digests"]["sha256"]}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
