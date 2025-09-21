from os import path
from setuptools import setup

from clipkit.version import __version__

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

CLASSIFIERS = [
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Topic :: Scientific/Engineering",
]

REQUIRES = [
    "biopython>=1.82; python_version < '3.13'",
    "biopython>=1.84; python_version >= '3.13'",
    "numpy>=1.24.0,<2.0; python_version < '3.12'",
    "numpy>=1.26.0,<2.1; python_version == '3.12'",
    "numpy>=2.0.0,<2.2; python_version >= '3.13'"
]

setup(
    name="clipkit",
    description="Alignment trimming software for phylogenetics.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/clipkit",
    packages=["clipkit"],
    classifiers=CLASSIFIERS,
    entry_points={"console_scripts": ["clipkit = clipkit.clipkit:main"]},
    version=__version__,
    include_package_data=True,
    install_requires=REQUIRES,
    python_requires=">=3.9",
)

## push new version to pypi
# rm -rf dist
# python3 setup.py sdist bdist_wheel --universal
# twine upload dist/* -r pypi
