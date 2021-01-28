from os import path
from setuptools import setup

from clipkit.version import __version__

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

CLASSIFIERS = [
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Topic :: Scientific/Engineering',
]

REQUIRES = ["biopython==1.76", "numpy==1.18.2", "tqdm==4.45.0"]

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
)

## push new version to pypi
# rm -rf dist
# python setup.py sdist bdist_wheel --universal
# twine upload dist/* -r pypi