from os import path
from setuptools import setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

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
    entry_points={"console_scripts": ["clipkit = clipkit.clipkit:main"]},
    version="0.1.0",
    include_package_data=True,
    install_requires=REQUIRES,
)
