from setuptools import setup

REQUIRES = ["biopython"]

setup(
    name = "clipkit",
    packages = ["clipkit"],
    entry_points = {
        "console_scripts": ["clipkit = clipkit.clipkit:main"]
    },
    version = "0.0.1",
    install_requires=REQUIRES,
)