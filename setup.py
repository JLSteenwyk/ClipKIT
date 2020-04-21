from setuptools import setup
from Cython.Build import cythonize

REQUIRES = ["biopython"]

setup(
    name = "clipkit",
    packages = ["clipkit"],
    entry_points = {
        "console_scripts": ["clipkit = clipkit.clipkit:main"]
    },
    ext_modules=cythonize("clipkit/*.pyx", compiler_directives={'language_level' : "3"}),
    version = "0.0.1",
    install_requires=REQUIRES,
)