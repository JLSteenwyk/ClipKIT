"""clipkit.__main__: executed when clipkit is called as script"""
import sys

from .cli import main

main(sys.argv[1:])
