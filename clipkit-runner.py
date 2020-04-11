#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Convenience wrapper for running clipkit directly from source tree."""
import sys

from clipkit.clipkit import main

if __name__ == '__main__':
    main(sys.argv[1:])
