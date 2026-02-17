# ClipKIT Documentation

This directory contains the Sphinx source for the ClipKIT documentation.

For user-facing documentation, use these canonical sources:

- Project README: `../README.md`
- Main docs index: `index.rst`
- Advanced CLI docs: `advanced/index.rst`
- Changelog: `change_log/index.rst`

Quick start:

```bash
pip install clipkit
clipkit input.fa
```

Notes:

- CLI behavior and option defaults are defined by `clipkit/parser.py`.
- If documentation text and CLI help ever differ, treat CLI help as authoritative.
