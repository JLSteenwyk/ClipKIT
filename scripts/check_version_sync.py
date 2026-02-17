#!/usr/bin/env python3
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def extract_package_version() -> str:
    text = read_text(ROOT / "clipkit" / "version.py")
    match = re.search(r'__version__\s*=\s*"([^"]+)"', text)
    if not match:
        raise ValueError("Could not parse __version__ from clipkit/version.py")
    return match.group(1)


def extract_changelog_version() -> str:
    text = read_text(ROOT / "CHANGELOG.md")
    match = re.search(r"^##\s+([0-9]+\.[0-9]+\.[0-9]+)\s*$", text, flags=re.MULTILINE)
    if not match:
        raise ValueError("Could not parse top version header from CHANGELOG.md")
    return match.group(1)


def extract_docs_changelog_version() -> str:
    text = read_text(ROOT / "docs" / "change_log" / "index.rst")
    match = re.search(r"^\*\*([0-9]+\.[0-9]+\.[0-9]+)\*\*\s*$", text, flags=re.MULTILINE)
    if not match:
        raise ValueError(
            "Could not parse top version header from docs/change_log/index.rst"
        )
    return match.group(1)


def main() -> int:
    package_version = extract_package_version()
    changelog_version = extract_changelog_version()
    docs_changelog_version = extract_docs_changelog_version()

    mismatches = []
    if package_version != changelog_version:
        mismatches.append(
            f"clipkit/version.py ({package_version}) != CHANGELOG.md ({changelog_version})"
        )
    if package_version != docs_changelog_version:
        mismatches.append(
            "clipkit/version.py "
            f"({package_version}) != docs/change_log/index.rst ({docs_changelog_version})"
        )

    if mismatches:
        print("Version consistency check failed:")
        for mismatch in mismatches:
            print(f"- {mismatch}")
        return 1

    print(f"Version consistency check passed: {package_version}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
