#!/usr/bin/env python

"""
Tracks the version number. If git is installed and file script
is located within a git repository, git describe is used to get
the version information. This version string is sanitized to
comply with PEP 386 and stored in the RELEASE-VERSION file.

If git describe can not be run, the RELEASE-VERSION file is used
for version information instead.

"""

__all__ = ("get_git_version")

from subprocess import check_output
from os import path
current_dir = path.dirname(path.abspath(__file__))
version_file = path.join(current_dir, "RELEASE-VERSION")

def call_git_describe(abbrev=4):
    try:
        return check_output(["git", "describe",  "--tags",
            "--abbrev=%d" % abbrev], dir=current_dir).strip()
    except:
        return None


def read_release_version():
    try:
        with open(version_file, "r") as infile:
            version = infile.read().strip()
        if len(version) == 0:
            version = None
        return version
    except:
        return None


def write_release_version(version):
    with open(version_file, "w") as outfile:
        outfile.write("%s\n" % version)


def get_git_version(abbrev=4):
    # Read in the version that's currently in RELEASE-VERSION.

    release_version = read_release_version()

    # First try to get the current version using "git describe".

    version = call_git_describe(abbrev)

    #adapt to PEP 386 compatible versioning scheme
    version = pep386adapt(version)

    # If that doesn't work, fall back on the value that's in
    # RELEASE-VERSION.

    if version is None:
        version = release_version

    # If we still don't have anything, that's an error.

    if version is None:
        raise ValueError("Cannot find the version number!")

    # If the current version is different from what's in the
    # RELEASE-VERSION file, update the file to be current.

    if version != release_version:
        write_release_version(version)

    # Finally, return the current version.

    return version


def pep386adapt(version):
    if version is None:
        return
    if '-' in version:
        # adapt git-describe version to be in line with PEP 386
        parts = version.split('-')
        parts[-2] = 'post'+parts[-2]
        version = '.'.join(parts[:-1])
    return version


if __name__ == "__main__":
    print get_git_version()
