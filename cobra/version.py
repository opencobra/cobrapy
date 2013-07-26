#!/usr/bin/env python

__all__ = ("get_version")

from subprocess import check_output
from os import path, name
current_dir = path.dirname(path.abspath(__file__))
version_file = path.join(current_dir, "RELEASE-VERSION")

git_command = "git"
if name == "nt":
    # TODO try to find better
    None

def call_git_describe(abbrev=4):
    try:
        return check_output(["git", "describe",  "--tags",
            "--abbrev=%d" % abbrev], cwd=current_dir).strip()
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


def get_version():
    """Tracks the version number.

    The file RELEASE-VERSION will contain the version of the last release. If
    this is not a git repository, it is safe to assume that the version is not
    being  incremented and the version returned will be the release version as
    read from the file.

    However, if the script is located within an active git repository,
    git-describe is used to get the version information.

    The file RELEASE-VERSION will need to be changed by manually. This only
    needs to occur twice per release:
      - Once right before running git tag (set to the same as the version in
        the tag).
      - Once right after running git tag (set to next_version.dev)
    """

    git_version = call_git_describe()
    if git_version is None:  # not a git repository
        return read_release_version()
    return git_version


if __name__ == "__main__":
    print get_version()
