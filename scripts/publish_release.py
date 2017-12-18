#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Insert a TOML header into the latest release note."""

from __future__ import absolute_import, print_function

import sys
from datetime import date
from glob import glob
from builtins import open
from os.path import join, basename
from shutil import copy


def insert_break(lines, break_pos=9):
    """
    Insert a <!--more--> tag for larger release notes.

    Parameters
    ----------
    lines : list of str
        The content of the release note.
    break_pos : int
        Line number before which a break should approximately be inserted.

    Returns
    -------
    list of str
        The text with the inserted tag or no modification if it was
        sufficiently short.
    """
    def line_filter(line):
        if len(line) == 0:
            return True
        return any(line.startswith(c) for c in "-*+")

    if len(lines) <= break_pos:
        return lines
    newlines = [
        i for i, line in enumerate(lines[break_pos:], start=break_pos)
        if line_filter(line.strip())]
    if len(newlines) > 0:
        break_pos = newlines[0]
    lines.insert(break_pos, "<!--more-->\n")
    return lines


def build_hugo_md(filename, tag, bump):
    """
    Build the markdown release notes for Hugo.

    Inserts the required TOML header with specific values and adds a break
    for long release notes.

    Parameters
    ----------
    filename : str, path
        The release notes file.
    tag : str
        The tag, following semantic versioning, of the current release.
    bump : {"major", "minor", "patch", "alpha", "beta"}
        The type of release.

    """
    header = [
        '+++\n',
        'date = "{}"\n'.format(date.today().isoformat()),
        'title = "{}"\n'.format(tag),
        'author = "The COBRApy Team"\n',
        'release = "{}"\n'.format(bump),
        '+++\n',
        '\n'
    ]
    with open(filename, "r") as file_h:
        content = insert_break(file_h.readlines())
    header.extend(content)
    with open(filename, "w") as file_h:
        file_h.writelines(header)


def intify(filename):
    """
    Turn a release note filename into something sortable.

    Parameters
    ----------
    filename : str
        A release note of expected filename format '<major>.<minor>.<patch>.md'.

    Returns
    -------
    tuple
        A pair of the major and minor versions as integers.

    """
    tmp = filename[:-3].split(".")
    return int(tmp[0]), int(tmp[1])


def find_bump(target, tag):
    """Identify the kind of release by comparing to existing ones."""
    tmp = tag.split(".")
    existing = [intify(basename(f)) for f in glob(join(target, "[0-9]*.md"))]
    latest = max(existing)
    if int(tmp[0]) > latest[0]:
        return "major"
    elif int(tmp[1]) > latest[1]:
        return "minor"
    else:
        return "patch"


def main(argv):
    """
    Identify the release type and create a new target file with TOML header.

    Requires three arguments.

    """
    source, target, tag = argv
    if "a" in tag:
        bump = "alpha"
    if "b" in tag:
        bump = "beta"
    else:
        bump = find_bump(target, tag)
    filename = "{}.md".format(tag)
    destination = copy(join(source, filename), target)
    build_hugo_md(destination, tag, bump)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage:\n{} <source dir> <target dir> <tag>"
              "".format(sys.argv[0]))
        sys.exit(2)
    sys.exit(main(sys.argv[1:]))
