# -*- coding: utf-8 -*-

"""Test functions of util.py"""

from __future__ import absolute_import

from cobra.util import show_versions


def test_show_versions(capsys):
    show_versions()
    captured = capsys.readouterr()
    lines = captured.out.split("\n")
    assert lines[1].startswith("System Information")
    assert lines[2].startswith("==================")
    assert lines[3].startswith("OS")
    assert lines[4].startswith("OS-release")
    assert lines[5].startswith("Python")

    assert lines[7].startswith("Package Versions")
    assert lines[8].startswith("================")
