# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
import sys
from contextlib import contextmanager

from six import StringIO


@contextmanager
def captured_output():
    """A context manager to test the IO summary methods."""
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr

    finally:
        sys.stdout, sys.stderr = old_out, old_err


def check_line(output, expected_entries,
               pattern=re.compile(r"\s")):
    """Ensure each expected entry is in the output."""
    output_set = set(
        pattern.sub("", line) for line in output.splitlines())
    for elem in expected_entries:
        assert pattern.sub("", elem) in output_set


def check_in_line(output, expected_entries,
                  pattern=re.compile(r"\s")):
    """Ensure each expected entry is contained in the output."""
    output_strip = [pattern.sub("", line) for line in
                    output.splitlines()]
    for elem in expected_entries:
        assert any(
            pattern.sub("", elem) in line for line in output_strip), \
            "Not found: {} in:\n{}".format(pattern.sub("", elem),
                                           "\n".join(output_strip))
