#!/usr/bin/env python


from sys import version_info
from warnings import warn

from setuptools import setup


if version_info[:2] < (3, 6):
    warn(
        "We only explicitly test Python 3.6 and later, however, earlier versions may "
        "be functional if you can manage the dependencies yourself."
    )


if __name__ == "__main__":
    setup()
