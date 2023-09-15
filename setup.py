#!/usr/bin/env python


"""Set up the cobra package."""


from sys import version_info
from warnings import warn

from setuptools import setup


if version_info[:2] < (3, 6):
    warn(
        "We only explicitly test Python 3.6 and later, however, earlier versions may "
        "be functional if you can manage the dependencies yourself."
    )


# All other arguments are defined in `setup.cfg`.
if __name__ == "__main__":
    setup(version="0.27.0")
