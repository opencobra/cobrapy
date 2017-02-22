# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function

# set the warning format to be on a single line
import warnings as _warnings
from os.path import (abspath as _abspath, dirname as _dirname)
from os import name as _name

from cobra.core import (Object, Metabolite, Gene, Reaction, Model, DictList,
    Species)
from cobra import (io, flux_analysis, design)

__version__ = "0.5.11"

# set the warning format to be prettier and fit on one line
_cobra_path = _dirname(_abspath(__file__))
if _name == "posix":
    _warning_base = "%s:%s \x1b[1;31m%s\x1b[0m: %s\n"  # colors
else:
    _warning_base = "%s:%s %s: %s\n"


def _warn_format(message, category, filename, lineno, file=None, line=None):
    shortname = filename.replace(_cobra_path, "cobra", 1)
    return _warning_base % (shortname, lineno, category.__name__, message)
_warnings.formatwarning = _warn_format
