# set the warning format to be on a single line
import warnings as _warnings
def _warn_format(message, category, filename, lineno, file=None, line=None):
    return "%s:%s %s: %s\n" % (filename, lineno, category.__name__, message)
_warnings.formatwarning = _warn_format

from .version import get_version
__version__ = get_version()
from .core import Object, Formula, Metabolite, Gene, Reaction, Model, DictList, Species
from . import io, flux_analysis

try:
    from .core import ArrayBasedModel
except ImportError:
    None

del get_version
