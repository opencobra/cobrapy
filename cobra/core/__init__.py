from .DictList import DictList
from .Object import Object
from .Gene import Gene
from .Metabolite import Metabolite
from .Reaction import Reaction
from .Solution import Solution, LazySolution
from .Model import Model
from .Species import Species

try:
    import scipy
except ImportError:
    scipy = None

if scipy:
    from .ArrayBasedModel import ArrayBasedModel
else:
    from warnings import warn
    warn("ArrayBasedModel requires scipy")
    del warn
del scipy
