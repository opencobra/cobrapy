from warnings import warn
warn("cobra.flux_analysis.objective is deprecated.")

#cobra.flux_analysis.objective.py
#functions for analyzing / creating objective functions
from ..core.Reaction import Reaction
from ..manipulation import initialize_growth_medium
import sys
if hasattr(sys, 'maxsize') and sys.maxsize > 2**32:
    try:
        from numpy import int64, int32
    except:
        int32 = int64 = int
else:
    int32 = int
    int64 = int


    
def assess_objective(model, objective=None,
                      objective_cutoff=0.001, growth_medium=None):
    """DEPRECATED
    """
    from warnings import warn
    warn("cobra.flux_analysis.objective.assess_objective is deprecated.  " +\
         "Please use cobra.flux_analysis.reaction.assess instead")
    return(False)
    
    
 

def update_objective(cobra_model, the_objectives):
    """Revised to take advantage of the new Reaction classes.

    cobra_model:  A cobra.Model

    the_objectives: A list or a dictionary.  If a list then
    a list of reactions for which the coefficient in the
    linear objective is set as 1.  If a dictionary then the
    key is the reaction and the value is the linear coefficient
    for the respective reaction.

    """
    return cobra_model.change_objective(the_objectives)

