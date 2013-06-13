#cobra.core.Solution.py
##########################
#BEGIN Class Solution
#
from .Object import Object
class Solution(Object):
    """Stores the solution from optimizing a cobra.Model.  This is
    used to provide a single interface to results from different
    solvers that store their values in different ways.

    NOTE: This class might be deprecated in favor of associating the
    values with the Reactions and Metabolites in the cobra.Model.

    f: The objective value
    
    the_time: Float.  Sometimes indicates how long it took to solve a
    problem.  As this is typically negligible and not used in cobra pie,
    it might be deprecated.
    
    the_solver: A string indicating which solver package was used.

    x: List or Array of the values from the primal.

    x_dict: A dictionary of reaction ids that maps to the primal values.

    y: List or Array of the values from the dual.

    y_dict: A dictionary of reaction ids that maps to the dual values.
    
    """
    def __init__(self, the_f, x=None,
                 x_dict=None, y=None, y_dict=None,
                 the_solver=None, the_time=0, status='NA'):
        Object.__init__(self, the_f)
        self.solver = the_solver
        self.f = the_f
        self.x = x
        self.x_dict = x_dict
        self.status = status
        self.y = y
        self.y_dict = y_dict
    def dress_results(self, model):
        """Attaches results from FBA simulations to the Model's Reactions and
        Metabolites.

        model: The model that matches the Solution.

        """
        [setattr(k, 'x', v) for k, v in zip(model.reactions, self.x)];
        [setattr(k, 'y', v) for k, v in zip(model.metabolites, self.y)];
#
#END Class Solution
#########################
