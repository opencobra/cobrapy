#cobra.core.Solution.py
##########################
#BEGIN Class Solution
#
from Object import Object
class Solution(Object):
    """Stores the solution from optimizing a cobra.Model
    """
    def __init__(self, the_f, the_x, the_stat=1,
                 the_solver=None, the_time=0, x_dict={},
                 status='NA'):
        Object.__init__(self, the_f)
        self.solver = the_solver
        self.f = the_f
        self.x = the_x
        self.x_dict = x_dict
        self.status = status
#
#END Class Solution
#########################
