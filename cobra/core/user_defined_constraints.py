
"""
Some documentation:

Reactions: v1; v2; v3;
Stochiometric Contraints (metabolic nework, SBML Reactions):
v1 - v2 = 0

UserDefinedConstraints:
v1 + v2 < 2

Non-constant parameters (x1: not reactions; )

expression =: v1 + v2 + 1.2 x1 < 2
v1 + v2 + 1.2 x1*x1 < 2

-> examples:
-> optimize: (v1, v2, x1); how to access the optimal solution?

Test Cases:
------------
1. add user constraint to simple model
A -v1-> B -v2-> C
optimize: v2; bounds v1[0; 10]; v2[0, 5]
optimal value: 5;

UserConstraint: v1 + v2 <= 4
optimal value: 2;

2. only user constraint;
empty cobra model

x1 <= 2
maximize: x1
flux distribution optimal: x1=2
"""

from warnings import warn
from cobra import DictList, Object


class UserDefinedConstraints(Object):
    """Class representation of constraints defined by
    user.

    Parameters
    ----------
    id : string
       An identifier for the chemical species
    name : string
       A human readable name.
    lower_bound : 'int', 'float'
        lower bound on constraint expression
    upper_bound : 'int', 'float'
        upper bound on constraint expression
    const_comp : list
        a list of constraint components
    """

    def __init__(self, id=None, name=None, lower_bound: ['int', 'float']=None,
                 upper_bound: ['int', 'float']=None, const_comps: 'list'=None):
        Object.__init__(self, id, name)
        self._model = None

        self._lower_bound = None
        self._upper_bound = None
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

        self._constraint_comps = DictList()
        self.add_constraint_comps(const_comps)

    @property
    def lower_bound(self):
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, value):
        if not isinstance(value, int) or \
                not isinstance(value, float):
            raise TypeError("The 'lower_bound' must be of "
                            "type 'number' (int, float):"
                            " {}".format(value))
        else:
            self._lower_bound = value

    @property
    def upper_bound(self):
        return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, value):
        if not isinstance(value, int) or \
                not isinstance(value, float):
            raise TypeError("The 'upper_bound' must be of "
                            "type 'number' (int, float):"
                            " {}".format(value))
        else:
            self._upper_bound = value

    @property
    def constraint_comps(self):
        return self._constraint_comps

    def add_constraint_comps(self, value):
        # COBRApy doesn't support parameters. So constraints
        # present on parameter will not be added
        if self._model is not None:
            raise ValueError("The constraint has already been "
                             "added to model. Can't add more "
                             "constraint components".format(value))

        if not isinstance(value, list):
            if isinstance(value, UserDefinedConstraintComponents):
                warn("Pass the Constraint Components inside a "
                     "list: {}".format(value))
                value = [value]
            else:
                raise TypeError("Pass the Constraint Components "
                                "inside a list: {}".format(value))

        for item in value:
            if not isinstance(item, UserDefinedConstraintComponents):
                raise TypeError("The constraint component should be of "
                                "type 'UserDefinedConstraintComponents'"
                                ": {}".format(item))

            # FIXME: what to do with non-constant parameters
            if item.ref_var not in self._model.reactions:
                warn("The referenced variable is not present "
                     "inside this model's reaction list: "
                     "{}".format(item.ref_var))
            self.constraint_comps.append(item)

    def remove_constraint_comps(self, value):
        if self._model is not None:
            raise ValueError("The constraint has already been "
                             "added to model. Can't remove any "
                             "constraint components".format(value))

        if not isinstance(value, list):
            if isinstance(value, UserDefinedConstraintComponents):
                warn("Pass the Constraint Components inside a "
                     "list: {}".format(value))
                value = [value]
            else:
                raise TypeError("Pass the Constraint Components "
                                "inside a list: {}".format(value))

        for item in value:
            if not isinstance(item, UserDefinedConstraintComponents):
                raise TypeError("The constraint component should be of "
                                "type 'UserDefinedConstraintComponents'"
                                ": {}".format(item))
            self.constraint_comps.remove(item)


class UserDefinedConstraintComponents(Object):

    variable_types = ('linear', 'quadratic')

    def __init__(self, id=None, name=None, ref_var=None,
                 coefficient=None, variable_type=None):
        Object.__init__(self, id, name)
        self._ref_var = None
        self._coefficient = None
        self._variable_type = None
        self.ref_var = ref_var
        self.coefficient = coefficient
        self.variable_type = variable_type

    @property
    def ref_var(self):
        return self._ref_var

    @ref_var.setter
    def ref_var(self, value):
        if not isinstance(value, str):
            raise TypeError("The 'ref_var' have to be an "
                            "COBRA object id and must be of"
                            " type string: {}".format(value))
        else:
            self._ref_var = value

    @property
    def coefficient(self):
        return self._coefficient

    @coefficient.setter
    def coefficient(self, value):
        if not isinstance(value, int) or \
                not isinstance(value, float):
            raise TypeError("The 'coefficient' must be of "
                            "type 'number' (int, float):"
                            " {}".format(value))
        else:
            self._coefficient = value

    @property
    def variable_type(self):
        return self._variable_type

    @variable_type.setter
    def variable_type(self, value):
        if value not in self.variable_types:
            raise ValueError("The 'variable_type' must be one"
                             "of 'linear' or 'quadratic': "
                             "{}".format(value))
        else:
            self._variable_type = value
