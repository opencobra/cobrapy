# -*- coding: utf-8 -*-

from warnings import warn

from cobra.core import DictList
from cobra.core.object import Object


class UserDefinedConstraint(Object):
    """Class representation of constraints defined by
    user.

    Parameters
    ----------
    id : string
       An identifier for the chemical species
    name : string
       A human readable name.
    lower_bound : int, float
        lower bound on constraint expression
    upper_bound : int, float
        upper bound on constraint expression
    const_comps : list
        a list of constraint components
    """

    def __init__(self, id=None, name=None, lower_bound: [int, float]=None,
                 upper_bound: [int, float]=None, const_comps: list=None):
        if id is not None:
            if id.startswith('_'):
                warn("Use of '_' before publicly set id is "
                     "discouraged. Ids with an underscore before "
                     "them are for private use only.")
        Object.__init__(self, id, name)
        self._model = None

        self._lower_bound = None
        self._upper_bound = None
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

        self._constraint_comps = DictList()
        self._const_comp_ids = set()
        if const_comps is not None:
            self.add_constraint_comps(const_comps)

    @staticmethod
    def constraint_from_expression(id=None, expression: 'str' = '',
                                   lower_bound=None, upper_bound=None):
        """
        Method to add a user defined constraint via an expression of
        type string. Following rules must be followed while making the
        expression:

            1. The coefficient must be added before the variable without
               parentheses and must have a multiplicative sign (*) only
               between itself and the variable (without any whitespaces
               around). The coefficient itself can be an integer, float
               or a number in fractional form (a/b). Coefficient with any
               other form is not allowed. For instance:

                    variable : Allowed
                    1*variable : Allowed
                    2*variable : Allowed
                    2 * variable : Not allowed
                    2.variable : Not Allowed
                    2.0*variable : Allowed
                    2/3*variable : Allowed
                    2*2*variable : Not allowed
                    (2/4)*variable : Not allowed

            2. The exponent for the variable must be set via the 'caret'
               operator (^) only. It must also not have whitespaces around
               it. Also, since only linear and quadratic variables are all-
               wed, the exponent can be 1 or 2 only. For instance:

                    variable : Allowed
                    variable^1 : Allowed
                    variable^2 : Allowed
                    variable^3 : Not allowed
                    variable ^ 2 : Not allowed
                    variable*variable : Not allowed

            3. The only possible sign between two or more variables are '+'
               and '-' (addition and subtraction). These signs must have a
               single whitespace character on both of its side, with only
               exception that if the sign comes before very first variable,
               it may drop the whitespace before it. For instance:

                    variable1 - variable2 : Allowed
                    variable1+variable2 : Not allowed
                    - variable1 + variable2 : Allowed
                     + variable1 - variable2 : Allowed
                    variable1 * variable2 : Not allowed


        Parameters
        ----------
        id : str, optional
            identifier to attach with this constraint
        expression : str
            the expression defining the constraint
        lower_bound : [int, float]
            the lower bound on constraint expression
        upper_bound : [int, float]
            the lower bound on constraint expression

        Returns
        -------
        A UserDefinedConstraint object
        """
        constraint = UserDefinedConstraint(id=id, lower_bound=lower_bound,
                                           upper_bound=upper_bound)
        if expression is None or expression == '':
            return constraint

        expression = expression.strip()

        expression = expression.replace('+ ', '+')
        expression = expression.replace('- ', '-')

        variables = expression.split(' ')
        for variable in variables:
            ind = variable.find('*')
            if ind == -1:
                if variable[0] == '-':
                    coefficient = -1
                else:
                    coefficient = 1
            else:
                parts = variable.split('*')

                def convert_to_float(frac_str):
                    """
                    Function to convert a numeric string value
                    to float.
                    """
                    try:
                        return float(frac_str)
                    except ValueError:
                        num, denom = frac_str.split('/')
                        try:
                            leading, num = num.split(' ')
                            whole = float(leading)
                        except ValueError:
                            whole = 0
                        frac = float(num) / float(denom)
                        return whole - frac if whole < 0 else whole + frac

                coefficient = convert_to_float(parts[0])
                variable = parts[1]
            ind = variable.find('^')
            if ind == -1:
                exponent = '1'
            else:
                parts = variable.split('^')
                exponent = parts[1]
                variable = parts[0]
            variable = variable.replace('+', '')
            variable = variable.replace('-', '')
            if exponent == '1':
                variable_type = 'linear'
            elif exponent == '2':
                variable_type = 'quadratic'
            else:
                raise ValueError("Only 'linear' or 'quadratic' "
                                 "variables are allowed. Variable"
                                 " {} is raised to exponent {}"
                                 ".".format(variable, exponent))
            const_comp = ConstraintComponent(coefficient=coefficient,
                                             variable=variable,
                                             variable_type=variable_type)
            constraint.add_constraint_comps([const_comp])

        return constraint

    @property
    def lower_bound(self):
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, value):
        if value is None:
            self._lower_bound = value
        elif not (isinstance(value, int) or
                  isinstance(value, float)):
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
        if value is None:
            self._upper_bound = value
        elif not (isinstance(value, int) or
                  isinstance(value, float)):
            raise TypeError("The 'upper_bound' must be of "
                            "type 'number' (int, float):"
                            " {}".format(value))
        else:
            self._upper_bound = value

    @property
    def constraint_comps(self):
        return self._constraint_comps

    def add_constraint_comps(self, value):
        """Adds a UserDefinedConstraintComponent in this constraint

        Parameters
        ----------
        value: list
            the constraint component to add in the model

        """
        if self._model is not None:
            raise ValueError("The constraint has already been "
                             "added to model. Can't add more "
                             "constraint components".format(value))

        if not isinstance(value, list):
            if isinstance(value, ConstraintComponent):
                warn("Pass the Constraint Components inside a "
                     "list: {}".format(value))
                value = [value]
            else:
                raise TypeError("Pass the Constraint Components "
                                "inside a list: {}".format(value))

        for item in value:
            if not isinstance(item, ConstraintComponent):
                raise TypeError("The constraint component should be of "
                                "type 'UserDefinedConstraintComponents'"
                                ": {}".format(item))
            if item.id is None or item.id == "":
                item.id = "_internal_comp_id" + str(len(self._const_comp_ids))
            self._const_comp_ids.add(item.id)
            self.constraint_comps.append(item)

    def remove_constraint_comps(self, value):
        """Removes a UserDefinedConstraintComponent from this constraint

        Parameters
        ----------
        value: ConstraintComponent
            the constraint component to br removed

        """
        if self._model is not None:
            raise ValueError("The constraint has already been "
                             "added to model. Can't remove any "
                             "constraint components".format(value))

        if not isinstance(value, list):
            if isinstance(value, ConstraintComponent):
                warn("Pass the Constraint Components inside a "
                     "list: {}".format(value))
                value = [value]
            else:
                raise TypeError("Pass the Constraint Components "
                                "inside a list: {}".format(value))

        for item in value:
            if not isinstance(item, ConstraintComponent):
                raise TypeError("The constraint component should be of "
                                "type 'UserDefinedConstraintComponents'"
                                ": {}".format(item))
            self._const_comp_ids.remove(item.id)
            self.constraint_comps.remove(item)


class ConstraintComponent(Object):
    """Class representation of component of a user-defined
    constraint.

    Parameters
    ----------
    id: str
        An identifier for the chemical species
    name : string
        A human readable name.
    variable: str
        the id of variable referenced by this component
    coefficient: int, float
        coefficient of the variable in constraint expression
    variable_type: str
        exponent of the variable in constraint expression

    """

    variable_types = ('linear', 'quadratic')

    def __init__(self, id=None, name=None, variable=None,
                 coefficient=1.0, variable_type='linear'):
        if id is not None:
            if id.startswith('_'):
                warn("Use of '_' before publicly set id is "
                     "discouraged. Ids with an underscore before "
                     "them are for private use only.")
        Object.__init__(self, id=id, name=name)
        self._variable = None
        self._coefficient = None
        self._variable_type = None
        self.variable = variable
        self.coefficient = coefficient
        self.variable_type = variable_type

    @property
    def variable(self):
        return self._variable

    @variable.setter
    def variable(self, value):
        if not isinstance(value, str):
            raise TypeError("The 'variable' have to be an "
                            "COBRA object id and must be of"
                            " type string: {}".format(value))
        else:
            self._variable = value

    @property
    def coefficient(self):
        return self._coefficient

    @coefficient.setter
    def coefficient(self, value):
        if value is None:
            self._coefficient = 1
        elif not (isinstance(value, int) or
                  isinstance(value, float)):
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
