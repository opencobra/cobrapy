import ast
from ast import (
    Add,
    BinOp,
    Div,
    Expression,
    Mod,
    Mult,
    Name,
    NodeTransformer,
    Num,
    Sub,
    UAdd,
    UnaryOp,
    USub,
    copy_location,
)
from typing import List, Union
from warnings import warn

from cobra.core import DictList
from cobra.core.object import Object


class UserDefinedConstraint(Object):
    """Class representation of constraints defined by
    user.

    The id attribute of UserDefinedConstraint is optional according
    to SBML specifications. But DictList requires the id of the object
    to be set. So temporary inner id's will be generated to store such
    objects.

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

    def __init__(
        self,
        id=None,
        name=None,
        lower_bound: [int, float] = None,
        upper_bound: [int, float] = None,
        const_comps: List = None,
    ):
        Object.__init__(self, id, name)
        self._model = None

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

        self._constraint_comps = DictList()
        self._const_comp_ids = set()
        if const_comps is not None:
            self.add_constraint_comps(const_comps)

    class ComputeNumericNodes(NodeTransformer):
        """ Compute the value of nodes which are solvable i.e
        node containing numeric value on both sides.
        """

        def visit_BinOp(self, node: BinOp):
            """This method visits every node of ast tree and solve the nodes
            which are computable i.e having numeric value on both sides and
            a defined BioOp (or UnaryOp)
            """
            node.left = self.visit(node.left)
            node.right = self.visit(node.right)
            if isinstance(node.left, Num) and isinstance(node.right, Num):
                if isinstance(node.op, Add):
                    result = Num(n=node.left.n + node.right.n)
                    return copy_location(result, node)
                elif isinstance(node.op, Sub):
                    result = Num(n=node.left.n - node.right.n)
                    return copy_location(result, node)
                elif isinstance(node.op, Mult):
                    result = Num(n=node.left.n * node.right.n)
                    return copy_location(result, node)
                elif isinstance(node.op, Div):
                    result = Num(n=node.left.n / node.right.n)
                    return copy_location(result, node)
                elif isinstance(node.op, Mod):
                    result = Num(n=node.left.n % node.right.n)
                    return copy_location(result, node)
            return node

    def add_single_node(
        self, ast_node: Union[Name, UnaryOp, BinOp], negate: bool = False
    ) -> None:
        """
        The final node to add constraint component inside constraint.

        Parameters
        ----------
        ast_node : AST Node classes
            Final node representation of a constraint component
        negate : bool
            Whether to add a negative sign before this component

        """
        # variable of type 'v1'
        if isinstance(ast_node, Name):
            coeff = 1
            if negate:
                coeff = -1
            self.add_constraint_comps(
                [
                    ConstraintComponent(
                        coefficient=coeff, variable=ast_node.id, variable_type="linear"
                    )
                ]
            )
        # variable of type '-v1'
        elif isinstance(ast_node, UnaryOp):
            if isinstance(ast_node.op, UAdd):
                coeff = 1
            elif isinstance(ast_node.op, USub):
                coeff = -1
            else:
                raise ValueError(f"Unsupported Unary Operand: {ast_node.op}")
            if negate:
                coeff = -1 * coeff
            self.add_constraint_comps(
                [
                    ConstraintComponent(
                        coefficient=coeff,
                        variable=ast_node.operand.id,
                        variable_type="linear",
                    )
                ]
            )
        # variable of type '2*v1', '3*v1*v1' etc.
        elif isinstance(ast_node, BinOp):
            if not isinstance(ast_node.op, Mult):
                raise ValueError(
                    f"Unsupported operand type between the variables:"
                    f"{ast_node.left} and {ast_node.right}"
                )
            if not isinstance(ast_node.right, Name):
                raise ValueError(
                    f"The second argument must be a single variable"
                    f": {ast_node.right}"
                )
            # variable of type '2*v1'
            if isinstance(ast_node.left, Num):
                coeff = ast_node.left.n
                if negate:
                    coeff = -1 * coeff
                var = ast_node.right.id
                self.add_constraint_comps(
                    [
                        ConstraintComponent(
                            coefficient=coeff, variable=var, variable_type="linear"
                        )
                    ]
                )
            # variables of type 'v1*v1'
            elif isinstance(ast_node.left, Name):
                if ast_node.left.id != ast_node.right.id:
                    raise ValueError(
                        f"Multiplication of two different variables is not "
                        f"allowed as per SBML FBC-V3: {ast_node.left.id} and "
                        f"{ast_node.right.id}"
                    )
                coeff = 1
                if negate:
                    coeff = -1
                self.add_constraint_comps(
                    [
                        ConstraintComponent(
                            coefficient=coeff,
                            variable=ast_node.right.id,
                            variable_type="quadratic",
                        )
                    ]
                )
            # variables of type '-v1*v1'
            elif isinstance(ast_node.left, UnaryOp):
                if isinstance(ast_node.left.operand, Num):
                    if isinstance(ast_node.left.op, UAdd):
                        coeff = 1
                    elif isinstance(ast_node.left.op, USub):
                        coeff = -1
                    else:
                        raise ValueError(
                            f"Unsupported Unary Operand: {ast_node.left.op}"
                        )
                    coeff *= ast_node.left.operand.n
                    if negate:
                        coeff = -1 * coeff
                    comp = ConstraintComponent(
                        coefficient=coeff,
                        variable=ast_node.right.id,
                        variable_type="linear",
                    )
                    self.add_constraint_comps([comp])
                elif isinstance(ast_node.left.operand, Name):
                    if ast_node.left.operand.id != ast_node.right.id:
                        raise ValueError(
                            f"Multiplication of two different variables is not allowed "
                            f"as per SBML FBC-V3: {ast_node.left.operand.id}"
                            f" and {ast_node.right.id}"
                        )
                    if isinstance(ast_node.left.op, UAdd):
                        coeff = 1
                    elif isinstance(ast_node.left.op, USub):
                        coeff = -1
                    else:
                        raise ValueError(
                            f"Unsupported Unary Operand: {ast_node.left.op}"
                        )
                    if negate:
                        coeff = -1 * coeff
                    comp = ConstraintComponent(
                        coefficient=coeff,
                        variable=ast_node.left.operand.id,
                        variable_type="quadratic",
                    )
                    self.add_constraint_comps([comp])
            # variables of type '2*v1*v1'
            elif isinstance(ast_node.left, BinOp):
                if ast_node.left.right.id != ast_node.right.id:
                    raise ValueError(
                        f"Multiplication of two different variables is not "
                        f"allowed as per SBML FBC-V3: {ast_node.left.name.id} "
                        f"and {ast_node.right.id}"
                    )
                # variables of type '-2*v1*v1'
                if isinstance(ast_node.left.left, UnaryOp):
                    if isinstance(ast_node.left.left.op, USub):
                        coeff = -1
                    elif isinstance(ast_node.left.left.op, UAdd):
                        coeff = 1
                    else:
                        raise ValueError("Invalid expression.")
                    coeff = coeff * ast_node.left.left.operand.n
                elif isinstance(ast_node.left.left, Num):
                    coeff = ast_node.left.left.n
                else:
                    raise ValueError("Invalid expression.")
                if negate:
                    coeff = -1 * coeff
                comp = ConstraintComponent(
                    coefficient=coeff,
                    variable=ast_node.right.id,
                    variable_type="quadratic",
                )
                self.add_constraint_comps([comp])

    def add_components_via_ast_nodes(
        self,
        ast_node: Union[Name, UnaryOp, BinOp, Expression, Num],
        negate: bool = False,
    ) -> None:
        """
        Add the constraint components to model via ast node
        representation of constraint expression.

        Parameters
        ----------
        ast_node : Expression
            The ast of the constraint expression
        negate : bool
            Whether to add a '-' sign before this node

        """

        if isinstance(ast_node, BinOp):
            if isinstance(ast_node.op, Add):
                self.add_components_via_ast_nodes(ast_node.left, negate)
                self.add_components_via_ast_nodes(ast_node.right, negate)
            elif isinstance(ast_node.op, Sub):
                self.add_components_via_ast_nodes(ast_node.left, negate)
                negate = not negate
                self.add_components_via_ast_nodes(ast_node.right, negate)
            elif isinstance(ast_node.op, Mult):
                self.add_single_node(ast_node, negate)
            else:
                raise ValueError(f"Unsupported operation between variables: {ast_node}")
        elif isinstance(ast_node, Name):
            self.add_single_node(ast_node, negate)
        elif isinstance(ast_node, UnaryOp):
            self.add_single_node(ast_node, negate)
        elif isinstance(ast_node, Expression):
            self.add_components_via_ast_nodes(ast_node.body, negate)
        else:
            raise ValueError(f"Unsupported variable type: {ast_node}")

    @staticmethod
    def constraint_from_expression(
        id: str = None,
        expression: str = "",
        lower_bound: Union[int, float] = None,
        upper_bound: Union[int, float] = None,
    ) -> "UserDefinedConstraint":
        """
        Method to add a user defined constraint via an expression of
        type string. Following rules must be followed while making the
        expression:

            1. The coefficient must be added before the variable and must
               have a multiplicative sign (*) only between itself and the
               variable. It can be an expression containing numbers:

                    2*variable : Allowed
                    2 * variable : Allowed
                    2.variable : Not Allowed
                    2.0*variable : Allowed
                    (4/2+5%2)*variable : Allowed
                    (2/4)*variable : Allowed
                    variable*2 : Not Allowed

            2. The 'caret' operator must not be used to set the exponent.
               Use multiplicative sign instead. Also, since only 'linear' and
               'quadratic' variables are allowed, we will need multiplication
               of variable at most one (in quadratic case).

                    variable : Allowed
                    variable^2 : Not allowed
                    variable * variable : Allowed

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
        constraint = UserDefinedConstraint(
            id=id, lower_bound=lower_bound, upper_bound=upper_bound
        )
        if expression is None or expression == "":
            return constraint

        expression = expression.strip()

        tree = ast.parse(source=expression, mode="eval")
        compute_nodes = UserDefinedConstraint.ComputeNumericNodes()
        tree = compute_nodes.visit(tree)
        print((ast.dump(tree)))
        constraint.add_components_via_ast_nodes(ast_node=tree)

        return constraint

    @property
    def lower_bound(self) -> Union[int, float]:
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, value: Union[int, float]) -> None:
        if value is None:
            self._lower_bound = value
        elif not (isinstance(value, int) or isinstance(value, float)):
            raise TypeError(
                f"The 'lower_bound' must be of type 'number' (int, float): {value}"
            )
        else:
            self._lower_bound = value

    @property
    def upper_bound(self) -> Union[int, float]:
        return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, value: Union[int, float]) -> None:
        if value is None:
            self._upper_bound = value
        elif not (isinstance(value, int) or isinstance(value, float)):
            raise TypeError(
                f"The 'upper_bound' must be of type 'number' (int, float): {value}"
            )
        else:
            self._upper_bound = value

    @property
    def constraint_comps(self) -> List:
        return self._constraint_comps

    def add_constraint_comps(self, value: List) -> None:
        """Adds a UserDefinedConstraintComponent in this constraint

        Parameters
        ----------
        value: list
            the constraint component to add in the model

        """
        if self._model is not None:
            raise ValueError(
                f"The constraint has already been added to model. Can't add more "
                f"constraint components: {value}"
            )

        if not isinstance(value, list):
            if isinstance(value, ConstraintComponent):
                warn(f"Pass the Constraint Components inside a list: {value}")
                value = [value]
            else:
                raise TypeError(
                    f"Pass the Constraint Components inside a list: {value}"
                )

        for item in value:
            if not isinstance(item, ConstraintComponent):
                raise TypeError(
                    f"The constraint component should be of type "
                    f"'UserDefinedConstraintComponents': {item}"
                )
            if item.id is None or item.id == "":
                item.id = "$_internal_comp_id" + str(len(self._const_comp_ids))
            self._const_comp_ids.add(item.id)
            self.constraint_comps.append(item)

    def remove_constraint_comps(self, value: "ConstraintComponent") -> None:
        """Removes a UserDefinedConstraintComponent from this constraint

        Parameters
        ----------
        value: ConstraintComponent
            the constraint component to br removed

        """
        if self._model is not None:
            raise ValueError(
                f"The constraint has already been added to model. Can't remove any "
                f"constraint components: {value}"
            )

        if not isinstance(value, list):
            if isinstance(value, ConstraintComponent):
                warn(f"Pass the Constraint Components inside a list: {value}")
                value = [value]
            else:
                raise TypeError(
                    f"Pass the Constraint Components inside a list: {value}"
                )

        for item in value:
            if not isinstance(item, ConstraintComponent):
                raise TypeError(
                    f"The constraint component should be of type "
                    f"'UserDefinedConstraintComponents': {item}"
                )
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

    variable_types = ("linear", "quadratic")

    def __init__(
        self,
        id: str = None,
        name: str = None,
        variable: str = None,
        coefficient: int = 1.0,
        variable_type: str = "linear",
    ):
        Object.__init__(self, id=id, name=name)
        self.variable = variable
        self.coefficient = coefficient
        self.variable_type = variable_type

    @property
    def variable(self) -> str:
        return self._variable

    @variable.setter
    def variable(self, value: str) -> None:
        if not isinstance(value, str):
            raise TypeError(
                f"The 'variable' have to be an COBRA object id and must be of"
                f" type string: {value}"
            )
        else:
            self._variable = value

    @property
    def coefficient(self) -> Union[int, float]:
        return self._coefficient

    @coefficient.setter
    def coefficient(self, value: Union[int, float]) -> None:
        if value is None:
            self._coefficient = 1
        elif not (isinstance(value, int) or isinstance(value, float)):
            raise TypeError(
                f"The 'coefficient' must be of type 'number' (int, float): {value}"
            )
        else:
            self._coefficient = value

    @property
    def variable_type(self) -> str:
        return self._variable_type

    @variable_type.setter
    def variable_type(self, value: str) -> None:
        if value not in self.variable_types:
            raise ValueError(
                f"The 'variable_type' must be one of 'linear' or 'quadratic': {value}"
            )
        else:
            self._variable_type = value
