
## Quadratic Programming

# This example is available as an IPython [notebook](http://nbviewer.ipython.or
# g/github/opencobra/cobrapy/blob/master/documentation_builder/qp.ipynb).

# Suppose we want to minimize the Euclidean distance of the solution to the
# origin while subject to linear constraints. This will require a quadratic
# objective function.
# 
# > **min** $\frac{1}{2}\left(x^2 + y^2 \right)$
# 
# > *subject to*
# 
# > $x + y = 2$
# 
# > $x \ge 0$
# 
# > $y \ge 0$
# 
# The objective can be rewritten as $\frac{1}{2} v^T \cdot \mathbf Q \cdot v$,
# where
# $v = \left(\begin{matrix} x \\ y\end{matrix} \right)$ and
# $\mathbf Q = \left(\begin{matrix} 1 & 0\\ 0 & 1 \end{matrix}\right)$
# 
# The matrix $\mathbf Q$ can be passed into a cobra model as the quadratic
# objective.

import scipy

from cobra import Reaction, Metabolite, Model, solvers


# The quadratic objective $\mathbf Q$ should be formatted as a scipy sparse
# matrix.

Q = scipy.sparse.eye(2).todok()
Q
# Output:
# <2x2 sparse matrix of type '<type 'numpy.float64'>'
# 	with 2 stored elements in Dictionary Of Keys format>

# In this case, the quadratic objective is simply the identity matrix

Q.todense()
# Output:
# matrix([[ 1.,  0.],
#         [ 0.,  1.]])

# We need to use a solver that supports quadratic programming, such as gurobi
# or cplex. If a solver which supports quadratic programming is installed, this
# function will return its name.

print(solvers.get_solver_name(qp=True))
# Prints:
# gurobi

c = Metabolite("c")
c._bound = 2
x = Reaction("x")
y = Reaction("y")
x.add_metabolites({c: 1})
y.add_metabolites({c: 1})
m = Model()
m.add_reactions([x, y])
sol = m.optimize(quadratic_component=Q, objective_sense="minimize")
sol.x_dict
# Output:
# {'x': 1.0, 'y': 1.0}

# Suppose we change the problem to have a mixed linear and quadratic objective.

# > **min** $\frac{1}{2}\left(x^2 + y^2 \right) - y$
# 
# > *subject to*
# 
# > $x + y = 2$
# 
# > $x \ge 0$
# 
# > $y \ge 0$

y.objective_coefficient = -1
sol = m.optimize(quadratic_component=Q, objective_sense="minimize")
sol.x_dict
# Output:
# {'x': 0.5, 'y': 1.5}
