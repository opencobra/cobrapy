For installation help, please use the [Google Group]
(http://groups.google.com/group/cobra-pie).

For usage instructions, please see the 
[documentation](https://cobrapy.readthedocs.org/en/latest/)

--------------------------------------------------------------------------------

All releases require Python 2.7 to be installed before proceeding.
Mac OS X (10.7+) and Ubuntu ship with Python. Windows users without python 
can download and install python from the [python 
website](http://www.python.org/download/releases/2.7.7/).

Python 3 support is still considered experimental.

#Installation of cobrapy

## Stable version installation

### Mac OS X
0. [install pip] (http://pip.readthedocs.org/en/latest/installing.html).
1. In a terminal, run ```sudo pip install cobra --pre```

### GNU/Linux
0. [install pip] (http://pip.readthedocs.org/en/latest/installing.html).
1. Install the glpk library. On debian-based systems (including Ubuntu
   and Mint), this can be done with ```sudo apt-get install libglpk-dev```
3. In a terminal, run ```sudo pip install cobra --pre```

### Microsoft Windows
Download and install the appropriate 32 bit or 64 bit installer,
both of which can be downloaded from the [python package
index](https://pypi.python.org/pypi/cobra/).


## Hacking version installation
Use pip to install [Cython](http://cython.org/). Install libglpk using your
package manger. This would be ```brew install glpk``` on a Mac and
```sudo apt-get install libglpk-dev``` on debian-based systems (including
Ubuntu and Mint). This can also be installed by compiling GLPK from source.

Clone the git repository using your preferred mothod. Cloning from your
own [github fork](https://help.github.com/articles/fork-a-repo) is recommended!
Afterwards, open a terminal, enter the cobrapy repository and run the following
command:

    python setup.py develop --user

# Installation of optional dependencies
## Optional Dependencies
On windows, these can downloaded from [this site]
(http://www.lfd.uci.edu/~gohlke/pythonlibs/). On Mac/Linux, they can be
installed using pip, from binary installers, or from package managers.

1. [libsbml](http://sbml.org) >= 5.10 to read/write SBML files
  * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#libsbml)
  * Use ```sudo pip install python-libsbml-experimental``` on Mac/Linux
2. [numpy](http://numpy.org) >= 1.6.1 for double_deletion_analysis
  * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy)
3. [scipy](http://scipy.org) >= 0.11 for ArrayBasedModel and saving to *.mat files.
  * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy)

## Other solvers
cobrapy comes with bindings to the GNU Linear Programming Kit ([glpk]
(http://www.gnu.org/software/glpk/)) using its own bindings called "cglpk" in
cobrapy. In addition, cobrapy currently cobrapy supports these linear
programming solvers:

 * ILOG/CPLEX
  [Academic](https://www.ibm.com/developerworks/university/academicinitiative/)
  [Commercial](http://www.ibm.com/software/integration/optimization/cplex-optimizer/)
 * [gurobi](http://gurobi.com)
 * GLPK through [pyGLPK](http://tfinley.net/software/pyglpk/)

ILOG/CPLEX and Gurobi are commercial software packages that, currently, 
provide free licenses for academics and support both linear and quadratic 
programming. GLPK is an opensource linear programming solver; however, it 
does not support quadratic programming and is not as robust as the 
commercial solvers when it comes to mixed-integer linear programming.


# Testing your installation
1. Start python
2. Type the following into the Python shell

```python
from cobra.test import test_all
test_all()
```

