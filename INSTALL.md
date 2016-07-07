#Installation of cobrapy

For installation help, please use the
[Google Group](http://groups.google.com/group/cobra-pie).
For usage instructions, please see the 
[documentation](https://cobrapy.readthedocs.org/en/latest/).

--------------------------------------------------------------------------------

All releases require Python 2.7+ or 3.4+ to be installed before proceeding.
Mac OS X (10.7+) and Ubuntu ship with Python. Windows users without python 
can download and install python from the [python 
website](https://www.python.org/ftp/python/2.7.9/python-2.7.9.amd64.msi).
Please note that though Anaconda and other python distributions may work with
cobrapy, they are not explicitly supported (yet!).

## Stable version installation

cobrapy can be installed with any recent installation of pip. Instructions
for several operating systems are below:

### Mac OS X or Linux
0. [install pip](http://pip.readthedocs.org/en/latest/installing.html).
1. In a terminal, run ```sudo pip install cobra```

### Microsoft Windows
The preferred installation method on Windows is also to use pip. The latest
Windows installers for Python 2.7 and 3.4 include pip, so if you use those you
will already have pip.

1. In a terminal, run ```C:\Python27\Scripts\pip.exe install cobra```
   (you may need to adjust the path accordingly).

To install without pip, you will need to download and use the appropriate
installer for your version of python from the [python package
index](https://pypi.python.org/pypi/cobra/).


## Hacking version installation
Use pip to install [Cython](http://cython.org/). Install libglpk 
using your package manger. This would be
```brew install homebrew/science/glpk``` on a Mac
and ```sudo apt-get install libglpk-dev``` on debian-based systems
(including Ubuntu and Mint). GLPK can also be compiled from the
released source.

Clone the git repository using your preferred mothod. Cloning from your
own [github fork](https://help.github.com/articles/fork-a-repo) is recommended!
Afterwards, open a terminal, enter the cobrapy repository and run the following
command:

    python setup.py develop --user

# Installation of optional dependencies
## Optional Dependencies
On windows, these can downloaded from [this site]
(http://www.lfd.uci.edu/~gohlke/pythonlibs/). On Mac/Linux, they can be
installed using pip, or from the OS package manager (e.g brew, apt, yum).

1. [libsbml](http://sbml.org) >= 5.10 to read/write SBML level 2 files
    * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#libsbml)
    * Use ```sudo pip install python-libsbml``` on Mac/Linux
2. [lxml](http://lxml.de/) to speed up read/write of SBML level 3 files.
3. [numpy](http://numpy.org) >= 1.6.1 for double deletions
    * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy)
4. [scipy](http://scipy.org) >= 0.11 for ArrayBasedModel and saving to *.mat files.
    * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy)

## Other solvers
cobrapy comes with bindings to the GNU Linear Programming Kit ([glpk]
(http://www.gnu.org/software/glpk/)) using its own bindings called "cglpk" in
cobrapy. In addition, cobrapy currently supports these linear programming
solvers:

 * ILOG/CPLEX (available with
   [Academic](https://www.ibm.com/developerworks/university/academicinitiative/)
   and
   [Commercial](http://www.ibm.com/software/integration/optimization/cplex-optimizer/)
   licenses).
 * [gurobi](http://gurobi.com)
 * [QSopt_ex esolver](http://www.dii.uchile.cl/~daespino/ESolver_doc/main.html)
 * [MOSEK](http://www.mosek.com/)
 * [coin-or clp and cbc](http://coin-or.org/) through
   [cylp](https://github.com/coin-or/CyLP).

ILOG/CPLEX, MOSEK, and Gurobi are commercial software packages that currently
provide free licenses for academics and support both linear and quadratic
programming. GLPK and clp are open source linear programming solvers; however,
they may not be as robsut as the commercial solvers for mixed-integer and
quadratic programming. QSopt_ex esolver is also open source, and can solve
linear programs using rational operations, giving exact solutions.


# Testing your installation
1. Start python
2. Type the following into the Python shell

```python
from cobra.test import test_all
test_all()
```

You should see some skipped tests and expected failures, and the function should return ```False```.

