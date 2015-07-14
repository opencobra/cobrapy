For installation help, please use the [Google Group]
(http://groups.google.com/group/cobra-pie).

For usage instructions, please see the 
[documentation](https://cobrapy.readthedocs.org/en/latest/)

--------------------------------------------------------------------------------

All releases require Python 2.7+ or 3.4+ to be installed before proceeding.
Mac OS X (10.7+) and Ubuntu ship with Python. Windows users without python 
can download and install python from the [python 
website](https://www.python.org/ftp/python/2.7.10/python-2.7.10.amd64.msi).
Please note that though Anaconda and other python distributions may work with
cobrapy, they are not explicitly supported at this time.


#Installation of cobrapy

## Stable version installation

### Mac OS X
0. [install pip](http://pip.readthedocs.org/en/latest/installing.html).
1. In a terminal, run ```sudo pip install cobra```

### Ubuntu Linux
0. [install pip](http://pip.readthedocs.org/en/latest/installing.html).
1. Install the python and glpk development libraries. On debian-based
   systems (including Ubuntu and Mint), this can be done with
   ```sudo apt-get install python-dev libglpk-dev```
3. In a terminal, run ```sudo pip install cobra```

### Microsoft Windows
The preferred installation method on Windows is also to use pip. The latest
Windows installers for Python 2.7 and 3.4 include pip, so if you use those you
will already have pip.

In order to use pip on Windows, you will need to use command prompt.
To open the command prompt, hit the start button and search for
```cmd```. Clicking on ```cmd.exe``` will launch the command prompt.
For more information about using the command line, many
[guides](http://dosprompt.info/basics.asp) are readily available.

0. Make sure you have pip installed already. To do this, run
   ```C:\Python27\Scripts\pip.exe --version```. On Python 3.4, you
   will need to adjust the path to ```C:\Python34\Scripts\pip.exe --version```
1. In a terminal, run ```C:\Python27\Scripts\pip.exe install cobra```
   (you will need to replace 27 with 34 if you are usuing Python 3.4).

To install without pip, you will need to download and use the appropriate
installer for your version of python from the [python package
index](https://pypi.python.org/pypi/cobra/).


## Hacking version installation
Use pip to install [Cython](http://cython.org/). Install libglpk 
using your package manger. This would be ```brew install glpk``` on a Mac
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

1. [libsbml](http://sbml.org) >= 5.10 to read/write older SBML files. This is
   not necessary for newer SBML level 3 files (such as those from
   http://bigg.ucsd.edu), which cobrapy can read/write natively.
  * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#libsbml)
  * Use ```sudo pip install python-libsbml-experimental``` on Mac/Linux
2. [numpy](http://numpy.org) >= 1.6.1 for double deletions
  * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy)
3. [scipy](http://scipy.org) >= 0.11 for ArrayBasedModel and saving to *.mat files.
  * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy)
4. [IPython notebook](http://ipython.org/) for a notebook interface. Cobrapy
   support tab-completion in the IPython notebook.
  * Use ```sudo pip install ipython[notebook]``` on Mac/Linux
  * Use ```C:\Python27\Scripts\pip.exe install ipython[notebook]``` on Windows with Python 2.7
  * Use ```C:\Python34\Scripts\pip.exe install ipython[notebook]``` on Windows with Python 3.4

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

