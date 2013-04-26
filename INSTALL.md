Installation instructions for Python and Jython are detailed below. The 
Python instructions will be approrpriate for most users. For installation 
help, please see the [sourceforge help install 
forum](https://sourceforge.net/projects/opencobra/forums/forum/4194275)

For usage instructions, please see the 
[documentation](https://cobrapy.readthedocs.org/en/latest/)

--------------------------------------------------------------------------------

# Python
All releases require Python 2.7 to be installed before proceeding. 
Mac OS X (10.7+) and Ubuntu ship with Python. Windows users without python 
can download and install python from the [python 
website](http://www.python.org/download/releases/2.7.3/).

Generally, installation should follow these steps:

1. Install cobrapy (either the stable, development, or hacking version).
2. Install an appropriate solver
3. Install optional dependencies
4. Test your installation


## Stable version installation in Python
On Mac OS X or GNU/Linux, use easy_install (or pip if you have it) by running
the following command in a terminal.

    sudo easy_install cobrapy

For Windows, download and install the appropriate 32 bit or 64 bit installer,
both of which can be downloaded from the [python package
index](https://pypi.python.org/pypi/cobra/).

## Development version installation in Python
Use easy_install (or pip if you have it)

    sudo easy_install https://github.com/opencobra/cobrapy/archive/master.zip

## Hacking version installation in Python
First, clone the git repository using your preferred mothod. Cloning from your
own github fork is recommended! Afterwards, open a terminal, enter the cobrapy
repository and run the following command:

    python setup.py develop --user

If the command fails with an error about the --user option not being recognized,
it means setuptools is not installed. Either install setuptools before
trying again, or instead run ```sudo python setup.py develop```

## Installation of a Solver in Python
Currently cobrapy supports three linear programming solvers: ILOG/CPLEX, 
Gurobi, and GLPK (through [pyGLPK](http://tfinley.net/software/pyglpk/)). 
ILOG/CPLEX and Gurobi are commercial software packages that, currently, 
provide free licenses for academics and support both linear and quadratic 
programming. GLPK is an opensource linear programming solver; however, it 
does not support quadratic programming and is not as robust as the 
commercial solvers when it comes to mixed-integer linear programming.

### Installation of a commercial solver in Python
See the linked instructions for each solver
* [ILOG/CPLEX Academic](https://www.ibm.com/developerworks/university/academicinitiative/)
* [ILOG/CPLEX Commercial](http://www.ibm.com/software/integration/optimization/cplex-optimizer/)
* [Gurobi Academic & Commercial](http://gurobi.com)

### Installation of pyGLPK in Python
Please note that pyGLPK is not the same as python-glpk.

#### GNU/Linux Installation of pyGLPK in Python
1. Install the glpk and gmp library packages. You will need the development 
versions if they are available. You will also need development headers for 
Python itself. For example, Ubuntu and Debian Wheezy (7.0) users would type the following into the 
command line:
```
    sudo apt-get install libglpk-dev libgmp-dev python-dev python-setuptools
```
Debian Squeeze (6.0) users will need to build libgmp from source.

2. install pyglpk with easy_install using the following command in the terminal:
```
    sudo easy_install glpk 
```

#### MAC OS X Installation of pyGLPK in Python
1. Install homebrew if you don't have it. This may require downloading Xcode 
from the AppStore and Command Line Tools for XCode from 
https://developer.apple.com/devcenter/mac/index.action. If you're already 
using macports then just use that to install glpk.
```
    brew install glpk
```

2. install pyglpk with easy_install using the following command in the terminal:
```
    sudo easy_install glpk
```

#### Windows Installation of pyGLPK in Python
Download and install the executable from [here](https://sourceforge.net/projects/opencobra/files/python/cobra/extras/pyGLPK/).

## Installation of Optional Dependencies
Installation instructions are not provided for these libraries. However, 
many of them can be easily built on GNU/Linux with easy_install. On windows, 
many can downloaded from [this site](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

1. [libsbml](http://sbml.org) >= 4.0 to read/write SBML files
  * [Windows installer](http://www.lfd.uci.edu/~gohlke/pythonlibs/#libsbml)
2. [numpy](http://numpy.org) >= 1.6.1 and [scipy](http://scipy.org) >= 0.11 for 
ArrayBasedModel, double_deletion analysis, and saving to MAT formats.
  * Windows installers for 
  [numpy](http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy) and 
  [scipy](http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy)
3. [Parallel Python](http://parallelpython.org) for parallel processing.
4. [MATLAB](http://mathworks.com) and 
[mlabwrap](http://mlabwrap.sourceforge.net) for connecting to the COBRA 
Toolbox for MATLAB.
  * Installation is tricky on most platforms.



## Testing your installation
1. Start python
2. Type the following into the Python shell

```python
from cobra.test import test_all
test_all()
```

--------------------------------------------------------------------------------


#INSTALLATION - Jython
All releases require Jython (2.5+).  JSBML (http://sbml.org) is required for reading / writing SBML files.
NOTE: ArrayBasedModel, double_deletion analysis, saving to MAT formats, parallel processing,  and connecting to the
COBRA Toolbox for MATLAB are currently unavailable when using Jython.

I. Releases - Hosted on http://opencobra.sourceforge.net
   A. Download the most recent bzipped archive.
   B. Unzip the archive and make sure that the toplevel directory is in your Java CLASSPATH.
   C. Install a linear programming solver (Section III).

II. Development Code - Hosted on github
  NOTE: Not intended for general users.  Some functions require advanced capabilities / settings.  Unless you're
  willing to deal with headaches and heartaches it's better to install one of the releases.
  1. git pull https://github.com/opencobra/cobrapy.git
  2. Add the cobrapy directory to your Java CLASSPATH.
  3. Install a linear programming solver (Section III).

III. Installation of linear programming solvers
    On Jython, cobrapy currently supports two linear programming solvers: ILOG/CPLEX and Gurobi.

    A. ILOG/CPLEX and Gurobi are commercial software packages that, currently, provide free licenses for academics and
    support both linear and quadratic programming.
        1. Please download the software from their respective sites, install according to their instructions, and make
        sure that their Java jars are in your Java CLASSPATH.
        2. Current links are listed below. If they don't work then search using google.
          * [ILOG/CPLEX Academic](https://www.ibm.com/developerworks/university/academicinitiative/)
          * [ILOG/CPLEX Commercial](http://www.ibm.com/software/integration/optimization/cplex-optimizer/)
          * [Gurobi Academic & Commercial](http://gurobi.com)

    B. GLPK:  We are exploring the possibility of using GLPK through [GLPK for Java](http://glpk-java.sourceforge.net);
    however, we've encountered irregular memory errors that we'll need to trace before listing it as a supported solver.
    Advanced users may try http://glpk-java.sourceforge.net at their own risk.

IV. Test your installation.
  A. Start Jython
  B. Enter: from cobra.test import test_all
  C. Enter: test_all()
