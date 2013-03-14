Installation instructions for Python and Jython are detailed below.  For installation help, please see the
sourceforge help install forum:
    https://sourceforge.net/projects/opencobra/forums/forum/4194275

Documentation is available here:
    http://opencobra.sourceforge.net/openCOBRA/opencobra_documentation/cobra_py/index.html

A tutorial is available here:


Example codes are available here:

------------------------------------------------------------------------------------------------------------------------

INSTALLATION - Python
All releases require Python 2.7 (or 2.6) to be installed before proceeding. Mac OS X (10.7+) ships with Python.
LibSBML (4.0+) should be installed from sbml.org to read / write SBML files.

I. Releases - Hosted on opencobra.sourceforge.net
  A. The easy way for GNU/Linux and Mac OS X (10.7+). If you've got easy_install or pip on your system.
    1. Open a Terminal
    2. Enter: easy_install cobra
    3. Install a linear programming solver (Section III) and any additional packages

  B. The slightly less easy way for GNU/Linux, Mac OS X, and Microsoft Windows 7:
    1. Download the most recent version for your operating system from here:
       http://sourceforge.net/projects/opencobra/files/python/cobra/
    2. Follow the installation instructions for your operating system.
       a. GNU/Linux
         i. If you have easy_install
           a. Open a Terminal
           b. Enter: easy_install the_path_to_the_file_and_the_file_that_you_downloaded
         ii. If you don't have easy_install on your machine then
           a. Unzip the file that you downloaded
           b. Open a Terminal
           c. Change to the top level directory of the file that you installed
           d. Enter: python setup.py install
       b. Mac OS/X
         i. Open a Terminal
         ii. Enter: easy_install the_path_to_the_file_and_the_file_that_you_downloaded
       c. Microsoft Windows 7 (or XP)
         i. Run the executable
    3. Install a linear programming solver (Section III) and any additional packages


II. Development Code - Hosted on github
  NOTE: Not intended for general users.  Some functions require advanced capabilities / settings.  Unless you're
  willing to deal with headaches and heartaches it's better to install one of the releases.
  1. git pull https://github.com/opencobra/cobrapy.git
  2. Use python setup.py install, easy_install, or add the cobrapy directory to your python path.
  3. Install a linear programming solver (Section III) and any additional packages

III. Installation of linear programming solvers
  cobrapy currently supports three linear programming solvers: ILOG/CPLEX, Gurobi, and GLPK (through pyGLPK).

  A. ILOG/CPLEX and Gurobi are commercial software packages that, currently, provide free licenses for academics and
  support both linear and quadratic programming.
    1. Please download the software from their respective sites, install according to their instructions, and make sure
   that their Python modules are installed or in your Python path.
    2. Current links are listed below.  If they don't work then search using google.
      ILOG/CPLEX Academic: https://www.ibm.com/developerworks/university/academicinitiative/
      ILOG/CPLEX Commercial: http://www.ibm.com/software/integration/optimization/cplex-optimizer/
      Gurobi Academic & Commercial: http://gurobi.com

  B. GLPK through PyGLPK.  GLPK is an opensource linear programming solver; however, it does not support quadratic
    programming and is not as robust as the commercial solvers when it comes to mixed-integer linear programming.
    NOTE: PyGLPK is not the same as python-glpk

    1. GNU/Linux
       a. Install GLPK using the package installer for your version of Linux.
          NOTE: If you get crashes when running cobrapy then you may need to compile GLPK from source on your
          system and perform the following steps because some distros have linkage problems.
       b. Download PyGLPK from here: http://tfinley.net/software/pyglpk/download.html
       c. If you have easy_install
         i. Open a Terminal
         ii. Enter: easy_install the_path_to_the_file_and_the_file_that_you_downloaded
       d. If you don't have easy_install on your machine then
         i. Unzip the file that you downloaded
         ii. Open a Terminal
         iii. Change to the top level directory of the file that you installed
         iv. Enter: python setup.py install

    2. Mac OS X
       a. Install homebrew if you don't have it (may require downloading Xcode from the AppStore and Command Line Tools
        for XCode from https://developer.apple.com/devcenter/mac/index.action)
         -if you're already using macports then just use that to install glpk
       b. Open a Terminal
       c. Enter: brew install glpk
       d. Download PyGLPK from here: http://tfinley.net/software/pyglpk/download.html
       e. In the Terminal, enter: easy_install the_path_to_the_file_and_the_file_that_you_downloaded


    3. Microsoft Windows 7 (or XP)
      a. Download and run the appropriate installer from here:
         https://sourceforge.net/projects/opencobra/files/python/cobra/extras/pyGLPK/

IV. Test your installation.
  A. Start Python
  B. Enter: from cobra.test import test_all
  C. Enter: test_all()
  NOTE: If you are using the GLPK Solver some tests may not work.


V. Optional Python packages for added features:
  A. Numpy >= 1.6.1 & Scipy >= 0.10 (http://scipy.org) for ArrayBasedModel, MoMA, double_deletion analysis, and
  saving to MAT formats.
  B. Parallel Python (http://parallelpython.org) for parallel processing.
  C. MATLAB (http://mathworks.com) and mlabwrap (http://mlabwrap.sourceforge.net) for connecting to the COBRA
  Toolbox for MATLAB.

------------------------------------------------------------------------------------------------------------------------


INSTALLATION - Jython
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
    On Jython, cobrapy currently supports tw linear programming solvers: ILOG/CPLEX and Gurobi.

    A. ILOG/CPLEX and Gurobi are commercial software packages that, currently, provide free licenses for academics and
    support both linear and quadratic programming.
      1. Please download the software from their respective sites, install according to their instructions, and make sure
     that their Java jars are in your Java CLASSPATH.
      2. Current links are listed below.  If they don't work then search using google.
        ILOG/CPLEX Academic: https://www.ibm.com/developerworks/university/academicinitiative/
        ILOG/CPLEX Commercial: http://www.ibm.com/software/integration/optimization/cplex-optimizer/
        Gurobi Academic & Commercial: http://gurobi.com

    B. GLPK:  We are exploring the possibility of using GLPK through GLPK for Java (http://glpk-java.sourceforge.net);
    however, we've encountered irregular memory errors that we'll need to trace before listing it as a supported solver.
    Advanced users may try http://glpk-java.sourceforge.net at their own risk.

IV. Test your installation.
  A. Start Jython
  B. Enter: from cobra.test import test_all
  C. Enter: test_all()