import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
setup(
    name = "cobra",
    version = '0.01a.dev',
    packages = find_packages(exclude=['internal', 'oven', 'db_tools', 'omics', 'general',]),
    #scripts = [''],
    #put in numpy, scipy, libsbml, and pyglpk
    setup_requires = [],
    install_requires = [],
    extras_require = {
        'parallel': ['pp>=1.6.0'],
        'matlab': ["mlabwrap>=1.1"],
        'R': ["rpy2>=2.2.2"]
        },

    package_data = {

         '': ['*.txt', '*.html','LICENSE','README','test/data/*','documentation/html/*',
              'examples/*py', 'examples/files/*']},

    author = "Daniel Robert Hyduke",
    author_email = "danielhyduke@gmail.com",
    description = "COBRApy is a package for constraints-based modeling of biological networks",
    license = "GPL V3.0",
    keywords = "metabolism biology linear programming optimization flux balance analysis fba",
    url = "http://opencobra.sourceforge.net"
    )
    
