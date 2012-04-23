import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
setup(
    name = "cobra",
    version = '0.01a.dev',
    packages = find_packages(exclude=['internal', 'oven', 'db_tools', 'omics', 'general','mlab']),
    #scripts = [''],
    setup_requires = ['numpy>=1.6.1',
                        'scipy>=0.10.1',
                        'rpy2>=2.2.2',
                      'libsbml>=5.1.0_b0',
                      'pp>=1.6.0'],
    install_requires = setup_requires,
    extras_require = {
        'matlab': ["mlabwrap>=1.1"]
        },

    ## package_data = {
    ##     '': ['*.txt', '*.html','gpl-3.0','README']},

    author = "Daniel Robert Hyduke",
    author_email = "danielhyduke@gmail.com",
    description = "COBRApy is a package for constraints-based modeling of biological networks",
    license = "GPL V3.0",
    keywords = "metabolism biology linear programming optimization flux balance analysis fba",
    url = "http://opencobra.sourceforge.net"
    )
    
