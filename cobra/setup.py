from setuptools import setup, find_packages
setup(
    name = "COBRApy",
    version = '0.01a.dev',
    packages = find_packages(exclude=['internal', 'oven', 'db_tools', 'omics', 'general']),
    #scripts = [''],
    setup_requires = ['numpy>=2.0.0.dev',
                        'scipy>=0.10.0.dev',
                        'rpy2>=2.2.2',
                        'libsbml=5.1.0_b0'],
    install_requires = ['numpy>=2.0.0.dev',
                        'scipy>=0.10.0.dev',
                        'rpy2>=2.2.2',
                        'libsbml=5.1.0_b0'],
    extras_require = {
        'parallel': ["pp>=1.6"],
        'matlab': ["mlabwrap>=1.1"]
        },

    package_data = {
        '': ['*.txt', '*.html','gpl-3.0','README']},

    author = "Daniel Robert Hyduke",
    author_email = "danielhyduke@gmail.com",
    description = "COBRApy is a package for constraints-based modeling of biological networks",
    license = "GPL V3.0",
    keywords = "metabolism biology linear programming optimization flux balance analysis fba",
    url = "http://opencobra.sourceforge.net"
    )
    
