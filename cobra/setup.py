from setuptools import setup, find_packages
setup(
    name = "COBRApy",
    version = '0.01a.dev',
    packages = ['collections',
                'core',
                'external',
                'flux_analysis',
                'general',
                'io',
                'manipulation',
                'mlab',
                'query',
                'stats',
                'tools',
                'topology'],
    scripts = [''],

    install_requires = ['numpy>=2.0.0',
                        'scipy>=0.10.0',
                        'rpy2>=2.2.2'],
    extras_require = {
        'parallel': ["pp>=1.6"],
        'matlab': ["mlabwrap>=1.1"]
        },

    package_data = {
        '': ['*.txt', '*.html']},

    author = "Daniel Robert Hyduke",
    author_email = "danielhyduke@gmail.com",
    description = "COBRApy is a package for constraints-based modeling of biological networks",
    license = "GPL V3.0",
    keywords = "metabolism biology linear programming optimization flux balance analysis fba",
    url = "http://opencobra.sourceforge.net"
    )
    
