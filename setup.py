# -*- coding: utf-8 -*-

from __future__ import absolute_import

from sys import argv, version_info
from warnings import warn

from setuptools import find_packages, setup


if version_info[:2] == (3, 4):
    warn("Support for Python 3.4 was dropped by pandas. Since cobrapy is a "
         "pure Python package you can still install it but will have to "
         "carefully manage your own pandas and numpy versions. We no longer "
         "include it in our automatic testing.")

setup_kwargs = dict()
setup_requirements = []
# prevent pytest-runner from being installed on every invocation
if {'pytest', 'test', 'ptr'}.intersection(argv):
    setup_requirements.append("pytest-runner")

extras = {
    'array': ["scipy"],
    'sbml': ["python-libsbml", "lxml"]
}
extras["all"] = sorted(extras.values())

try:
    with open('README.rst') as handle:
        readme = handle.read()
    with open('INSTALL.rst') as handle:
        install = handle.read()
    setup_kwargs["long_description"] = readme + "\n\n" + install
except IOError:
    setup_kwargs["long_description"] = ''


if __name__ == "__main__":
    setup(
        name="cobra",
        version="0.15.1",
        packages=find_packages(),
        setup_requires=setup_requirements,
        install_requires=[
            "six",
            "future",
            "swiglpk",
            "ruamel.yaml>=0.15",
            "numpy>=1.13",
            "pandas>=0.17.0",
            "optlang>=1.4.2",
            "tabulate",
            "depinfo",
            "python-libsbml-experimental>=5.17.2",
        ],
        tests_require=[
            "jsonschema > 2.5",
            "pytest",
            "pytest-benchmark"
        ],
        extras_require=extras,
        package_data={
             '': [
                 'test/data/*',
                 'mlab/matlab_scripts/*m'
             ]
        },
        author="The cobrapy core team",
        author_email="cobra-pie@googlegroups.com",
        maintainer="Moritz E. Beber",
        maintainer_email="morbeb@biosustain.dtu.dk",
        description="COBRApy is a package for constraints-based modeling of "
        "biological networks",
        license="LGPL/GPL v2+",
        keywords=("metabolism biology linear programming optimization flux"
                  " balance analysis fba"),
        url="https://opencobra.github.io/cobrapy",
        test_suite="cobra.test.suite",
        download_url='https://pypi.python.org/pypi/cobra',
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU Lesser General Public License v2'
                ' or later (LGPLv2+)',
            'License :: OSI Approved :: GNU General Public License v2'
                ' or later (GPLv2+)',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: Implementation :: CPython',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
        platforms="GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7",
        **setup_kwargs
    )
