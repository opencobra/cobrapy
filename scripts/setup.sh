#!/usr/bin/env bash

if [ "${TOX_ENV}" == "pep8" | "${TOX_ENV}" == "rstcheck" ];then
    exit 0
fi


git clone https://github.com/matthew-brett/multibuild.git
cd multibuild
git checkout e6ebbfa
# matplotlib non-compatible as testing runs in venv (non-framework)
TEST_DEPENDS="swiglpk optlang sympy decorator cython codecov coverage numpy scipy python-libsbml jsonschema six pytest pytest-cov pytest-benchmark tabulate"
BUILD_DEPENDS="swiglpk optlang sympy cython numpy scipy pandas"
source multibuild/common_utils.sh
source multibuild/travis_steps.sh
before_install
