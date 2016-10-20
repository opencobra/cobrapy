#!/bin/bash

set -e

for PYBIN in /Library/Frameworks/Python.framework/Versions/*/bin; do
	${PYBIN}/pip install wheel delocate
    ${PYBIN}/python setup.py bdist_wheel
done

# Bundle external shared libraries into the wheels
for whl in dist/cobra*.whl; do
    delocate-wheel $whl -w dist/wheelhouse/
done

