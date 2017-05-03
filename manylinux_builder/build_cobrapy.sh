#!/bin/bash

for PYBIN in /opt/python/*/bin; do
    ${PYBIN}/pip wheel cobra --pre
done

# Bundle external shared libraries into the wheels
for whl in cobra*.whl; do
    auditwheel repair $whl -w /io/wheelhouse/
done

