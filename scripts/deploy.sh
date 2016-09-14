#!/bin/bash
echo -e " ... running twine to deploy ... "
pip install twine

## twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${PYPI_PASSWORD}" ${TRAVIS_BUILD_DIR}/wheelhouse/*
twine upload --r "https://testpypi.python.org/pypi" --skip-existing --username testcobra --password cobra_metabolic_engineering ${TRAVIS_BUILD_DIR}/wheelhouse/*
