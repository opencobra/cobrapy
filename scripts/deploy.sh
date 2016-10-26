#!/bin/bash
echo -e " ... running twine to deploy ... "
pip install twine

twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${PYPI_PASSWORD}" ${TRAVIS_BUILD_DIR}/wheelhouse/*
