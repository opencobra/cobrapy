#!/bin/bash


echo -e " starting deploy for branch ${TRAVIS_BRANCH} .."

pip install twine
twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${PYPI_PASSWORD}" ${TRAVIS_BUILD_DIR}/wheelhouse/*
