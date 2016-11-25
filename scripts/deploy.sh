#!/bin/bash

echo -e "
[distutils]
index-servers=
    pypirepository

[pypirepository]
repository: https://testpypi.python.org/pypi
username: henred
password: pippi_langstrump
" > ~/.pypirc

echo -e " starting deploy for branch ${TRAVIS_BRANCH} .."

if [[ "$TRAVIS_BRANCH" == "fix/deploy-problem" ]]; then
	echo -e " ... running twine to TEST deploy ... "
	pip install twine
	twine upload --skip-existing ${TRAVIS_BUILD_DIR}/wheelhouse/* -r pypirepository
elif [[ "$TRAVIS_BRANCH" == "master" ]]; then
	echo -e " ... running twine to deploy proper ... "
	# pip install twine
	# twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${PYPI_PASSWORD}" ${TRAVIS_BUILD_DIR}/wheelhouse/*
fi
