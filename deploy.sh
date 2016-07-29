#!/bin/bash
set -e

# seemingly bug in travis: TRAVIS_BRANCH doesn't get set to the branch but to the tag

if [[ -n "$TRAVIS_TAG" ]]; then
	echo -e " ... running twine to deploy ... "
	pip install twine
	twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${PYPI_PASSWORD}" ${TRAVIS_BUILD_DIR}/wheelhouse/*
else
	echo -e " ... skipping deploy as no tag detected: $TRAVIS_TAG ... "
fi
exit 0;
