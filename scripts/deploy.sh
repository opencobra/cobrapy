#!/bin/bash
echo -e " ... running twine to deploy ... "
echo -e "
[pypirepository]
repository = ${PYPI_REPOSITORY}
" > ~/.pypirc

pip install twine
twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${PYPI_PASSWORD}" \
	  ${TRAVIS_BUILD_DIR}/wheelhouse/* -r pypirepository
