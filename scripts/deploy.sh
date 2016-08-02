#!/bin/bash
echo -e " ... running twine to deploy ... "
echo -e "
[distutils]
index-servers=
    pypirepository

[pypirepository]
repository: ${PYPI_REPOSITORY}
" > ~/.pypirc

cat ~/.pypirc

pip install twine
twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${TEST_PYPI_PASSWORD}" \
	  ${TRAVIS_BUILD_DIR}/wheelhouse/* -r pypirepository
