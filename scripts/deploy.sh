#!/bin/bash
echo -e " ... running twine to deploy ... "
pip install twine

# only upload mac wheels for now, linux wheels are broken due to https://github.com/pypa/manylinux/issues/80
twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${PYPI_PASSWORD}" ${TRAVIS_BUILD_DIR}/wheelhouse/*macosx*.whl
