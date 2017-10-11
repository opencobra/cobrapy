#!/usr/bin/env bash

set -e

target=/tmp/cobrapy-website

git clone "https://github.com/opencobra/cobrapy-website.git" "${target}"
cd "${target}"
git checkout master
python "${TRAVIS_BUILD_DIR}"/scripts/publish_releases.py "${TRAVIS_BUILD_DIR}"/release-notes "${target}"/content/releases "${TRAVIS_TAG}"
git add .
git commit -m "feat: publish ${TRAVIS_TAG} on $(date +'%F %T')"
git push origin master
cd "${TRAVIS_BUILD_DIR}"

