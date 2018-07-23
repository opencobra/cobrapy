#!/usr/bin/env bash

set -eu

target=/tmp/cobrapy-website

git clone "https://github.com/opencobra/cobrapy-website.git" "${target}"
cd "${target}"

git config user.name "Deployment Bot"
git config user.email "deploy@travis-ci.org"
git remote rm origin
git remote add origin "https://user:${GITHUB_TOKEN}@github.com/opencobra/cobrapy-website.git" &> /dev/null

git checkout master
python "${TRAVIS_BUILD_DIR}"/scripts/publish_release.py "${TRAVIS_BUILD_DIR}"/release-notes "${target}"/content/releases "${TRAVIS_TAG}"
git add .
git commit -m "feat: publish ${TRAVIS_TAG} on $(date +'%F %T')"
git push origin master
cd "${TRAVIS_BUILD_DIR}"

