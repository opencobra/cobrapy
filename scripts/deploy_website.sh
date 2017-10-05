#!/usr/bin/env bash

set -e

target=/tmp/cobrapy-website

git clone "https://github.com/opencobra/cobrapy-website.git" "${target}"
cp -f "${TRAVIS_BUILD_DIR}"/release-notes/[0-9]*.md "${target}/content/releases/"
cd "${target}"
git add .
git commit -m "feat: publish ${TRAVIS_TAG} on $(date +'%F %T')"
cd "${TRAVIS_BUILD_DIR}"

