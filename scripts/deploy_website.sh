#!/usr/bin/env bash

set -eu

target="${HOME}/cobrapy-website"

git clone "https://github.com/opencobra/cobrapy-website.git" "${target}"
pushd "${target}"

git config user.name "COBRApy Bot"
git config user.email "opencobra.cobrapy@gmail.com"
git remote rm origin
git remote add origin "https://user:${WEBSITE_DEPLOY_TOKEN}@github.com/opencobra/cobrapy-website.git" &> /dev/null

git checkout master
python "${WORKSPACE}/scripts/publish_release.py" "${WORKSPACE}/release-notes" "${target}/content/releases" "${TAG}"
git add .
git commit -m "feat: publish ${TAG} on $(date --utc --iso-8601=seconds)"
git push origin master

popd

