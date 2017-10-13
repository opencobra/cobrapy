#!/usr/bin/env bash

set -e

if [[ -n "${TRAVIS_TAG}" && "${TRAVIS_OS_NAME}" == "linux" && "${MB_PYTHON_VERSION}" == "3.6" ]]; then
    echo "Parsing ${TRAVIS_BUILD_DIR}/release-notes/${TRAVIS_TAG}.md"
    export TAG_NOTES=$(cat "${TRAVIS_BUILD_DIR}/release-notes/${TRAVIS_TAG}.md")
fi

set +e
