#!/bin/bash -

# Clone the black-box test-suite (swipe-tests) at the branch matching
# the swipe branch being tested:
# - a swipe-tests branch with the same name (twin branches), if any,
# - otherwise main for master (releases), and dev for other branches.
#
# usage: bash .github/scripts/fetch_tests.sh BRANCH [DIRECTORY]

set -euo pipefail

# TESTS_URL can be overridden (forks, local rehearsals)
readonly TESTS_URL="${TESTS_URL:-https://github.com/frederic-mahe/swipe-tests.git}"
readonly BRANCH="${1:?usage: fetch_tests.sh BRANCH [DIRECTORY]}"
readonly DIRECTORY="${2:-swipe-tests}"

if git ls-remote --exit-code --heads "${TESTS_URL}" "${BRANCH}" > /dev/null ; then
    TESTS_BRANCH="${BRANCH}"
elif [[ "${BRANCH}" == "master" ]] ; then
    TESTS_BRANCH="main"
else
    TESTS_BRANCH="dev"
fi

printf "swipe branch: %s, swipe-tests branch: %s\n" "${BRANCH}" "${TESTS_BRANCH}"
git clone --quiet --depth 1 --branch "${TESTS_BRANCH}" "${TESTS_URL}" "${DIRECTORY}"
