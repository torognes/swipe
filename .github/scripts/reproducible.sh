#!/bin/bash -
#
# Reproducible builds: build the committed sources (HEAD) twice, in two
# directories with different names, with a different time zone, locale
# and umask, and check that the two binaries are byte-identical.
#
# usage: bash .github/scripts/reproducible.sh [make arguments...]
#
# e.g. bash .github/scripts/reproducible.sh RELEASE=1 LDFLAGS=-static
#
# Only the committed files are built (git archive): uncommitted changes
# are not checked.

set -o errexit -o nounset -o pipefail

TOP=$(git rev-parse --show-toplevel)
WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

build() {
    local directory="${1}"
    shift
    mkdir -p "${directory}"
    git -C "${TOP}" archive --format=tar HEAD | tar -x -C "${directory}"
    make -C "${directory}" -j "$(getconf _NPROCESSORS_ONLN)" "${@}" > /dev/null
}

build "${WORK}/first" "${@}"
(
    umask 077
    export TZ="Pacific/Kiritimati" LC_ALL=C
    build "${WORK}/second_build_in_a_longer_directory_name" "${@}"
)

FIRST="${WORK}/first/swipe"
SECOND="${WORK}/second_build_in_a_longer_directory_name/swipe"
sha256sum "${FIRST}" "${SECOND}" | sed "s|${WORK}/||"
if cmp -s "${FIRST}" "${SECOND}" ; then
    echo "REPRODUCIBLE: the two builds are byte-identical"
else
    echo "NOT REPRODUCIBLE: the two builds differ"
    cmp "${FIRST}" "${SECOND}" | head -n 1 || true
    exit 1
fi
