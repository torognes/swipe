#!/bin/bash -
#
# The portable vector operations (arch/generic/intrinsics_to_functions.h)
# must give the same results as the x86-64 intrinsics: build
# generic_operations.cc with both, run both, compare the outputs (one
# line per operation and input: every pair of 8-bit lanes, edge and
# pseudo-random 16-bit lanes). On x86-64 only.
#
# usage: bash .github/scripts/generic_operations.sh [compiler]

set -o errexit -o nounset -o pipefail
export LC_ALL=C

CXX="${1:-${CXX:-g++}}"
SOURCES="$(cd "$(dirname "${0}")/../.." && pwd)"
WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

FLAGS=(-std=c++11 -O2 -Wall -Wextra -Wpedantic -Werror -I "${SOURCES}")
"${CXX}" "${FLAGS[@]}" -o "${WORK}/x86_64" "${SOURCES}/.github/scripts/generic_operations.cc"
"${CXX}" "${FLAGS[@]}" -DSWIPE_ARCH_GENERIC -o "${WORK}/generic" "${SOURCES}/.github/scripts/generic_operations.cc"
"${WORK}/x86_64" > "${WORK}/x86_64.txt"
"${WORK}/generic" > "${WORK}/generic.txt"

printf "%d results compared\n" "$(wc -l < "${WORK}/x86_64.txt")"
if ! cmp -s "${WORK}/x86_64.txt" "${WORK}/generic.txt" ; then
    echo "DIFFER: the portable operations do not give the x86-64 results"
    diff "${WORK}/x86_64.txt" "${WORK}/generic.txt" | head -n 20
    exit 1
fi
echo "OK: the portable operations give the x86-64 results"
