#!/bin/bash -
#
# The AVX2 kernels (*_avx2.cc) are compiled with -mavx2, the rest of
# swipe without: the AVX2 code must only run when the processor has it
# (cpu_features.avx2). An inline function defined in an AVX2 object
# and in another object (a header function: View, std:: templates,
# swipe.h helpers) is emitted once in the binary, and the linker may
# keep the copy compiled with -mavx2: the other objects would then run
# AVX instructions on any processor. This script checks, after a build,
# that none of these shared functions contains an AVX instruction (VEX
# encoding: mnemonics starting with v, or ymm registers).
#
# usage: bash .github/scripts/avx2_isolation.sh [build directory]

set -o errexit -o nounset -o pipefail
export LC_ALL=C

BUILD="${1:-.}"
cd "${BUILD}"
[[ -x swipe ]] || { echo "no swipe binary in ${BUILD}" ; exit 2 ; }

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

# the weak symbols (inline functions, templates) of the AVX2 objects,
# and of the other objects
weak_symbols() { nm "${@}" | awk '$2 ~ /^[WV]$/ { print $3 }' | sort -u ; }
others=()
for object in ./*.o ; do
    [[ "${object}" == *_avx2.o ]] || others+=("${object}")
done
weak_symbols ./*_avx2.o > "${WORK}/avx2"
weak_symbols "${others[@]}" > "${WORK}/other"
comm -12 "${WORK}/avx2" "${WORK}/other" > "${WORK}/shared"

# the functions of the binary that contain an AVX instruction
objdump -d --no-show-raw-insn swipe | \
    awk '/^[0-9a-f]+ <.*>:$/ { name = $2 ; gsub(/^<|>:$/, "", name) ; next }
         /\tv[a-z0-9]+[ \t]/ || /%ymm/ { print name }' | \
    sort -u > "${WORK}/with_avx"

comm -12 "${WORK}/shared" "${WORK}/with_avx" > "${WORK}/leaks"

printf "%d inline functions shared by the AVX2 objects and the others\n" \
       "$(wc -l < "${WORK}/shared")"
if [[ -s "${WORK}/leaks" ]] ; then
    echo "AVX instructions in shared functions (run on any processor):"
    c++filt < "${WORK}/leaks" | sed 's/^/  /'
    exit 1
fi
echo "OK: none contains an AVX instruction"
