#!/bin/bash -
#
# The x86-64 kernels of arch/x86_64/AVX2/ are compiled with -mavx2,
# those of arch/x86_64/SSSE3/ with -mssse3, the rest of swipe without:
# that code must only run when the processor has the instruction set
# (arch/x86_64/search_dispatch.cc). An inline function defined in one
# of these objects and in another object (a header function: View,
# align_cells(), std:: templates, swipe.h helpers) is emitted once in
# the binary, and the linker may keep the copy compiled with -mavx2 or
# -mssse3: the other objects would then run these instructions on any
# processor. This script checks, after a build, that none of these
# shared functions contains an AVX instruction (VEX encoding:
# mnemonics starting with v, or ymm registers), nor, for the SSSE3
# objects, an SSE3 or SSSE3 instruction.
# [formerly avx2_isolation.sh, AVX2 objects only]
#
# The check is meant for release builds. At -O0 (make DEBUG=1) nothing
# is inlined: every template an AVX2 kernel uses (std::array members,
# the _GLIBCXX_DEBUG error formatter) gets an out-of-line AVX2 copy,
# and which copy the linker keeps varies from one linker to another.
#
# usage: bash .github/scripts/isa_isolation.sh [build directory]

set -o errexit -o nounset -o pipefail
export LC_ALL=C

BUILD="${1:-.}"
cd "${BUILD}"
[[ -x swipe ]] || { echo "no swipe binary in ${BUILD}" ; exit 2 ; }

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

# the weak symbols (inline functions, templates) of some objects
weak_symbols() { nm "${@}" | awk '$2 ~ /^[WV]$/ { print $3 }' | sort -u ; }

# the functions of the binary with an instruction matching a pattern
functions_with() {
    objdump -d --no-show-raw-insn swipe | \
        awk -v pattern="${1}" \
            '/^[0-9a-f]+ <.*>:$/ { name = $2 ; gsub(/^<|>:$/, "", name) ; next }
             $0 ~ pattern { print name }' | \
        sort -u
}

shopt -s nullglob globstar
avx2_objects=(arch/x86_64/AVX2/*.o)
ssse3_objects=(arch/x86_64/SSSE3/*.o)
not_avx2=()
baseline=()
for object in ./**/*.o ; do
    case "${object}" in
        ./arch/x86_64/AVX2/*) ;;
        ./arch/x86_64/SSSE3/*) not_avx2+=("${object}") ;;
        *) not_avx2+=("${object}") ; baseline+=("${object}") ;;
    esac
done
if (( ${#avx2_objects[@]} == 0 || ${#ssse3_objects[@]} == 0 )) ; then
    echo "no object in arch/x86_64/AVX2/ or arch/x86_64/SSSE3/ (ARCH=x86_64 build?)"
    exit 2
fi

status=0

# $1: name of the instruction set, $2: pattern of its instructions,
# $3: file of the weak symbols of its objects, $4: file of the weak
# symbols of the objects that run on processors without it
check() {
    comm -12 "${3}" "${4}" > "${WORK}/shared"
    functions_with "${2}" > "${WORK}/with_isa"
    comm -12 "${WORK}/shared" "${WORK}/with_isa" > "${WORK}/leaks"
    printf "%d inline functions shared by the %s objects and the others\n" \
           "$(wc -l < "${WORK}/shared")" "${1}"
    if [[ -s "${WORK}/leaks" ]] ; then
        echo "${1} instructions in shared functions (run on any processor):"
        c++filt < "${WORK}/leaks" | sed 's/^/  /'
        status=1
        return
    fi
    echo "OK: none contains an ${1} instruction"
}

weak_symbols "${avx2_objects[@]}" > "${WORK}/avx2"
weak_symbols "${not_avx2[@]}" > "${WORK}/not_avx2"
check AVX '\tv[a-z0-9]+[ \t]|%ymm' "${WORK}/avx2" "${WORK}/not_avx2"

weak_symbols "${ssse3_objects[@]}" > "${WORK}/ssse3"
weak_symbols "${baseline[@]}" > "${WORK}/baseline"
check SSSE3 '\t(pshufb|palignr|phadd|phsub|pabs|pmaddubsw|pmulhrsw|psign|lddqu|movddup|movshdup|movsldup|addsubp|haddp|hsubp|fisttp)[a-z]*[ \t]' \
      "${WORK}/ssse3" "${WORK}/baseline"

exit "${status}"
