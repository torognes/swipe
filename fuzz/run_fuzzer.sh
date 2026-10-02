#!/bin/bash -
#
# Fuzz the ASN.1 header parser (asnparse.cc) with libFuzzer.
#
# usage: bash fuzz/run_fuzzer.sh [seconds] [libFuzzer options...]
#
# Builds fuzz/build/fuzz_asnparse (clang, AddressSanitizer and
# UndefinedBehaviorSanitizer), writes seed headers made by makeblastdb
# (BLAST+) into fuzz/build/seeds, then fuzzes for the given number of
# seconds (default: 60), with the corpus kept in fuzz/build/corpus
# between runs. A crashing input is saved as fuzz/build/crash-<sha1>;
# replay it with:
#
#   fuzz/build/fuzz_asnparse fuzz/build/crash-<sha1>
#
# Environment: CXX (default: the first of clang++, clang++-20 to
# clang++-14 found in the PATH).

set -o errexit -o nounset -o pipefail

SECONDS_TO_RUN="${1:-60}"
shift || true

FUZZ_DIR=$(cd "$(dirname "${0}")" && pwd)
BUILD_DIR="${FUZZ_DIR}/build"
SEEDS="${BUILD_DIR}/seeds"
CORPUS="${BUILD_DIR}/corpus"
mkdir -p "${BUILD_DIR}" "${SEEDS}" "${CORPUS}"

if [[ -z "${CXX:-}" ]] ; then
    for candidate in clang++ clang++-{20..14} ; do
        if command -v "${candidate}" > /dev/null ; then
            CXX="${candidate}"
            break
        fi
    done
fi
[[ -n "${CXX:-}" ]] || { echo "clang++ not found (set CXX)" ; exit 1 ; }

## build ---------------------------------------------------------------------

# fatal() jumps back to the harness (no exception), and the harness
# replaces operator new and delete (see fuzz_asnparse.cc)
FLAGS=(-std=c++11 -fsized-deallocation -g -O1 -fno-omit-frame-pointer -pthread
       "-fsanitize=fuzzer,address,undefined" -fno-sanitize-recover=undefined
       -Wall -Wextra)
"${CXX}" "${FLAGS[@]}" \
         -o "${BUILD_DIR}/fuzz_asnparse" \
         "${FUZZ_DIR}/fuzz_asnparse.cc" "${FUZZ_DIR}/../asnparse.cc"

## seeds ---------------------------------------------------------------------

# one header per seed file: the slices of the .phr file given by the
# header offsets of the .pin file (BLAST database version 4)
split_headers() {
    python3 - "${1}" "${2}" << 'END_OF_PYTHON'
import struct, sys
base, prefix = sys.argv[1], sys.argv[2]
pin = open(base + ".pin", "rb").read()
phr = open(base + ".phr", "rb").read()
pos = 8
title_length, = struct.unpack(">I", pin[pos:pos + 4]); pos += 4 + title_length
date_length, = struct.unpack(">I", pin[pos:pos + 4]); pos += 4 + date_length
pos += (8 - pos % 8) % 8
seqcount, = struct.unpack(">I", pin[pos:pos + 4]); pos += 4 + 8 + 4
offsets = struct.unpack(">%dI" % (seqcount + 1), pin[pos:pos + 4 * (seqcount + 1)])
for i in range(seqcount):
    open("%s_%d" % (prefix, i), "wb").write(phr[offsets[i]:offsets[i + 1]])
END_OF_PYTHON
}

if command -v makeblastdb > /dev/null ; then
    WORK=$(mktemp -d)
    trap 'rm -rf "${WORK}"' EXIT
    # every kind of seq-id that asnparse.cc decodes, titles with and
    # without spaces, several deflines per header (^A), taxids
    DEFLINES=(
        'gi|123|sp|P12345.2|NAME_HUMAN some protein'
        'sp|Q12345|NAME_MOUSE'
        'tr|A0A000|A0A000_9ZZZZ unreviewed entry'
        'ref|NP_000001.1| reference protein'
        'gb|AAA00001.1| genbank entry'
        'emb|CAA00001.1| embl entry'
        'dbj|BAA00001.1| ddbj entry'
        'pir||A00001 pir entry'
        'prf||1234567A prf entry'
        'tpg|DAA00001.1| tpg entry'
        'tpe|CBB00001.1| tpe entry'
        'tpd|FAA00001.1| tpd entry'
        'gpp|XX_000001.1| gpipe entry'
        'pdb|1ABC|A upper case chain'
        'pdb|1ABC|a lower case chain'
        'pdb|7XYZ|AAA long chain name'
        'gnl|mydb|abc123 general string id'
        'gnl|mydb|42 general numeric id'
        'lcl|local_id local'
        'lcl|17'
        'pat|US|5543331|1 granted patent'
        'pgp|EP|0238993|7 patent application'
        'bbs|1234'
        'bbm|5678'
        "gi|1|sp|P1|A_HUMAN first"$'\x01'"sp|P2|B_MOUSE second"$'\x01'"pdb|2XYZ|b third"
        'no_seqid_title_only with <xml> & "quotes"'
    )
    TAXIDS=(9606 10090 0 2147483647)
    number=0
    for defline in "${DEFLINES[@]}" ; do
        taxid="${TAXIDS[$(( number % ${#TAXIDS[@]} ))]}"
        printf ">%s\nMKVLAAGIVALLLAAGCSSA\n" "${defline}" > "${WORK}/in.fa"
        # some deflines are not parsed as seq-ids: try both ways
        for parse in -parse_seqids "" ; do
            if makeblastdb -in "${WORK}/in.fa" -dbtype prot ${parse} \
                           -blastdb_version 4 -taxid "${taxid}" \
                           -out "${WORK}/db" > /dev/null 2>&1 ; then
                split_headers "${WORK}/db" "${SEEDS}/seed_${number}${parse}"
            fi
        done
        number=$(( number + 1 ))
    done
    rm -rf "${WORK}"
    trap - EXIT
else
    echo "makeblastdb not found: fuzzing from the existing corpus only"
fi

## fuzz ----------------------------------------------------------------------

# the parser's own messages go to stderr with libFuzzer's; no symbol
# download from debuginfod
DEBUGINFOD_URLS="" \
UBSAN_OPTIONS="print_stacktrace=1:${UBSAN_OPTIONS:-}" \
"${BUILD_DIR}/fuzz_asnparse" \
    -max_total_time="${SECONDS_TO_RUN}" \
    -dict="${FUZZ_DIR}/asnparse.dict" \
    -artifact_prefix="${BUILD_DIR}/" \
    "${@}" \
    "${CORPUS}" "${SEEDS}"
