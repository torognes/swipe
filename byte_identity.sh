#!/bin/bash -
#
# Byte-identity matrix between two swipe binaries.
#
# usage: bash byte_identity.sh <binary_a> <binary_b>
#
# Runs the same command matrix with both binaries and compares their
# outputs (stdout, stderr, exit status, and the -o file): a refactoring
# must not change a single byte. Ported from swarm's and vsearch's
# byte_identity.sh. The databases and queries are generated here
# (makeblastdb, NCBI BLAST+), small and deliberately varied: every
# kind of sequence id, several deflines per sequence, taxids, long
# titles, lowercase residues, ambiguities, stop codons, ties, and
# more hits than shown. For the large datasets, use the conformity
# harness (perf-results/hyperfine/conformity/, not tracked).
#
# Left out of the comparison: the lines that change from run to run
# (dates, times, speed), and the version number (so that a release
# can be compared with the next one).

BIN_A="${1}"
BIN_B="${2}"

[[ -x "${BIN_A}" ]] || { echo "not executable: ${BIN_A}" ; exit 1 ; }
[[ -x "${BIN_B}" ]] || { echo "not executable: ${BIN_B}" ; exit 1 ; }
command -v makeblastdb > /dev/null || { echo "makeblastdb not found" ; exit 1 ; }

# each run happens inside the data directory, so a relative binary
# path would no longer resolve there
BIN_A=$(cd "$(dirname "${BIN_A}")" && pwd)/$(basename "${BIN_A}")
BIN_B=$(cd "$(dirname "${BIN_B}")" && pwd)/$(basename "${BIN_B}")

# Compare like with like: a DEBUG build (sanitizers) reports sanitizer
# messages on stderr that have nothing to do with the code under test.
# (The size of the binaries is no guide: with or without debugging
# information, a release build can be twice as large as another.)
has_sanitizer() { grep -q -a "__asan_init" "${1}" ; }
if has_sanitizer "${BIN_A}" && ! has_sanitizer "${BIN_B}" || \
   ! has_sanitizer "${BIN_A}" && has_sanitizer "${BIN_B}" ; then
    echo "REFUSING: one binary has the sanitizers (DEBUG=1), the other not"
    echo "  (${BIN_A})"
    echo "  (${BIN_B})"
    echo "Rebuild both the same way (both DEBUG=1, or both RELEASE=1) and re-run."
    exit 2
fi

WORKDIR=$(mktemp -d)
DATA="${WORKDIR}/data"
mkdir -p "${DATA}"
# keep the directory when something differed, to look at the files
# shellcheck disable=SC2329  # invoked indirectly, by the EXIT trap below
cleanup() { [[ "${differences:-0}" -eq 0 ]] && rm -rf "${WORKDIR}" ; }
trap cleanup EXIT

# --- inputs -----------------------------------------------------------------

# random sequences, fixed seed: $1 = alphabet, $2 = count, $3 = length,
# $4 = id prefix, $5 = seed
random_sequences() {
    awk -v alphabet="${1}" -v count="${2}" -v seqlength="${3}" \
        -v prefix="${4}" -v seed="${5}" '
        BEGIN {
            srand(seed)
            for (i = 0 ; i < count ; i++) {
                s = ""
                for (j = 0 ; j < seqlength ; j++) {
                    s = s substr(alphabet, int(rand() * length(alphabet)) + 1, 1)
                }
                printf ">%s%d random sequence %d\n%s\n", prefix, i, i, s
            }
        }'
}

AMINO_ACIDS="ACDEFGHIKLMNPQRSTVWY"
NUCLEOTIDES="ACGT"

# protein database: hand-written entries (every seq-id kind that
# asnparse.cc decodes, several deflines, a long title, XML special
# characters), then random sequences, ten of them twice (ties)
{
    printf ">gi|123|sp|P12345.2|NAME_HUMAN protein one\x01pdb|1ABC|a chain a\x01pdb|7XYZ|AAA long chain\n"
    printf "MKVLAAGIVALLLAAGCSSAKEETPVQVEPEQAAPAEEKKAEPKPKKAEKK\n"
    printf ">ref|NP_000001.1| reference <protein> & \"quotes\"\x01gnl|mydb|abc123 general\x01lcl|local7\n"
    printf "MDSRPAQSEGSWVMHIMFKNWCVYQYKTPLELVSLLQPLKQGMFTVIRKMF\n"
    printf ">tr|A0A000|A0A000_9ZZZZ unreviewed %s\n" "$(printf 'long title %.0s' {1..40})"
    printf "MKVLAAGIVALLLAAGCSSAKEETPVQVEPEQAAPAEEKKAEPKPKKAEKKXBZ*\n"
    printf ">pat|US|5543331|1 patent\x01gb|AAA00001.1| genbank\x01emb|CAA00001.1| embl\n"
    printf "HRYMYWQDYKGAAFSNISSETFFDLAEDRSQSWWSINDPKVQAAYGPEDN\n"
    random_sequences "${AMINO_ACIDS}" 200 80 "sp|R" 1
    random_sequences "${AMINO_ACIDS}" 10 80 "sp|R" 1 | sed "s/^>sp|R/>sp|D/"
} > "${DATA}/prot.fasta"

# nucleotide database: the same layout
{
    printf ">gi|42|gb|AB000001.1| nucleotide one\x01emb|X00001.1| embl copy\n"
    printf "ATGAAAGTTCTGGCTTGGCCGGATGAAAGTTCTGGCTTGGCCATGAAATAA\n"
    printf ">ref|NM_000001.1| ambiguities RYKMSWBDHVN\n"
    printf "ATGRYKMSWBDHVNACGTACGTNNNNNNNNNNACGTACGTGGATCCAAGC\n"
    random_sequences "${NUCLEOTIDES}" 200 150 "gb|R" 2
    random_sequences "${NUCLEOTIDES}" 10 150 "gb|R" 2 | sed "s/^>gb|R/>gb|D/"
} > "${DATA}/nucl.fasta"

# taxids for some sequences (-H, -x)
{
    printf "P12345.2 9606\nNP_000001.1 10090\nAB000001.1 9606\n"
    for (( i = 0 ; i < 200 ; i += 3 )) ; do printf "R%d %d\n" "${i}" $(( i % 7 + 1 )) ; done
} > "${DATA}/taxid.map"
printf "9606\n3\n5\n" > "${DATA}/taxids.txt"

for TYPE in prot nucl ; do
    makeblastdb -in "${DATA}/${TYPE}.fasta" -dbtype "${TYPE}" \
                -blastdb_version 4 -parse_seqids -title "${TYPE} ids" \
                -taxid_map "${DATA}/taxid.map" \
                -out "${DATA}/${TYPE}" > /dev/null 2>&1 || {
        echo "makeblastdb failed (${TYPE})" ; exit 2 ; }
    # the same sequences, ids not parsed (BL_ORD_ID), three volumes
    makeblastdb -in "${DATA}/${TYPE}.fasta" -dbtype "${TYPE}" \
                -blastdb_version 4 -title "${TYPE} split" \
                -max_file_sz 8000 \
                -out "${DATA}/${TYPE}_split" > /dev/null 2>&1 || {
        echo "makeblastdb failed (${TYPE}, split)" ; exit 2 ; }
done

# queries: database sequences, random ones, lowercase, odd residues
{
    sed -n 1,2p "${DATA}/prot.fasta"
    sed -n 7,8p "${DATA}/prot.fasta"
    random_sequences "${AMINO_ACIDS}" 3 60 "q" 3
    printf ">lower case\nmkvlaagivalllaagcssakeetpvqvep\n"
    printf ">odd residues\nMKVLAXBZJUO*-KEETPVQ\n"
} > "${DATA}/prot_query.fasta"
{
    sed -n 1,2p "${DATA}/nucl.fasta"
    sed -n 5,6p "${DATA}/nucl.fasta"
    random_sequences "${NUCLEOTIDES}" 3 120 "q" 4
    printf ">lower case ambiguous\natgaaagttctggcttggccnnnryatgaaagttctggcttgg\n"
} > "${DATA}/nucl_query.fasta"

# a score matrix file (-M FILE): BLOSUM62 as printed by swipe is not
# available, so a small identity-like matrix
{
    printf "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *\n"
    for row in A R N D C Q E G H I L K M F P S T W Y V '*' ; do
        printf "%s" "${row}"
        for col in A R N D C Q E G H I L K M F P S T W Y V '*' ; do
            if [[ "${row}" == "${col}" ]] ; then printf " %2d" 5 ; else printf " %2d" -2 ; fi
        done
        printf "\n"
    done
} > "${DATA}/matrix.txt"

# --- running and comparing --------------------------------------------------

without_volatile_lines() {
    sed -E \
        -e "s/SWIPE [0-9]+\.[0-9]+\.[0-9]+/SWIPE x.y.z/" \
        -e "s|<version>[^<]*</version>|<version>x.y.z</version>|" \
        "${1}" | \
    grep -v -E \
         -e "^(Search started|Search completed|Elapsed|Speed): " \
         -e "<(searchStarted|searchCompleted|searchElapsedTime|searchSpeed)>"
}

run_one() {
    # $1 = binary, $2 = output directory, rest = options
    local binary="${1}" outdir="${2}"
    shift 2
    rm -rf "${outdir}"
    mkdir -p "${outdir}"
    (
        cd "${DATA}" || exit 1
        # both binaries named swipe (argv[0] is shown by the help)
        (exec -a swipe "${binary}" "${@}") < /dev/null \
            > "${outdir}/raw.stdout" 2> "${outdir}/raw.stderr"
        echo "$?" > "${outdir}/out.status"
    )
    without_volatile_lines "${outdir}/raw.stdout" > "${outdir}/out.stdout"
    without_volatile_lines "${outdir}/raw.stderr" > "${outdir}/out.stderr"
    # the output file of -o out.file, written in the data directory
    if [[ -e "${DATA}/out.file" ]] ; then
        without_volatile_lines "${DATA}/out.file" > "${outdir}/out.file"
        rm -f "${DATA}/out.file"
    fi
}

COMPARED=(stdout stderr status file)

differences=0
comparisons=0
runs=0
failed_runs=0
empty_runs=0

compare_run() {
    run_one "${BIN_A}" "${WORKDIR}/a" "${@}"
    run_one "${BIN_B}" "${WORKDIR}/b" "${@}"
    runs=$(( runs + 1 ))
    [[ "$(cat "${WORKDIR}/a/out.status")" -ne 0 ]] && failed_runs=$(( failed_runs + 1 ))
    [[ -s "${WORKDIR}/a/out.stdout" ]] || empty_runs=$(( empty_runs + 1 ))
    [[ -n "${VERBOSE:-}" ]] && \
        echo "status $(cat "${WORKDIR}/a/out.status"), $(wc -l < "${WORKDIR}/a/out.stdout") lines: swipe ${*}"
    local suffix
    for suffix in "${COMPARED[@]}" ; do
        # no output file in either run: nothing to compare
        if [[ ! -e "${WORKDIR}/a/out.${suffix}" && ! -e "${WORKDIR}/b/out.${suffix}" ]] ; then
            continue
        fi
        comparisons=$(( comparisons + 1 ))
        if ! cmp -s "${WORKDIR}/a/out.${suffix}" "${WORKDIR}/b/out.${suffix}" ; then
            differences=$(( differences + 1 ))
            echo "DIFF: swipe ${*} [${suffix}]"
            diff "${WORKDIR}/a/out.${suffix}" "${WORKDIR}/b/out.${suffix}" | head -6
            cp -r "${WORKDIR}/a" "${WORKDIR}/diff_${differences}_a"
            cp -r "${WORKDIR}/b" "${WORKDIR}/diff_${differences}_b"
        fi
    done
}

# a database or query that swipe cannot read would make every run
# fail the same way in both binaries, and the matrix would report
# BYTE-IDENTICAL without comparing any search: check first
smoke_test() {
    run_one "${BIN_A}" "${WORKDIR}/smoke" "${@}" -m 8
    if [[ "$(cat "${WORKDIR}/smoke/out.status")" -ne 0 ]] || \
       [[ $(wc -l < "${WORKDIR}/smoke/out.stdout") -lt 20 ]] ; then
        echo "REFUSING: no hits from swipe ${*} -m 8"
        tail -n 2 "${WORKDIR}/smoke/out.stderr"
        exit 2
    fi
}

PROT=(-d prot -i prot_query.fasta)
NUCL=(-d nucl -i nucl_query.fasta)
smoke_test "${PROT[@]}"
smoke_test "${NUCL[@]}" -p 0

# output formats, symbol types, threads
for outfmt in 0 7 8 9 99 ; do
    compare_run "${PROT[@]}" -m "${outfmt}"
    compare_run "${PROT[@]}" -m "${outfmt}" -I -H -e 1000 -a 4
    compare_run "${NUCL[@]}" -m "${outfmt}" -p 0
    compare_run -d prot -i nucl_query.fasta -m "${outfmt}" -p 2
    compare_run -d nucl -i prot_query.fasta -m "${outfmt}" -p 3
    compare_run "${NUCL[@]}" -m "${outfmt}" -p 4
    compare_run "${PROT[@]}" -m "${outfmt}" -p 5
done
for threads in 1 3 16 ; do
    compare_run "${PROT[@]}" -m 0 -e 1000 -a "${threads}"
    compare_run "${NUCL[@]}" -m 0 -e 1000 -p 4 -a "${threads}"
done

# multi-volume databases, ordinal ids
compare_run -d prot_split -i prot_query.fasta -m 0 -e 1000 -a 3
compare_run -d nucl_split -i nucl_query.fasta -m 8 -p 0 -e 1000 -a 3
compare_run -d nucl_split -i prot_query.fasta -m 99 -p 3 -v 5 -b 3

# hit limits and thresholds
for limits in "-v 0 -b 0" "-v 1 -b 1" "-v 5 -b 2" "-v 300 -b 300" \
              "-e 1e-5" "-e 1000" "-k 1e-3 -e 1000" "-c 30" "-u 40" ; do
    # shellcheck disable=SC2086  # split on purpose
    compare_run "${PROT[@]}" -m 0 ${limits}
    # shellcheck disable=SC2086
    compare_run "${NUCL[@]}" -m 8 -p 0 ${limits}
done
compare_run "${PROT[@]}" -m 8 -z 1000000000

# scoring: matrices, gap penalties, nucleotide scores
for matrix in BLOSUM45 BLOSUM50 BLOSUM62 BLOSUM80 BLOSUM90 PAM30 PAM70 PAM250 matrix.txt ; do
    compare_run "${PROT[@]}" -m 8 -e 1000 -M "${matrix}" -G 11 -E 1
done
compare_run "${PROT[@]}" -m 0 -G 9 -E 2
for scores in "1 -3" "1 -2" "2 -3" "5 -4" "1 -1" ; do
    # shellcheck disable=SC2086
    set -- ${scores}
    compare_run "${NUCL[@]}" -m 8 -p 0 -e 1000 -r "${1}" -q "${2}"
done

# strands and genetic codes
for strand in 1 2 3 ; do
    compare_run "${NUCL[@]}" -m 0 -p 0 -S "${strand}"
    compare_run -d prot -i nucl_query.fasta -m 8 -p 2 -S "${strand}" -e 1000
done
for code in 1 2 4 11 23 ; do
    compare_run -d prot -i nucl_query.fasta -m 8 -p 2 -Q "${code}" -e 1000
    compare_run -d nucl -i prot_query.fasta -m 8 -p 3 -D "${code}" -e 1000
done

# taxids, dumps, output file, error paths
compare_run "${PROT[@]}" -m 0 -x taxids.txt -H
compare_run "${NUCL[@]}" -m 8 -p 0 -x taxids.txt
for dump in 1 2 ; do
    compare_run -d prot -N "${dump}"
    compare_run -d nucl -p 0 -N "${dump}"
    compare_run -d prot_split -N "${dump}"
done
compare_run "${PROT[@]}" -m 8 -o out.file
compare_run "${PROT[@]}" -m 0 -o out.file
compare_run -d missing -i prot_query.fasta
compare_run "${PROT[@]}" -M NO_SUCH_MATRIX
compare_run "${PROT[@]}" -G 1 -E 0 -m 8
compare_run -h
compare_run --version

echo "---"
echo "${comparisons} comparisons, ${differences} differences"
echo "${runs} runs: ${failed_runs} exited non-zero (2 are expected: the missing database and matrix), ${empty_runs} with an empty stdout"
[[ "${differences}" -eq 0 ]] && echo "BYTE-IDENTICAL" || echo "NOT IDENTICAL (files kept in ${WORKDIR})"
exit "${differences}"
