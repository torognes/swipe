#!/bin/sh
# Print swipe's version, after checking that every place carrying it agrees.
#
#   usage: version.sh [tag]
#
# swipe has no single source for its version (phase 6a of the
# maintenance plan), so a release has to keep these in step:
#
#   swipe.h        SWIPE_VERSION, the version "swipe -h" prints
#   CHANGES        the first entry ("* Version X.Y.Z"), which must no
#                  longer be marked "(in development)"
#   CITATION.cff   what Zenodo and citation managers read (checked when
#                  present: it arrives with GitHub PR #39)
#   the git tag    what the release page is built from
#
# Run from the top of the source tree. The tag is only checked when
# given; both "v2.1.2" and "2.1.2" are accepted. The version is written
# to GITHUB_OUTPUT when that variable is set, so that the workflow can
# hand it to the other jobs. Modeled on swarm's .github/scripts/version.sh.
set -eu

fail() { echo "::error::$*" >&2 ; exit 1 ; }

header=$(sed -n 's/^#define SWIPE_VERSION "\([^"]*\)".*/\1/p' swipe.h)
changes_line=$(grep -m 1 '^[[:space:]]*\* Version ' CHANGES || true)
changes=$(printf '%s\n' "${changes_line}" | sed -n 's/^[[:space:]]*\* Version \([^ ]*\).*/\1/p')

test -n "${header}" || fail "cannot read SWIPE_VERSION from swipe.h"
test -n "${changes}" || fail "cannot read the first version entry of CHANGES"

echo "swipe.h:      ${header}" >&2
echo "CHANGES:      ${changes}" >&2

test "${header}" = "${changes}" || fail "CHANGES says ${changes}, swipe.h says ${header}"
case "${changes_line}" in
  *"in development"*) fail "the CHANGES entry of ${changes} is still marked (in development)" ;;
esac

if [ -f CITATION.cff ] ; then
  cff=$(sed -n 's/^version: *\([^ ]*\).*/\1/p' CITATION.cff | tr -d '\r"')
  test -n "${cff}" || fail "cannot read the version from CITATION.cff"
  echo "CITATION.cff: ${cff}" >&2
  test "${header}" = "${cff}" || fail "CITATION.cff says ${cff}, swipe.h says ${header}"
else
  echo "::warning::no CITATION.cff (GitHub PR #39): not checked" >&2
fi

# Accept the tag with or without its leading "v": releases are tagged
# v2.1.2, but a manual run may well be given the bare version.
tag="${1:-}"
if [ -n "${tag}" ] ; then
  echo "tag:          ${tag}" >&2
  test "${tag#v}" = "${header}" || fail "tag ${tag} does not match version ${header}"
fi

test -z "${GITHUB_OUTPUT:-}" || echo "version=${header}" >> "${GITHUB_OUTPUT}"
echo "${header}"
