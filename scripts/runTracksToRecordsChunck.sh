#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
macrosdir="${scriptdir}/../macros"

#valgrind --tool=memcheck \
root -l -b -q "${macrosdir}/runTracksToRecords.C($1, $2)"
