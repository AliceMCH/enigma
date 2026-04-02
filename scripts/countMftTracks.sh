#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
macrosdir="${scriptdir}/../macros"

root -l -b -q "${macrosdir}/countMftTracks.C(\"$1\")"
