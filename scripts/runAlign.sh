#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

echo "${scriptdir}/../macros/runAlign.C"
root -l -b -q ${scriptdir}/../macros/runAlign.C 2>&1 | stdbuf -oL grep -v "Failed to solve locals"
