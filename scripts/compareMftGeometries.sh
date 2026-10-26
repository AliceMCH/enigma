#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

GEOM1="$1"
GEOM2="$2"

echo "root -l -b -q \"${scriptdir}/../macros/compareMftGeometries.C(\\\"$GEOM1\\\", \\\"$GEOM2\\\")\""
root -l -b -q "${scriptdir}/../macros/compareMftGeometries.C(\"$GEOM1\", \"$GEOM2\")"
