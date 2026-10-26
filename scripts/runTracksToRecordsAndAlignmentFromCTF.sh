#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

CTFLIST=$(realpath "$1")

if [ ! -e mft_mille_records.root ]; then
    echo "${scriptdir}/runTracksToRecords.sh \"$CTFLIST\""
    ${scriptdir}/runTracksToRecordsFromCTF.sh "$CTFLIST"

    #hadd -ff mft_mille_records.root ?/mft_mille_records-*.root
    #hadd -ff mft_align_point.root ?/mft_align_point-*.root
fi

if [ ! -e mft_mille_records.root ]; then
    echo "Missing mft_mille_records.root file, exiting"
    exit 1
fi


if [ ! -e mft_alignment.root ]; then
    echo "${scriptdir}/runAlign.sh"
    ${scriptdir}/runAlign.sh >& align.txt
fi

if [ ! -e mft_alignment.root ]; then
    echo "Missing mft_alignment.root file, exiting"
    exit 1
fi


echo "${scriptdir}/../macros/exportGeomWrtIdeal.C"
root -l -b -q ${scriptdir}/../macros/exportGeomWrtIdeal.C

if [ ! -e mft_alignment_wrt_ideal.root ]; then
    echo "Missing mft_alignment_wrt_ideal.root file, exiting"
    exit 1
fi

rm -rf temp
mkdir -p temp/MFT/Calib/Align
cp mft_alignment_wrt_ideal.root temp/MFT/Calib/Align/snapshot.root

(cd temp && o2-create-aligned-geometry-workflow  --condition-remap file://.=MFT/Calib/Align  --configKeyValues "HBFUtils.startTime=1770037200000" -b >& create-aligned-geometry.txt)
