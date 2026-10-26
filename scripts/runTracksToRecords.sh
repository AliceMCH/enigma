#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

DIRLIST=$(realpath "$1")

NJOBS=10
NCHUNKS=$(cat "$DIRLIST" | wc -l)

NCHUNKSPERJOB=0
if [ $NCHUNKS -le $NJOBS ]; then
    NCHUNKSPERJOB=1
else
    NCHUNKSPERJOB=$((NCHUNKS/NJOBS))
    if [ $((NCHUNKSPERJOB*NJOBS)) -lt $NCHUNKS ]; then
        NCHUNKSPERJOB=$((NCHUNKSPERJOB+1))
    fi
fi

echo "NJOBS:         $NJOBS"
echo "NCHUNKS:       $NCHUNKS"
echo "NCHUNKSPERJOB: $NCHUNKSPERJOB"

# create input folders if not exsiting yet

while IFS= read -r line; do
    #echo "$line"
    run_number=$(basename $(dirname $(dirname $(dirname "$line"))))
    chunk_folder=$(basename $(dirname "$line"))
    tf_folder=$(basename "$line")
    local_folder=input/${run_number}/${chunk_folder}/${tf_folder}

    mkdir -p "${local_folder}"
done < "$DIRLIST"

rm -f stop

# process the tracks

ID=0
CHUNKID=0
while [ $CHUNKID -lt $NCHUNKS ]; do

    mkdir -p $ID && cd $ID
    ln -s ../input ../o2sim_geometry.root ../o2sim_geometry-aligned.root ../o2_mft_dictionary.root .
    if [ -e ../realign ]; then ln -s ../realign .; fi
    
    CHUNKID2=$((CHUNKID+NCHUNKSPERJOB-1))
    if [ $CHUNKID2 -ge $NCHUNKS ]; then
        CHUNKID2=$((NCHUNKS-1))
    fi

    echo "${scriptdir}/runTracksToRecordsChunck.sh \"$DIRLIST\" $CHUNKID $CHUNKID2"
    ${scriptdir}/runTracksToRecordsChunck.sh "$DIRLIST" $CHUNKID $CHUNKID2 >& tracksToRecords.txt&

    cd ..

    ID=$((ID+1))
    CHUNKID=$((CHUNKID+NCHUNKSPERJOB))
    
    #break
    
done

jobs -p
wait $(jobs -p)

hadd -ff mft_mille_records.root ?/mft_mille_records_full.root
#hadd -ff mft_align_point.root ?/mft_align_point_full.root
