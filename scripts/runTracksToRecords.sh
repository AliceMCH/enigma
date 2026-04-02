#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

NJOBS=5
NCHUNKS=$(find input -name "mfttracks.root" | wc -l)

NCHUNKSPERJOB=$((NCHUNKS/NJOBS+1))

echo "NJOBS:         $NJOBS"
echo "NCHUNKS:       $NCHUNKS"
echo "NCHUNKSPERJOB: $NCHUNKSPERJOB"

ID=0
CHUNKID=0
while [ $CHUNKID -lt $NCHUNKS ]; do

    (rm -rf $ID && \
     mkdir -p $ID && \
     cd $ID && \
     ln -s ../input ../o2sim_geometry.root ../o2sim_geometry-aligned.root ../o2_mft_dictionary.root . && \
     ${scriptdir}/runTracksToRecordsChunck.sh $CHUNKID $((CHUNKID+NCHUNKSPERJOB-1)) >& tracksToRecords.txt)&

    ID=$((ID+1))
    CHUNKID=$((CHUNKID+NCHUNKSPERJOB))
    
    #break
    
done

jobs -p
wait $(jobs -p)
