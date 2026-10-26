#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
macrosdir="${scriptdir}/../macros"

CTFLIST="$1"
CHUNKMIN=$2
CHUNKMAX=$3

USE_VERTEX=${USE_VERTEX:-0}


CHUNKID=$CHUNKMIN
while [ $CHUNKID -le $CHUNKMAX ]; do

    if [ -e ../stop ]; then break; fi

    rm -f mft_mille_records.root mft_align_point.root
    
    CTF=$(cat "$CTFLIST" | sed -n $((CHUNKID+1))p)

    alien_path="$CTF"

    alien_path=$(cat "$CTFLIST" | sed -n $((CHUNKID+1))p)
    run_number=$(basename $(dirname $(dirname $(dirname "${alien_path}"))))
    chunk_folder=$(basename $(dirname "${alien_path}"))
    tf_folder=$(basename "${alien_path}")
    local_folder=input/${run_number}/${chunk_folder}/${tf_folder}
    mkdir -p ${local_folder}

    delete_vertex=0
    delete_mftclusters=0
    delete_tfidinfo=0
    delete_mfttracks=0
    if [ -e realign ]; then
        # if the realignment folder exists, then run the MFT reconstruction
        echo "Reconstructing MFT tracks with new alignment"
        echo "${scriptdir}/run-mft-reco-ctf.sh ${CTF} ${local_folder} realign"
        ${scriptdir}/run-mft-reco-ctf.sh ${CTF} ${local_folder} realign
        delete_mfttracks=1
    else
        echo "Reconstructing MFT tracks with standard alignment"
        echo "${scriptdir}/run-mft-reco-ctf.sh ${CTF} ${local_folder}"
        ${scriptdir}/run-mft-reco-ctf.sh ${CTF} ${local_folder}
        delete_mfttracks=1
    fi

    if [ -e ../stop ]; then break; fi

    root -l -b -q "${macrosdir}/runTracksToRecords.C($CHUNKID, $CHUNKID, ${USE_VERTEX})"
    rm -f mft_mille_records-$CHUNKID.root
    mv mft_mille_records.root mft_mille_records-$CHUNKID.root
    #cp mft_align_point.root mft_align_point-$CHUNKID.root
    rm mft_align_point.root

    if [ x"${delete_vertex}" = "x1" ]; then
        rm -f ${local_folder}/o2_primary_vertex.root
    fi

    if [ x"${delete_tfidinfo}" = "x1" ]; then
        rm -f ${local_folder}/o2_tfidinfo.root
    fi

    if [ x"${delete_mftclusters}" = "x1" ]; then
        rm -f ${local_folder}/mftclusters.root
    fi

    if [ x"${delete_mfttracks}" = "x1" ]; then
        rm -f ${local_folder}/mfttracks.root
    fi

    CHUNKID=$((CHUNKID+1))

done

hadd -ff mft_mille_records_full.root mft_mille_records-*.root
#hadd -ff mft_align_point_full.root mft_align_point-*.root
