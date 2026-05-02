#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
macrosdir="${scriptdir}/../macros"

DIRLIST="$1"
CHUNKMIN=$2
CHUNKMAX=$3

USE_VERTEX=${USE_VERTEX:-0}


CHUNKID=$CHUNKMIN
while [ $CHUNKID -le $CHUNKMAX ]; do

    if [ -e ../stop ]; then break; fi

    rm -f mft_mille_records.root mft_align_point.root

    alien_folder=$(cat "$DIRLIST" | sed -n $((CHUNKID+1))p)
    run_number=$(basename $(dirname $(dirname $(dirname "${alien_folder}"))))
    chunk_folder=$(basename $(dirname "${alien_folder}"))
    tf_folder=$(basename "${alien_folder}")
    local_folder=input/${run_number}/${chunk_folder}/${tf_folder}

    delete_vertex=0
    if [ x"${USE_VERTEX}" = "x1" ]; then
        if [ ! -e "${local_folder}/o2_primary_vertex.root" ]; then
            rootfile=o2_primary_vertex.root
            alien_cp alien://${alien_folder}/${rootfile} file://./${local_folder}/${rootfile}
            delete_vertex=1
        fi
    fi

    delete_mftclusters=0
    if [ ! -e "${local_folder}/mftclusters.root" ]; then
        rootfile=mftclusters.root
        alien_cp alien://${alien_folder}/${rootfile} file://./${local_folder}/${rootfile}
        delete_mftclusters=1
    fi

    delete_tfidinfo=0
    if [ ! -e "${local_folder}/o2_tfidinfo.root" ]; then
        rootfile=o2_tfidinfo.root
        alien_cp alien://${alien_folder}/${rootfile} file://./${local_folder}/${rootfile}
        delete_tfidinfo=1
    fi

    delete_mfttracks=0
    if [ -e realign ]; then
        # if the realignment folder exists, then run the MFT reconstruction
        echo "Reconstructing MFT tracks with new alignment"
        echo "${scriptdir}/run-mft-reco.sh ${local_folder} realign"
        ${scriptdir}/run-mft-reco.sh ${local_folder} realign
        delete_mfttracks=1
    else
        if [ ! -e "${local_folder}/mfttracks.root" ]; then
            rootfile=mfttracks.root
            alien_cp alien://${alien_folder}/${rootfile} file://./${local_folder}/${rootfile}
            delete_mfttracks=1
        fi
    fi

    if [ -e ../stop ]; then break; fi

    root -l -b -q "${macrosdir}/runTracksToRecords.C($CHUNKID, $CHUNKID, ${USE_VERTEX})"
    cp mft_mille_records.root mft_mille_records-$CHUNKID.root
    cp mft_align_point.root mft_align_point-$CHUNKID.root

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
