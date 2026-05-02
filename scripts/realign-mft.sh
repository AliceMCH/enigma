#! /bin/bash

export SCRIPTDIR=$(readlink -f $(dirname $0))
echo "SCRIPTDIR: ${SCRIPTDIR}"

rm -f stop-realign

ALIGN_FOLDER=$2

nfolders=$(cat "$1" | wc -l)

index=1
while [ $index -le $nfolders ]; do

    line=$(cat "$1" | sed -n ${index}p)
    
    #echo "$line"
    #index=$((index+1))
    #continue;
    
    alien_folder="$line"
    run_number=$(basename $(dirname $(dirname $(dirname "$line"))))
    chunk_folder=$(basename $(dirname "$line"))
    tf_folder=$(basename "$line")
    local_folder=input/${run_number}/${chunk_folder}/${tf_folder}
    echo "${index}/${nfolders} -> ${local_folder}"
    ls "${local_folder}"

    TFIDINFO=${local_folder}/o2_tfidinfo.root
    MFTCLUSTERS=${local_folder}/mftclusters-filtered.root
    MFTTRACKSOUT=${local_folder}/mfttracks-realign2.root

    if [ -e stop-realign ]; then
        break;
    fi

    if [ ! -e ${MFTCLUSTERS} ]; then
        index=$((index+1))
        continue;
    fi
    if [ ! -e ${TFIDINFO} ]; then
        index=$((index+1))
        continue;
    fi
    if [ -e ${MFTTRACKSOUT} ]; then
        index=$((index+1))
        continue;
    fi

    echo "$SCRIPTDIR/run-mft-reco.sh \"${local_folder}\" ${ALIGN_FOLDER}"
    $SCRIPTDIR/run-mft-reco.sh "${local_folder}" ${ALIGN_FOLDER} >& ${local_folder}/realign.txt

    echo "Realignment done -> ${local_folder}"
    echo ""
    index=$((index+1))

    #break
    
done

echo "Realignment of ${index}/${nfolders} folders finished"
