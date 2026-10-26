#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

DIRLIST="$1"
VZSHIFTOPT=""
if [ -n "$2" ]; then VZSHIFTOPT=", $2"; fi
OUTDIR=DCA

mkdir -p $OUTDIR
rm -f stop

NDIRS=$(cat "${DIRLIST}" | wc -l)
I=1
while [ $I -le $NDIRS ]; do

    if [ -e stop ]; then
        break;
    fi

    if [ -e "$OUTDIR"/DCA-${I}.root ]; then
        I=$((I+1))
        continue
    fi
    if [ -e "$OUTDIR"/done-${I} ]; then
        I=$((I+1))
        continue
    fi

    DIR=$(cat "${DIRLIST}" | head -n $I | tail -n 1)
    echo ""
    echo "I=${I} => processing \"${DIR}\" ..."

    alien_folder="$DIR"
    run_number=$(basename $(dirname $(dirname $(dirname "${alien_folder}"))))
    chunk_folder=$(basename $(dirname "${alien_folder}"))
    tf_folder=$(basename "${alien_folder}")
    local_folder=input/${run_number}/${chunk_folder}/${tf_folder}
    mkdir -p ${local_folder}

    delete_mfttracks=0
    delete_mftclusters=0
    delete_tfidinfo=0
    if [ -e realign ]; then
        if [ ! -e "${local_folder}/mftclusters.root" ]; then
            rootfile=mftclusters.root
            alien_cp alien://${alien_folder}/${rootfile} file://./${local_folder}/${rootfile}
            echo "Copied ${alien_folder}/${rootfile} to ${local_folder}/${rootfile}"
            delete_mftclusters=1
        fi
        
        if [ ! -e "${local_folder}/o2_tfidinfo.root" ]; then
            rootfile=o2_tfidinfo.root
            alien_cp alien://${alien_folder}/${rootfile} file://./${local_folder}/${rootfile}
            echo "Copied ${alien_folder}/${rootfile} to ${local_folder}/${rootfile}"
            delete_tfidinfo=1
        fi

        # if the realignment folder exists, then run the MFT reconstruction
        echo "Reconstructing MFT tracks with new alignment"
        echo "${scriptdir}/run-mft-reco.sh ${local_folder} realign"
        ${scriptdir}/run-mft-reco.sh ${local_folder} realign >& "$OUTDIR"/mft-reco-${I}.txt
        gzip -f "$OUTDIR"/mft-reco-${I}.txt
        delete_mfttracks=1
    fi

    #echo -n "Pausing before processing..."
    #read dummy

    rm -f DCA.root
    echo "root -l -b -q \"${scriptdir}/../macros/runDCA.C(\\\"$DIR\\\" ${VZSHIFTOPT})\""
    root -l -b -q "${scriptdir}/../macros/runDCA.C(\"$DIR\" ${VZSHIFTOPT})" >& "$OUTDIR"/log-${I}.txt
    gzip -f "$OUTDIR"/log-${I}.txt
    echo "... done."


    if [ -e DCA.root ]; then
        touch "$OUTDIR"/done-${I}
        cp DCA.root $OUTDIR/DCA-$I.root
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


    I=$((I+1))

    #break
    
done

echo ""
echo "Processing finished"
