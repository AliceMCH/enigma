#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

CTFLIST="$1"
OUTDIR=DCA

mkdir -p $OUTDIR
rm -f stop

NCTF=$(cat "${CTFLIST}" | wc -l)
I=1
while [ $I -le $NCTF ]; do

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

    CTF=$(cat "${CTFLIST}" | head -n $I | tail -n 1)

    alien_path="$CTF"
    run_number=$(basename $(dirname $(dirname $(dirname "${alien_path}"))))
    chunk_folder=$(basename $(dirname "${alien_path}"))
    tf_folder=$(basename "${alien_path}")
    local_folder=input/${run_number}/${chunk_folder}/${tf_folder}
    mkdir -p ${local_folder}

    delete_mfttracks=0
    delete_mftclusters=0
    delete_tfidinfo=0
    if [ -e realign ]; then
        # if the realignment folder exists, then run the MFT reconstruction
        echo "Reconstructing MFT tracks with new alignment"
        echo "${scriptdir}/run-mft-reco-ctf.sh ${CTF} ${local_folder} realign"
        ${scriptdir}/run-mft-reco-ctf.sh ${CTF} ${local_folder} realign >& "$OUTDIR"/mft-reco-${I}.txt
        delete_mfttracks=1
    fi

    #echo -n "Pausing before processing..."
    #read dummy

    echo "I=${I} => processing \"${DIR}\" ..."
    rm -f DCA.root
    echo "root -l -b -q \"${scriptdir}/../macros/runDCA.C(\\\"$DIR\\\")\""
    root -l -b -q "${scriptdir}/../macros/runDCA.C(\"$DIR\")" >& "$OUTDIR"/log-${I}.txt
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

    break
    
done

echo ""
echo "Processing finished"
