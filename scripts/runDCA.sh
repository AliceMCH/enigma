#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

DIRLIST="$1"
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

    echo "I=${I} => processing \"${DIR}\" ..."
    rm -f DCA.root
    root -l -b -q "${scriptdir}/../macros/runDCA.C(\"$DIR\")" >& "$OUTDIR"/log-${I}.txt
    gzip -f "$OUTDIR"/log-${I}.txt
    echo "... done."


    if [ -e DCA.root ]; then
        touch "$OUTDIR"/done-${I}
        cp DCA.root $OUTDIR/DCA-$I.root
    fi

    I=$((I+1))
    
done

echo ""
echo "Processing finished"
