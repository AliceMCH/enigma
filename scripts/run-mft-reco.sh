#! /bin/bash

export INFOLOGGER_MODE=stdout
export SCRIPTDIR=$(readlink -f $(dirname $0))
echo "SCRIPTDIR: ${SCRIPTDIR}"

export XRD_REQUESTTIMEOUT=1200


local_folder="$1"
condition_remapping_folder="$2"

CONDITION_REMAP=""
#CONDITION_REMAP="--condition-remap \"http://alice-ccdb.cern.ch/TestReco/MFTAlign2024/=MFT/Config/Geometry,GLO/Config/GeometryAligned\""
#CONDITION_REMAP="--condition-remap \"http://alice-ccdb.cern.ch/TestReco/MFTAlign2023/=MFT/Config/Geometry,GLO/Config/GeometryAligned\""
#CONDITION_REMAP="--condition-remap \"file://new-align/=MFT/Config/Geometry,GLO/Config/GeometryAligned\""

if [ -n "${condition_remapping_folder}" ]; then
    CONDITION_REMAP="--condition-remap \"file://${condition_remapping_folder}=MFT/Config/Geometry,GLO/Config/GeometryAligned\""
fi

#ROOTOUTPUT_OPT="--disable-root-output"
ROOTOUTPUT_OPT=""

#MFTRECOCONF="MFTTracking.RBins=30;MFTTracking.PhiBins=120;MFTTracking.ZVtxMin=-13;MFTTracking.ZVtxMax=13;MFTTracking.MFTRadLength=0.084;MFTClustererParam.maxBCDiffToMaskBias=-10;MFTClustererParam.maxBCDiffToSquashBias=10;MFTTracking.FullClusterScan=true;MFTTracking.irFramesOnly=0"
MFTRECOCONF="MFTTracking.RBins=30;MFTTracking.PhiBins=120;MFTTracking.ZVtxMin=-13;MFTTracking.ZVtxMax=13;MFTTracking.MFTRadLength=0.084;MFTClustererParam.maxBCDiffToMaskBias=-10;MFTClustererParam.maxBCDiffToSquashBias=10;MFTTracking.irFramesOnly=0"

TFIDINFO=${local_folder}/o2_tfidinfo.root
MFTCLUSTERS=${local_folder}/mftclusters.root
MFTTRACKSOUT=${local_folder}/mfttracks-realign.root

#ARGS_ALL="--session default --shm-segment-size 16000000000 --timeframes-rate-limit 10 --timeframes-rate-limit-ipcid 0  --hbfutils-config ${TFIDINFO},upstream --max-tf -1"
ARGS_ALL="--session default --shm-segment-size 16000000000 --hbfutils-config ${TFIDINFO},upstream"

WORKFLOW=""
WORKFLOW+="o2-reader-driver-workflow ${ARGS_ALL} | "
WORKFLOW+="o2-mft-cluster-reader-workflow ${ARGS_ALL} --mft-cluster-infile ${MFTCLUSTERS} | "
WORKFLOW+="o2-mft-reco-workflow ${ARGS_ALL} ${CONDITION_REMAP} --clusters-from-upstream --disable-mc --mft-cluster-writer \"--outfile /dev/null\" --mft-track-writer \"--outfile ${MFTTRACKSOUT}\" --nThreads 4 --pipeline mft-tracker:4 ${ROOTOUTPUT_OPT} --configKeyValues \"${MFTRECOCONF}\" | "

WORKFLOW+="o2-dpl-run --session default --shm-segment-size 16000000000 ${CONDITION_REMAP} --batch --run"

echo $WORKFLOW > ${local_folder}/workflow.sh
eval $WORKFLOW >& ${local_folder}/realign.txt

exit 0
