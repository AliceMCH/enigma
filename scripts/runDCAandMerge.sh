#! /bin/bash

scriptdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")

${scriptdir}/runDCA.sh "$@"

echo ""
echo "Merging outputs"
hadd -f DCAFull.root DCA/DCA-*.root
echo ""
echo "Merging finished"

exit

echo ""
echo "Running DCA analysis"
$HOME/Workflows/MuonMatching/AliceGlobalMuonAnalysis/alignment/muonGlobalAlignmentVertexShift.sh DCAFull.root
echo ""
echo "DCA analysis finished"
