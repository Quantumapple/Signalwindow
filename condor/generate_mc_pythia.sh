#!/bin/bash

delphesDIR=/cms/ldap_home/jongho/L1PixelTrigger/Delphes/Delphes-3.4.1 

source ${delphesDIR}/setenv.sh  
cd ${delphesDIR}

mkdir -p Results/job_$1
./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_my200PU.tcl examples/Pythia8/generatePileUp.cmnd Results/job_$1/SingleEle200PU.root

cp Example9.C Example9_$1.C
sed -i "718s/results.root/results_$1.root/g" Example9_$1.C

root -l -b << EOF
.L Example9_$1.C
Example9("Results/job_$1/SingleEle200PU.root")
.q
EOF

mv results_$1.root Results/job_$1/results_$1.root
rm Results/job_$1/SingleEle200PU.root
rm Example9_$1.C
rm job_$1.err job_$1.out
