#!/bin/bash

delphesDIR=/cms/ldap_home/jongho/L1PixelTrigger/Delphes/Delphes-3.4.1 

source ${delphesDIR}/setenv.sh  
cd ${delphesDIR}

mkdir -p Results/job_$1
./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_my200PU.tcl examples/Pythia8/generatePileUp.cmnd Results/job_$1/SingleEle200PU.root

cp Example10.C Example10_$1.C
sed -i "718s/results.root/results_$1.root/g" Example10_$1.C

root -l -b << EOF
.L Example10_$1.C
Example10("Results/job_$1/SingleEle200PU.root")
.q
EOF

xrdcp file:///$PWD/results_$1.root root://cluster142.knu.ac.kr//store/user/jongho/MinBias200PU_Delphes/Test1/results_$1.root
#mv results_$1.root Results/job_$1/results_$1.root
rm Results/job_$1/SingleEle200PU.root
rm Example10_$1.C
rm results_$1.root
