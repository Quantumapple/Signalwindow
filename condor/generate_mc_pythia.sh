#!/bin/bash

delphesDIR=/cms/scratch/jongho/Delphes/delphes

mkdir -p Results err out log

source ${delphesDIR}/setenv.sh  
cd ${delphesDIR}

tempNum=$1
sleepTime=$(echo "${tempNum}*7" | bc)

mkdir -p Results/job_${tempNum}

sleep ${sleepTime}
#./DelphesPythia8 cards/CMS_PhaseII/test.tcl examples/Pythia8/configParticleGun.cmnd Results/job_${tempNum}/SingleEle200PU.root 
./DelphesPythia8 cards/CMS_PhaseII/test2.tcl examples/Pythia8/generatePileUpCMS.cmnd Results/job_${tempNum}/SingleEle200PU.root 
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_Final200PU.tcl examples/Pythia8/configParticleGun.cmnd Results/job_${tempNum}/SingleEle200PU.root 
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_Final200PU.tcl examples/Pythia8/generatePileUpCMS.cmnd Results/job_${tempNum}/SingleEle200PU.root 

cp Example25.C Example25_${tempNum}.C
sed -i "34s/results.root/results_${tempNum}.root/g" Example25_${tempNum}.C

root -l -b << EOF
.L Example25_${tempNum}.C
Example25("Results/job_${tempNum}/SingleEle200PU.root")
.q
EOF

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0/
eval `scramv1 runtime -sh`
cd -

rm Results/job_${tempNum}/SingleEle200PU.root
xrdcp -f results_${tempNum}.root root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Minbias/0029
#xrdcp -f results_${tempNum}.root root://cms-xrdr.private.lo:2094//xrd/store/user/jongho/SingleEle200PU/0000
#mv results_${tempNum}.root Results/job_${tempNum}/results_${tempNum}.root
rm Example25_${tempNum}.C
rm results_${tempNum}.root

