#!/bin/bash

delphesDIR=/cms/scratch/jongho/Delphes/Delphes-3.4.1

source ${delphesDIR}/setenv.sh  
cd ${delphesDIR}

#tar -axvf submit.tar 

tempNum=$1
mkdir -p Results/job_${tempNum}
./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_Final200PU.tcl examples/Pythia8/configParticleGun.cmnd Results/job_${tempNum}/SingleEle200PU.root 
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_Final200PU.tcl examples/Pythia8/generatePileUpCMS.cmnd Results/job_${tempNum}/SingleEle200PU.root 

sleep 3

cp Example21.C Example21_${tempNum}.C
sed -i "698s/results.root/results_${tempNum}.root/g" Example21_${tempNum}.C

sleep 3

root -l -b << EOF
.L Example21_${tempNum}.C
Example21("Results/job_${tempNum}/SingleEle200PU.root")
.q
EOF

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0/
eval `scramv1 runtime -sh`
cd -

#cp x509up_u* /tmp
#voms-proxy-info -all

#xrdcp -f results_${tempNum}.root root://cms-xrdr.private.lo:2094//xrd/store/user/jongho/Minbias/0000
xrdcp -f results_${tempNum}.root root://cms-xrdr.private.lo:2094//xrd/store/user/jongho/SingleEle200PU/0000
#mv results_${tempNum}.root Results/job_${tempNum}/results_${tempNum}.root
rm Results/job_${tempNum}/SingleEle200PU.root
rm Example21_${tempNum}.C
rm results_${tempNum}.root

