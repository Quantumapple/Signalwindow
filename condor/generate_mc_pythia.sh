#!/bin/bash

delphesDIR=/cms/scratch/jongho/Delphes/Delphes-3.4.1

source ${delphesDIR}/setenv.sh  
cd ${delphesDIR}

#tar -axvf submit.tar 

tempNum=$1
mkdir -p Results/job_${tempNum}
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_myV3.tcl examples/Pythia8/configParticleGun.cmnd Results/job_$1/SingleEle200PU.root ## Single Ele NOPU
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_my200PU.tcl examples/Pythia8/configParticleGun.cmnd Results/job_$1/SingleEle200PU.root ## Single Ele 200PU
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_EcalUpdate200PU.tcl examples/Pythia8/configNoLHE.cmnd Results/job_${tempNum}/SingleEle200PU.root ## Minbias 200PU
./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_190727_200PU.tcl examples/Pythia8/generatePileUpCMS.cmnd Results/job_${tempNum}/SingleEle200PU.root 
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_myV3.tcl examples/Pythia8/configLHE.cmnd Results/job_${tempNum}/SingleEle200PU.root ## ZToNuNu 200PU

cp Example18.C Example18_${tempNum}.C
sed -i "587s/results.root/results_${tempNum}.root/g" Example18_${tempNum}.C

root -l -b << EOF
.L Example18_${tempNum}.C
Example18("Results/job_${tempNum}/SingleEle200PU.root")
.q
EOF

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0/
eval `scramv1 runtime -sh`
cd -

#cp x509up_u* /tmp
#voms-proxy-info -all

#xrdcp results_${tempNum}.root root://cluster142.knu.ac.kr//store/user/jongho/MinBias200PU_Delphes/0000
xrdcp -f results_${tempNum}.root root://cms-xrdr.private.lo:2094//xrd/store/user/jongho/Minbias/0000
#xrdcp -f results_${tempNum}.root root://cms-xrdr.private.lo:2094//xrd/store/user/jongho/SingleEle/0000
#mv results_${tempNum}.root Results/job_${tempNum}/results_${tempNum}.root
rm Results/job_${tempNum}/SingleEle200PU.root
rm Example18_${tempNum}.C
rm results_${tempNum}.root

