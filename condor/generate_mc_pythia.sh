#!/bin/bash

delphesDIR=/cms/ldap_home/jongho/L1PixelTrigger/Delphes/Delphes-3.4.1 

source ${delphesDIR}/setenv.sh  
cd ${delphesDIR}

mkdir -p Results/job_$1
./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_my.tcl examples/Pythia8/configParticleGun.cmnd Results/job_$1/SingleEle200PU.root ## Single Ele NOPU
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_my200PU.tcl examples/Pythia8/configParticleGun.cmnd Results/job_$1/SingleEle200PU.root ## Single Ele 200PU
#./DelphesPythia8 cards/CMS_PhaseII/CMS_PhaseII_PixelStudies_my200PU.tcl examples/Pythia8/generatePileUp.cmnd Results/job_$1/SingleEle200PU.root ## Minbias 200PU

cp Example11.C Example11_$1.C
sed -i "647s/results.root/results_$1.root/g" Example11_$1.C

root -l -b << EOF
.L Example11_$1.C
Example11("Results/job_$1/SingleEle200PU.root")
.q
EOF

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0/
eval `scramv1 runtime -sh`
cd -

xrdcp file:///$PWD/results_$1.root root://cluster142.knu.ac.kr//store/user/jongho/MinBias200PU_Delphes/Test1/results_$1.root
#mv results_$1.root Results/job_$1/results_$1.root
rm Results/job_$1/SingleEle200PU.root
rm Example11_$1.C
rm results_$1.root

