#!/bin/bash

### Clear previous logs
#rm /cms/ldap_home/jongho/Delphes-3.4.2/condor_logs/out/*
#rm /cms/ldap_home/jongho/Delphes-3.4.2/condor_logs/err/*
#rm /cms/ldap_home/jongho/Delphes-3.4.2/condor_logs/log/*

### ROOT environment
source /cms/ldap_home/jongho/root/bin/thisroot.sh

tempNum=$1

mkdir -p $_CONDOR_SCRATCH_DIR/job_${tempNum}

## Give the rest on job to avoid duplication in random seed
sleepTime=$(echo "${tempNum}*2" | bc)
sleep ${sleepTime}

cd /cms/ldap_home/jongho/Delphes-3.4.2 
source DelphesEnv.sh

mkdir -p /cms_scratch/jongho/job_${tempNum}

###############################################################################
## Running DelphesPythia8, in this step, you need detector card(.tcl) and Pythia configuration card(.cmnd)

## Generate QCD with 200 pile-up sample
#./DelphesPythia8 cards/CMS_PhaseII/New_CMS_PhaseII_PixelStudies_200PU.tcl examples/Pythia8/configNoLHE.cmnd Results2/job_${tempNum}/Delphes4trkiso.root 

## Generate SingleElectron with no pile-up sample
#./DelphesPythia8 cards/CMS_PhaseII/New_CMS_PhaseII_PixelStudies_noPU_my.tcl examples/Pythia8/configParticleGun.cmnd Results2/job_${tempNum}/Delphes4trkiso.root 
#./DelphesPythia8 cards/CMS_PhaseII/New_CMS_PhaseII_PixelStudies_noPU_my.tcl examples/Pythia8/configParticleGun_noPU.cmnd Results2/job_${tempNum}/Delphes4trkiso.root 

## Generate Minimum-bias with no pile-up sample, this sample is used to make pile-up merging file
#./DelphesPythia8 cards/CMS_PhaseII/New_CMS_PhaseII_PixelStudies_noPU_my.tcl examples/Pythia8/generatePileUpCMS.cmnd Results2/job_${tempNum}/Delphes4trkiso.root

## Generate SingleElectron with 200 pile-up smaple
#./DelphesPythia8 cards/CMS_PhaseII/New_CMS_PhaseII_PixelStudies_200PU_my.tcl examples/Pythia8/configParticleGun.cmnd $_CONDOR_SCRATCH_DIR/job_${tempNum}/Delphes4trkiso.root 

## Generate Minimum-bias with 200 pile-up sample
#./DelphesPythia8 cards/CMS_PhaseII/New_CMS_PhaseII_PixelStudies_200PU_my.tcl examples/Pythia8/generatePileUpCMS.cmnd $_CONDOR_SCRATCH_DIR/job_${tempNum}/Delphes4trkiso.root 
./DelphesPythia8 cards/CMS_PhaseII/New_CMS_PhaseII_PixelStudies_200PU_my.tcl examples/Pythia8/generatePileUpCMS.cmnd /cms_scratch/jongho/job_${tempNum}/Delphes4trkiso.root 

###############################################################################

#cd $_CONDOR_SCRATCH_DIR
cp Example35.C disksmear_Minbias.C disksmear.h /cms_scratch/jongho
cd /cms_scratch/jongho

#fileEx=$(ls $_CONDOR_SCRATCH_DIR/job_${tempNum} | grep "Delphes4trkiso.root")
fileEx=$(ls ${PWD}/job_${tempNum} | grep "Delphes4trkiso.root")
echo "\n Successfully produced ${fileEx}"

## Smear Generator level info, Tower object, Pixel barrel layer
cp Example35.C Example35_${tempNum}.C
sed -i "62s/result.root/result_${tempNum}.root/g" Example35_${tempNum}.C

root -l -b << EOF
.L Example35_${tempNum}.C
Example35("job_${tempNum}/Delphes4trkiso.root")
.q
EOF

#echo ""
#ls -l $_CONDOR_SCRATCH_DIR

rm job_${tempNum}/Delphes4trkiso.root
rm Example35_${tempNum}.C

## Smear Pixel foward disk
#cp disksmear_SE.C disksmear_${tempNum}.C  # for SingleElectron sample
cp disksmear_Minbias.C disksmear_${tempNum}.C  # for Minimum-bias sample
sed -i "2s/disksmear.h/disksmear_${tempNum}.h/g" disksmear_${tempNum}.C

cp disksmear.h disksmear_${tempNum}.h
sed -i "128,132s/result.root/result_${tempNum}.root/g" disksmear_${tempNum}.h
sed -i "138s/results.root/final_${tempNum}.root/g" disksmear_${tempNum}.h

#echo ""
#ls -l $_CONDOR_SCRATCH_DIR

root -l -b << EOF
.L disksmear_${tempNum}.C
disksmear a
a.Loop()
.q
EOF

#echo ""
#ls -l $_CONDOR_SCRATCH_DIR

rm disksmear_${tempNum}.C
rm disksmear_${tempNum}.h
#rm result_${tempNum}.root

# Copy the sample to xrootd directory using xrdcp.
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0/
eval `scramv1 runtime -sh`
cd -

#xrdcp -f final_${tempNum}.root root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/SE/0009
#xrdcp -f final_${tempNum}.root root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/Bkg/0019
xrdcp -f final_${tempNum}.root root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/Bkg/0035
#rm final_${tempNum}.root
