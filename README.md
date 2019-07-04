Level-1 pixel track trigger based on Delphes
============================================

Incomplete part will be filled soon.  
Standalone root is available.  

# Submit condor job at KISTI server

## STEP 0
Connect ui20 machine  
```
ssh -p 4280 username@ui20.sdfarm.kr
```

## STEP 1
Install [root v6.06.00](https://root.cern.ch/content/release-60600), [Delphes-3.4.1](https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/QuickTour) and [pythia8-v2.3.5](https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/Pythia8)

## STEP 2
Copy MinBias_100k.pileup from lxplus to KISTI Delphes directory     
`scp <username>@lxplus6.cern.ch:/eos/cms/store/group/upgrade/delphes/PhaseII/MinBias_100k.pileup .`

## STEP 3
Fork the git repository and do `git clone -b jongho git@github.com:<username>/Signalwindow.git`

## STEP 4
Edit the paths in DelphesCards/\*.tcl files, condor/setenv.sh

## STEP 5
`condor_submit submit.jds`

# Measure signal windows with delphes sample

## STEP 0
You can pass this step if you did git clone   
Fork the git repository and do `git clone -b jongho git@github.com:<username>/Signalwindow.git`

## STEP 1
Prepare root files with quantization using "Example9.C"  

## STEP 2
Let's say that the root file name is "results.root" (Strongly recommended!)  
```
cd SignalWindows   
root mySw.C   
mySw a   
a.Loop()   
```

## STEP 3 
```
cp sw.root fit_median
root Make2Dplots.C
Make2Dplots a
a.Loop(1)
```
You need to change the number in a.Loop(1) from 1 to 6 and repeatedly run the code.    
I recommend you do `mkdir -p plots` and move all png files into plots directory.  
And follow the step 9-1 to the end in master branch README.  



