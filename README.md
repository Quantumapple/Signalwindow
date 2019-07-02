Level-1 pixel track trigger based on Delphes
============================================

Incomplete part will be filled soon.  
Standalone root is available.  

## STEP 0
Fork the git repository and do `git clone -b jongho git@github.com:<username>/Signalwindow.git`

## STEP 1
Prepare root files with quantization using "Example7.C"  

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



