#!/bin/bash

for i in 0000 0001 0002 0003 0004 0005 0006 0007 0008 0009 0010 0011 0012 0013 0014 0015 0016 0017 0018 0019 0020 0021 0022 0023 0024
do
    #num=$(ls /xrootd/store/user/jongho/Delphes/mergeSE/${i}/*.root | wc -l)
    #echo ${i}, $((${num}/4 + 4)) 
    echo python submit_control.py ${i} 
    python submit_control.py ${i}
    echo ""
done
