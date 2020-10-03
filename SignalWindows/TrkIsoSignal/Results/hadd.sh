#!/bin/bash

basedir=/cms/ldap_home/jongho/TrkIso/TrkIsoSignal/Results

cd ${HOME}/CMSSW_9_3_7/src
eval `scramv1 runtime -sh`
cd -

rm sig_final.root

for d in SE* ; 
do
    echo "$d/output/Tree"
    cd $d/output/Tree
    hadd Tree_SE_PU200.root *.root
    mv Tree_SE_PU200.root ../../
    cd ${basedir}
done

hadd -f sig_final.root */*.root
