#!/bin/bash

#baseDIR=/cms/ldap_home/jongho/L1PixelTrigger/Delphes

source /cms/ldap_home/jongho/L1PixelTrigger/Delphes/root/bin/thisroot.sh
export PYTHONPATH=`pwd`/python:$PYTHONPATH
export LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH
#echo "This is for standalone ROOT 6.06.00"
#source ${baseDIR}/root/bin/thisroot.sh

#rootvr=`root-config --version`
#echo "Root version ${rootvr} is sourced!"
