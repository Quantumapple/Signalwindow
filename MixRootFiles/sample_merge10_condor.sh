#!/bin/bash

dirn=0000
maxfile=$(ls /xrootd/store/user/jongho/Delphes/SE/${dirn}/*.root | wc -l)
initfile=0

echo "Proceeding with ${maxfile} root files."
echo ""

tmp1=Results
tmp2=process
tmp3=output
tmp4=condorLogs/out
tmp5=condorLogs/err
tmp6=condorLogs/log

if [ -d ${tmp1} ]; then rm ${tmp1}/*; fi
if [ -d ${tmp2} ]; then rm ${tmp2}/*.log; fi
if [ -d ${tmp3} ]; then rm ${tmp3}/*.root; fi
if [ ! -d ${tmp1} ]; then mkdir -p $tmp1; fi
if [ ! -d ${tmp2} ]; then mkdir -p $tmp2; fi
if [ ! -d ${tmp3} ]; then mkdir -p $tmp3; fi
if [ ! -d ${tmp4} ]; then mkdir -p $tmp4; fi
if [ ! -d ${tmp5} ]; then mkdir -p $tmp5; fi
if [ ! -d ${tmp6} ]; then mkdir -p $tmp6; fi

for tempNum in $(seq $initfile $maxfile); do
cat << EOF > merge_${tempNum}.sh
#!/bin/bash

cd /cms/ldap_home/jongho/CMSSW_9_3_7/src
eval \`scramv1 runtime -sh\`

cd /cms_scratch/jongho/mix 
xrdcp -f root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/SE/${dirn}/final_${tempNum}.root . 

cp mix13.cc mix13_${tempNum}.cc
sed -i "24s/mix13/mix13_${tempNum}/g" mix13_${tempNum}.cc
sed -i "27s/input/final_${tempNum}/g" mix13_${tempNum}.cc
sed -i "33s/merged_0.root/merged_${tempNum}.root/g" mix13_${tempNum}.cc

root -l -q -b mix13_${tempNum}.cc > process/process_${tempNum}.log

xrdcp -f output/merged_${tempNum}.root root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/mergeSE/${dirn}/merged_${tempNum}.root

mv mix13_${tempNum}.cc ${tmp1}/ 
mv job_${tempNum}.jds merge_${tempNum}.sh ${tmp1}/ 
rm final_${tempNum}.root
EOF

chmod 777 merge_${tempNum}.sh

cat << EOF > job_${tempNum}.jds
executable = merge_${tempNum}.sh
universe   = vanilla
output     = condorLogs/out/condorOut_${tempNum}.out
error      = condorLogs/err/condorErr_${tempNum}.err
log        = condorLogs/log/condorLog_${tempNum}.log
getenv     = True
transfer_input_files = mix13.cc 
accounting_group=group_cms
when_to_transfer_output = ON_EXIT
requirements = (Machine=!="cms-t3-wn3018.sdfarm.kr") && (Machine=!="cms-t3-wn3021.sdfarm.kr") && (Machine=!="bio-wn3006.sdfarm.kr") && (Machine=!="bio-wn3015.sdfarm.kr") && (Machine=!="bio-wn3013.sdfarm.kr")
request_Cpus=2
JobBatchName = Delphes_merging
queue 1
EOF

#condor_submit job_${tempNum}.jds

done
