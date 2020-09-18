#!/bin/bash

maxfile=499
initfile=0
fileName=0

tmp1=Results
tmp2=process
tmp3=output
tmp4=condorLogs/out
tmp5=condorLogs/err
tmp6=condorLogs/log

if [ -d ${tmp3} ]; then rm ${tmp3}/*.root; fi
if [ ! -d ${tmp1} ]; then mkdir -p $tmp1; fi
if [ ! -d ${tmp2} ]; then mkdir -p $tmp2; fi
if [ ! -d ${tmp3} ]; then mkdir -p $tmp3; fi
if [ ! -d ${tmp4} ]; then mkdir -p $tmp4; fi
if [ ! -d ${tmp5} ]; then mkdir -p $tmp5; fi
if [ ! -d ${tmp6} ]; then mkdir -p $tmp6; fi

for tempNum in $(seq $initfile $maxfile); do
cat << EOF > merge_${fileName}_${tempNum}.sh
#!/bin/bash

cd /cms/ldap_home/jongho/CMSSW_9_3_7/src
eval \`scramv1 runtime -sh\`

cd /cms_scratch/jongho/merging
xrdcp -f root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/SE/0002/final_${tempNum}.root . 

cp mix10.cc mix10_${fileName}_${tempNum}.cc
sed -i "543s/mix10/mix10_${fileName}_${tempNum}/g" mix10_${fileName}_${tempNum}.cc
sed -i "547s/input/final_${tempNum}/g" mix10_${fileName}_${tempNum}.cc
sed -i "561s/merged_0.root/merged_${tempNum}.root/g" mix10_${fileName}_${tempNum}.cc

root -l -q -b mix10_${fileName}_${tempNum}.cc > process/process_${fileName}_${tempNum}.log
#.L mix10_${fileName}_${tempNum}.cc
#mix10_${fileName}_${tempNum}()
#.q

mv mix10_${fileName}_${tempNum}.cc Results
mv job_${fileName}_${tempNum}.jds merge_${fileName}_${tempNum}.sh Results
rm final_${tempNum}.root
xrdcp -f output/merged_${tempNum}.root root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/mergeSE/0002/merged_${tempNum}.root
EOF

chmod 777 merge_${fileName}_${tempNum}.sh

cat << EOF > job_${fileName}_${tempNum}.jds
executable = merge_${fileName}_${tempNum}.sh
universe   = vanilla
output     = condorLogs/out/condorOut_${fileName}_${tempNum}.out
error      = condorLogs/err/condorErr_${fileName}_${tempNum}.err
log        = condorLogs/log/condorLog_${fileName}_${tempNum}.log
getenv     = True
transfer_input_files = mix10.cc 
accounting_group=group_cms
when_to_transfer_output = ON_EXIT
requirements = (Machine=!="cms-t3-wn3005.sdfarm.kr") && (Machine=!="cms-t3-wn3008.sdfarm.kr") && (Machine =!= "cms-t3-wn3013.sdfarm.kr") && (Machine =!= "cms-t3-wn3010.sdfarm.kr") && (Machine =!= "cms-t3-wn3017.sdfarm.kr") && (Machine =!= "kiaf-wn1001.sdfarm.kr")
request_Cpus=2
JobBatchName = Delphes_merging
queue 1
EOF

condor_submit job_${fileName}_${tempNum}.jds
#nohup ./sample_merge10_${fileName}_${tempNum}.sh >/dev/null 2>&1 &

done
