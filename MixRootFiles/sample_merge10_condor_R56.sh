#!/bin/bash

dirn=0020
maxfile=999
initfile=0

echo "Proceeding with $((${maxfile}+1)) root files."
echo ""

storage=/xrootd/store/user/jongho/Delphes/SE
workdir=`pwd`

tmp4=condorLogs/out
tmp5=condorLogs/err
tmp6=condorLogs/log

mkdir -p ./${tmp4}
mkdir -p ./${tmp5}
mkdir -p ./${tmp6}

for tempNum in $(seq $initfile $maxfile); do

    mkdir -p ./job_${tempNum}
    cp mix13_r56.cc ./job_${tempNum}/mix13_r56_${tempNum}.cc

cat << EOF > ./job_${tempNum}/merge_${tempNum}.sh
#!/bin/bash

cd /cms/ldap_home/jongho/CMSSW_9_3_7/src
eval \`scramv1 runtime -sh\`

cd ${workdir}/job_${tempNum} 
mkdir -p ./dummy

##### Check the root file exists or not #####
xrdfs cms-xrdr.private.lo:2094 stat /xrd/store/user/jongho/Delphes/SE/${dirn}/final_${tempNum}.root 2>&1
XRDEXIT=\$?

if [[ \${XRDEXIT} -eq 0 ]]; then 
    
    xrdcp -f root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/SE/${dirn}/final_${tempNum}.root ./final_${tempNum}.root
    
    sed -i "24s/mix13_r56/mix13_r56_${tempNum}/g" mix13_r56_${tempNum}.cc
    sed -i "27s/input/final_${tempNum}/g" mix13_r56_${tempNum}.cc
    sed -i "33s/merged_0.root/merged_${tempNum}.root/g" mix13_r56_${tempNum}.cc
    
    root -l -q -b mix13_r56_${tempNum}.cc >& process_${tempNum}.log

    xrdcp -f merged_${tempNum}.root root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/mergeSE/${dirn}/merged_${tempNum}.root

    mv merge_${tempNum}.sh merged_${tempNum}.root ./dummy
    mv ../job_${tempNum}.jds ./dummy
    rm final_${tempNum}.root
else
    echo "File final_${tempNum}.root doesn't exist."
    echo "Exit the job"
    mv merge_${tempNum}.sh ./dummy
    mv ../job_${tempNum}.jds ./dummy

    exit 1
fi
EOF

chmod 777 ./job_${tempNum}/merge_${tempNum}.sh 

cat << EOF > job_${tempNum}.jds
executable = job_${tempNum}/merge_${tempNum}.sh
universe   = vanilla
output     = ${workdir}/condorLogs/out/condorOut_${tempNum}.out
error      = ${workdir}/condorLogs/err/condorErr_${tempNum}.err
log        = ${workdir}/condorLogs/log/condorLog_${tempNum}.log
getenv     = True
transfer_input_files = mix13_r56.cc 
accounting_group=group_cms
when_to_transfer_output = ON_EXIT
requirements = (Machine=!="cms-t3-wn3018.sdfarm.kr") && (Machine=!="cms-t3-wn3021.sdfarm.kr") && (Machine=!="bio-wn3006.sdfarm.kr") && (Machine=!="bio-wn3015.sdfarm.kr") && (Machine=!="bio-wn3013.sdfarm.kr")
request_Cpus=2
JobBatchName = Delphes_merging
queue 1
EOF

condor_submit job_${tempNum}.jds

done
