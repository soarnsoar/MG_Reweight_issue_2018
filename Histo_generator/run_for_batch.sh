#!/bin/bash
#touch /cms/ldap_home/jhchoi/generator_group/reweight_issue/Histo_generator/jobstart

StartTime=$(date +%s)

SECTION=`printf %03d $1`
WORKDIR=`pwd`    
echo "#### Extracting cmssw ####" 
tar -zxvf INPUT.tar.gz
echo "#### cmsenv ####" 
export CMS_PATH=/cvmfs/cms.cern.ch
source $CMS_PATH/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
#touch /cms/ldap_home/jhchoi/generator_group/reweight_issue/Histo_generator/tarend


cd CMSSW_10_3_0/src                                                                                                                           
scram build ProjectRename
eval `scramv1 runtime -sh`
cd ../../
#touch /cms/ldap_home/jhchoi/generator_group/reweight_issue/Histo_generator/scram




sample=__SAMPLE__
idx=__IDX__


#root -b -l -q run.C
#root -b -l .L histo_generator.cc+ test() 
#gSystem->Exit(0);
#root -b -l -q 'run.C("__SAMPLE__",__IDX__)'
#touch /cms/ldap_home/jhchoi/generator_group/reweight_issue/Histo_generator/runroot

root -b -l <<EOF
.L histo_generator.cc+
//run_job_i_batch("$sample",$idx)
run_jobs("$sample",$idx)
//test()
.q
EOF
#root -b -l -q run.C 
#touch /cms/ldap_home/jhchoi/generator_group/reweight_issue/Histo_generator/endroot

#mv $sample/OUTPUT.root /cms/ldap_home/jhchoi/generator_group/reweight_issue/Histo_generator/batch.root
#touch /cms/ldap_home/jhchoi/generator_group/reweight_issue/Histo_generator/jobend
EndTime=$(date +%s)
echo $EndTime
echo "runtime : $(($EndTime - $StartTime)) sec"
