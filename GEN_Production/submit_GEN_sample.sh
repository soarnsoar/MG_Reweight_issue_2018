#!/bin/bash
##This is for KISTI job
function run_template(){
scramarch=$1
cmsrel=$2
echo '#!/bin/bash' > run_template.sh                           
echo 'SECTION=`printf %03d $1`' >> run_template.sh
echo 'WORKDIR=`pwd`'>> run_template.sh
echo 'echo "#### Extracting cmssw ####"'>> run_template.sh
echo 'tar -zxvf INPUT.tar.gz'>> run_template.sh
echo 'echo "#### cmsenv ####"'>> run_template.sh
echo 'export CMS_PATH=/cvmfs/cms.cern.ch'>> run_template.sh
echo 'source $CMS_PATH/cmsset_default.sh'>> run_template.sh
echo "export SCRAM_ARCH=$scramarch">> run_template.sh

echo "cd ${cmsrel}/src">> run_template.sh
echo 'scram build ProjectRename'>> run_template.sh
echo 'eval `scramv1 runtime -sh`'>> run_template.sh
echo 'cd ../../'>> run_template.sh
echo 'cmsRun __SCRIPT__.py'>> run_template.sh

}
function batch_creater(){
    echo "===batch_creater_jhchoi.sh==="
#############Set variable##########
    #CURDIR=`pwd`
    JOBNAME="JOB_"$1
    NJOBS=$2
    GRIDPACK=$3
    echo "JOB DIR ="$JOBNAME
##input tar's name => INPUT.tar.gz
    
###################################
    
    
    if [ -z $2 ];then
	echo "default NJOBs = 1"
	NJOBS=1
    fi
    
    
    
    
    
#########Check input argument######
    if [ -z $1 ];then
	echo "Need argument"
	
	
#########Check alreay job env######
    elif [ -d $JOBNAME ];then
	echo "The job directory alreay exists"
	
    else
	
##Make tar input
	
#tar -cvzf INPUT.tar.gz *
	echo "===Make INPUT.tar.gz==="
	
	
	tar -czf INPUT_${1}.tar.gz CMSSW* *.py $GRIDPACK
	
#$JOBNAME
	
	mkdir $JOBNAME
	pushd $JOBNAME
	mv ../INPUT_${1}.tar.gz INPUT.tar.gz
	
	echo "===Make submit.jds==="
	
	echo "executable = run_${1}.sh" >> submit.jds 
	echo "universe   = vanilla" >> submit.jds
	echo "arguments  = \$(Process)" >> submit.jds
	echo "requirements = OpSysMajorVer == 6" >> submit.jds
	echo "log = condor.log" >> submit.jds
	echo "getenv     = True" >> submit.jds
	echo "should_transfer_files = YES" >> submit.jds
	echo "when_to_transfer_output = ON_EXIT" >> submit.jds
	echo "output = job_\$(Process).log" >> submit.jds
	echo "error = job_\$(Process).err" >> submit.jds
	echo "transfer_input_files = INPUT.tar.gz" >> submit.jds
#echo "use_x509userproxy = true" >> submit.jds
#	echo "transfer_output_remaps = \"*inDQM.root = OUTPUT_\$(Process).root\"" >> submit.jds
	echo "transfer_output_remaps = \"OUTPUT.root = OUTPUT_\$(Process).root\"" >> submit.jds
	echo "queue $NJOBS" >> submit.jds
	
	
	echo "===Make submit.jds DONE.==="
	
    fi
###################################
    
    #cd $CURDIR
    popd
##then move to $JOBNAME directory and
##make run.sh script
##make output file name to OUTPUT.root
#condor_submit submit.jds
    
    
    echo "batch_creater_jhchoi.sh DONE."
    
}

function submit_batch(){
    JOB=$1
    NJOB=$2
    PYTHON=$3
    GRIDPACK=$4
    batch_creater $JOB $NJOB $GRIDPACK
    pushd JOB_$JOB
    cp ../run_template.sh run_$JOB.sh
    
    find . -name run_$JOB.sh | xargs perl -pi -e s/__SCRIPT__/$PYTHON/g
    condor_submit submit.jds
    popd
}


##This is main###
### settings to modify
# specify batch system 
BATCH=LSF # SGE LSF 
# number of jobs 
NJOBS=50
# number of events per job 
NEVTS=50000  
# path to submit jobs 
WORKDIR=`pwd -P`
# path for private fragments not yet in cmssw
FRAGMENTDIR=${WORKDIR}
# release setup 
#
#export SCRAM_ARCH=slc6_amd64_gcc630
export SCRAM_ARCH=slc6_amd64_gcc700
#SCRAM_ARCH=slc6_amd64_gcc481
#RELEASE=CMSSW_9_3_8
RELEASE=CMSSW_10_3_0
#RELEASE=CMSSW_7_1_30
# path to store output files
ODIR=${WORKDIR}/samples

# define output tags as well as corresponding gridpacks and shower fragments 
OTAGLIST=()
GRIDPACKLIST=() 
GENFRAGMENTLIST=()

#OTAGLIST+=( dyee012j_261_false_pdfwgt )
OTAGLIST+=( dyee012j_5f_LO_242_false_pdfwgt dyee012j_5f_LO_242_true_pdfwgt dyee012j_5f_LO_26x_false_pdfwgt dyee012j_5f_LO_26x_true_pdfwgt )

#GRIDPACKLIST+=( ${WORKDIR}/dyellell01234j_5f_LO_MLM_VMG5_26x_false_pdfwgt_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz ${WORKDIR}/dyellell01234j_5f_LO_MLM_VMG5_26x_true_pdfwgt_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz ${WORKDIR}/dyellell01234j_5f_LO_MLM_VMG5_261_false_pdfwgt_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz ${WORKDIR}/dyellell01234j_5f_LO_MLM_VMG5_261_true_pdfwgt_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz)
#GENFRAGMENTLIST+=( Hadronizer_TuneCUETP8M1_13TeV_MLM_5f_max2j_LHE_pythia8_cff ) 
GRIDPACKLIST+=( dyee012j_5f_LO_242_false_pdfwgt_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz dyee012j_5f_LO_242_true_pdfwgt_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz dyee012j_5f_LO_26x_false_pdfwgt_slc6_amd64_gcc630_CMSSW_9_3_8_tarball.tar.xz dyee012j_5f_LO_26x_true_pdfwgt_slc6_amd64_gcc630_CMSSW_9_3_8_tarball.tar.xz )

GENFRAGMENTLIST+=( Hadronizer_TuneCP5_13TeV_MLM_5f_max2j_LHE_pythia8_cff Hadronizer_TuneCP5_13TeV_MLM_5f_max2j_LHE_pythia8_cff Hadronizer_TuneCP5_13TeV_MLM_5f_max2j_LHE_pythia8_cff Hadronizer_TuneCP5_13TeV_MLM_5f_max2j_LHE_pythia8_cff )
#GENFRAGMENT=Hadronizer_TuneCUETP8M1_13TeV_MLM_5f_max4j_LHE_pythia8_cff # wjets/zjets
#GENFRAGMENT=Hadronizer_TuneCUETP8M1_13TeV_aMCatNLO_FXFX_5f_max2j_max0p_LHE_pythia8_cff # zjets fxfx
#GENFRAGMENT=Hadronizer_TuneCUETP8M1_13TeV_aMCatNLO_FXFX_5f_max2j_max1p_LHE_pythia8_cff # ttbar fxfx 
### done with settings 



### setup release 
if [ -r ${WORKDIR}/${RELEASE}/src ] ; then 
    echo release ${RELEASE} already exists
else
    scram p CMSSW ${RELEASE}
fi
cd ${WORKDIR}/${RELEASE}/src
eval `scram runtime -sh`


### checkout generator configs 
git-cms-addpkg --quiet Configuration/Generator


### copy additional fragments if needed 
if [ -d "${FRAGMENTDIR}" ]; then 
    cp ${FRAGMENTDIR}/*.py ${CMSSW_BASE}/src/Configuration/Generator/python/. 
fi


### scram release 
scram b 


### start tag loop for setups to be validated  
NTAG=`echo "scale=0; ${#OTAGLIST[@]} -1 " | bc` 

for ITAG in `seq 0 ${NTAG}`; do
    OTAG=${OTAGLIST[${ITAG}]}
    GRIDPACK=${GRIDPACKLIST[${ITAG}]}
    GENFRAGMENT=${GENFRAGMENTLIST[${ITAG}]}
    
    ### move to python path 
    cd ${CMSSW_BASE}/src/Configuration/Generator/python/
    
    ### check that fragments are available 
    echo "Check that fragments are available ..."
    if [ ! -s ${GENFRAGMENT}.py ] ; then 
	echo "... cannot find ${GENFRAGMENT}.py"
	exit 0;
    else
	echo "... found required fragments!"
    fi
    

    
    ### create generator fragment 
    CONFIG=${OTAG}_cff.py
    if [ -f ${CONFIG} ] ; then 
	rm ${CONFIG} 
    fi
    
    cat > ${CONFIG} <<EOF

import FWCore.ParameterSet.Config as cms
externalLHEProducer = cms.EDProducer('ExternalLHEProducer', 
args = cms.vstring('${WORKDIR}/${GRIDPACK}'),
nEvents = cms.untracked.uint32(5000),
numberOfParameters = cms.uint32(1),  
outputFile = cms.string('cmsgrid_final.lhe'),
scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)
EOF
    cat ${GENFRAGMENT}.py >> ${CONFIG}
    
       
    ### make validation fragment 
    echo "make validation fragment" 
    cmsDriver.py Configuration/Generator/python/${CONFIG} \
	-n ${NEVTS} --mc --no_exec --python_filename cmsrun_${OTAG}.py \
	-s LHE,GEN --datatier GEN,LHE --eventcontent RAWSIM,LHE \
	--conditions auto:run2_mc_FULL --beamspot Realistic8TeVCollision \
	--fileout file:OUTPUT.root
    echo "move to submission directory"
    ### move to submission directory 
    cd ${WORKDIR}


    ### prepare submission script 



#    cp ${CMSSW_BASE}/src/Configuration/Generator/python/cmsrun_${OTAG}.py . 
    echo "adjust random numbers "
### adjust random numbers 
   
    echo '#####Set Random number####' > cmsrun_${OTAG}.py 
    echo 'from datetime import datetime' >> cmsrun_${OTAG}.py
    echo 'dt = datetime.now()' >> cmsrun_${OTAG}.py
    echo 'myseed=dt.microsecond' >> cmsrun_${OTAG}.py
    echo 'print myseed' >> cmsrun_${OTAG}.py
    echo 'dt2 = datetime.now()' >> cmsrun_${OTAG}.py
    echo 'myseed2=dt2.microsecond+1' >> cmsrun_${OTAG}.py
    echo 'print myseed2' >> cmsrun_${OTAG}.py
    echo 'dt3 = datetime.now()' >> cmsrun_${OTAG}.py
    echo 'myseed3=dt3.microsecond+2' >> cmsrun_${OTAG}.py
    echo 'print myseed3' >> cmsrun_${OTAG}.py
    
########################## 
    echo "cp the config"

    cat ${CMSSW_BASE}/src/Configuration/Generator/python/cmsrun_${OTAG}.py >> cmsrun_${OTAG}.py
    
    LINE=`egrep -n Configuration.StandardSequences.Services_cff cmsrun_${OTAG}.py | cut -d: -f1 `
   # SEED=`echo "5267+${OFFSET}" | bc`
    echo "add set seed"
    sed -i "${LINE}"aprocess.RandomNumberGeneratorService.generator.initialSeed=myseed cmsrun_${OTAG}.py  
    #SEED=`echo "289634+\${OFFSET}" | bc`
    sed -i "${LINE}"aprocess.RandomNumberGeneratorService.externalLHEProducer.initialSeed=myseed2 cmsrun_${OTAG}.py  
#JOB=$1
#NJOB=$2
#PYTHON=$3
    run_template $SCRAM_ARCH $RELEASE
    submit_batch $OTAG $NJOBS cmsrun_$OTAG $GRIDPACK

done # end of tag loop 






