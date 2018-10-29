#!/bin/bash

G_RELEASE="CMSSW_10_3_0"
G_JOBDIR=`pwd`
G_PYTHON='temp.py'
G_RUNFILE='runtemp.sh'
function run_template(){
_pytorun=$1
_runfile=$2
echo '#!/bin/bash' > $_runfile                           
echo 'SECTION=`printf %03d $1`' >> $_runfile
echo 'WORKDIR=`pwd`'>> $_runfile
echo 'echo "#### Extracting cmssw ####"'>> $_runfile
echo 'tar -zxvf INPUT.tar.gz'>> $_runfile
echo 'echo "#### cmsenv ####"'>> $_runfile
echo 'export CMS_PATH=/cvmfs/cms.cern.ch'>> $_runfile
echo 'source $CMS_PATH/cmsset_default.sh'>> $_runfile
echo 'export SCRAM_ARCH=slc6_amd64_gcc630'>> $_runfile

echo "cd $G_RELEASE/src">> $_runfile
echo 'scram build ProjectRename'>> $_runfile
echo 'eval `scramv1 runtime -sh`'>> $_runfile
echo 'cd ../../'>> $_runfile
echo "cmsRun $_pytorun">> $_runfile
G_RUNFILE=$_runfile
}


function batch_creater(){
    echo "===batch_creater_jhchoi.sh==="
#############Set variable##########
    #CURDIR=`pwd`
    JOBNAME="JOB_"$1
    NJOBS=1
    echo "JOB DIR ="$JOBNAME
##input tar's name => INPUT.tar.gz
    
###################################
    G_JOBDIR=$JOBNAME
      
    
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
	
	
	tar -czf INPUT_${1}.tar.gz CMSSW* *.py
	
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
	echo "transfer_output_remaps = \"Ntuple.root = OUTPUT_\$(Process).root\"" >> submit.jds
	echo "queue $NJOBS" >> submit.jds

	popd
	echo "===Make submit.jds DONE.==="

    fi
###################################
    
    #cd $CURDIR
  
##then move to $JOBNAME directory and
##make run.sh script
##make output file name to OUTPUT.root
#condor_submit submit.jds
    
    
    echo "batch_creater_jhchoi.sh DONE."
    
}

function LOpython_script(){
_version=$1 #26xLO
_PDFWGT=$2 #false_pdfwgt
_IJOB=$3
cp run_weight_assign_${_version}_template.py run_ntupler_${_version}_${_PDFWGT}.py
find . -name run_ntupler_${_version}_${_PDFWGT}.py | xargs perl -pi -e s/PDFWGT/${_PDFWGT}/g
find . -name run_ntupler_${_version}_${_PDFWGT}.py | xargs perl -pi -e s/IJOB/${_IJOB}/g
G_PYTHON="run_ntupler_${_version}_${_PDFWGT}.py"
}

function submit_LOjob(){
    version=$1
    pdfwgt=$2
    NTOTAL=$3
#    LASTi=`expr $NTOTAL - 1`
    i=0
    while [ $i -lt $NTOTAL ]
    do
	echo $i
	LOpython_script $version $pdfwgt $i
	batch_creater ${version}_${pdfwgt}_$i 
	mv $G_PYTHON $G_JOBDIR/ 
	run_template $G_PYTHON run_${version}_${pdfwgt}_$i.sh
	mv $G_RUNFILE $G_JOBDIR/
	pushd $G_JOBDIR 
	condor_submit submit.jds
	popd                                                                                                                                 
	echo "end submit"
	pwd
        i=`expr $i + 1`
	echo $i

    done
}


#Version="242LO"
#pdfwgt="false_pdfwgt"
#NTOTAL=50
#i=0

#LOpython_script $version $pdfwgt $i
#batch_creater ${version}_${pdfwgt}_$i
#mv $G_PYTHON $G_JOBDIR/
#run_template $G_PYTHON run_${version}_${pdfwgt}_$i.sh
#mv $G_RUNFILE $G_JOBDIR/
#pushd $G_JOBDIR
#condor_submit submit.jds                                                                                                                   
#popd
#i=0   
#i=`expr $i + 1`
#echo $i

#submit_LOjob 242LO false_pdfwgt 50
submit_LOjob 242LO true_pdfwgt 50
submit_LOjob 26xLO false_pdfwgt 50
submit_LOjob 26xLO true_pdfwgt 50

#N=100
#echo `expr $N - 1` 