
function submit_tag_i(){

tag=$1
i=$2

source batch_creator.sh ${tag}_${i} 
cp run_for_batch.sh JOB_${tag}_${i}/run_${tag}_${i}.sh
pushd JOB_${tag}_${i}
find . -name run_${tag}_${i}.sh | xargs perl -pi -e s/__IDX__/$i/g
find . -name run_${tag}_${i}.sh | xargs perl -pi -e s/__SAMPLE__/$tag/g
condor_submit submit.jds
popd
}


submit_tag_i 242LO_false_pdfwgt 50
submit_tag_i 242LO_true_pdfwgt 50
submit_tag_i 26xLO_false_pdfwgt 50
submit_tag_i 26xLO_true_pdfwgt 50