StartTime=$(date +%s)

sample=$1 #242LO_false_pdfwgt
idx=$2
#root -b -l -q run.C
#root -b -l .L histo_generator.cc+ test() 
#gSystem->Exit(0);

root -b -l <<EOF
.L histo_generator.cc+
//test("242LO_false_pdfwgt")
run_jobs("$sample",$idx)
EOF

EndTime=$(date +%s)
echo $EndTime
echo "runtime : $(($EndTime - $StartTime)) sec"