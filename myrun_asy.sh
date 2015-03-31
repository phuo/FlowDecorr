#!/bin/bash

echo "Starting the run"
#cat input.txt | sed s/,/\\n/g > input2.txt; 
echo "List is made"

#root -b -l <<EOF
#.L flowAsyTwist.C+
#.q
#EOF
#1 calib 2 scalar-product 3 zbin 4 addstat 5 nevts
root -b -l <<EOF
    gSystem->Load("flowAsyTwist_C.so");
    flowAsyTwist* flow = new flowAsyTwist();
    char *input="input2.txt";
    char *out = "Analysis_calib$1_wtq$2_addstat$4.root";
    flow->from = 0;
    flow->to = 20;
    flow->filelist = input;
    flow->outfile  = out;
    flow->calib =$1 ;
    flow->perpartQ = $2;
    flow->Zbin_Size = $3;
    flow->addstat=$4;
    flow->nevts=$5;
    flow->Init_histos();
    flow->run();
    flow->Save_histos();
.q; 
EOF

echo "DONE RUNNING THE SCRIPT"
