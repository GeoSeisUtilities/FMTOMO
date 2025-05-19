#!/bin/bash

rm src.dat invert_*/src.dat
cd invert_p/mkmodel/
echo "--------------------------------"
echo "preparing P"
./makemodel_copy.sh
cd ..
./prepare_synthetic.sh mkmodel/sourceslocal.in
cd ../invert_s/mkmodel/
echo "--------------------------------"
echo "preparing S"
./makemodel_copy.sh
cd ../
./prepare_synthetic.sh mkmodel/sourceslocal.in
cd ../
echo "--------------------------------"
echo "running the tomography"
jointimes
tomo3dj
echo "--------------------------------"
echo "calculating the norm"
cd invert_p
mprop > mprop.txt
cd ../invert_s
mprop > mprop.txt
echo "--------------------------------"
echo "inverting vp/vs"
cd ../vpvs
./run_tomography.sh
echo "end of the tomography process"
