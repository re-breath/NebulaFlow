#!/usr/bin/bash
init_address=$PWD
startime=$(date +%s)
for i in $(find $init_address/ -type d -name "train-*")
do
  cd $i
  mpirun -np 1 vasp
  cd $init_address
done
endtime=$(date +%s)
alltime=$((endtime-startime))
min_time=$((alltime/60))
echo -e "--------------------->run_train-file.sh 完美完成工作,总共用时$min_time min OVO"