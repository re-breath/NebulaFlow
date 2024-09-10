#!/bin/bash
#进入当前文件夹所有的（job-*）进入后提交vasp作业
for i in  $(find . -type d -name "job-00*" )
do
     cd $i
     dos2unix * || true
     mpirun -np 1 vasp
     cd ..
done

echo "sbatch.sh 已完成对所有文件夹的计算。"
