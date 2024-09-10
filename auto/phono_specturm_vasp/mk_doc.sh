#!/bin/bash
#寻找生成的POSCAR-[0-9]+文件的数量并分发到其他文件夹中
m=0
m=`find ./ -type f -regex ".*POSCAR-[0-9]+" |wc -l`
for i in $(seq -f "-%03g" 1 $m)
do
    mkdir -p job$i
    mv ./POSCAR$i ./job$i
    cp ./INCAR   ./job$i
    cp ./POTCAR  ./job$i
    cp ./KPOINTS  ./job$i || true

done
echo "mk_doc.sh的任务完成啦  >_<"
