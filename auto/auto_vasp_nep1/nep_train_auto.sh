#!/bin/bash
#from rebreath
#-------------------------------->6<-----------------------------------
#本脚本为自动化的其中一环，执行序列为6号（第6个执行该脚本）
#本脚本的作用为简化操作的脚本，使用之前训练的训练集，进行训练nep

set -e
if [ -n "$1" ]
then
    need_perturbfile=$1
elif [ -n "$pertu_poscar" ]
then
    need_perturbfile=$pertu_poscar
else
    need_perturbfile=POSCAR
fi
startime=$(date +%s)
init_address=$PWD
if [ -n "$library_address" ]
then
    library_address=$library_address
else
    library_address=~/.rebreath/auto/datesets_Construction/auto_vasp_set
fi
mkdir -p $init_address/nep_train
cp $init_address/train.xyz $init_address/nep_train/test.xyz
cp $init_address/train.xyz $init_address/nep_train/train.xyz
cp $init_address/nep.in $init_address/nep_train/nep.in
cd $init_address/nep_train
nep
cd $init_address
endtime=$(date +%s)
alltime=$((endtime-startime))
min_time=$((alltime/60))
echo -e "--------------------->nep_train_auto.sh 完美完成工作,总共用时$min_time min OVO"
