#!/bin/bash
#form rebreath
#本脚本的目的：为了在计算完md后自动化预测与挑选合适的构型到训练集与测试集，而后进行新一轮的nep势训练
#本脚本的用法：bash predict_elect_configs_to_newnep.sh
#脚本的使用场景：需要在md的文件夹下运行
#-----------------------------settings------------------------
elect_configs_num=18

#-------------------------------------------------------------
set -e
init_address=$PWD   #该位置应当为md的文件夹位置
startime=$(date +%s)
#进行预测
cd ..
if [ ! -d "prediction_md" ];then
  mkdir prediction_md
else
  rm -r prediction_md
  mkdir prediction_md
fi
cd $init_address
#大召唤术
files="nep.in hplt_nep_results.py nep.txt"
for i in $files
do
  cp ../../prediction/$i  ../prediction_md/
done
cp $init_address/train.xyz ../prediction_md/train.xyz
cp $init_address/test.xyz ../prediction_md/test.xyz

prestr=$(grep "prediction" ../prediction_md/nep.in)
if [ -n "$prestr" ]
then
  sed -i '/prediction/c\prediction 1' ../prediction_md/nep.in
else
  echo 'prediction 1 ' >> nep.in
fi
cd ../prediction_md
nep
python3 hplt_nep_results.py

#进行挑选构型
cd $init_address/..
if [ ! -d "elect" ];then
  mkdir elect
else
  rm -r elect
  mkdir elect
fi
cp ~/.rebreath/auto/modify_nep/get_max_rmse_xyz.py elect/
cp $init_address/train.xyz elect/train_md.xyz
cp ../test.xyz elect/test.xyz
cp ../train.xyz elect/train.xyz
cp prediction_md/*.out elect/
cd elect
python3 get_max_rmse_xyz.py train_md.xyz virial_train.out $elect_configs_num
cat find_out.xyz >> train.xyz
cat reserve.xyz >>test.xyz
if [ ! -d "../../nep_modfy" ];then
  mkdir ../../nep_modfy
else
  rm -r ../../nep_modfy
  mkdir ../../nep_modfy
fi
cp train.xyz ../../nep_modfy/train.xyz
cp test.xyz ../../nep_modfy/test.xyz
cp ../../nep.in ../../nep_modfy/nep.in
cd ../../

cd nep_modfy
nep
endtime=$(date +%s)
alltime=$((endtime-startime))
min_time=$((alltime/60))
echo -e "---------------->predict_elect_configs_to_newnep.sh 完成工作，用时$min_time min OVO"
