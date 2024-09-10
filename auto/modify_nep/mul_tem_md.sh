#!/usr/bin/bash
# 本脚本的作用为对当前文件夹下的run.in文件进行修改，\
# 将其中的replace置换为不同温度，并放入不同的文件夹进行MD模拟，\
# 之后每个MD提取20个构型进行单点能计算，之后收集outcar信息生成一个训练集
set -e
init_address=$PWD
allstartime=$(date +%s)
for i in {400..1000..100};do
  startime=$(date +%s)
  if [ ! -d "md_$i" ];then
    mkdir md_$i
  else
    rm -r md_$i
    mkdir md_$i
  fi
  cp run.in md_$i
  cp model.xyz md_$i
  cp nep.txt md_$i
  cp train-dumpxyz_to_POSCAR.py md_$i
  cp KPOINTS POTCAR "md_$i"
  cp INCAR_single_point_energy md_$i/INCAR
  cd md_$i
  sed -i "s/replace/$i/g" run.in
  gpumd
  if [ -d "train_folders" ];then
    rm -r train_folders
  fi
  python3 train-dumpxyz_to_POSCAR.py
  cd $init_address
  endtime=$(date +%s)
  alltime=$((endtime-startime))
  min_time=$((alltime/60))
  echo -e "---------------->md_$i的MD任务已完成，用时$min_time min"
done
cd $init_address
startime_alltrain=$(date +%s)
bash run_alltrain-file.sh
endtime_alltrain=$(date +%s)
alltime_alltrain=$((endtime_alltrain-startime_alltrain))
min_time_alltrain=$((alltime_alltrain/60))
echo -e "---------------->所有的单点能计算已完成，用时$min_time_alltrain min"
cd $init_address
bash deal_outcar_to_train_auto.sh
allendtime=$(date +%s)
alltime=$((allendtime-allstartime))
min_time=$((alltime/60))
echo -e "---------------->所有任务已完成，总用时$min_time min OVO"