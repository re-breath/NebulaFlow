#!/bin/bash
#form rebreath
#本脚的作用：自动进行nep的优化工作，本工作流将会运行3次，将nep优化三次
#优化的方案：先进行一次md的运算，生成一系列的构型，对每个构型进行一次单点能计算（本脚本只能使用vasp来计算）,将这些md计算的构型用nep进行预测一下，使用pynep选出最大差别的数据放到训练集重新进行训练
set -e


lib_address=~/.rebreath/auto/modify_nep
group_size=66 #原子个数加二
core_num=1
incar_single_energy='INCAR_single_point_energy'
elect_file=virial_train.out   #挑选依据的文件
elect_num=20     #挑选时候挑选的数目
vasp_exe=vasp_std_gpu
export CUDA_VISIBLE_DEVICES=2
dump_configs=20  #一个温度下GPUMD跑出来的构型数量
max_tem=1200   #md模拟的最大温度
nep_step=100000  #nep训练的步数
kspacing=0.2

check_gpumd_complete(){
  if [ ! -d "train_folders" ];then
     if [ $(grep "^" thermo.out |wc -l)  -ne 100 ];then
       gpumd
     fi
  fi
}
check_vasp_complete() {
  if [ ! -f  "OUTCAR" ];then
    echo "当前目录下没有找到OUTCAR文件，正在生成..."
    mpirun -np $core_num $vasp_exe
  elif [ -z "$(grep "General timing and accounting informations for this job" OUTCAR)" ];then
    echo "$PWD 目录下OUTCAR不完整，开始重新计算"
    mpirun -np $core_num $vasp_exe
  fi
}

add_kspacing_to_incar(){
  mv KPOINTS KPOINTS_backup || true
  if [ -z "$(grep "KSAPCING" INCAR)" ];then
     sed -i "\$a\KSPACING = $kspacing "  INCAR
  else
     sed "/KSPACING/c\KSPACING = $kspacing/" INCAR
  fi
}

modify_initdir=$PWD
#大召唤术
modify_nep_sh="mul_tem_md_vasp.sh train-dumpxyz_to_POSCAR.py deal_outcar_to_train_auto.sh hplt_nep_results_prediction.py run_alltrain-file.sh hplt_nep_results.py get_max_rmse_xyz.py"
for file in $modify_nep_sh;do
  if [ ! -f "$PWD/$file" ];then
    echo "当前文件夹没有找到$file，正在复制脚本文件..."
    cp $lib_address/$file  $modify_initdir
  fi
done

if [ ! -f "INCAR_single_point_energy" ];then
  mv  $incar_single_energy INCAR_single_point_energy
fi

for loop in {1..3};do
  if [ ! -d "nep_$loop" ];then
      mkdir -p "nep_$loop"
  fi
  cd $modify_initdir/nep_$loop
  mkdir -p md
  tomdfiles="nep.txt run.in model.xyz mul_tem_md_vasp.sh train-dumpxyz_to_POSCAR.py KPOINTS POTCAR run_alltrain-file.sh INCAR_single_point_energy deal_outcar_to_train_auto.sh"
  for i in $tomdfiles;do
    if [ ! -f "$i" ];then
      echo "当前文件夹没有找到$i，正在复制脚本文件..."
      cp $modify_initdir/$i md
    fi
  done
  cd md
  source mul_tem_md_vasp.sh
  mv $modify_initdir/nep_$loop/md/train.xyz $modify_initdir/nep_$loop/md/train_md_nep_$loop.xyz
  cd $modify_initdir/nep_$loop
  mkdir -p prediction_elect
  cd prediction_elect
  toprediction_electfiles="hplt_nep_results_prediction.py nep.txt nep.in get_max_rmse_xyz.py "
  for i in $toprediction_electfiles;do
    cp $modify_initdir/$i .
  done
  if [ -n "$(grep "prediction" $modify_initdir/nep_$loop/prediction_elect/nep.in)" ]
  then
     sed -i '/prediction/c\prediction 1' $modify_initdir/nep_$loop/prediction_elect/nep.in
  else
     echo -e "\nprediction 1 " >> $modify_initdir/nep_$loop/prediction_elect/nep.in
  fi
  cp  $modify_initdir/nep_$loop/md/train_md_nep_$loop.xyz train.xyz
  nep
  mv train.xyz train_md_nep_$loop.xyz
  python3 hplt_nep_results_prediction.py
  python3 get_max_rmse_xyz.py  train_md_nep_$loop.xyz $elect_file $elect_num
  cp $modify_initdir/train.xyz $modify_initdir/test.xyz .
  cat find_out.xyz >> train.xyz
  cat reserve.xyz >>test.xyz
  mkdir -p $modify_initdir/nep_$loop/neptrain
  if [ -z  "$(grep "100000" $modify_initdir/nep_$loop/neptrain/loss.out)" ];then
    cp $modify_initdir/nep_$loop/prediction_elect/train.xyz $modify_initdir/nep_$loop/neptrain/train.xyz
    cp $modify_initdir/nep_$loop/prediction_elect/test.xyz $modify_initdir/nep_$loop/neptrain/test.xyz
    cp $modify_initdir/nep.in $modify_initdir/nep_$loop/neptrain/nep.in
    cd $modify_initdir/nep_$loop/neptrain
    nep
  else
    cd $modify_initdir/nep_$loop/neptrain
  fi
  cp $modify_initdir/hplt_nep_results.py  .
  python3 hplt_nep_results.py
  cp nep.txt $modify_initdir/nep.txt
  m=$((loop-1))
  mv $modify_initdir/train.xyz  $modify_initdir/train_$m.xyz
  mv $modify_initdir/test.xyz $modify_initdir/test_$m.xyz
  cp train.xyz $modify_initdir/train.xyz
  cp test.xyz $modify_initdir/test.xyz
  cd $modify_initdir
done





