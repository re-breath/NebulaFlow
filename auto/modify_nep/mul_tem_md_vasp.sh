#!/usr/bin/bash
# 本脚本的作用为对当前文件夹下的run.in文件进行修改，\
# 将其中的replace置换为不同温度，并放入不同的文件夹进行MD模拟，\
# 之后每个MD提取20个构型进行单点能计算，之后收集outcar信息生成一个训练集
set -e
group_size=${group_size:=56}   #一个构型的行数
dump_configs=${dump_configs:=20}
kspacing=${kspacing:=0}
core_num=${core_num:=1}

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

sed  -i "/group_size =/c\    group_size = $group_size" train-dumpxyz_to_POSCAR.py
init_address=$PWD
allstartime=$(date +%s)
for i in {400..1200..100};do
  startime=$(date +%s)
  goto=$(grep "Lattice" md_$i/dump.xyz | wc -l)
  if [ $goto -eq $dump_configs ];then
    break
  fi
  mkdir -p md_$i
  cp run.in md_$i
  cp model.xyz md_$i
  cp nep.txt md_$i
  cp train-dumpxyz_to_POSCAR.py md_$i
  cp KPOINTS POTCAR "md_$i"
  cp INCAR_single_point_energy md_$i/INCAR
  cd md_$i
  sed -i "s/replace/$i/g" run.in
  check_gpumd_complete
  if [ -d "train_folders" ];then
    echo "正在进行之前的工作，如果之前的工作有问题请删除文件夹后再重新开始任务吧"
    else
    python3 train-dumpxyz_to_POSCAR.py
  fi

  cd $init_address
  endtime=$(date +%s)
  alltime=$((endtime-startime))
  min_time=$((alltime/60))
  echo -e "---------------->md_$i的MD任务已完成，用时$min_time min"
done
cd $init_address
startime_alltrain=$(date +%s)
for i in $(seq 3);do
  source run_alltrain-file.sh || true
done
endtime_alltrain=$(date +%s)
alltime_alltrain=$((endtime_alltrain-startime_alltrain))
min_time_alltrain=$((alltime_alltrain/60))
echo -e "---------------->所有的单点能计算已完成，用时$min_time_alltrain min"
cd $init_address
if [ -f "train.xyz" ];then
  mv train.xyz train_pro_md.xyz
fi
bash deal_outcar_to_train_auto.sh
allendtime=$(date +%s)
alltime=$((allendtime-allstartime))
min_time=$((alltime/60))
echo -e "---------------->所有任务已完成，总用时$min_time min OVO"