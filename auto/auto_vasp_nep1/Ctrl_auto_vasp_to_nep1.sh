#!/bin/source
#from rebreath
#本脚本的作用为使用vasp生成gpumd的训练集,并使用其生成初步的nep势
#使用方法为bash Ctrl_auto_vasp_nep1.sh，使用后即可一键生成nep势
#--------------------->workstation

Temp="800 1000"      #设置想要跑AIMD的温度有哪些，举例Tem="1100 1200 1300"
core_num=1       #mpirun使用的cpu或gpu个数
pertu_poscar=aln_wuzitble_poscar16_atoms   #设置需要进行微扰的poscar文件
aimd_poscar=aln_wuzitble_64         #设置需要进行aimd的poscar文件
single_point_energy_incar=INCAR_single_point_energy          #设置计算单点能的INCAR文件
aimd_incar=INCAR_AIMD                 #设置进行aimd的incar输入文件
library_address=~/.rebreath/auto/auto_vasp_nep1        #.rebreath的自动生成nep1文件库的所在地址
inter_configuration_step=50      #设置xdatcar文件每隔多少个构型提取一个
all_configuration=5000           #xdatcar文件总的构型数
vasp_exe=vasp_std_gpu
kspacing=0.2
export CUDA_VISIBLE_DEVICES=2
#-------------------------------------

caltime() {
  #计算程序运行的时间
  startime=$(date +%s)
  "$@"
  endtiem=$(date +%s)
  alltime=$((endtiem-startime))
  min_time=$((alltime/60))
  echo "执行命令 $* 用时 $min_time min"
}

check_vasp_complete() {
  #检查vasp是否计算完成，如果没有完成计算则会重新计算
  if [ ! -f  "OUTCAR" ];then
    echo "当前目录下没有找到OUTCAR文件，正在生成..."
    mpirun -np $core_num $vasp_exe
  elif [ -z "$(grep "General timing and accounting informations for this job" OUTCAR)" ];then
    echo "$PWD 目录下OUTCAR不完整，开始重新计算"
    mpirun -np $core_num $vasp_exe
  fi
}

add_kspacing_to_incar(){
  #如果kspacing的值不为零，则添加kspacing
  mv KPOINTS KPOINTS_backup || true
  if [ -z "$(grep "KSAPCING" INCAR)" ];then
     sed -i "\$a\KSPACING = $kspacing "  INCAR
  else
     sed "/KSPACING/c\KSPACING = $kspacing/" INCAR
  fi
}

sing_point_energy_outcar_to_train(){
  #本函数的目的是辅助生成单点能的训练集
  for i in $(find $PWD -type d -name "*aimd*" | xargs -n 1 basename);do
     if [ -f $i/OUTCAR ];then
       echo "将会将 ${i} 的OUTCAR文件进行改名"
       mv $i/OUTCAR $i/outcar_back || true
     fi
  done
  bash deal_outcar_to_train_auto.sh
  for i in $(find $PWD -type d -name "*aimd*" | xargs -n 1 basename);do
     if [ -f $i/outcar_back ];then
       echo "将会将 ${i} 的OUTCAR文件进行改名"
       mv $i/outcar_back $i/OUTCAR || true
     fi
  done
}
startime_all=$(date +%s)
set -e
Initial_location=$PWD
if [ -n "$(which figlet)" ]
then
    figlet Rebreath
fi
#dos2unix *
init_address=$PWD

#大召唤术
all_sh_pyfile="vasp_aimd_auto.sh deal_xdatcar_to_poscar_auto.sh single_point_energy_vaspset_auto.sh deal_outcar_to_train_auto.sh perturb_posfile_sets_auto.sh nep_train_auto.sh perturb.py"
for i in $all_sh_pyfile
do
  if [ ! -f "$PWD/$i" ];then
    echo "当前文件夹没有找到$i，正在复制脚本文件..."
    cp $library_address/$i  $init_address
  fi
done

sed -i "/NSW/ c\   NSW=$all_configuration" $aimd_incar
#1.进行AIMD模拟
caltime source vasp_aimd_auto.sh
cd $Initial_location
echo -e "\n--------------------->Ctrl_auto_vasp_nep1.sh向您汇报：控制第一环节——AIMD的 vasp_aimd_auto.sh已完成任务" 'O.<'

#2.对AIMD的轨迹问件——XDATCAR文件的处理
caltime source deal_xdatcar_to_poscar_auto.sh
cd $Initial_location
echo -e "\n--------------------->Ctrl_auto_vasp_nep1.sh向您汇报： 控制第二环节——对AIMD的轨迹问件——XDATCAR文件的处理的 deal_xdatcar_to_poscar_auto.sh已完成任务" 'O.<'

#3.计算单点能
caltime source single_point_energy_vaspset_auto.sh
cd $Initial_location
echo -e "\n--------------------->Ctrl_auto_vasp_nep1.sh向您汇报：控制第三环节——计算单点能的 single_point_energy_vaspset_auto.sh已完成任务" 'O.<'

if [ -f  "train.xyz" ]
then
    mv train.xyz train_back.xyz
fi
#4.处理单点能的计算结果outcar并生成训练集(本步骤已被优化)
#caltime bash deal_outcar_to_train_auto.sh

cd $Initial_location
echo -e "\n--------------------->Ctrl_auto_vasp_nep1.sh向您汇报：控制第四环节——处理单点能的计算结果outcar并生成训练集的 deal_outcar_to_train_auto.sh已完成任务" 'O.<'


#5.对少数原子进行微扰后将其产生的构型进行单点能计算生成训练集，与前面的训练集进行合并
#注意这里进行微扰的poscar文件与之前的不一样，原子数目需要少一些，AIN的计算中取了8个原子进行微扰

caltime source perturb_posfile_sets_auto.sh
cd $Initial_location
echo -e "\n--------------------->Ctrl_auto_vasp_nep1.sh向您汇报：控制第五环节——对少数原子进行微扰后将其产生的构型进行单点能计算生成训练集的 perturb_posfile_sets_auto.sh已完成任务，用时$min_time5 min"'O.<'

#形成训练集（包括了微扰与aimd的单点能）
sing_point_energy_outcar_to_train


#6.使用上面生成的训练集训练出一个nep势
#注意，这里产生的训练集为初步的nep势能，后续的需要你自己主动进行进一步的优化

caltime source nep_train_auto.sh
cd $Initial_location
echo -e "\n--------------------->Ctrl_auto_vasp_nep1.sh向您汇报：控制第五环节——使用上面生成的训练集训练出一个nep势的 nep_train_auto.sh已完成任务" 'O.<'


endtime_all=$(date +%s)
alltime_all=$((endtime_all-startime_all))
min_time_all=$((alltime_all/60))
echo -e "\n---------------->已完成nep势的初步构建，之后需要继续对nep进行改良"
echo -e "\n--------------------->控制台已完美完成所有的工作，总共用时$min_time_all min OVO"