#!/bin/bash
#from rebreath
#本脚本的作用为使用vasp生成微扰后的gpumd的训练集
#本脚本可以指定文件当POSCAR，也可以使用默认的POSACR
#-------------------------------->5<-----------------------------------
#本脚本为自动化的其中一环，执行序列为5号（第5个执行该脚本）
#本脚本的作用：建立工作文件夹对相应的文件夹进行微扰并计算单点能，生成一个训练集

set -e
startime=$(date +%s)
init_address=$PWD
#下面为对于自动化控制台参数进行判定的判断
#-----------------------独立使用该脚本下面的这部分无需太注意-------------------------
vasp_exe=${vasp_exe:="vasp"}
core_num=${core_num:=1}
if [ -n "$single_point_energy_incar" ]
then
    incarfile=$single_point_energy_incar
else
    incarfile=INCAR
fi

if [ -n "$pertu_poscar" ];then
    poscarfile=$pertu_poscar
else
    poscarfile=POSCAR
fi
library_address=${library_address:=~/.rebreath/auto/auto_vasp_nep1}
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

#-------------------------------------------------------------------------------
need_py='perturb.py'
if [ -n "$1" ];then
    need_perturbfile=$poscarfile
else
    need_perturbfile=$1
fi
if [ ! -f  "$need_py" ]
then
    echo "当前文件夹没有找到$need_py，正在复制脚本文件..."
    cp $library_address/$need_py  $init_address
fi

mkdir -p  $init_address/perturb_posfile_sets
mvfile="POTCAR KPOINTS $need_py"
for j in $mvfile
do
  cp $init_address/$j $init_address/perturb_posfile_sets
done
cp $init_address/$incarfile $init_address/perturb_posfile_sets/INCAR
cp $init_address/$poscarfile $init_address/perturb_posfile_sets/POSCAR
cd $init_address/perturb_posfile_sets
if [ -z "$(grep "complete all pertubecompute" $init_address/pertubeset.log)" ];then
python3 $need_py
fi
echo "微扰需要的文件已经组装完毕，即将对微扰后的poscar进行单点能计算"
for r in $(find $init_address/perturb_posfile_sets -type d -name "train-*")
do
  cd $r
  add_kspacing_to_incar
  check_vasp_complete
  cd ..
done
echo "complete all pertubecompute" > $init_address/pertubeset.log
cd $init_address
endtime=$(date +%s)
alltime=$((endtime-startime))
min_time=$((alltime/60))
echo -e "--------------------->perturb_posfile_sets.sh 完美完成工作,总共用时$min_time min OVO"

