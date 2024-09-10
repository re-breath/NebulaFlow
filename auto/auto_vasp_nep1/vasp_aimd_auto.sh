#!/bin/bash
#from rebreath
#-------------------------------->1<-----------------------------------
#本脚本为自动化的其中一环，执行序列为1号（第1个执行该脚本）
#本脚本的作用：建立工作文件夹自动进行多次不同温度下的AIMD
#处理的方式：建立多个文件夹将vasp的输入文件们放到其中并进行vasp计算
set -e
#下面为对于自动化控制台参数进行判定的判断
startime=$(date +%s)
init_adress=$PWD
vasp_exe=${vasp_exe:=vasp}
#-----------------------独立使用该脚本下面的这部分无需太注意-------------------------
core_num=${core_num:=1}
Temp=${Temp:=1000}
template_file=${aimd_incar:=INCAR}
pos_file=${aimd_poscar:=POSCAR}
#---------------------------------------------------------------------------

 #--------------------->workstation
template_file=$template_file
pos_file=$pos_file
kpoint_file=KPOINTS
pseudopotential_file=POTCAR

#-----------------------------------------
check_vasp_complete() {
  if [ ! -f  "OUTCAR" ];then
    echo "当前目录下没有找到OUTCAR文件，正在生成..."
    mpirun -np $core_num $vasp_exe
  elif [ -z "$(grep "General timing and accounting informations for this job" OUTCAR)" ];then
    echo "$PWD 目录下OUTCAR不完整，开始重新计算"
    mpirun -np $core_num $vasp_exe
  fi
}
first=${aimd_poscar:-train}
cd $init_adress
for ii in $Temp ; do
    work_dir=${first}_aimd_${ii}
        mkdir -p $work_dir
    sed  "s/replace/${ii}/g" $template_file > $work_dir/INCAR

    cp $kpoint_file $work_dir/KPOINTS
    cp $pseudopotential_file $work_dir/POTCAR
    cp $pos_file $work_dir/POSCAR
    cd $work_dir
    if [ ! -f "xdatcar" ];then
    check_vasp_complete
    fi
    cp OUTCAR outcar || true
    cd $init_adress
done
endtime=$(date +%s)
alltime=$((endtime-startime))
min_time=$((alltime/60))
echo -e "--------------------->vasp_aimd_auto.sh 完美完成工作,总共用时$min_time min OVO"
