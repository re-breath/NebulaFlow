#!/bin/bash
#from rebreath
#-------------------------------->3<-----------------------------------
#本脚本为自动化的其中一环，执行序列为3号（第3个执行该脚本）
#本脚本的作用：将序列2的流程产生的构型都进行一次单点能
#处理的方式：进入当前文件夹的所有的文件夹，构建工作文件夹（一般为single_point_energy_*），序列2产生的构型库中的每个POSCAR-*将会在单点能库建立一个工作文件夹来计算每个构型的单点能
set -e
#-----------------------独立使用该脚本下面的这部分无需太注意-------------------------
core_num=${core_num:=1}
vasp_exe=${vasp_exe:="vasp"}
if [ -n "$single_point_energy_incar" ]
then
    singlepointenergy_inp_file=$single_point_energy_incar
else
    singlepointenergy_inp_file=INCAR
fi
#-------------------------------------------------------------------------------

#--------------------->workstation
core_num=$core_num
singlepointenergy_inp_file=$singlepointenergy_inp_file
kpoint_file=KPOINTS
pseudopotential_file=POTCAR
vasp_exe=${vasp_exe:="vasp_std"}
kspacing=${kspacing:=0}
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

add_kspacing_to_incar(){
  mv KPOINTS KPOINTS_backup
  if [ -z "$(grep "KSAPCING" INCAR)" ];then
     sed -i "\$a\KSPACING = $kspacing "  INCAR
  else
     sed "/KSPACING/c\KSPACING = $kspacing/" INCAR
  fi
}

init_adress=$PWD
startime=$(date +%s)
for i in $(find $init_adress/ -mindepth 1 -maxdepth 1 -type d -name "*aimd*"); do
         cd $i
         
         spe_dir=single_point_energy_$(basename $i)
         mkdir -p $spe_dir
         for ii in $(find $i  -mindepth 1 -name "POSCAR-*" | sort -t '-' -k2 -n);do
                
                filename=singlepointenergy_config-$(basename $ii | cut -c 8- )
                mkdir -p  $i/$spe_dir/$filename
                cp $ii $i/$spe_dir/$filename/POSCAR
                cp $init_adress/$singlepointenergy_inp_file  $i/$spe_dir/$filename/INCAR
                cp $init_adress/$kpoint_file  $i/$spe_dir/$filename/KPOINTS
                cp $init_adress/$pseudopotential_file  $i/$spe_dir/$filename/POTCAR
                cd $i/$spe_dir/$filename
                if [ ! -z $kspacing ];then
                    add_kspacing_to_incar
                fi
                check_vasp_complete
        done
        cd $init_adress
done
endtime=$(date +%s)
alltime=$((endtime-startime))
min_time=$((alltime/60))
echo "---------------------------->single_point_energy_vaspset_auto.sh的任务已完成,总共用时$min_time min"

