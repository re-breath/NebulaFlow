#!/usr/bin/bash
init_address=$PWD
core_num=${core_num:=1}
vasp_exe=${vasp_exe:=vasp_std_gpu}
startime=$(date +%s)


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


for i in $(find $init_address/ -type d -name "train-*" |sort)
do
  cd $i
  add_kspacing_to_incar
  if [ ! -f "OUTCAR" ];then
    echo "文件夹$i下没有找到OUTCAR文件，正在对其进行计算"
    check_vasp_complete
    echo "$i的退出码为$?" >> $init_address/run_train-file.log
  fi
 complete=$(grep "General timing and accounting informations for this job" OUTCAR)
 if [ -z "$complete" ];then
   echo "文件夹$i下的OUTCAR文件不完整，正在重新计算"
   echo 'SYMPREC=1E-03' >> $i/INCAR
   check_vasp_complete
   echo "$i的退出码为$?" >> $init_address/run_train-file.log
fi
  cd $init_address
done
time_log=$(date)
echo "完成的时间为:  $time_log"
endtime=$(date +%s)
alltime=$((endtime-startime))
min_time=$((alltime/60))
echo -e "--------------------->run_train-file.sh 完美完成工作,总共用时$min_time min OVO"