#!/bin/bash
#from rebreath
#-------------------------------->4<-----------------------------------
#本脚本为自动化的其中一环，执行序列为4号（第4个执行该脚本）
#本脚本的作用：提取当前文件夹的所有的OUTCAR信息，转化为训练集（train.xyz）
#处理的方式：找出当前文件夹内所有的OUTCAR文件，将需要的信息按nep需要的格式放到train.xyz中，（注意本脚本每个OUTCAR只提取一个构型到训练集）
#注意：本脚本的使用格式为bash 本脚本名，因为其中有的语法dash无法解释，很多机器的sh是指向dash的，因此sh会出错
#本脚本的使用方式bash deal_outcar_to_train_auto.sh [  指定地址或默认地址 ]（地址不填为默认地址，填写后则处理指定的地址）
set -e
startime=$(date +%s)
writ_file="train.xyz"
init_address=$PWD
if [ -n "$1" ]
then
       aim_address=$1
else
      aim_address=$init_address
fi
N_case=$(find $aim_address -type f -name "OUTCAR" | wc -l)
echo "当前文件夹内有$N_case个OUTCAR文件"
if [ $N_case -eq 0 ]
then
  echo "当前文件夹内有$N_case个OUTCAR文件"
  exit 1
fi
N_cout=0
for outcar in $(find $aim_address -type f -name "OUTCAR");do
   all_atom=$(grep "number of ions" $outcar | tail -n 1 | awk '{ print $NF }')
   if [ -z "$all_atom" ];then
     continue
   fi

   if [ $? -ne 0 ]; then
        echo "Error processing OUTCAR file: $outcar" >> error_log.txt
        continue # 跳到下一个文件
    fi
   config_type=$(basename $(dirname $outcar))
   weight=1.0
   lattice=$(grep -A 7 "VOLUME and BASIS-vectors are now" $outcar | tail -n 3 | awk '{ print $1,$2,$3 }' | xargs)
   if [ -z "$lattice" ];then
     continue
   fi
   energy=$(grep  'free  energy   TOTEN' $outcar | awk -F "=" '{ print $NF }' |tail -n 1 | sed 's/ //g')
   virial=$(grep -A 20 "FORCE on cell =-STRESS" $outcar | grep "Total" | tail -n 1 | awk '{ print $2,$5,$7,$5,$3,$6,$7,$6,$4 }')
   echo "$all_atom" >> $writ_file
   if   [ -n "$virial" ]
   then
     echo "Config_type=$config_type Weight=$weight Lattice=\"$lattice\" Energy=$energy Virial=\"$virial\" Properties=species:S:1:pos:R:3:force:R:3" >> $writ_file
   else
     echo "Config_type=$config_type Weight=$weight Lattice=\"$lattice\" Energy=$energy Properties=species:S:1:pos:R:3:force:R:3" >> $writ_file
   fi
   element_str=$(grep "VRHFIN" $outcar | awk -F"=" '{print $2}' |awk -F":" '{print $1}')
   ion_element_array=($element_str)
   ion_num_array=($(grep "ions per type"  $outcar | tail -n 1 | awk -F"=" '{print $2}'))
   
   for((i=0;i<${#ion_element_array[@]};i++))
   do
     for((j=0;j<${ion_num_array[i]};j++))
     do
        echo "${ion_element_array[$i]}" >> atom_name
     done
   done
   grep -A $((all_atom+1)) "TOTAL-FORCE (eV/Angst)" $outcar | tail -n $all_atom >> atom_pos
   paste atom_name atom_pos >> $writ_file
   rm atom_name atom_pos
   N_cout=$((N_cout+1))
done
echo "当前工作目录为$init_address"
echo "训练集新增构型数目为 $N_cout" '让我们欢迎新增的构型>o<'
echo "完成数目任务数与总任务之比 $N_cout/$N_case"
endtime=$(date +%s)
alltime=$((endtime-startime))
min_time=$((alltime/60))
echo -e "--------------------->deal_outcar_to_train_auto.sh 完美完成工作,总共用时$min_time min OVO"
