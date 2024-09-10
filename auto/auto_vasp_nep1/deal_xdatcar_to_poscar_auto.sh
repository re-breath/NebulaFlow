#!/bin/bash
#from rebreath
#-------------------------------->2<-----------------------------------
#本脚本为自动化的其中一环，执行序列为2号（第2个执行该脚本）
#本脚本的作用：进入当前文件夹下所有的工作文件夹内寻找XDATCAR文件进行处理（放到工作文件夹下，运行该脚本）
#处理的方式：遍历工作的文件夹，进入每一个文件夹后将在当前文件夹创立一个文件夹（一般为configuration_files)，将从XDATCAR拆分出的构型放入其中
set -e
#-----------------------独立使用该脚本下面的这部分无需太注意-------------------------
inter_configuration_step=${inter_configuration_step:=50}
all_configuration=${all_configuration:=5000}

#-------------------------------------------------------------------------------

#-----------------------------下面为控制台部分----------------------
inter_configuration_step=$inter_configuration_step
all_configuration=$all_configuration
#-------------------------------------------------------------------

init_address=$PWD

for i in $(find $init_address/ -mindepth 1 -maxdepth 1 -type d -name "*aimd*"); do
    prefix_name=$(basename $(dirname $i/XDATCAR))
    cd $i
    if [ ! -f XDATCAR ];then
       echo 'Error : XDATCAR不存在(轨迹文件不存在），你再想想你跑的是不是AIMD'
       exit 1
    fi

    workdir=${prefix_name}_configurations
    mkdir -p $workdir
    sed '/Direct configuration=/{ s/ //g }'  $i/XDATCAR >xdatcar
    for ii in $(seq 1 $inter_configuration_step $all_configuration); do
        head -n 7 $i/XDATCAR > $workdir/POSCAR-$ii
        start_config="Directconfiguration=$ii$"
        #$表示末尾，这样可以使得其严格匹配Directconfiguration=1，而不是11等
        end_config="Directconfiguration=$((ii+1))$"
        echo "start_config: $start_config"
        echo "end_config: $end_config"
        sed -n "/$start_config/,/$end_config/p" $i/xdatcar >> $workdir/POSCAR-$ii
        sed  -i "/$end_config/d" $workdir/POSCAR-$ii
        echo "complete POSCAR-$ii"
    done
    cd $init_adress
done

echo '---------------------------->deal_xdatcar_to_poscar_auto.sh的任务已完成'





