#!/bin/bash

# 该脚本的运行方法 source run.sh >nohup.log 2>&1 &

# 该函数用于重复运行gpumd,将会创建n个文件夹，将model.xyz,nep.txt,run.in复制到每个文件夹中
repeat_buildgpumd(){
    local n=$1
    for ((i=1;i<=$n;i++))
    do
        mkdir -p run$i
        cp model.xyz nep.txt run.in run$i/
    done
}

# 该函数用于重复运行gpumd,将会创建n个文件夹，将model.xyz,nep.txt,run.in复制到每个文件夹中
repeat_rungpumd(){
    local n=$1
    for ((i=1;i<=$n;i++))
    do
        cd run$i
        free_time_run 'nohup gpumd 2>&1 &'
        sleep 5
        cd ..
    done
}

# 进入当前所有的文件夹中，构建run1，run2，run3，run4，run5五个文件夹，将model.xyz,nep.txt,run.in复制到每个文件夹中，而后排队运行gpumd
for i in $PWD/*/
do
    cd $i
    repeat_buildgpumd 5
    repeat_rungpumd 5
done