#!/bin/sh
#本脚本将会进入将当前文件夹的所有的vasprun.xml*的文件进行改名，改名为所属文件夹的一部分,并且将会将其全部拖到当前文件夹
init_address=$(pwd)
for file in $(find $init_address -type f -name "vasprun.xml*");do
    addressname=$( dirname $file )
    cd $addressname
    folder_name=$(basename $(pwd) | cut -c 5-7)
    #cut将会截取当前文件夹的字符
    mv $file $addressname/vasprun.xml-$folder_name
    cp $addressname/vasprun.xml-$folder_name  $init_address
    cd $init_address
done
echo " name_mul.sh 任务完成 >_O"
