#!/bin/sh
for file in `find . -type f -name "POSCAR-*"`;do
    dirname=`dirname $file`
    mv $file $dirname/POSCAR
done

echo "name.sh完成任务。"
