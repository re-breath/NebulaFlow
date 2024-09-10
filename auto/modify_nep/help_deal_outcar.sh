#!/bin/bash
#from rebreath
#用于辅助提取outcar的信息，在md模拟后计算单点能可能会遇到一个神秘的bug导致几乎无法提取信息，该脚本可辅助提取train-*文件夹的构型
# 假设deal_outcar_to_train_auto.sh脚本在当前目录下
DEAL_SCRIPT="deal_outcar_to_train_auto.sh"

init_address=$PWD
# 遍历所有以train-开头的子目录
n=1
for dir in $(find $init_address -type d  -name "train-*"); do
    if [ -d "$dir" ]; then
        # 复制deal_outcar_to_train_auto.sh到子目录
        cp $DEAL_SCRIPT $dir

        # 进入子目录并执行脚本
        cd $dir
        bash "$DEAL_SCRIPT"

        # 检查是否生成了train.xyz文件，如果是，则将其内容追加到最终的train.xyz文件中
        if [ -f "train.xyz" ]; then
            cat train.xyz >> $init_address/train.xyz
            # 可选：删除子目录中的train.xyz文件，如果你不想保留它们的话
            rm train.xyz
            rm deal_outcar_to_train_auto.sh
        fi

        # 返回原始目录
        cd $init_address
    fi
    n=$((n+1))
done

echo "All train.xyz files have been collected into the current directory."
echo -e "共有$n个构型添加"