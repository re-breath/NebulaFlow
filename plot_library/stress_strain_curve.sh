#!/bin/bash
#该脚本用于提取thermo.out文件中的应力-应变曲线，自动判断是计算什么样的应力应变曲线（对文件名进行关键词检索）
#from rebreath

#注意，该脚本需要保证gpumd的thermo.out文件的格式不变，否则会导致提取结果错误，目前gpumd的thermo.out文件的格式为：
#column   1 2 3 4  5  6  7   8   9   10 11 12 13 14 15 16 17 18
#quantity T K U Pxx Pyy Pzz Pyz Pxz Pxy ax ay az bx by bz cx cy cz
#之前还有一个版本为12列的，目前已经删除，现在只有18列

date_file=${1:-"thermo.out"}

# 检查data_file文件有多少列，如果不是18列那就报错
if [ $(head -n 1 $date_file | wc -w) -ne 18 ]; then
    echo "Error: $date_file GPUMD新版本要求thermo.out必须是18列的文件，无法提取应力-应变曲线/n请检查文件格式是否正确，或者是否是旧版本的thermo.out文件"
    exit 1
fi

flag=0

check_directory_name() {
    dir_name=$(basename "$(pwd)")

    # 使用正则表达式检查名称下划线后的第一个字母是否为x、y或z
    if [[ $dir_name =~ _([xyz]) ]]; then
        flag=${BASH_REMATCH[1]}
        echo '找到 >.O ！'"检索的结果为计算 ${BASH_REMATCH[1]} 方向的 stress-strain 曲线"
    else
        echo "文件夹的格式不正确哦，其中必须有_x、_y或_z "
    fi
}

process_data_x() {
    Lx_0=$(awk 'NR==1 {print $10}' $date_file)
    output_file="stress_strain_curve.txt"
    if [ -z "$Lx_0" ]; then
        echo "Error: Lx_0 is not set."
        return 1
    fi

    awk -v Lx_0="$Lx_0" ' 
        function abs(x) { return x < 0 ? -x : x; }
        {
            # 计算应变，ε = (Lx - Lx_0) / Lx_0
            strain = ($10 - Lx_0) / Lx_0
            # 应力σ已在Px中以GPa给出
            stress = $4
            stress_abs = abs(stress)
            print strain, stress_abs
        }
    ' "$date_file" > "$output_file"
}

process_data_y() {
    Ly_0=$(awk 'NR==1 {print $15}' $date_file)
    output_file="stress_strain_curve.txt"
    if [ -z "$Ly_0" ]; then
        echo "Error: Ly_0 is not set."
        return 1
    fi

    awk -v Ly_0="$Ly_0" ' 
        function abs(x) { return x < 0 ? -x : x; }
        {
            strain = ($15 - Ly_0) / Ly_0
            stress = $5
            stress_abs = abs(stress)
            print strain, stress_abs
        }
    ' "$date_file" > "$output_file"
}

process_data_z() {
    Lz_0=$(awk 'NR==1 {print $18}' $date_file)
    output_file="stress_strain_curve.txt"
    if [ -z "$Lz_0" ]; then
        echo "Error: Ly_0 is not set."
        return 1
    fi

    awk -v Lz_0="$Lz_0" ' 
        function abs(x) { return x < 0 ? -x : x; }
        {
            strain = ($18 - Lz_0) / Lz_0
            stress = $6
            stress_abs = abs(stress)
            print strain, stress_abs
        }
    ' "$date_file" > "$output_file"
}

check_variable_and_process_data() {
    case $flag in
        x)
            echo "下面开始计算x轴的应力应变曲线"
            process_data_x
            ;;
        y)
            echo "下面开始计算y轴的应力应变曲线"
            process_data_y
            ;;
        z)
            echo "下面开始计算z轴的应力应变曲线"
            process_data_z
            ;;
        *)
            echo "网络出错，请待会重试 ~"
            ;;
    esac
}

replot() {
#该函数将指定文件的前两列进行画图
    if [ "$#" -ne 1 ]; then
        echo "用法: replot <数据文件>"
        return 1
    fi

    python3 <<EOF
import matplotlib.pyplot as plt
import sys

datafile = "$1"

x, y = [], []
with open(datafile, 'r') as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) < 2:
            continue  # 跳过不正确的行
        x_val, y_val = map(float, parts[:2])
        x.append(x_val)
        y.append(y_val)

plt.figure(figsize=(10,6))
plt.plot(x, y, marker='o', linestyle='-',markersize=3)
plt.xlabel('strain')
plt.ylabel('stress(Gpa)')
plt.savefig('plot.png', format='png')
EOF
}

check_directory_name
check_variable_and_process_data
#replot stress_strain_curve.txt
cp ~/.rebreath/plot_library/plot_stress_strain_curve.py .
python3 plot_stress_strain_curve.py
