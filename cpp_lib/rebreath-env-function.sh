#rebreath的函数库
#环境变量部分的函数
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

export PATH=$PATH:$HOME/.rebreath
alias loadenv='source ~/.bashrc'
alias modenv='vim ~/.bashrc'
export PYTHONPATH=$PYTHONPATH:$HOME/.rebreath/nebula_pylib
#处理数据集数据的快捷键，可以使用来处理dump文件的数据


get_energy(){
#获得xyz文件中的所有的能量
#grep "attice" $xyz_file |  grep -oE '\b\w+="[^"]+"|\b\w+=[^ ]+\b' |grep 'energy'
    local xyz_file=${1:-'dump.xyz'}
    grep "attice" $xyz_file |  grep -oE '\b\w+="[^"]+"|\b\w+=[^ ]+\b' |grep -oE '[Ee]nergy.*' |grep -oE '[-]?[0-9]+[\.]?[0-9]+'
}

get_Lattice(){
#得到xyz文件中所有的晶格参数
    local xyz_file=${1:-'dump.xyz'}
    grep "attice" $xyz_file |  grep -oE '\b\w+="[^"]+"|\b\w+=[^ ]+\b' | grep "attice"

}
get_virial(){
    local xyz_file=${1:-'dump.xyz'}
    grep "attice" $xyz_file |  grep -oE '\b\w+="[^"]+"|\b\w+=[^ ]+\b' | grep "irial"
}
get_configs_num(){ 
#获得构型的数量
    local xyz_file=${1:-'train.xyz'}
    grep "attice" $1 |wc -l
}

get_V(){
#对xyz文件的体积进行计算，并输出为1列
    local xyz_file=${1:-'dump.xyz'}
    get_Lattice $xyz_file |sed 's/"/ /g' | awk '{printf "%.15f\n",$2*$6*$10}'
}
average_file(){
#将输入文件进行平均,使用C++实现，平均后输入到average.out文件中
    local cpplib="$HOME/.rebreath/cpp_lib"
    first_file="$1"
    avg_name=${1:0:4}
    g++ $cpplib/averagefiles.cpp -o average_file
    ./average_file "$@" 
    rm -f average_file   
}

average_file_s(){
#将输入文件进行平均,使用shell实现
    files_num="$#"
    tmp_file='temp'
    first_file="$1"
    avg_name=${1:0:4}
    paste "$@" >$tmp_file
    # 逐行读取处理，无论多少行都能处理
    while IFS= read -r line
    do
        echo "$line" | awk -v files_num="$files_num" '{
            sum=0;
            for(i=1; i<=NF/files_num; i++){
                sum=0;
                for(j=0; j<files_num; j++){
                    sum += $(i+j*NF/files_num)
                }
                printf "%.15f\t", sum/files_num;
            }
            printf "\n";
        }'
    done < "$tmp_file" >average_${avg_name}.out
    rm -f "$tmp_file"
}
average_file_c(){
#将输入文件进行平均,使用C++实现，平均后输入到average.out文件中,后将文件改名
    local cpplib="$HOME/.rebreath/cpp_lib"
    first_file="$1"
    avg_name=${1:0:4}
    g++ $cpplib/averagefiles.cpp -o average_file
    ./average_file "$@" 
    rm -f average_file
    mv average.out average_${avg_name}.out
}


check_hnemd_thermo(){
#检查hnemd的thermo输出文件，并输出平均值,输出最后的晶格参数
    average_file_s thermo_*
    echo $PWD
    echo Lattice :
    tail -n 1 average_ther.out | awk '{print $10,$11,$12}'
}

use_agent(){
    export http_proxy="http://192.168.1.10:3128"
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

plt.plot(x, y, marker='o', linestyle='-')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig('plot.png', format='png')
EOF
}

proxy_download() {
#该函数的作用为使用代理节点下载文件
    local remote_server="dhk@10.14.0.15" 
    local remote_path="/home/dhk/rebreath/proxy_downfiles" 
    local remote_temp_path="${remote_path}/temp" 

    # 用ssh在远程服务器的临时目录创建目录并下载文件
    ssh -i ~/.ssh/id_rsa ${remote_server} "mkdir -p ${remote_temp_path} && cd ${remote_temp_path} && $1"
    
    # 检查远程下载命令的返回值
    if [[ "$?" -ne 0 ]]; then
        echo "文件下载失败，请检查命令并重试。"
        return 1
    fi
    
    # 使用scp从临时文件夹复制所有文件到本地当前目录下
    scp -rv -i ~/.ssh/id_rsa ${remote_server}:${remote_temp_path}/* .
    
    # 检查scp命令的返回值
    if [[ "$?" -ne 0 ]]; then
        echo "文件复制失败，请检查远程服务器上的路径和文件名。"
        return 1
    fi

    # 清空远程服务器中的临时目录
    ssh -i ~/.ssh/id_rsa ${remote_server} "rm -rf ${remote_temp_path}/*"

    echo "文件已成功下载到本地并清理了临时文件夹。"
}

#关闭vasp组件未使用警告
export OMPI_MCA_btl_base_warn_component_unused=0

free_time_run() {
#监测到空闲的gpu后进行任务
  while true; do
    for gpu_id in $(nvidia-smi --query-gpu=index --format=csv,noheader,nounits); do
      mem_used=$(nvidia-smi --id=$gpu_id --query-gpu=memory.used --format=csv,noheader,nounits)
      if [ $mem_used -lt 200 ]; then
          export CUDA_VISIBLE_DEVICES=$gpu_id
          echo "Running task on GPU $gpu_id"
          eval $1
          break 2
      fi
      done
      echo "No free GPU found, waiting..."
      sleep 60
  done
  date '+%Y-%m-%d %H:%M:%S' >> run_train-file.log
  echo -e "执行 $1 \n" >> run_train-file.log    
}

source $HOME/.rebreath/re_function.sh

#调用ovito转换文件的格式
xyz_to_poscar() {
    python3.10 -c "from ovito.io import import_file, export_file; pipeline = import_file('$1'); export_file(pipeline, 'POSCAR_convered', 'vasp')"
}

poscar_to_xyz() {     
       	python3.10 -c "from ovito.io import import_file, export_file; pipeline = import_file('$1'); export_file(pipeline, 'model_conversed.xyz', 'xyz',columns=['Particle Type', 'Position.X', 'Position.Y', 'Position.Z'])"
}


plot_stress_strain_curve(){
#自动检测应变的轴，输出数据到文件中，并进行简单的画图
    source $HOME/.rebreath/plot_library/stress_strain_curve.sh  
}

plot_mul_stress_strain_curve(){
#该函数为自动检测多个deform_{xyz}+文件并自动将其数据进行处理画出图来
    source $HOME/.rebreath/plot_library/auto_plot_xyz_strain_stress_curve.sh
}

xyz_to_cssr() {
# 内嵌的Python脚本
    python3 - "$1" "$2" << 'EOF'
from ase.io import read, write

# 转换函数
def convert(file_in, file_out):
    # 读取XYZ文件
    atoms = read(file_in)
    
    # 写出CSSR文件
    write(file_out, atoms, format='cssr')

# 主函数
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print("用法: convert_xyz_to_cssr <input.xyz> <output.cssr>")
        sys.return(1)

    # 执行转换
    convert(sys.argv[1], sys.argv[2])

EOF
}

xyz_to_cif () 
#该函数可能会出问题
{ 
    python3 - "$1" "$2" <<'EOF'
from ase.io import read, write

def convert(file_in, file_out):
    # 读取XYZ文件
    atoms = read(file_in)
    
    # 输出为CIF文件
    write(file_out, atoms)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print("Usage: python $0 <input.xyz> <output.cif>")
        sys.return(1)
    
    file_in = sys.argv[1]
    file_out = sys.argv[2]
    
    convert(file_in, file_out)
EOF

}

plot_ultimate_nep(){
#画出调色后的nep图
    cp $HOME/.rebreath/plot_library/plot_nep_results_ultimate.py .
    python3 plot_nep_results_ultimate.py
    rm -f plot_nep_results_ultimate.py
}

compute_elastic_moduli(){
#使用calorine计算弹性模量 使用方法  compute_elastic_moduli   nepfile
    loc compute_lib=/home/dhk/.rebreath/compute_lib
    loc nep_file=${1:}
    sed "s/nepfile/$1/g" $compute_lib/calorine_compute_elastic.py > calorine_compute_elastic.py
    python3 calorine_compute_elastic.py > elastic_calorine.txt
    rm -f calorine_compute_elastic.py 
    cat elastic_calorine.txt
}
plot_nep(){
#画出结果的图
    python3 ~/.rebreath/plot_library/hplt_nep_results.py
}

plot_hnemd(){
#该函数用来画hnemd的图像，自动识别路径中的的_{xyz},以此来对特定方向的hnemd画图
    lib_address="$HOME/.rebreath/plot_library/"
    
    find $PWD -type f -regex '.*kappa.out' | while read -r file;do
      local hnemd_direct=0
      if [[ $file =~ _([xyz]) ]];then 
         hnemd_direct=${BASH_REMATCH[-1]} 
         echo -e '\n——>'"找到hnemd的方向了，该处的hnemd计算的为 $hnemd_direct 方向的热导率" 'O.<'
         else
         echo "没有找到hnemd的驱动力方向，请检查是否按约定命名方式命名 oooooo"
         return 1
      fi
    file_address=$(dirname $file)
    cd $file_address
    echo $PWD
    python3 "$lib_address/plot_hnemd_${hnemd_direct}.py" 2>&1
    cd - >/dev/null
    done
}

plot_mul_hnemd(){
#找到average文件夹，自动识别该文件路径中的方向，对其中的文件进行画图
#注意格式需要为kappa_[0-9]+.out
    local lib_address="$HOME/.rebreath/plot_library/"
    local hnemd_direct=0
    local kappa_num=$(find $PWD -type f -regex '.*kappa_[0-9]+.out' |wc -l)
    find $PWD -type d -regex '.*average_hnemd/kappa' | while read -r workdir;do
        cd $workdir
        if [ -f "average.out" ];then
            cp average.out kappa.out
        fi

        if [[ $workdir =~ _([xyz]) ]];then 
            hnemd_direct=${BASH_REMATCH[-1]} 
            echo -e '\n——>'"找到hnemd的方向了，该处的hnemd计算的为 $hnemd_direct 方向的热导率" 'O.<'
            else
            echo "没有找到hnemd的驱动力方向，请检查是否按约定命名方式命名 oooooo"
            return 1
        fi
        cat $lib_address/plot_hnemd_mul_${hnemd_direct}.py |sed "/for i in range(/c\for i in range(${kappa_num}):" > plot_hnemd_mul_${hnemd_direct}.py
        python3 "plot_hnemd_mul_${hnemd_direct}.py" 2>&1
        cd - >/dev/null
    done
}

deal_hnemd_data(){
#处理数据hnmed计算完成的数据，注意需要符合一定的命名规则
#该函数能够将hnemd_[0-9]+类型的文件夹进行整理，对每个将文件夹中的kappa.out进行平均并将其整合到一张图中
#场地要求：使用该函数的地方需要有很多的hnemd_[0-9]+类型的文件夹
    local lib_address="$HOME/.rebreath/deal_data/"
    source ${lib_address}/deal_hnemd_data.sh 
    #cd /average_hnemd/kappa
    plot_mul_hnemd
}

find_column_max(){
#检查文件指定列的最大值 使用方式为 find_column_max 4 thermo.out 4 (检查文件第四列的最大值) 
    if [[ $# -ne 2 ]]; then
    echo "错误：需要提供文件名和列号作为参数。"
    return 1
    fi

    filename=$1
    column=$2

    # 检查文件是否存在且可读
    if [[ ! -f "$filename" || ! -r "$filename" ]]; then
        echo "错误：文件 '$filename' 不存在或不可读。"
        return 1
    fi
    # 验证列号是否为正整数
    if [[ ! $column =~ ^[0-9]+$ ]]; then
        echo "错误：列号 '$column' 必须是正整数。"
        return 1
    fi
    # 使用双引号包含变量，并在 awk 脚本中安全地使用变量
    awk -v col="$column" 'BEGIN {max=0} 
        NR > 0 && ($col > max) {max=$col} 
        END { printf("第 %d 列的最大值为 ：%f\n", col, max) }' "$filename"
}

find_column_abs_max(){
#检查文件指定列的绝对值最大值 使用方式为 find_column_max thermo.out 4 (检查文件第四列的最大值) 
    if [[ $# -ne 2 ]]; then
    echo "错误：需要提供文件名和列号作为参数。"
    return 1
    fi

    filename=$1
    column=$2

    # 检查文件是否存在且可读
    if [[ ! -f "$filename" || ! -r "$filename" ]]; then
        echo "错误：文件 '$filename' 不存在或不可读。"
        return 1
    fi
    # 验证列号是否为正整数
    if [[ ! $column =~ ^[0-9]+$ ]]; then
        echo "错误：列号 '$column' 必须是正整数。"
        return 1
    fi
    # 使用双引号包含变量，并在 awk 脚本中安全地使用变量
    awk -v col="$column" 'BEGIN {max=0} 
        NR > 0 && (sqrt(($col)^2) > sqrt((max)^2)) {max=$col} 
        END { printf("第 %d 列的绝对值最大值为：%f\n", col, max) }' "$filename"
}

cell_expansion(){
#使用gpumd进行扩胞 使用的方式类似 cell_expansion 10 10 1
    local nepfile=${nepfile:='nep.txt'}
    local xyzfile=${xyzfile:='model.xyz'}
    tempfile="temp_$(date +%s%3N)"
    mkdir -p $tempfile
    cp $nepfile $xyzfile  $tempfile/
    cd $tempfile/
    cat >run.in << EOF
replicate $1 $2 $3 
potential       $nepfile
ensemble        nve
time_step       0
dump_exyz       1
run             1
EOF
gpumd > /dev/null 2>&1
if [ -f  "dump.xyz" ];then
   cp dump.xyz  ../expanded.xyz
fi
cd - >/dev/null
rm -rf $tempfile
echo "已经完成 $1  $2  $3 扩胞"
}


relib() {
# rebreath库的管理员   
#示例使用方式 relib -v -m hp pre sh

    lib="$HOME/.rebreath"
    view=0
    multi=0
    args=()

    # 处理参数
    while [ $# -gt 0 ]; do
        case "$1" in
            -v)
                view=1
                ;;
            -m)
                multi=1
                ;;
            *)
                args+=("$1")
                ;;
        esac
        shift
    done

    # 检查是否有实际的查找关键词
    if [ ${#args[@]} -eq 0 ]; then
        echo "请注意提供至少一个查找关键词。"
        exit 521
    fi

    # 用于存储中间结果的文件
    tmpfile=$(mktemp)

    # 初始查找
    find "$lib" -type f -name "*${args[0]}*" > "$tmpfile"

    # 逐层次查找
    for keyword in "${args[@]:1}"; do
        new_tmpfile=$(mktemp)
        while IFS= read -r line; do
            find "$line" -type f -name "*$keyword*" >> "$new_tmpfile"
        done < "$tmpfile"
        mv "$new_tmpfile" "$tmpfile"
    done

    # 执行查找和相应操作
    if [ $view -eq 1 ]; then
        cat "$tmpfile"
    else
        while IFS= read -r line; do
            cp "$line" ./
        done < "$tmpfile"
    fi

    # 移除临时文件
    rm -f "$tmpfile"
}



gpumdstart(){
#本程序建立是为了简化在曙光上使用gpumd的流程，该程序默认任务名为执行命令的文件夹的名字，默认使用1dcu与1cpu,也可以使用-n参数指定使用dcu的数量
dcu_num=1
for arg in "$@"; do 
   if [[ $arg = "-n" ]]; then
     dcu_num=${2:-1}
     break
   fi
done

job_name=$(basename "$PWD")

echo "#!/bin/bash
#SBATCH -p xahdnormal
#SBATCH -N  1
#SBATCH --ntasks-per-node=$dcu_num
#SBATCH --gres=dcu:$dcu_num
#SBATCH --time 144:00:00
#SBATCH --comment=GPUMD
#SBATCH -o $PWD/std.out.%j
#SBATCH -e $PWD/std.err.%j

# MARK_CMD
source /work/home/rebreath/sbatch_need/gpumd_env.sh
/work/share/acmtrwrxv5/GPUMD/src/gpumd" | sbatch -J "$job_name"
EOF
}

#+++++++++++++++++++++++++++nep训练集相关+++++++++++++++++++++++++++++++++

screening_reasonable_forces(){
#筛选nep的训练集，将训练集的合理的构型提取出
#使用的方法为 screening_reasonable_forces xyzfile min_force max_force
    local deal_lib="$HOME/.rebreath/deal_data/"
    python3 $deal_lib/elect_rely_force.py $1  $2  $3
}
screening_reasonable_energy(){
#筛选nep的训练集，将训练集的合理的能量提取出
#使用的方法为 screening_reasonable_energy xyzfile min_energy max_energy
    local deal_lib="$HOME/.rebreath/deal_data/"
    python3 $deal_lib/elect_rely_energy.py $1  $2  $3
}

screening_reasonable_virial(){
#筛选nep的训练集，将训练集的合理的位力提取出
#使用的方法为 screening_reasonable_virial xyzfile min_strain max_strain
    local deal_lib="$HOME/.rebreath/deal_data/"
    python3 $deal_lib/elect_rely_virial.py $1  $2  $3
}



#+++++++++++++++++++++++++++vasp计算相关++++++++++++++++++++++++++++++++++

add_kspacing_to_incar(){
#如果kspacing的值不为零，则添加kspacing
#使用的方法为 add_kspacing_to_incar 0.2
    kspacing=${1:=0.2}
    if [ -f "KPOINTS" ];then
    mv KPOINTS KPOINTS_backup
    fi
    if [ -z "$(grep "KSPACING" INCAR)" ];then
        echo "KSPACING = $kspacing " >> INCAR
    else
        sed -i -E "/^KSPACING.*=.*/c\KSPACING = $kspacing" INCAR
    fi
}

check_vasp_complete() {
#检查档期那目录的vasp是否完成了计算，如果没有完成就会运行vasp
    vasp_exe=${1:-'vasp'}
    core_num=${2:-'1'}
    if [ ! -f  "OUTCAR" ];then
        echo "当前目录下没有找到OUTCAR文件，正在生成..."
        mpirun -np $core_num $vasp_exe
    elif [ -z "$(grep "General timing and accounting informations for this job" OUTCAR)" ];then
        echo "$PWD 目录下OUTCAR不完整，开始重新计算"
        mpirun -np $core_num $vasp_exe
    fi
}

run_all_vasp_job(){
#检查目录下所有的train-*文件夹，并运行vasp
    starttime=$(date +%s)
    init_address=$PWD
    core_num=${core_num:=1}
    vasp_exe=${vasp_exe:=vasp}
    startime=$(date +%s)
    kspacing=${kspacing:=0.2}

    for i in $(find $init_address/ -type d -name "train-*" |sort)
    do
    cd $i
    add_kspacing_to_incar $kspacing
    if [ ! -f "OUTCAR" ];then
        echo "文件夹$i下没有找到OUTCAR文件，正在对其进行计算"
        check_vasp_complete
        echo "$i的退出码为$?" >> $init_address/run_train-file.log
    fi
    complete=$(grep "General timing and accounting informations for this job" OUTCAR)
    if [ -z "$complete" ];then
    echo "文件夹$i下的OUTCAR文件不完整，正在重新计算"
    #echo 'SYMPREC=1E-03' >> $i/INCAR
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
}

load_single_point_energy_dir(){
#检查该目录下所有的POSCAR类的文件，计算其单点能
#要求是poscar文件的目录下面有incar与potcar，另外incar
    orig_dir=$PWD
    find $orig_dir -type f -regex ".*POSCAR.*" |while IFS= read -r i; do
    dname=$(dirname $i)
    cd $dname
    filename=$(basename $i)
    suffix=${filename#POSCAR}
    mkdir -p single_energy/train-$suffix
    cp $dname/POTCAR $dname/INCAR single_energy/train-$suffix
    if [ -f "$dname/KPOINTS" ];then
    cp $dname/KPOINTS  single_energy/train-$suffix
    fi
    cp $i single_energy/train-$suffix/POSCAR
    cd $orig_dir
done
}

deal_outcar_to_train(){
#寻找当前所有的OUTCAR文件，将其整理成train.xyz
#注意其寻找的是所有的OUTCAR文件，因此需要确保当前目录下所有OUTCAR文件都为单点能计算的OUTCAR
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
    for outcar in $(find $aim_address -type f -name "OUTCAR" |sort);do
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
    echo -e "--------------------->完美完成工作,总共用时$min_time min OVO"
}

echo "NebulaFlow library loaded" 'O.<'

#-----------------------------模板-----------------------------------

module load compiler/gcc/12.2.0 || true

compute_high_of_tree() {
# 创建一个临时文件夹来保存生成的文件
    local tmpdir=$(mktemp -d -t high_of_tree-XXXXXX)
    
    # 使用指定的临时文件夹
    pushd "$tmpdir"
    
    # 使用shell调用c++的模板
    cat > high_of_tree.cpp <<EOF
    #include <iostream>
    #include <cmath>
    #include <numbers>

    using namespace std;

    void high_of_tree(){
        const double pi  {std::numbers::pi};
        double high {};
        double h {};
        double d {};
        double angle {};
        cout << "计算树的高度" << endl;
        cout << "please enter h ,d ,angle：";
        
        cin >> h >> d >> angle;
        cout << "h = " << h << "  d = " << d << "  angle = " << angle << endl;
        cout << "树的高度是" << h + d * std::tan(angle * pi / 180) << endl;
    }

    int main(){
        high_of_tree();
        return 0;
    }
EOF

# 编译C++代码
g++ -std=c++20 high_of_tree.cpp -o high_of_tree

# 运行程序
./high_of_tree

# 从临时文件夹中返回到原始目录
popd

# 删除临时文件夹及其内容
rm -rf "$tmpdir"
}
