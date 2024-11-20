#该部分的函数使用来处理gpumd，nep软件相关的函数


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


plot_nep(){
#画出结果的图
    python3 ~/.rebreath/plot_library/hplt_nep_results.py
}

plot_E_F_Vir_distribution(){
#画出类似nep的数据集中能量、力、位力分布的图
#使用方法：plot_E_F_Vir_distribution train.xyz

    local dump_file=${1:-train.xyz}
    > plot_E_F_Vir_distribution.py
    cat >> plot_E_F_Vir_distribution.py << EOF
import nebula
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

font = {'weight' : 'bold',  'size' : 7}
rc('font', **font)

configs = nebula.read_xyz('train_all.xyz')
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
E_avg = [float(config.energy / config.atom_num) for config in configs]
plt.hist(E_avg, bins=300, alpha=0.5,color='#835593FF')
plt.xlabel('Average Energy (eV/atom)',fontweight='bold')
plt.ylabel('Frequency',fontweight='bold')

# 画出力分布
plt.subplot(1, 3, 2)
color_force = ['#2CAC15FF','#2C3460FF','#1A87ABFF']
F_avg = np.concatenate([config.force for config in configs])
plt.hist(F_avg, bins=100, alpha=0.5,color=['#2CAC15FF','#2C3460FF','#1A87ABFF'])
plt.xlabel('Interatomic Force (eV/$\AA$)',fontweight='bold')

plt.ylabel('Frequency',fontweight='bold')

color_viral = ['#00CED1FF','#FFA500FF','#FF4500FF']
plt.subplot(1, 3, 3)
vir_avg = np.concatenate([np.array(config.virial[::3]) / config.atom_num for config in configs])
vir_avg = np.column_stack((vir_avg[::3],vir_avg[1::3],vir_avg[2::3]))

plt.hist(vir_avg, bins=300, alpha=0.5)
plt.xlabel('Average Viral (eV/atom)',fontweight='bold')
plt.ylabel('Frequency',fontweight='bold')   
plt.tight_layout()
plt.savefig('E_F_vir.png',dpi=600)
plt.close()
EOF
python3 plot_E_F_Vir_distribution.py
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

start_mul_hnemd(){
#该函数用来启动gpumd的hnemd方法来计算热导率，可以进行制定次数的热导率的计算
#使用方法：start_mul_hnemd nep.txt 6 1 将会使用gpumd进行6次热导率的计算，每个计算使用1个核
    local nepfile=${1:-nep.txt}
    local times=${2:-6}
    local times=$((times-1))
    local core_num=${3:-1}
    local init_address=$PWD
    for i in $(seq 1 $times);do
        mkdir -p hnemd_${i}
        cp $nepfile model.xyz run.in hnemd_${i}/
        cd hnemd_${i}
        #gpumd_dcu_394 -n $core_num
        cd $init_address
    done
}

start_mul_hnemd() {
    # 该函数用来启动gpumd的hnemd方法来计算热导率，可以进行制定次数的热导率的计算
    # 使用方法：start_mul_hnemd nep.txt 6 1 将会使用gpumd进行6次热导率的计算，每个计算使用1个核
    local nepfile=${1:-nep.txt}
    local times=${2:-6}
    local core_num=${3:-1}
    local init_address=$PWD

    # 生成大于10000的质数列表
    local prime_seeds=($(generate_large_primes $times 10001))

    for i in $(seq 1 $((times-1))); do
        mkdir -p hnemd_${i}
        cp $nepfile model.xyz run.in hnemd_${i}/
        cd hnemd_${i}

        # 替换 run.in 中的 seed 值
        sed -i "s/seed.*$/seed ${prime_seeds[$i]}/" run.in

        # 运行 gpumd (示例命令，根据实际情况调整)
        gpumd_dcu_394 -n $core_num

        cd $init_address
    done
}

deal_strain_fluctuation_to_elastic(){
#该函数用来处理应变波动法计算弹性模量得到的thermo文件，计算出其中的弹性模量
#使用方法：deal_strain_fluctuation_to_elastic 1200  # 1200为温度
    cp ~/.rebreath/deal_data/deal_strain_fluctuation.sh .
    bash deal_strain_fluctuation.sh $1 > elastics.txt
    rm -f deal_strain_fluctuation.sh
    cat elastics.txt
}

get_hnemd_data(){
    local init_address=$(pwd)
    local xyz="x y z"
    for i in $xyz;
    do 
      cd hnemd_${i}/average_hnemd/kappa/
      python3 plot_hnemd_mul_${i}.py 
      pwd
      cd $init_address

    done
}

verify_gpumd_result(){
#重新计算验算当前的gpumd算例
    local prediction_nep=$(ls |grep -oE "nep.*\.txt")
    local nepfile=${1:-$prediction_nep}
    local dir=$(basename $PWD)
    local verify_dir="verify_${dir}"
    if [ -d "$verify_dir" ];then
        rm -rf $verify_dir
    fi
    mkdir  $verify_dir
    cp "$nepfile" "$verify_dir"
    cp model.xyz run.in $verify_dir/
}


plot_ultimate_nep(){
#画出调色后的nep图
    cp $HOME/.rebreath/plot_library/plot_nep_results_ultimate.py .
    python3 plot_nep_results_ultimate.py
    rm -f plot_nep_results_ultimate.py
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

compute_elastic_moduli(){
#使用calorine计算弹性模量 使用方法  compute_elastic_moduli   nepfile
    local compute_lib=/home/dhk/.rebreath/compute_lib
    local nep_file=${1:-nep.txt}
    sed "s/nepfile/$nep_file/g" $compute_lib/calorine_compute_elastic.py > calorine_compute_elastic.py
    python3 calorine_compute_elastic.py > elastic_calorine.txt
    rm -f calorine_compute_elastic.py 
    cat elastic_calorine.txt
}



start_gpumd(){
#该函数用来启动gpumd的计算，对于dcu可以进行制定核数的计算,可以使用来代替gpumd的启动
    if [ $gpumd_exe -eq "gpumd" ]; then
        gpumd

    elif [ $gpumd_exe -eq "gpumdstart_dcu" ]; then
        dcu_num=${1:-1}
        gpumdstart_dcu -n $dcu_num
    else
        echo "错误：未知的 gpumd 启动方式，清修改 ~/.rebreath/.config 文件中的配置。"
        exit 521
    fi
}

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
#能够识别两种类型的文件，主要是特定的xyz文件与thermo.out文件
    local xyz_file=${1:-'dump.xyz'}
    if [[ $xyz_file =~ ".xyz" ]];then
    get_Lattice $xyz_file |sed 's/"/ /g' | awk '{printf "%.15f\n",$2*$6*$10}'
    elif [[ $xyz_file =~ ".out" ]];then
        python3 << EOF
import numpy as np
filename='$xyz_file'
data = np.loadtxt(filename)
def get_volume(data):
    lx = data[:,-3]
    ly = data[:,-2]
    lz = data[:,-1]
    return lx*ly*lz
V = get_volume(data)
for i in V:
    print(i)
EOF
    else
        echo "错误：未知的输入文件类型，请检查文件名后重试。"
fi
}

get_area_of_xy_and_volume() {
# 计算xy的晶面积,并顺便的计算出斜胞的体积(注意文件对象为gpumd风格的xyz文件)
# 输出一个文件 Volume_area_xy.txt, 第一列为斜胞体积，第二列为xy面的面积
# 注意该函数为临时函数，只计算xy平面的面积，为计算CaF2使用ase进行切片后的固定晶面所需，提交时候需要谨慎考虑
    local file=$1
    get_Lattice $1 |awk -F "=" '{print $2}' | sed 's/"//g'  > lattic.log
    python3 << EOF
import numpy as np
filename = 'lattic.log'
data = np.loadtxt(filename)
la = data[:,:3]
lb = data[:,3:6]
lc = data[:,6:9]
def get_area_xy(la,lb,lc):
    """
    计算xy面的面积(使用叉乘法)
    """
    areas = np.zeros(len(la))
    vols = np.zeros(len(la))
    for i in range(len(la)):
        temp = np.cross(la[i],lb[i])
        areas[i] = np.linalg.norm(temp)
        vols[i] = np.dot(lc[i],temp)

    return areas,vols

areas ,vols = get_area_xy(la,lb,lc)

np.savetxt('${file}_Volume_area_xy.txt',np.column_stack((vols,areas)))
EOF
    rm -f lattic.log
}


plot_E-frame(){
    # 该函数使用来直接画出能量的变化图
    get_energy $1 > E.log
    seq $(wc -l E.log | awk '{print $1}') > frame.log
    paste frame.log E.log > E-frame.txt
    replot E-frame.txt
    rm -f frame.log E.log
}

check_hnemd_thermo(){
#检查hnemd的thermo输出文件，并输出平均值,输出最后的晶格参数
    average_file_s thermo_*
    echo $PWD
    echo Lattice :
    tail -n 1 average_ther.out | awk '{print $10,$11,$12}'
}



# 特殊函数，与vasp联合部分，该部分主要使用来构建nep的数据集

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
    element_str=$(grep "VRHFIN" $outcar | awk -F"=" '{print $2}' |awk -F ":" '{print $1}')
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


plot_stress_strain_curve(){
#自动检测应变的轴，输出数据到文件中，并进行简单的画图
    source $HOME/.rebreath/plot_library/stress_strain_curve.sh  
}

plot_mul_stress_strain_curve(){
#该函数为自动检测多个deform_{xyz}+文件并自动将其数据进行处理画出图来
    source $HOME/.rebreath/plot_library/auto_plot_xyz_strain_stress_curve.sh
}
