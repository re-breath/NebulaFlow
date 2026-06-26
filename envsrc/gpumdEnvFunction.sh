# ======================================================================
# File:         gpumdEnvFunction.sh
# Project:      NebulaFlow
# Description:  GPUMD/NEP训练/热导率相关函数库
#               Functions for GPUMD simulation, NEP potential training,
#               HNEMD thermal conductivity, elastic moduli, and phonon.
# Author:       rebreath
# Dependencies: gpumd, python3, numpy, ase, calorine, nvidia-smi
# ======================================================================


# =============================================================================
# SECTION 1: NEP Training Set Management / NEP训练集管理
# =============================================================================

# ---------------------------------------------------------------------------
# Function: screening_reasonable_forces
# 功能: 筛选NEP训练集中力在指定范围内的构型
#       从xyz文件中提取原子间力在[min, max]范围内的构型，用于构建高质量训练集。
# 场景: 当你从AIMD或DFT计算生成了大量构型后，需要剔除受力异常（过大或过小）
#       的构型，保证训练集物理合理性。
# Usage: screening_reasonable_forces <xyzfile> <min_force> <max_force>
# Example:
#   screening_reasonable_forces train.xyz -10 10
#   # 筛选train.xyz中每个原子受力在[-10, 10] eV/A范围内的构型
# ---------------------------------------------------------------------------
screening_reasonable_forces(){
    local deal_lib="$HOME/.rebreath/deal_data/"
    python3 $deal_lib/elect_rely_force.py $1  $2  $3
}


# ---------------------------------------------------------------------------
# Function: screening_reasonable_energy
# 功能: 筛选NEP训练集中能量在指定范围内的构型
#       从xyz文件中提取能量在[min, max]范围内的构型。
# 场景: 数据集中可能存在DFT计算未收敛的构型（能量异常高或低），
#       使用此函数快速剔除这些离群点。
# Usage: screening_reasonable_energy <xyzfile> <min_energy> <max_energy>
# Example:
#   screening_reasonable_energy train.xyz -800 -780
#   # 筛选train.xyz中总能量在[-800, -780] eV范围内的构型
# ---------------------------------------------------------------------------
screening_reasonable_energy(){
    local deal_lib="$HOME/.rebreath/deal_data/"
    python3 $deal_lib/elect_rely_energy.py $1  $2  $3
}


# ---------------------------------------------------------------------------
# Function: screening_reasonable_virial
# 功能: 筛选NEP训练集中位力(virial)在指定范围内的构型
#       位力反映了体系的应力状态，筛选合理位力范围可排除异常应力构型。
# 场景: 当训练集包含极端压缩或拉伸的构型时，位力值会异常，
#       使用此函数剔除这些非物理构型。
# Usage: screening_reasonable_virial <xyzfile> <min_virial> <max_virial>
# Example:
#   screening_reasonable_virial train.xyz -100 100
#   # 筛选train.xyz中位力在[-100, 100] eV范围内的构型
# ---------------------------------------------------------------------------
screening_reasonable_virial(){
    local deal_lib="$HOME/.rebreath/deal_data/"
    python3 $deal_lib/elect_rely_virial.py $1  $2  $3
}


# =============================================================================
# SECTION 2: NEP Visualization / NEP 可视化
# =============================================================================

# ---------------------------------------------------------------------------
# Function: plot_nep
# 功能: 使用hplt_nep_results.py绘制NEP训练结果的标准图
#       画loss.out中能量/力/位力的训练和测试RMSE。
# 场景: 每次NEP训练完成后，在训练目录下执行此命令快速查看训练效果。
# Usage: plot_nep
#   # 在包含loss.out等文件的NEP训练目录下直接运行
#   cd my_nep_train/
#   plot_nep
# ---------------------------------------------------------------------------
plot_nep(){
    python3 ~/.rebreath/plot_library/hplt_nep_results.py
}


# ---------------------------------------------------------------------------
# Function: plot_ultimate_nep
# 功能: 使用优化配色方案绘制NEP训练结果图
#       比plot_nep使用更美观的配色和排版，适合论文/报告使用。
# 场景: 当你需要生成出版级别的NEP训练结果图时使用此函数。
# Usage: plot_ultimate_nep
#   # 在包含loss.out的NEP训练目录下运行
#   cd my_nep_train/
#   plot_ultimate_nep
# ---------------------------------------------------------------------------
plot_ultimate_nep(){
    cp $HOME/.rebreath/plot_library/plot_nep_results_ultimate.py .
    python3 plot_nep_results_ultimate.py
    rm -f plot_nep_results_ultimate.py
}


# ---------------------------------------------------------------------------
# Function: plot_E_F_Vir_distribution
# 功能: 绘制训练集中能量、力、位力的分布直方图
#       一张图三个子面板，直观展示数据集标签的统计分布。
# 场景: 在构建NEP训练集后，用于快速检查数据集的标签分布是否合理。
#       例如：力分布是否对称、能量分布是否连续、是否有离群值。
# Usage: plot_E_F_Vir_distribution [train.xyz]
# Example:
#   plot_E_F_Vir_distribution train.xyz
#   # 生成 E_F_vir.png 文件，包含三个直方图
# ---------------------------------------------------------------------------
plot_E_F_Vir_distribution(){
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

plt.subplot(1, 3, 2)
color_force = ['#2CAC15FF','#2C3460FF','#1A87ABFF']
F_avg = np.concatenate([config.force for config in configs])
plt.hist(F_avg, bins=100, alpha=0.5,color=['#2CAC15FF','#2C3460FF','#1A87ABFF'])
plt.xlabel('Interatomic Force (eV/\AA)',fontweight='bold')
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


# =============================================================================
# SECTION 3: HNEMD Thermal Conductivity / HNEMD 热导率
# =============================================================================

# ---------------------------------------------------------------------------
# Function: start_mul_hnemd
# 功能: 批量启动多次独立的HNEMD热导率计算
#       创建hnemd_0, hnemd_1...等文件夹，将输入文件复制到每个文件夹，
#       然后在每个文件夹中使用free_time_run在空闲GPU上启动gpumd。
# 场景: HNEMD模拟需要多次独立运行求统计平均（通常6-10次）。
#       此函数可一键创建所有独立运行目录并排队提交任务。
# Usage: start_mul_hnemd <nepfile> <times> <core_num>
# Example:
#   start_mul_hnemd nep.txt 6 1
#   # 创建6个独立HNEMD运行目录(nemd_0~5)，各使用1个CPU核
# ---------------------------------------------------------------------------
start_mul_hnemd() {
    local nepfile=${1:-nep.txt}
    local times=${2:-6}
    local core_num=${3:-1}
    local init_address=$PWD

    for i in $(seq 1 $((times-1))); do
        mkdir -p hnemd_${i}
        cp $nepfile model.xyz run.in hnemd_${i}/
        cd hnemd_${i}
        free_time_run 'nohup gpumd > nohup.out 2>&1 &'
        sleep 1
        cd $init_address
    done
}


# ---------------------------------------------------------------------------
# Function: start_mul_hnemd_shuguang
# 功能: 在曙光超算上批量启动HNEMD计算（使用质数种子避免随机数冲突）
# 场景: 在曙光DCU集群上运行HNEMD时使用，自动生成不同的随机数种子。
# Usage: start_mul_hnemd_shuguang <nepfile> <times> <core_num>
# Example:
#   start_mul_hnemd_shuguang nep.txt 8 4
#   # 创建8个目录，每个使用4个DCU核
# ---------------------------------------------------------------------------
start_mul_hnemd_shuguang() {
    local nepfile=${1:-nep.txt}
    local times=${2:-6}
    local core_num=${3:-1}
    local init_address=$PWD

    local prime_seeds=($(generate_large_primes $times 10001))

    for i in $(seq 1 $((times-1))); do
        mkdir -p hnemd_${i}
        cp $nepfile model.xyz run.in hnemd_${i}/
        cd hnemd_${i}
        sed -i "s/seed.*$/seed ${prime_seeds[$i]}/" run.in
        gpumd_dcu_394 -n $core_num
        cd $init_address
    done
}


# ---------------------------------------------------------------------------
# Function: plot_hnemd
# 功能: 自动识别目录路径中的热导方向(_x/_y/_z)并绘制对应方向的热导率图
# 场景: 完成HNEMD模拟后，在包含kappa.out的目录的上级目录执行，
#       自动遍历所有子目录找到kappa.out并画图。
# Usage: plot_hnemd
#   # 在包含多个hnemd_X目录的父目录下运行
#   plot_hnemd
# ---------------------------------------------------------------------------
plot_hnemd(){
    lib_address="$HOME/.rebreath/plot_library/"

    find $PWD -type f -regex '.*kappa.out' | while read -r file;do
      local hnemd_direct=0
      if [[ $file =~ _([xyz]) ]];then
         hnemd_direct=${BASH_REMATCH[-1]}
         echo -e '\n——> Found HNEMD direction: $hnemd_direct'
         else
         echo "HNEMD direction not found! Please follow naming convention with _{x|y|z}"
         return 1
      fi
    file_address=$(dirname $file)
    cd $file_address
    echo $PWD
    python3 "$lib_address/plot_hnemd_${hnemd_direct}.py" 2>&1
    cd - >/dev/null
    done
}


# ---------------------------------------------------------------------------
# Function: plot_hnemd_para
# 功能: 并行绘制所有HNEMD热导率图（自动检测CPU核数并分配任务）
# 场景: 当有大量HNEMD计算结果需要同时画图时使用，比plot_hnemd更快。
# Usage: plot_hnemd_para [root_dir]
# Example:
#   plot_hnemd_para           # 在当前目录搜索并并行画图
#   plot_hnemd_para /path/to/simulations
# ---------------------------------------------------------------------------
plot_hnemd_para() {
  local lib_address="$HOME/.rebreath/plot_library"
  local root="${1:-$PWD}"

  local total load1 used free jobs
  total=$(nproc --all)
  load1=$(awk '{print $1}' /proc/loadavg)
  used=${load1%.*}
  [[ -z "$used" ]] && used=0
  free=$(( total - used ))
  (( free < 1 )) && free=1
  jobs=$(( (free * 8) / 10 ))
  (( jobs < 1 )) && jobs=1
  (( jobs > total )) && jobs=$total

  echo "CPU total=$total, load1=$load1, est_free=$free, parallel_jobs=$jobs"

  find "$root" -type f -name "kappa.out" -print0 |
    while IFS= read -r -d '' file; do
      if [[ "$file" =~ _([xyz])(/|_) ]]; then
        hnemd_direct="${BASH_REMATCH[1]}"
        dir="$(dirname "$file")"
        printf '%s\t%s\n' "$dir" "$hnemd_direct"
      elif [[ "$file" =~ _([xyz]) ]]; then
        hnemd_direct="${BASH_REMATCH[1]}"
        dir="$(dirname "$file")"
        printf '%s\t%s\n' "$dir" "$hnemd_direct"
      else
        echo "[Skip] Direction marker not found: $file" >&2
      fi
    done |
    xargs -P "$jobs" -n 1 -I {} bash -c '
      set -euo pipefail
      line="{}"
      dir="${line%%$'\''\t'\''*}"
      direct="${line##*$'\''\t'\''}"
      lib_address="'"$lib_address"'"
      echo -e "\n——> Processing: $dir (direction=$direct)"
      cd "$dir"
      python3 "$lib_address/plot_hnemd_${direct}.py"
    '
}


# ---------------------------------------------------------------------------
# Function: plot_mul_hnemd
# 功能: 绘制多次HNEMD模拟的平均热导率图
#       自动检测average_hnemd/kappa目录，将所有独立运行的平均结果画在一张图上。
# 场景: 运行start_mul_hnemd启动多次独立的HNEMD模拟后，
#       使用此函数将多次结果汇总并画图，输出统计平均和误差。
# Usage: plot_mul_hnemd
#   # 在包含多个hnemd_1, hnemd_2...目录的父目录下运行
#   plot_mul_hnemd
# ---------------------------------------------------------------------------
plot_mul_hnemd(){
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
            echo -e '\n——> HNEMD direction: $hnemd_direct'
            else
            echo "HNEMD direction not found! Please follow naming convention."
            return 1
        fi
        cat $lib_address/plot_hnemd_mul_${hnemd_direct}.py |sed "/for i in range(/c\for i in range(${kappa_num}):" > plot_hnemd_mul_${hnemd_direct}.py
        python3 "plot_hnemd_mul_${hnemd_direct}.py" 2>&1
        cd - >/dev/null
    done
}


# ---------------------------------------------------------------------------
# Function: deal_hnemd_data
# 功能: 一键处理HNEMD数据：自动整理hnemd_[0-9]+文件夹，平均kappa.out，
#       并调用plot_mul_hnemd绘制最终热导率图。
# 场景: 完成start_mul_hnemd所有任务后，执行此函数一键完成数据处理和画图。
#       HNEMD模拟 → 等待所有任务完成 → deal_hnemd_data → 得到热导率图
# 前置条件: 当前目录下存在多个hnemd_[0-9]+文件夹，每个包含kappa.out
# Usage: deal_hnemd_data
# ---------------------------------------------------------------------------
deal_hnemd_data(){
    local lib_address="$HOME/.rebreath/deal_data/"
    source ${lib_address}/deal_hnemd_data.sh
    plot_mul_hnemd
}


# ---------------------------------------------------------------------------
# Function: check_hnemd_thermo
# 功能: 检查HNEMD的thermo.out输出文件，计算平均值和最终晶格参数
# 场景: HNEMD模拟运行中或运行后，快速检查thermo量是否正常收敛。
# Usage: check_hnemd_thermo
#   # 在包含thermo_*文件的目录下运行
#   check_hnemd_thermo
# ---------------------------------------------------------------------------
check_hnemd_thermo(){
    average_file_s thermo_*
    echo $PWD
    echo Lattice :
    tail -n 1 average_ther.out | awk '{print $10,$11,$12}'
}


# ---------------------------------------------------------------------------
# Function: get_hnemd_data
# 功能: 遍历x/y/z三个方向，依次进入各方向目录绘制HNEMD结果图
# 场景: 当你按hnemd_x/, hnemd_y/, hnemd_z/三个独立目录组织模拟时使用。
# Usage: get_hnemd_data
#   # 当前目录下需有 hnemd_x/ hnemd_y/ hnemd_z/ 三个子目录
#   get_hnemd_data
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Function: deal_strain_fluctuation_to_elastic
# 功能: 从应变波动法模拟的thermo.out中提取弹性模量
#       使用应变波动理论（stress-strain fluctuation formula）计算弹性常数。
# 场景: 运行GPUMD NVT模拟后，从能量和应力的涨落计算弹性模量。
#       适用于小体系（无法直接做拉伸模拟）的弹性性质估算。
# Usage: deal_strain_fluctuation_to_elastic <temperature_K>
# Example:
#   deal_strain_fluctuation_to_elastic 1200
#   # 从当前目录的thermo.out计算1200K下的弹性模量，输出到elastics.txt
# Dependencies: $HOME/.rebreath/deal_data/deal_strain_fluctuation.sh
# ---------------------------------------------------------------------------
deal_strain_fluctuation_to_elastic(){
    cp ~/.rebreath/deal_data/deal_strain_fluctuation.sh .
    bash deal_strain_fluctuation.sh $1 > elastics.txt
    rm -f deal_strain_fluctuation.sh
    cat elastics.txt
}


# =============================================================================
# SECTION 4: GPUMD Utilities / GPUMD 辅助工具
# =============================================================================

# ---------------------------------------------------------------------------
# Function: cell_expansion
# 功能: 使用GPUMD对xyz结构进行扩胞
#       通过生成临时run.in文件，调用gpumd的replicate功能完成扩胞。
# 场景: 构建大体系模型时需要将单胞扩展为超胞。
#       例如：将4原子的AlN原胞扩展为10x10x1的400原子超胞。
# Usage: cell_expansion <nx> <ny> <nz>
# Example:
#   cell_expansion 10 10 1
#   # 将当前目录的model.xyz在x和y方向各扩10倍，输出expanded.xyz
# Requires: model.xyz, nep.txt in current directory
# ---------------------------------------------------------------------------
cell_expansion(){
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
    echo "Cell expanded: $1 x $2 x $3"
}


# ---------------------------------------------------------------------------
# Function: verify_gpumd_result
# 功能: 创建verify_目录来重新验算当前GPUMD算例的结果
# 场景: 当怀疑GPUMD计算结果有误时，使用此函数创建独立目录重新计算。
# Usage: verify_gpumd_result [nepfile]
# Example:
#   verify_gpumd_result nep.txt
#   # 创建verify_<dirname>目录，复制输入文件用于重算
# ---------------------------------------------------------------------------
verify_gpumd_result(){
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


# ---------------------------------------------------------------------------
# Function: start_gpumd
# 功能: 根据配置文件(~/.rebreath/.config)自动选择GPUMD启动方式
#       支持普通GPU(gpumd)和DCU加速器(gpumdstart_dcu)两种模式。
# 场景: 在不同计算平台上使用统一的命令启动GPUMD，无需手动切换。
#       配置文件中设置gpumd_exe决定启动方式。
# Usage: start_gpumd [dcu_num]
# Example:
#   start_gpumd        # 普通GPU模式直接运行gpumd
#   start_gpumd 4      # DCU模式使用4个DCU核
# ---------------------------------------------------------------------------
start_gpumd(){
    if [ $gpumd_exe -eq "gpumd" ]; then
        gpumd
    elif [ $gpumd_exe -eq "gpumdstart_dcu" ]; then
        dcu_num=${1:-1}
        gpumdstart_dcu -n $dcu_num
    else
        echo "Error: Unknown gpumd startup method. Check ~/.rebreath/.config"
        exit 521
    fi
}


# ---------------------------------------------------------------------------
# Function: check_dataset_quality
# 功能: 对NEP训练集进行多维度质量诊断
#       包括：元素组成多样性、几何合理性(原子距离)、标签合理性(能量/力统计)、
#       描述符空间覆盖率(PCA)、孤立构型检测等。
# 场景: 每次构建/更新NEP训练集后，运行此诊断确保数据集质量。
#       - 发现原子重叠(距离过近)问题
#       - 检测能量/力标签中的离群值
#       - 评估训练集的化学空间覆盖度
# Usage: check_dataset_quality [--descriptor descriptor.out] [--outdir report_dir]
# Example:
#   check_dataset_quality train.xyz
#   check_dataset_quality train.xyz --descriptor descriptor.out --outdir diagnosis_report
# Dependencies: numpy, matplotlib, ase, scikit-learn
# ---------------------------------------------------------------------------
check_dataset_quality(){
    python3 $HOME/.rebreath/deal_data/dataset_quality_diagnosis.py $@
    echo "Dataset quality check complete"
}


# ---------------------------------------------------------------------------
# Function: compare_nepdata
# 功能: 比较训练集和测试集的基本信息（构型数量、元素组成等）
# 场景: 当你需要确认训练/测试集划分是否合理时使用。
# Usage: compare_nepdata [train.xyz] [test.xyz]
# Example:
#   compare_nepdata train.xyz test.xyz
#   # 输出: Training set has 2500 configs, Test set has 500 configs
# ---------------------------------------------------------------------------
compare_nepdata() {
    local train_file="${1:-train.xyz}"
    local test_file="${2:-test.xyz}"
    local num_train=""
    local num_test=""

    [[ -e "$train_file" ]] || echo "Warning: training set not found: $train_file" >&2
    [[ -e "$test_file"  ]] || echo "Warning: test set not found: $test_file" >&2

    if declare -F get_configs_num >/dev/null; then
        [[ -r "$train_file" ]] && num_train=$(get_configs_num "$train_file" 2>/dev/null)
        [[ -r "$test_file"  ]] && num_test=$(get_configs_num "$test_file" 2>/dev/null)
    else
        echo "Warning: get_configs_num not defined" >&2
    fi

    echo "Training set: ${num_train:-unknown} configs"
    echo "Test set:     ${num_test:-unknown} configs"
}


# =============================================================================
# SECTION 5: XYZ File Parsing Utilities / XYZ文件解析工具
# =============================================================================

# ---------------------------------------------------------------------------
# Function: get_energy
# 功能: 从GPUMD风格的xyz文件中提取所有构型的总能量
#       解析第二行注释中的Energy=字段。
# 场景: 需要批量查看/统计训练集中能量的分布或变化趋势时使用。
# Usage: get_energy [xyz_file]
# Example:
#   get_energy dump.xyz > energy_list.txt
#   get_energy train.xyz | head -5    # 查看前5个构型的能量
# ---------------------------------------------------------------------------
get_energy(){
    local xyz_file=${1:-'dump.xyz'}
    grep "attice" $xyz_file |  grep -oE '\b\w+="[^"]+"|\b\w+=[^ ]+\b' |grep -oE '[Ee]nergy.*' |grep -oE '[-]?[0-9]+[\.]?[0-9]+'
}


# ---------------------------------------------------------------------------
# Function: get_Lattice
# 功能: 从GPUMD风格的xyz文件中提取所有构型的晶格参数
#       解析第二行注释中的Lattice=字段，每行输出9个晶格常数。
# 场景: 需要检查轨迹文件的晶格参数变化，或为后续处理提取晶格信息。
# Usage: get_Lattice [xyz_file]
# Example:
#   get_Lattice dump.xyz | head -3   # 查看前3帧的晶格
#   get_Lattice train.xyz > lattice.log
# ---------------------------------------------------------------------------
get_Lattice(){
    local xyz_file=${1:-'dump.xyz'}
    grep "attice" $xyz_file |  grep -oE '\b\w+="[^"]+"|\b\w+=[^ ]+\b' | grep "attice"
}


# ---------------------------------------------------------------------------
# Function: get_virial
# 功能: 从GPUMD风格的xyz文件中提取所有构型的位力(virial)
#       解析第二行注释中的Virial=字段。
# 场景: 需要分析训练集中应力状态的分布时使用。
# Usage: get_virial [xyz_file]
# Example:
#   get_virial train.xyz | head -5
# ---------------------------------------------------------------------------
get_virial(){
    local xyz_file=${1:-'dump.xyz'}
    grep "attice" $xyz_file |  grep -oE '\b\w+="[^"]+"|\b\w+=[^ ]+\b' | grep "irial"
}


# ---------------------------------------------------------------------------
# Function: get_configs_num
# 功能: 统计xyz文件中构型的数量
#       通过统计包含"lattice"的行数来确定构型数。
# 场景: 快速确认训练集/测试集的大小。
# Usage: get_configs_num [xyz_file]
# Example:
#   get_configs_num train.xyz    # 输出构型总数
# ---------------------------------------------------------------------------
get_configs_num(){
    local xyz_file=${1:-'train.xyz'}
    grep "attice" $1 |wc -l
}


# ---------------------------------------------------------------------------
# Function: get_V
# 功能: 计算xyz轨迹文件或thermo.out文件中每帧的晶胞体积
#       对于xyz文件：从Lattice参数计算。对于thermo.out：从盒子向量计算。
# 场景: 分析MD模拟中体积随时间的变化（相变、热膨胀等）。
# Usage: get_V [file]
# Example:
#   get_V dump.xyz > volume_evolution.txt
#   get_V thermo.out > volume_evolution.txt
# ---------------------------------------------------------------------------
get_V(){
    local xyz_file=${1:-'dump.xyz'}
    if [[ $xyz_file =~ ".xyz" ]];then
    get_Lattice $xyz_file |sed 's/"/ /g' | awk '{printf "%.15f\n",$2*$6*$10}'
    elif [[ $xyz_file =~ ".out" ]];then
        python3 << EOF
import numpy as np
filename='$xyz_file'

def read_thermo(filename):
    """Read GPUMD thermo.out output file"""
    data = np.loadtxt(filename)
    thermo = dict()
    if data.shape[1] == 12:
        labels = ['T','K','U','Px','Py','Pz','Pyz','Pxz','Pxy','Lx','Ly','Lz']
        for i in range(12):
            thermo[labels[i]] = data[:,i]
    elif data.shape[1] == 18:
        labels = ['T','K','U','Px','Py','Pz','Pyz','Pxz','Pxy','ax','ay','az','bx','by','bz','cx','cy','cz']
        for i in range(18):
            thermo[labels[i]] = data[:,i]
    else:
        raise ValueError("Invalid number of columns in thermo.out file")
    return thermo

def get_volume(filename):
    """Compute cell volume for each frame"""
    thermo = read_thermo(filename)
    line_num = len(thermo['T'])
    V = np.zeros(line_num)
    if (len(thermo) == 18):
        for i in range(line_num):
            a1 = [thermo['ax'][i],thermo['ay'][i],thermo['az'][i]]
            a2 = [thermo['bx'][i],thermo['by'][i],thermo['bz'][i]]
            a3 = [thermo['cx'][i],thermo['cy'][i],thermo['cz'][i]]
            V[i] = np.dot(np.cross(a1,a2),a3)
        return V
    else:
        return thermo['Lx']*thermo['Ly']*thermo['Lz']

V = get_volume(filename)
for i in V:
    print(i)
EOF
    else
        echo "Error: Unknown file type. Supported: .xyz, .out"
fi
}


# ---------------------------------------------------------------------------
# Function: get_area_of_xy_and_volume
# 功能: 计算GPUMD风格xyz文件中每帧的xy面面积和斜胞体积
#       输出文件Volume_area_xy.txt，第一列为体积，第二列为xy面积。
# 场景: 为CaF2等材料使用ase切片后计算固定晶面所需。
#       注意：此为特定用途的临时函数。
# Usage: get_area_of_xy_and_volume <xyz_file>
# ---------------------------------------------------------------------------
get_area_of_xy_and_volume() {
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
    """Compute xy plane area using cross product and cell volume"""
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


# ---------------------------------------------------------------------------
# Function: plot_E_frame
# 功能: 绘制能量随帧数/MD步数的变化曲线
#       支持两种输入：dump.xyz（从xyz第二行提取能量）和thermo.out（GPUMD热力学输出）。
# 场景: MD模拟完成后，快速可视化能量随时间的变化，判断平衡状态。
#       - 对于dump.xyz：直接从注释行解析能量
#       - 对于thermo.out：计算K+U总能量，并做滑动平均
# Usage: plot_E_frame [file] [dump_interval]
# Example:
#   plot_E_frame dump.xyz          # 从xyz文件绘制能量曲线
#   plot_E_frame thermo.out 100    # 从thermo.out绘制，dump间隔为100步
#   # 输出: E-frame.png
# ---------------------------------------------------------------------------
plot_E_frame() {
    input_file=${1:-thermo.out}
    dump_interval=${2:-100}

    if [ ! -f "$input_file" ]; then
        echo "Error: file not found: $input_file"
        return 1
    fi

    case "$input_file" in
        *.xyz)
            echo "Detected xyz file: $input_file"
            energy_log=$(mktemp)
            frame_log=$(mktemp)
            get_energy "$input_file" > "$energy_log"
            seq "$(wc -l < "$energy_log")" > "$frame_log"
            paste "$frame_log" "$energy_log" > E-frame.txt
            replot E-frame.txt
            rm -f "$energy_log" "$frame_log"
            echo "Saved data: E-frame.txt"
            ;;

        *.out)
            echo "Detected thermo output file: $input_file"
            echo "dump_interval = $dump_interval"
            python3 - "$input_file" "$dump_interval" <<'PY'
import sys
import numpy as np
import matplotlib.pyplot as plt

thermo_file = sys.argv[1]
dump_interval = int(sys.argv[2])

data = np.loadtxt(thermo_file)
if data.ndim == 1:
    data = data.reshape(1, -1)
if data.shape[1] < 3:
    raise ValueError("thermo.out should contain at least 3 columns: T K U")

T = data[:, 0]
K = data[:, 1]
U = data[:, 2]
E = K + U

frame = np.arange(1, len(E) + 1)
step = frame * dump_interval

def moving_average(x, window):
    if window <= 1:
        return x
    kernel = np.ones(window) / window
    return np.convolve(x, kernel, mode="valid")

window = max(5, len(E) // 100)
if len(E) < window:
    window = max(1, len(E))

E_avg = moving_average(E, window)
U_avg = moving_average(U, window)
step_avg = step[window - 1:]

plt.figure(figsize=(8, 5))
plt.plot(step, K, alpha=0.4, label="Kinetic energy")
plt.plot(step, U, alpha=0.4, label="Potential energy")
plt.plot(step, E, alpha=0.4, label="Total energy K+U")
plt.plot(step_avg, U_avg, linewidth=2, label=f"Potential energy moving avg, window={window}")
plt.plot(step_avg, E_avg, linewidth=2, label=f"Total energy moving avg, window={window}")
plt.xlabel("MD step")
plt.ylabel("Energy (eV)")
plt.title("Energy evolution from thermo.out")
plt.legend()
plt.tight_layout()
plt.savefig("E-frame.png", dpi=300)
print("Saved figure: E-frame.png")
print(f"Number of thermo frames: {len(E)}")
print(f"Moving average window: {window} frames")
PY
            ;;

        *)
            echo "Error: unsupported file suffix: $input_file"
            echo "Supported: *.xyz (dump.xyz), *.out (thermo.out)"
            return 1
            ;;
    esac
}


# =============================================================================
# SECTION 6: Elastic Properties & Phonon / 弹性性质与声子
# =============================================================================

# ---------------------------------------------------------------------------
# Function: compute_elastic_moduli
# 功能: 使用calorine库计算弹性模量(Cij矩阵)
#       通过NEP模型对给定结构进行弹性常数计算。
# 场景: 训练好NEP势函数后，需要快速评估材料的弹性性质。
#       可计算弹性常数张量、体模量、剪切模量、杨氏模量等。
# Usage: compute_elastic_moduli [nepfile]
# Example:
#   compute_elastic_moduli nep.txt
#   # 输出弹性常数到elastic_calorine.txt
# Dependencies: calorine (pip install calorine)
# ---------------------------------------------------------------------------
compute_elastic_moduli(){
    local compute_lib=$HOME/.rebreath/compute_lib
    local nep_file=${1:-nep.txt}
    sed "s/nepfile/$nep_file/g" $compute_lib/calorine_compute_elastic.py > calorine_compute_elastic.py
    python3 calorine_compute_elastic.py > elastic_calorine.txt
    rm -f calorine_compute_elastic.py
    cat elastic_calorine.txt
}


# ---------------------------------------------------------------------------
# Function: compute_phonon_spectrum
# 功能: 使用GPUMD计算声子谱（通过phonopy接口）
#       调用gpumd_compute_phonon_spectrum.py进行声子色散计算。
# 场景: 训练好NEP势函数后，需要验证势函数能否正确描述晶格动力学。
#       声子谱是NEP势函数质量的重要验证手段。
# Usage: compute_phonon_spectrum
#   # 在当前目录下运行，需要model.xyz和nep.txt
#   compute_phonon_spectrum
# Dependencies: phonopy, gpumd_compute_phonon_spectrum.py
# ---------------------------------------------------------------------------
compute_phonon_spectrum() {
  local script="$HOME/.rebreath/compute_lib/gpumd_compute_phonon_spectrum.py"
  python3 - < "$script"
}


# =============================================================================
# SECTION 7: VASP-NEP Bridge / VASP与NEP桥接
# =============================================================================

# ---------------------------------------------------------------------------
# Function: add_kspacing_to_incar
# 功能: 为VASP的INCAR文件添加或修改KSPACING参数
#       用于使用单一k点间距替代KPOINTS文件（适用于大体系）。
# 场景: 在处理大量VASP单点能计算时，使用KSPACING比手写KPOINTS更方便。
#       此函数会自动备份原有的KPOINTS文件。
# Usage: add_kspacing_to_incar <kspacing_value>
# Example:
#   add_kspacing_to_incar 0.2
#   # 在INCAR中设置KSPACING = 0.2，并备份原KPOINTS
# ---------------------------------------------------------------------------
add_kspacing_to_incar(){
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


# ---------------------------------------------------------------------------
# Function: check_vasp_complete
# 功能: 检查当前目录的VASP计算是否完成，如未完成则提交运行
#       通过检查OUTCAR中是否包含完成标志来判断。
# 场景: 批量处理VASP计算时，自动跳过已完成的任务，只运行未完成的。
#       VASP单点能工作流中确保所有计算都完成的关键环节。
# Usage: check_vasp_complete [vasp_exe] [core_num]
# Example:
#   check_vasp_complete vasp_std 4
#   # 如果OUTCAR不存在或不完整，使用4核运行vasp_std
# ---------------------------------------------------------------------------
check_vasp_complete() {
    vasp_exe=${1:-'vasp'}
    core_num=${2:-'1'}
    if [ ! -f  "OUTCAR" ];then
        echo "OUTCAR not found, starting calculation..."
        mpirun -np $core_num $vasp_exe
    elif [ -z "$(grep "General timing and accounting informations for this job" OUTCAR)" ];then
        echo "OUTCAR incomplete in $PWD, restarting..."
        mpirun -np $core_num $vasp_exe
    fi
}


# ---------------------------------------------------------------------------
# Function: run_all_vasp_job
# 功能: 遍历所有train-*目录，依次确保每个目录的VASP任务完成
#       自动添加KSPACING、检查计算完整性、记录日志。
# 场景: VASP→NEP训练集构建流程中的核心步骤。
#       当你在根目录下创建了大量train-XXX文件夹后（每个包含一个POSCAR扰动构型），
#       使用此函数自动完成所有单点能计算。
# Usage: run_all_vasp_job
#   # 当前目录下需要有多个train-*子目录，每个包含INCAR, POTCAR, POSCAR
#   run_all_vasp_job
# Output: run_train-file.log (记录每个文件夹的计算状态)
# ---------------------------------------------------------------------------
run_all_vasp_job(){
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
        echo "Folder $i: OUTCAR not found, starting calculation"
        check_vasp_complete
        echo "Exit code: $? for $i" >> $init_address/run_train-file.log
    fi
    complete=$(grep "General timing and accounting informations for this job" OUTCAR)
    if [ -z "$complete" ];then
    echo "Folder $i: OUTCAR incomplete, restarting..."
    check_vasp_complete
    echo "Exit code: $? for $i" >> $init_address/run_train-file.log
    fi
    cd $init_address
    done
    time_log=$(date)
    echo "Completion time: $time_log"
    endtime=$(date +%s)
    alltime=$((endtime-startime))
    min_time=$((alltime/60))
    echo "run_train-file.sh completed! Total time: $min_time min"
}


# ---------------------------------------------------------------------------
# Function: load_single_point_energy_dir
# 功能: 在当前目录下查找所有POSCAR文件，为每个创建单点能计算目录
#       自动构建single_energy/train-XXX目录结构，复制INCAR/POTCAR/POSCAR。
# 场景: 当你有一组POSCAR构型文件（如从AIMD轨迹提取的），需要批量为它们
#       创建VASP单点能计算目录。此函数配合run_all_vasp_job使用。
# Usage: load_single_point_energy_dir
#   # 当前目录下需要包含多个POSCAR*文件和INCAR, POTCAR
#   load_single_point_energy_dir
#   # 然后可以运行 run_all_vasp_job 来批量计算
# ---------------------------------------------------------------------------
load_single_point_energy_dir(){
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


# ---------------------------------------------------------------------------
# Function: deal_outcar_to_train
# 功能: 将当前目录下所有OUTCAR文件（VASP单点能输出）转换为NEP训练集格式
#       遍历所有OUTCAR，提取晶格、能量、位力和原子受力，写入train.xyz。
# 场景: VASP→NEP训练集构建流程的最后一步。
#       完成所有单点能计算后，执行此函数一键生成NEP可用的train.xyz文件。
# Usage: deal_outcar_to_train [target_directory]
# Example:
#   deal_outcar_to_train              # 处理当前目录所有OUTCAR
#   deal_outcar_to_train /path/to/calcs  # 处理指定目录
# Output: train.xyz (NEP格式训练集文件)
# ---------------------------------------------------------------------------
deal_outcar_to_train(){
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
    echo "Found $N_case OUTCAR files"
    if [ $N_case -eq 0 ]
    then
    echo "No OUTCAR files found!"
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
            continue
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
    echo "Working directory: $init_address"
    echo "New configs added: $N_cout"
    echo "Completed: $N_cout / $N_case"
    endtime=$(date +%s)
    alltime=$((endtime-startime))
    min_time=$((alltime/60))
    echo "Total time: $min_time min"
}


# =============================================================================
# SECTION 8: Stress-Strain Analysis / 应力-应变分析
# =============================================================================

# ---------------------------------------------------------------------------
# Function: plot_stress_strain_curve
# 功能: 自动检测变形轴并绘制应力-应变曲线
#       识别deform_{xyz}+目录结构，提取应力数据并画图。
# 场景: 完成GPUMD的deform模拟（拉伸/压缩）后，使用此函数快速得到
#       应力-应变曲线，用于计算杨氏模量、屈服强度等力学性质。
# Usage: plot_stress_strain_curve
#   # 在当前目录下运行，自动检测deform_*目录
#   plot_stress_strain_curve
#   # 输出应力-应变数据和曲线图
# ---------------------------------------------------------------------------
plot_stress_strain_curve(){
    source $HOME/.rebreath/plot_library/stress_strain_curve.sh
}


# ---------------------------------------------------------------------------
# Function: plot_mul_stress_strain_curve
# 功能: 自动检测多个方向的deform_{xyz}+文件并绘制多方向应力-应变曲线
# 场景: 当你对同一材料在不同方向(x/y/z)进行了拉伸模拟后，
#       使用此函数一次性将所有方向的应力-应变曲线画在同一张图上对比。
# Usage: plot_mul_stress_strain_curve
#   # 当前目录下需有deform_x/, deform_y/, deform_z/等子目录
#   plot_mul_stress_strain_curve
# ---------------------------------------------------------------------------
plot_mul_stress_strain_curve(){
    source $HOME/.rebreath/plot_library/auto_plot_xyz_strain_stress_curve.sh
}


# =============================================================================
# SECTION 9: Cluster-Specific Startup / 超算集群启动器
# =============================================================================

# ---------------------------------------------------------------------------
# Function: gpumdstart_rebreath / gpumdstart_zwj
# 功能: 分别在曙光超算不同的DCU分区上提交GPUMD任务
#       自动生成SLURM提交脚本，根据配置文件选择GPUMD版本。
# 场景: 在曙光超算上运行GPUMD时，使用对应的启动器简化SLURM脚本编写。
# Usage:
#   gpumdstart_rebreath [-n dcu_num]   # rebreath分区
#   gpumdstart_zwj [-n dcu_num]        # zwj分区
# Example:
#   gpumdstart_rebreath -n 4
#   # 使用4个DCU核在rebreath分区提交GPUMD任务
# ---------------------------------------------------------------------------
gpumdstart_rebreath(){
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
    #SBATCH --time 240:00:00
    #SBATCH --comment=GPUMD
    #SBATCH -o $PWD/std.out.%j
    #SBATCH -e $PWD/std.err.%j

    # MARK_CMD
    source /work/home/rebreath/sbatch_need/gpumd_env.sh
    /work/share/acmtrwrxv5/GPUMD/src/gpumd" | sbatch -J "$job_name"
}

gpumdstart_zwj(){
    dcu_num=1
    for arg in "$@"; do
        if [[ $arg = "-n" ]]; then
        dcu_num=${2:-1}
        break
        fi
    done

    job_name=$(basename "$PWD")
    echo "
#!/bin/bash
#SBATCH -N  1
#SBATCH --ntasks-per-node=$dcu_num
#SBATCH --gres=dcu:$dcu_num
#SBATCH --comment=GPUMD
#SBATCH -p wzhdtest

# MARK_CMD
module purge
module use /public/software/modules /opt/hpc/software/modules
source /work/home/acyrhkta3v/apprepo/gpumd/3.9.5-dtk24.04/scripts/env.sh
/work/home/acyrhkta3v/apprepo/gpumd/3.9.5-dtk24.04/app/bin/gpumd" | sbatch -J "$job_name"

}
