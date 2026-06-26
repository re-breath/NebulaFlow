# ======================================================================
# File:         cp2kEnvFunction.sh
# Project:      NebulaFlow
# Description:  CP2K计算相关函数库 — restart转换、能量提取、任务提交
#               CP2K restart file conversion, energy extraction, job submission.
# Author:       rebreath
# Dependencies: Multiwfn, cp2k
# ======================================================================


# ---------------------------------------------------------------------------
# Function: cp2krestart2cif / cp2krestart2xyz
# 功能: 使用Multiwfn将CP2K restart文件转换为cif或xyz格式
# 场景: CP2K几何优化完成后，restart文件包含最终结构，需要提取为
#       标准格式用于后续分析或VASP/GPUMD计算。
# Usage:
#   cp2krestart2cif [cif_filename]     # 默认输出opt.cif
#   cp2krestart2xyz [xyz_filename]     # 默认输出opt.xyz
# Example:
#   cp2krestart2cif result.cif
#   cp2krestart2xyz optimized.xyz
# Dependencies: Multiwfn
# ---------------------------------------------------------------------------
cp2krestart2cif() {
    local cifilename=${1:-"opt.cif"}
    Multiwfn *.restart << EOF > /dev/null
    100
    2
    33
    $cifilename
    0
    q
EOF
}

cp2krestart2xyz() {
    local xyzfilename=${1:-"opt.xyz"}
    Multiwfn *.restart << EOF  > /dev/null
    100
    2
    2
    $xyzfilename
    0
    q
EOF
}


# ---------------------------------------------------------------------------
# Function: cp2kstart
# 功能: 在曙光超算上生成SLURM脚本并提交CP2K任务
# 场景: 在集群上运行CP2K几何优化或单点能计算时使用。
# Usage: cp2kstart [cp2k.inp] [cpu_num]
# Example:
#   cp2kstart cp2k.inp 32
#   # 使用32核提交CP2K任务
# ---------------------------------------------------------------------------
cp2kstart() {
    local inpfile=${1:-"cp2k.inp"}
    local cpu_num=${2:-32}
    local job_name=$(basename "$PWD")
    cat > temp.slurm <<EOF
#!/bin/bash
#SBATCH -J $job_name
#SBATCH -N 1
#SBATCH --ntasks-per-node=$cpu_num
#SBATCH -p xahcnormal

module purge
module load cp2k/2023.1-intelmpi-2018

srun --mpi=pmi2 cp2k.popt -i $inpfile -o cp2k.log
EOF
     sbatch  temp.slurm
}


# ---------------------------------------------------------------------------
# Function: get_cp2k_energy / update_cp2k_inp_cell_from_xyz
# 功能: 从CP2K log文件提取能量 / 从xyz更新CP2K输入文件的晶胞参数
# Usage:
#   get_cp2k_energy [cp2k.log]
#   update_cp2k_inp_cell_from_xyz model.xyz cp2k.inp
# ---------------------------------------------------------------------------
get_cp2k_energy() {
    local logfile=${1:-"cp2k.log"}
    grep "ENERGY|" $logfile |tail -n 5
}

update_cp2k_inp_cell_from_xyz() {
    local xyz_file="$1"; local cp2k_inp="$2"
    Lattice=$(get_Lattice $1 | grep -oP '(?<=Lattice=").*(?=")')
    cell_A=$(echo $Lattice | awk '{print $1,$2,$3}')
    cell_B=$(echo $Lattice | awk '{print $4,$5,$6}')
    cell_C=$(echo $Lattice | awk '{print $7,$8,$9}')
    sed -E "/^\s*A\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      A   $cell_A/" $cp2k_inp |sed -E "/^\s*B\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      B   $cell_B/" |sed -E "/^\s*C\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      C   $cell_C/" > ${cp2k_inp%.*}_up.inp
    sed -i "/@SET XYZFILE/s/.*/@SET XYZFILE    $1/" ${cp2k_inp%.*}_up.inp
    sed -i "/@SET PROJECT_NAME/s/.*/@SET PROJECT_NAME    ${xyz_file%.*}/" ${cp2k_inp%.*}_up.inp
}
