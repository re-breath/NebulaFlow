# ======================================================================
# File:         lammpsEnvFunction.sh
# Project:      NebulaFlow
# Description:  LAMMPS模拟相关函数库 — data/setting合并、移位、排序、任务提交
#               LAMMPS data/settings merging, shifting, sorting, job submission.
# Author:       rebreath
# Dependencies: lmp, python3, ovito
# ======================================================================


# ---------------------------------------------------------------------------
# Function: merge_lmp
# 功能: 合并LAMMPS的data文件和setting文件为一个完整的data文件
#       使用LAMMPS的write_data功能，自动处理力场参数合并。
# 场景: LigParGen等工具分开输出data和setting文件时，合并为单一文件方便管理。
# Usage: merge_lmp <datafile> <settings>
# Example:
#   merge_lmp system.data frag_GMX.set
#   # 输出: mergedata.lmp
# ---------------------------------------------------------------------------
merge_lmp() {
    local datafile="$1"; local settings="$2"
    dpmd || true
    cat > mergelmp.in << EOF
    units real
    atom_style      full
    bond_style      harmonic
    angle_style     harmonic
    dihedral_style  opls
    improper_style  harmonic
    pair_style      lj/cut/coul/long   10
    kspace_style    pppm 0.0001
    special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes dihedral yes
    read_data ${datafile}
    include ${settings}
    write_data mergedata.lmp
EOF
    mpirun -np 1 lmp -i mergelmp.in > merge.log
}


# ---------------------------------------------------------------------------
# Function: shift_lmp_data
# 功能: 使用OVITO平移LAMMPS data文件中的原子并wrap回周期性盒子
# 场景: 分子在盒子边界处被切割时，需要整体平移使分子完整。
# Usage: shift_lmp_data <in.data> <out.data> [dx] [dy] [dz]
# Example:
#   shift_lmp_data system.data shifted.data 10 0 0
#   # 在x方向平移10埃
# ---------------------------------------------------------------------------
shift_lmp_data() {
    local in_file="$1"; local out_file="$2"; local shiftx="${3:-10}"; local shifty="${4:-10}"; local shiftz="${5:-10}"
    python3 << EOF
import sys,numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import AffineTransformationModifier, WrapPeriodicImagesModifier
pipeline=import_file('$in_file')
pipeline.modifiers.append(AffineTransformationModifier(translation=np.array([$shiftx,$shifty,$shiftz])))
pipeline.modifiers.append(WrapPeriodicImagesModifier())
export_file(pipeline,'$out_file',"lammps/data",atom_style="full")
print(f"Written: $out_file")
EOF
}


# ---------------------------------------------------------------------------
# Function: lammpsrun
# 功能: 在曙光超算CPU分区提交LAMMPS任务
# Usage: lammpsrun <ncore> <in_file>
# Example:
#   lammpsrun 64 lmp.in
#   # 使用64核提交LAMMPS任务
# ---------------------------------------------------------------------------
lammpsrun(){
    local ncore="$1"; local in_file="$2"
    cat > lmp_cpu.slurm << EOF
#!/bin/bash
#SBATCH -J lammps
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=$ncore

module purge
source /work/home/rebreath/apprepo/lammps/12Dec18-intelmpi2017/scripts/env.sh
mpirun -np $ncore lmp -i  $in_file
EOF
    sbatch lmp_cpu.slurm
    echo "LAMMPS running $in_file with $ncore cores"
}


# ---------------------------------------------------------------------------
# Function: sortlmpdata
# 功能: 对LAMMPS data文件中的原子ID进行重新排序
# Usage: sortlmpdata <input.data>
# Example:
#   sortlmpdata system.data
#   # 输出: system_sort.data
# ---------------------------------------------------------------------------
sortlmpdata(){
    local in_file="$1"; local out_file="${in_file%.*}_sort.data"
    cp ~/.rebreath/deal_data/sortlmpdata.py .
    python3 sortlmpdata.py $in_file $out_file
}
