# 该函数库存放lammps相关的函数

# 合并LAMMPS的data和setting文件
# datafile: 输入的data文件
# settings: 输入的setting文件
merge_lmp() {
    local datafile="$1"
    local settings="$2"

    dpmd || true

    cat > mergelmp.in << EOF
    units real
    atom_style      full
    bond_style      harmonic
    angle_style     harmonic
    dihedral_style  opls
    improper_style  harmonic #cvff #有时候LigParGen输出的文件中improper_style的格式为cvff
    pair_style      lj/cut/coul/long   10
    
    kspace_style    pppm 0.0001
    special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes dihedral yes
       
    read_data ${datafile}        
    include ${settings}

    # 输出数据文件
    write_data mergedata.lmp
EOF
    mpirun -np 1  lmp -i mergelmp.in > merge.log

}

# 该函数使用ovitos来处理LAMMPS的data文件
# in_file: 输入的data文件
# out_file: 输出的data文件
# shift: 平移的距离，默认值为10.0 10.0 10.0
shift_lmp_data() {
    local in_file="$1"
    local out_file="$2"
    local shiftx="${3:-10}"
    local shifty="${4:-10}"
    local shiftz="${5:-10}"

python3 << EOF
import sys
import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import AffineTransformationModifier, WrapPeriodicImagesModifier

# ---------------- 参数 ----------------
in_file  = '$in_file'
out_file = '$out_file'
shift    = np.array([$shiftx, $shifty, $shiftz])  # Δx, Δy, Δz (Å)

# ---------------- 读入 ----------------
pipeline = import_file(in_file)

# 1) 整体平移
pipeline.modifiers.append(
    AffineTransformationModifier(translation=shift)
)

# 2) 周期性 wrap
pipeline.modifiers.append(
    WrapPeriodicImagesModifier()
)

# ---------------- 导出 ----------------
export_file(
    pipeline,
    out_file,
    "lammps/data",
    atom_style="full"   # 与 LAMMPS atom_style 保持一致，可改
)
print(f"Shifted & wrapped file written to {out_file}")
EOF
}

# 该函数使用lammps来运行模拟
# 使用方式： lammpsrun 64 lmp.inp 表示将会使用64个核心来运行lmp.inp文件
lammpsrun(){
    local ncore="$1"
    local in_file="$2"

    cat > lmp_cpu.slurm << EOF
#!/bin/bash
#SBATCH -J lammps #作业名称
#SBATCH -p xahcnormal #队列名称
#SBATCH -N 1 #节点数量
#SBATCH --ntasks-per-node=$ncore 
 
module purge
source /work/home/rebreath/apprepo/lammps/12Dec18-intelmpi2017/scripts/env.sh
 
mpirun -np $ncore lmp -i  $in_file
EOF

    sbatch lmp_cpu.slurm
    echo "Lammps runing $in_file with $ncore cores"
}