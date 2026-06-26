# ======================================================================
# File:         vaspEnvFunction.sh
# Project:      NebulaFlow
# Description:  VASP计算相关函数库 — SLURM提交、POTCAR管理、声子谱、弹性模量
#               VASP job submission, POTCAR/k-points management, and phonon.
# Author:       rebreath
# Dependencies: vasp, vaspkit, python3, ase, yaml, matplotlib
# ======================================================================


# =============================================================================
# SECTION 1: VASP Job Submission / VASP 任务提交
# =============================================================================

# ---------------------------------------------------------------------------
# Function: vasprun
# 功能: 生成SLURM脚本并提交VASP CPU计算任务
#       自动检测INCAR/POSCAR/POTCAR文件，生成sbatch脚本提交到CPU队列。
# 场景: 在集群上运行VASP结构优化、单点能等标准计算时使用。
#       适合CPU集群(xahcnormal分区)，默认64核。
# Usage: vasprun [core_num]
# Example:
#   vasprun        # 使用64核提交VASP任务
#   vasprun 32     # 使用32核提交VASP任务
# Requires: INCAR, POSCAR, POTCAR in current directory
# ---------------------------------------------------------------------------
vasprun(){
    local taskname=$(basename $PWD)
    local corenum=${1:-64}

    for file in INCAR POSCAR POTCAR; do
        if [ ! -f "$file" ]; then
            echo "Error: Required VASP input file '$file' does not exist."
            return 1
        fi
    done

    local slurm_script="vasp.slurm"
    cat > "$slurm_script" <<EOF
#!/bin/bash
#SBATCH -J $taskname
#SBATCH -N 1
#SBATCH --ntasks-per-node=$corenum
#SBATCH -p xahcnormal
module purge
source $HOME/apprepo/vasp/6.4.3-ioptcell_intelmpi2017_hdf5_libxc/scripts/env.sh
export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=AVX2
export I_MPI_PIN_DOMAIN=numa
echo "Starting VASP job at \$(date)"
srun --mpi=pmi2 vasp_std
echo "VASP job completed at \$(date)"
EOF
    echo "Submitting VASP job $taskname with $corenum cores..."
    sbatch "$slurm_script"
    local sbatch_status=$?
    if [ $sbatch_status -eq 0 ]; then
        echo "VASP job $taskname submitted successfully."
    else
        echo "Error: Failed to submit VASP job '$taskname'."
        return 1
    fi
}


# ---------------------------------------------------------------------------
# Function: vasprun_dcu
# 功能: 生成SLURM脚本并提交VASP DCU加速计算任务
#       与vasprun类似，但使用DCU加速卡(xahdnormal分区)。每个DCU配8个CPU核。
# 场景: 在曙光DCU集群上运行VASP计算时使用。
# Usage: vasprun_dcu [dcu_num]
# Example:
#   vasprun_dcu       # 使用1个DCU提交VASP任务
#   vasprun_dcu 4     # 使用4个DCU提交VASP任务
# Requires: INCAR, POSCAR, POTCAR in current directory
# ---------------------------------------------------------------------------
vasprun_dcu() {
    local taskname=$(basename "$PWD")
    local corenum_dcu=${1:-1}
    local core_cpu=$((corenum_dcu * 8))

    for file in INCAR POSCAR POTCAR; do
        if [ ! -f "$file" ]; then
            echo "Error: Required VASP input file '$file' does not exist."
            return 1
        fi
    done
    local slurm_script="vasp.slurm"
    cat > "$slurm_script" <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -n $core_cpu
#SBATCH -J $taskname
#SBATCH --gres=dcu:$corenum_dcu
#SBATCH -p xahdnormal

module purge
source \$HOME/apprepo/vasp/6.4.2-dtk23.10_cpvasp/scripts/env.sh
VASP_EXE=vasp_std
config=config.\${SLURM_JOB_ID}
echo -e "-genv OMP_NUM_THREADS 6 \c" > \$config
for i in \$(scontrol show hostnames \$SLURM_NODELIST)
do
  for ((j=0; j<$corenum_dcu; j++))
  do
    echo "-host \$i -env HIP_VISIBLE_DEVICES \$j -n 1 numactl --cpunodebind=\$j --membind=\$j \$VASP_EXE" >> \$config
  done
done

ulimit -s unlimited
export NCCL_IB_HCA="mlx5_0"
export HSA_FORCE_FINE_GRAIN_PCIE=1
mpirun -configfile \$config
EOF
    echo "Submitting VASP job $taskname with $corenum_dcu DCU(s)..."
    sbatch "$slurm_script"
    local sbatch_status=$?
    if [ $sbatch_status -eq 0 ]; then
        echo "VASP job $taskname submitted successfully."
    else
        echo "Error: Failed to submit VASP job '$taskname'."
        return 1
    fi
}


# =============================================================================
# SECTION 2: VASP Input File Management / VASP 输入文件管理
# =============================================================================

# ---------------------------------------------------------------------------
# Function: add_potcar
# 功能: 使用vaspkit自动生成POTCAR文件（根据POSCAR中的元素）
#       调用vaspkit的103功能，根据POSCAR自动匹配推荐赝势。
# 场景: 每次创建新的VASP计算目录时，不需要手动拼接POTCAR，
#       直接在包含POSCAR的目录下运行此命令即可自动生成POTCAR。
# Usage: add_potcar
# Requires: POSCAR, vaspkit
# ---------------------------------------------------------------------------
add_potcar() {
    echo "103" | vaspkit > /dev/null
}


# ---------------------------------------------------------------------------
# Function: add_kpoints
# 功能: 使用vaspkit生成KPOINTS文件（Gamma-centered MP网格）
#       指定k点间距(单位: 2π/Å)，自动根据倒格矢计算k点网格。
# 场景: 当你需要根据统一的k点密度标准为不同体系生成KPOINTS时使用。
#       对周期性体系，0.03是常用精度标准。
# Usage: add_kpoints [kpoint_spacing]
# Example:
#   add_kpoints 0.03         # 标准精度的k点
#   add_kpoints 0.02         # 更高精度的k点
#   add_kpoints 0.05         # 快速计算用
# Requires: POSCAR, vaspkit
# ---------------------------------------------------------------------------
add_kpoints() {
    local level=${1:-0.03}
    echo -e "102 \n 1 \n $level" | vaspkit > /dev/null
}


# ---------------------------------------------------------------------------
# Function: add_kspacing_to_incar
# 功能: 为INCAR添加或修改KSPACING参数，使用单一k点间距代替KPOINTS文件
#       适用于大体系（几百原子以上），KSPACING比KPOINTS网格更方便。
# 场景: 处理大体系的VASP单点能计算时，使用KSPACING简化输入。
#       此函数会自动备份原有KPOINTS文件。
# Usage: add_kspacing_to_incar <kspacing_value>
# Example:
#   add_kspacing_to_incar 0.2
#   # 在INCAR中设置KSPACING = 0.2，备份原KPOINTS为KPOINTS_backup
# ---------------------------------------------------------------------------
add_kspacing_to_incar(){
    kspacing=${1:-0.2}
    mv KPOINTS KPOINTS_backup || true
    if [ -z "$(grep "KSPACING" INCAR)" ];then
        sed -i "\$a\KSPACING = $kspacing "  INCAR
    else
        sed "/KSPACING/c\KSPACING = $kspacing/" INCAR
    fi
}


# ---------------------------------------------------------------------------
# Function: clean_vasp_out
# 功能: 删除当前目录下所有VASP输出文件，保留输入文件
#       清理WAVECAR, CHGCAR, DOSCAR等大文件，释放磁盘空间。
# 场景: VASP计算完成后，输出文件占用大量磁盘空间（尤其是WAVECAR），
#       确认不需要重启计算后使用此命令清理。
# Usage: clean_vasp_out
#   cd my_vasp_calc/
#   clean_vasp_out
# Note: 此操作不可逆，确保你不再需要这些输出文件
# ---------------------------------------------------------------------------
clean_vasp_out(){
    rm nohup.out CONTCAR OSZICAR WAVECAR PCDAT REPORT slurm-* vasp* config.* DOSCAR C* wait.log XDATCAR IBZKPT EIGENVAL REPORT OUTCAR PROCAR > /dev/null
}


# =============================================================================
# SECTION 3: Structure Manipulation / 结构操作
# =============================================================================

# ---------------------------------------------------------------------------
# Function: fix_poscar_zFrc
# 功能: 使用vaspkit和ASE按分数固定POSCAR中z轴下半部分的原子
#       计算z轴坐标范围，固定指定分数以下的原子。
# 场景: 做表面/界面计算时，需要固定底部几层原子模拟体相环境。
#       例如：固定底部70%的原子，只让表面30%的原子自由弛豫。
# Usage: fix_poscar_zFrc [posfile] [frac]
# Example:
#   fix_poscar_zFrc POSCAR 0.7
#   # 固定z轴下半部分70%的原子，输出POSCAR_FIX.vasp
# Dependencies: ase, vaspkit
# ---------------------------------------------------------------------------
fix_poscar_zFrc() {
    local posfile=${1:-POSCAR}
    local frac=${2:-0.7}
    read z_min z_max < <(python3 << EOF
def find_min_max_pos_z(frac):
    import ase.io as ai
    atoms = ai.read('$posfile')
    pos = atoms.get_positions()
    z_pos = pos[:, 2]
    z_min = min(z_pos)
    z_max = max(z_pos)
    z_interval = z_max - z_min
    z_max_new = z_min + z_interval * frac
    print(f"{z_min} {z_max_new}")
find_min_max_pos_z($frac)
EOF
)
    echo "z_min: $z_min"
    echo "z_max: $z_max"
    if [ ! -f "POSCAR" ]; then
        cp "$posfile" POSCAR > /dev/null
    fi
    echo -e "402 \n 1 \n 3 \n $z_min $z_max\n 2\n all" |vaspkit > /dev/null
    unfix_num=$(grep "T" POSCAR_FIX.vasp |wc -l)
    unfix_num=$((unfix_num - 1))
    echo "Unfixed atom numbers: $unfix_num"
}


# ---------------------------------------------------------------------------
# Function: fix_poscar_zCar
# 功能: 使用vaspkit按绝对坐标固定POSCAR中指定z范围内的原子
#       与fix_poscar_zFrc不同，此函数使用绝对z值而非分数。
# 场景: 当你明确知道要固定的z范围时（如固定z=0到z=18埃的原子）。
# Usage: fix_poscar_zCar [posfile] [z_max] [z_min]
# Example:
#   fix_poscar_zCar POSCAR 18 0
#   # 固定z坐标在[0, 18]范围内的原子
# ---------------------------------------------------------------------------
fix_poscar_zCar() {
    local posfile=${1:-POSCAR}
    local z_max=${2:-18}
    local z_min=${3:-0}

echo "z_min: $z_min"
echo "z_max: $z_max"
if [ ! -f "POSCAR" ]; then
    cp "$posfile" POSCAR > /dev/null
fi
echo -e "402 \n 1 \n 3 \n $z_min $z_max\n 2\n all" |vaspkit > /dev/null
unfix_num=$(grep "T" POSCAR_FIX.vasp |wc -l)
unfix_num=$((unfix_num - 1))
echo "Unfixed atom numbers: $unfix_num"
}


# =============================================================================
# SECTION 4: Analysis / 分析
# =============================================================================

# ---------------------------------------------------------------------------
# Function: outcar_get_virial
# 功能: 从VASP的OUTCAR文件中提取总位力(virial/stress)
#       提取FORCE on cell =-STRESS部分的Total行，输出9个分量。
# 场景: 需要从完成的VASP计算中提取应力张量时使用。
#       例如后续用于构建包含位力标签的NEP训练集。
# Usage: outcar_get_virial [outcar_file]
# Example:
#   outcar_get_virial OUTCAR
#   outcar_get_virial OUTCAR > virial.txt
# ---------------------------------------------------------------------------
outcar_get_virial() {
  local outcar=${1:-OUTCAR}
  grep -A 20 "FORCE on cell =-STRESS" $outcar | grep "Total" | tail -n 1 | awk '{ print $2,$5,$7,$5,$3,$6,$7,$6,$4 }'
}


# ---------------------------------------------------------------------------
# Function: generate_band_plot
# 功能: 从phonopy输出的band.yaml文件绘制声子色散图
#       使用Python的yaml和matplotlib读取band.yaml并画图。
# 场景: 完成VASP+phonopy声子谱计算后，使用此函数快速可视化声子色散。
# Usage: generate_band_plot <band.yaml>
# Example:
#   generate_band_plot band.yaml
#   # 输出文件: phonon_band_structure.png
# Dependencies: pyyaml, numpy, matplotlib
# ---------------------------------------------------------------------------
generate_band_plot() {
    if [ $# -ne 1 ]; then
        echo "Usage: generate_band_plot <band.yaml>"
        return 1
    fi

    band_file=$1

    if [ ! -f "$band_file" ]; then
        echo "File not found: $band_file"
        return 1
    fi

    python3 <<EOF
import yaml
import numpy as np
import matplotlib.pyplot as plt

with open("$band_file", 'r') as file:
    data = yaml.safe_load(file)

q_points = [entry['distance'] for entry in data['phonon']]

frequencies = []
for entry in data['phonon']:
    freqs = [band['frequency'] for band in entry['band']]
    frequencies.append(freqs)

frequencies = list(map(list, zip(*frequencies)))

for freq in frequencies:
    plt.plot(q_points, freq)

plt.xlabel('Wave Vector')
plt.ylabel('Frequency (THz)')
plt.title('Phonon Band Structure')
plt.savefig('phonon_band_structure.png')
plt.show()
EOF

    echo "Phonon band structure plot generated."
}
