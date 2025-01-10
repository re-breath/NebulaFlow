# 该函数用来存放VASP环境变量设置和函数

vasprun(){
# 该函数用来提交VASP任务到SLURM队列中
# 调用方式：vasprun [corenum]

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
export MKL_DEBUG_CPU_TYPE=5 # 加速代码
export MKL_CBWR=AVX2 # 使cpu默认支持avx2
export I_MPI_PIN_DOMAIN=numa # 内存位置与cpu位置绑定，加速内存读取。对于内存带宽要求高的计算提速明显
echo "Starting VASP job at $(date)"
srun --mpi=pmi2 vasp_std
echo "VASP job completed at $(date)"
EOF
    echo "Submitting VASP job $taskname with $corenum cores..."
    sbatch "$slurm_script"
    local sbatch_status=$?
    if [ $sbatch_status -eq 0 ]; then
        echo "VASP job $taskname submitted successfully."
        # rm "$slurm_script"
    else
        echo "Error: Failed to submit VASP job '$taskname'."
        return 1
    fi
}


vasprun_dcu() {
    # 该函数用来提交VASP任务到SLURM队列中
    # 调用方式：vasprun [corenum]
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
source $HOME/apprepo/vasp/6.4.2-dtk23.10_cpvasp/scripts/env.sh
# 指定版本, std, gam 或 ncl 版
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
    echo "Submitting VASP job $taskname with $corenum_dcu DCU cores..."
    sbatch "$slurm_script"
    local sbatch_status=$?
    if [ $sbatch_status -eq 0 ]; then
        echo "VASP job $taskname submitted successfully."
        # rm "$slurm_script"
    else
        echo "Error: Failed to submit VASP job '$taskname'."
        return 1
    fi
}

add_potcar() {
# 该函数用来添加POTCAR文件到当前目录
    echo "103" | vaspkit > /dev/null
}


add_kspacing_to_incar(){
#如果kspacing的值不为零，则添加kspacing
    kspacing=${1:-0.2}
    mv KPOINTS KPOINTS_backup || true
    if [ -z "$(grep "KSAPCING" INCAR)" ];then
        sed -i "\$a\KSPACING = $kspacing "  INCAR
    else
        sed "/KSPACING/c\KSPACING = $kspacing/" INCAR
    fi
}

fix_poscar_zFrc() {
# 该函数将会固定z轴的下半部分，可以指定一个分数来，固定该分数内的原子
# 调用方式：fix_poscar_z [frac] file
# 使用ase计算出z轴的下半部分的位置，然后调动vaspkit来固定
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
echo "unfix atom numbers : $unfix_num"
}


fix_poscar_zCar() {
# 该函数将会固定z轴的下半部分，指定最大和最小的z值
# 调用方式：fix_poscar_zCar file [z_max] [z_min] 

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
echo "unfix atom numbers : $unfix_num"
}

clean_vasp_out(){
# 该函数用来删除VASP输出文件
    rm nohup.out CONTCAR OSZICAR WAVECAR PCDAT REPORT slurm-* vasp* config.* DOSCAR C* wait.log XDATCAR IBZKPT EIGENVAL REPORT OUTCAR PROCAR > /dev/null
}

add_kpoints() {
# 该函数用来添加KPOINTS文件到当前目录
    local level=${1:-0.03}
    echo -e "102 \n 1 \n $level" | vaspkit > /dev/null
}

generate_band_plot() {
    #函数说明：该函数用来指定DFPT的band.yaml文件拿来画图
    # 检查是否传入了 band.yaml 文件路径
    if [ $# -ne 1 ]; then
        echo "用法: generate_band_plot <band.yaml文件路径>"
        return 1
    fi

    band_file=$1

    # 检查文件是否存在
    if [ ! -f "$band_file" ]; then
        echo "文件不存在: $band_file"
        return 1
    fi

    # 调用Python脚本处理band.yaml并使用matplotlib绘图
    python3 <<EOF
import yaml
import numpy as np
import matplotlib.pyplot as plt

# 读取YAML文件
with open("$band_file", 'r') as file:
    data = yaml.safe_load(file)

# 提取q点的路径
q_points = [entry['distance'] for entry in data['phonon']]

# 提取所有频率
frequencies = []
for entry in data['phonon']:
    freqs = [band['frequency'] for band in entry['band']]
    frequencies.append(freqs)

# 转置列表以将解包频率与其对应的q点对齐
frequencies = list(map(list, zip(*frequencies)))

# 绘制每个频率
for freq in frequencies:
    plt.plot(q_points, freq)

# 设置坐标轴标签和标题
plt.xlabel('Wave Vector')
plt.ylabel('Frequency (THz)')
plt.title('Phonon Band Structure')
plt.savefig('phonon_band_structure.png')
plt.show()
EOF

    echo "声子谱图生成完毕。"
}

