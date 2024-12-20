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
source /work/home/rebreath/apprepo/vasp/6.4.3-ioptcell_intelmpi2017_hdf5_libxc/scripts/env.sh
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

