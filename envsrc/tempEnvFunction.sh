# 该库使用来放置临时函数，多为针对性极强的函数，泛用性一般没有


del_caf2xyz_F(){
# 删除caf2文件中的F原子,用于处理xyz文件
filename=${1:-"dump.xyz"}
file2name=${2:-"del_${filename}"}
> caf2_tmpfile.cpp
cat >> caf2_tmpfile.cpp <<EOF
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
using namespace std;
pair<string, string> get_needReplace(const string& filename) {
    ifstream file(filename);
    string firstLine;
    getline(file, firstLine);  // 读取第一行
    int needReplace = stoi(firstLine);
    string replace2 = to_string(needReplace / 3);
    string needReplaceStr = to_string(needReplace);
    file.close();
    return make_pair(needReplaceStr, replace2);
}

void del_F(const string& filename) {
    ifstream input_file(filename);
    ofstream output_file("$file2name");

    pair<string, string> result = get_needReplace(filename);
    string needReplace = result.first;
    string replace2 = result.second;

    string line;
    while (getline(input_file, line)) {
        if (line.empty()) {
            continue;  // 跳过空行
        }

        stringstream ss(line);
        string token;
        vector<string> tokens;
        while (getline(ss, token, ' ')) {
            tokens.push_back(token);
        }

        if (tokens[0] == needReplace) {
            tokens[0] = replace2;
        } else if (tokens[0] == "F") {
            continue;
        }

        for (size_t i = 0; i < tokens.size(); ++i) {
            output_file << tokens[i];
            if (i != tokens.size() - 1) {
                output_file << " ";
            }
        }
        output_file << "\n";
    }

    input_file.close();
    output_file.close();
}

int main() {
    string filename = "$filename";

    auto start_time = std::chrono::high_resolution_clock::now();
    del_F(filename);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    cout << "Using Time: " << elapsed_time << " ms" << endl;

    return 0;
}
EOF

g++  -O3 caf2_tmpfile.cpp -o caf2_tmpfile
./caf2_tmpfile
rm caf2_tmpfile.cpp caf2_tmpfile
}


verify_deform(){
#临时函数，用来验算拉伸曲线的结果
    verify_gpumd_result
    cd verify_*
    sed -i  "/dump_thermo/c\dump_thermo  100" run.in
    sed -i  "/run/c\run 1000000" run.in
    gpumdstart
}

deal_box(){
    #临时函数，用来处理box.raw文件
    filename=${1:-"box.raw"}
    python3 <<EOF
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("$filename")
Lx ,Ly ,Lz = data[:,0] ,data[:,4],data[:,8]
Volume = Lx*Ly*Lz
plt.plot(Volume)
#print(data)
plt.savefig("box.png")
EOF
}


gpumdstart_dcu(){
#本程序建立是为了简化在曙光上使用gpumd的流程，该程序默认任务名为执行命令的文件夹的名字，默认使用1dcu与1cpu,也可以使用-n参数指定使用dcu的数量
dcu_num=1
for arg in "$@"; do 
   if [[ $arg = "-n" ]]; then
     dcu_num=${2:-1}
     break
   fi
done
gpumd_version=${gpumd_version:"GPUMD"}
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
/work/share/acmtrwrxv5/${gpumd_version}/src/gpumd" | sbatch -J "$job_name"

}

update_cp2k_inp_cell_from_xyz() {
    # 辅助函数：更新CP2K输入文件中CELL_PARAMETERS部分的晶胞参数
    # 使用方式：update_cp2k_inp_cell_from_xyz model.xyz cp2k.inp
    # 达成目的：修改cp2k输入文件中的CELL部分

    local xyz_file="$1"
    local cp2k_inp="$2"
    Lattice=$(get_Lattice $1 | grep -oP '(?<=Lattice=").*(?=")')
    cell_A=$(echo $Lattice | awk '{print $1,$2,$3}')
    cell_B=$(echo $Lattice | awk '{print $4,$5,$6}')
    cell_C=$(echo $Lattice | awk '{print $7,$8,$9}')
    sed -E "/^\s*A\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      A   $cell_A/"       $cp2k_inp |sed -E "/^\s*B\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      B   $cell_B/ " |sed -E "/^\s*C\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      C   $cell_C/" >      ${cp2k_inp%.*}_up.inp

    sed -i "/@SET XYZFILE/s/.*/@SET XYZFILE    $1/"  ${cp2k_inp%.*}_up.inp
}


cp2kstart() {
# 临时函数，用来在曙光上启动 CP2K
# 使用方式：cp2kstart cp2k.inp 32
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

find_dcu_speed() {
    find $PWD -type f -name "std.out*" -exec sh -c '
        for file; do
            grep -A 10 "GPU information" "$file"
            grep "Speed of this run =" "$file" | awk -F "=" "{print \$2}"
        done
    ' sh {} + > dcu.txt
}

build_stages(){
# 用于臭氧的项目
# 该函数使用来构建多个stage文件夹，并且将输入文件放进去
    local init=$PWD
    for i in {1..4};do
        mkdir -p stage$i
        cp  *.inp *stage${i}.xyz stage$i
        cd stage$i
        echo "----->stage${i} : $PWD" 
        update_cp2k_inp_cell_from_xyz *stage${i}.xyz *.inp
        cp2kstart Temp_up_up.inp
        cd $init
    done
}

check_stage_energy(){
# 该函数是用来检查各个阶段的能量
    echo "stage1: "
    grep "ENERGY|" stage1/cp2k.log | head -n 3
    echo "stage2: "
    grep "ENERGY|" stage2/cp2k.log | head -n 3
    echo "stage3: "
    grep "ENERGY|" stage3/cp2k.log | head -n 3
    echo "stage4: "
    grep "ENERGY|" stage4/cp2k.log | head -n 3
}



vaspstart_geo_optstage1(){
# 该函数为临时函数，使用来快速的构建vasp的结构优化的任务
    add_potcar
    cat > INCAR <<EOF
 ISTART =  0            (Read existing wavefunction, if there)
 ISPIN  =  1            (Non-Spin polarised DFT)
 LREAL  = Auto       (Projection operators: automatic)
 ENCUT  =  300        (Cut-off energy for plane wave basis set, in eV)
 PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
 LWAVE  = .FALSE.        (Write WAVECAR or not)
 LCHARG = .FALSE.        (Write CHGCAR or not)
 ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
 ISIF  =  2
 ALGO  =  Fast
 IBRION = 2
 #EDIFFG = -0.02
 
 Static Calculation
 ISMEAR =  0            (gaussian smearing method)
 SIGMA  =  0.05        (please check the width of the smearing)
 LORBIT =  11           (PAW radii for projected DOS)
 NEDOS  =  2001         (DOSCAR points)
 #NELM   =  260           (Max electronic SCF steps)
 NSW    =   30
 EDIFF  =  1E-02      (SCF energy convergence, in eV)
 
 LUSE_VDW = .TRUE.
 IVDW = 11
 KSPACING = 0.3
EOF
    vasprun_dcu 4

}

vaspstart_geo_optstage2(){
# 该函数为临时函数，使用来快速的构建vasp的结构优化的任务
    add_potcar
    cat > INCAR <<EOF
 ISTART =  0            (Read existing wavefunction, if there)
 ISPIN  =  1            (Non-Spin polarised DFT)
 LREAL  = Auto       (Projection operators: automatic)
 ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
 PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
 LWAVE  = .FALSE.        (Write WAVECAR or not)
 LCHARG = .FALSE.        (Write CHGCAR or not)
 ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
 ISIF  =  2
 ALGO  =  Fast
 IBRION = 2
 #EDIFFG = -0.02
 
 Static Calculation
 ISMEAR =  0            (gaussian smearing method)
 SIGMA  =  0.05        (please check the width of the smearing)
 LORBIT =  11           (PAW radii for projected DOS)
 NEDOS  =  2001         (DOSCAR points)
 #NELM   =  260           (Max electronic SCF steps)
 NSW    =   100
 EDIFF  =  1E-04     (SCF energy convergence, in eV)
 LUSE_VDW = .TRUE.
 IVDW = 11
 KSPACING = 0.3
EOF
    vasprun_dcu 4

}

vaspstart_geo_optstage3(){
# 该函数为临时函数，使用来快速的构建vasp的结构优化的任务
    add_potcar
    cat > INCAR <<EOF
 ISTART =  0            (Read existing wavefunction, if there)
 ISPIN  =  1            (Non-Spin polarised DFT)
 LREAL  = Auto       (Projection operators: automatic)
 ENCUT  =  500        (Cut-off energy for plane wave basis set, in eV)
 PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
 LWAVE  = .FALSE.        (Write WAVECAR or not)
 LCHARG = .FALSE.        (Write CHGCAR or not)
 ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
 ISIF  =  2
 ALGO  =  Fast
 IBRION = 2
 EDIFFG = -0.02
 
 Static Calculation
 ISMEAR =  0            (gaussian smearing method)
 SIGMA  =  0.05        (please check the width of the smearing)
 LORBIT =  11           (PAW radii for projected DOS)
 NEDOS  =  2001         (DOSCAR points)
 #NELM   =  260           (Max electronic SCF steps)
 NSW    =   100
 EDIFF  =  1E-04     (SCF energy convergence, in eV)
 LUSE_VDW = .TRUE.
 IVDW = 11
 KSPACING = 0.2
EOF
    vasprun_dcu 4

}



vaspstart_single_energy_after_relaxation(){
    # 该函数为临时函数，使用来快速的构建vasp的单点能量计算的任务
    mkdir -p calc_SE
    cd calc_SE
    cp ../CONTCAR POSCAR
    add_potcar
    cat > INCAR <<EOF
ISTART =  0            
ISPIN  =  1           
LREAL  = Auto      
ENCUT  =  500     
PREC   =  Accurate   
LWAVE  = .FALSE.      
LCHARG = .FALSE.      
ADDGRID= .TRUE.      
ISIF  =  2
ALGO  =  Fast
IBRION = 0

Static Calculation
ISMEAR =  0            
SIGMA  =  0.05        
LORBIT =  11           
NEDOS  =  2001        
NELM   =  120         
EDIFF  =  1E-05        
LUSE_VDW = .TRUE.
IVDW = 11
KSPACING = 0.2 
EOF
    echo "task adress : $PWD"
    vasprun
    echo "Single-point energy calculation"
}


vaspstart_single_energy(){
    # 该函数为临时函数，使用来快速的构建vasp的单点能量计算的任务
    add_potcar
    cat > INCAR <<EOF
ISTART =  0            
ISPIN  =  1           
LREAL  = Auto      
ENCUT  =  500     
PREC   =  Accurate   
LWAVE  = .FALSE.      
LCHARG = .FALSE.      
ADDGRID= .TRUE.      
ISIF  =  2
ALGO  =  Fast
IBRION = 0

Static Calculation
ISMEAR =  0            
SIGMA  =  0.05         
LORBIT =  11           
NEDOS  =  2001        
NELM   =  120         
EDIFF  =  1E-05        
LUSE_VDW = .TRUE.
IVDW = 11 
KSPACING = 0.2 
EOF
    echo "task adress : $PWD"
    vasprun
    echo "Single-point energy calculation"
}

wait_geo_run_SE(){
    # 该函数用于检查当前文件夹中的 VASP 任务是否完成，如果已完成，则自动提交单点能计算任务
    > wait.log
    while true; do
        if grep -q "General timing and accounting informations for this job" OUTCAR; then
            echo -e "$(date) - Geometry optimization completed, starting single-point energy calculation." >> wait.log
            vaspstart_single_energy_after_relaxation
            break
        else
            echo -e "$(date) - Geometry optimization is running, waiting 100s\n" | tee -a wait.log
            sleep 100
        fi
    done
}

wait_complete_vasp_run(){
    # 该函数用于检查当前文件夹中的 VASP 任务是否完成，如果已完成，进行下一个任务
    local log_file="wait.log"
    > "$log_file" 

    local max_retries=1000  
    local retry_interval=100  # 每次重试间隔时间（秒）
    local retries=0

    while [ $retries -lt $max_retries ]; do
        if [ -f "OUTCAR" ] && grep -q "General timing and accounting informations for this job" OUTCAR; then
            echo -e "$(date) - Geometry optimization completed, starting single-point energy calculation." | tee -a "$log_file"
            eval "$1"  # 执行传递的命令
            break
        else
            echo -e "$(date) - Geometry optimization is running, waiting ${retry_interval}s\n" | tee -a "$log_file"
            sleep $retry_interval
            ((retries++))
        fi
    done

    if [ $retries -ge $max_retries ]; then
        echo -e "$(date) - Maximum retries reached, geometry optimization may have failed or stalled." | tee -a "$log_file"
        exit 1
    fi
}

fix_adsorption_model(){
# 该脚本使用来创建固定的文件夹并对指定的文件进行固定
    local model_name=$1
    local max_high=$2
    mkdir fix_$max_high
    cp $model_name fix_$max_high
    cd fix_$max_high
    fix_poscar_zCar $model_name $max_high
}

clean_poscar(){
# 该函数使用来转化poscar的格式，ase生成的poscar可能会有很多重复的元素，使用该函数将其进行清理
    poscar_to_xyz $1
    xyz_to_poscar model_conversed.xyz
}

use_pbc_run(){
    mkdir -p use_pbc
    cp del_poisson/elctron.inp del_poisson/*.xyz use_pbc/
    cd use_pbc/

    sed -i 's/PERIODIC NONE/PERIODIC XYZ/' elctron.inp 
    cp2kstart elctron.inp 

}

use_sccs_run(){
    mkdir -p use_sccs
    cp elctron.inp *xyz use_sccs/
    cd use_sccs/

    sed -i '/&END SCF/a\    &SCCS\n      ALPHA [N*m^-1] 0.0\n      BETA [kbar] 0.0\n      GAMMA [mN/m] 0.0\n      DIELECTRIC_CONSTANT 78.36 #Water\n      EPS_SCCS 1E-6 #Default. Requested accuracy for the SCCS iteration cycle\n      EPS_SCF 0.5 #Default. SCCS iteration is activated only if SCF iteration is converged to this threshold\n      MAX_ITER 100 #Default. Maximum number of SCCS iteration steps\n      DERIVATIVE_METHOD FFT #Default. Method for calculation of numerical derivatives. Can also be CD3, CD5, CD7\n      &ANDREUSSI\n        RHO_MAX 0.0035 #Default\n        RHO_MIN 0.0001 #Default\n      &END ANDREUSSI\n    &END SCCS' elctron.inp 

    cp2kstart elctron.inp
}

# 该函数用于重复运行gpumd,将会创建n个文件夹，将model.xyz,nep.txt,run.in复制到每个文件夹中
repeat_buildgpumd(){
    local n=$1
    for ((i=1;i<=$n;i++))
    do
        mkdir -p run$i
        cp model.xyz nep.txt run.in run$i/
    done
}

# 该函数用于重复运行gpumd,将会创建n个文件夹，将model.xyz,nep.txt,run.in复制到每个文件夹中
repeat_rungpumd(){
    local n=$1
    for ((i=1;i<=$n;i++))
    do
        cd run$i
        free_time_run 'nohup gpumd 2>&1 &'
        sleep 5
        cd ..
    done
}

cp2desktop(){
    cp -r $1 /mnt/c/Users/rebreath.REBREATH-TX4/Desktop
}

mv2desktop(){
    mv $1 /mnt/c/Users/rebreath.REBREATH-TX4/Desktop
}

# AlN临时函数，对所有的dump.xyz分析结晶原子数
# 该函数将会查找当前文件夹下所有的dump.xyz文件，并对其进行并行化晶型分析
AlN_analyze_phase() {
    while IFS= read -r -d '' file; do
        dir=$(dirname "$file")
        base=$(basename "$file")

        echo "Processing: $file"

        (
            cd "$dir" || exit 1
            analyze_crystallinity_fraction SC "$base"
            analyze_crystallinity_fraction HEX_DIAMOND "$base"
        ) &

        # (
        #     cd "$dir" || exit 1
        #     analyze_crystallinity_fraction HEX_DIAMOND "$base"
        # ) &
    done < <(find . -name "dump.xyz" -print0)

    wait
}

# 该函数用于检查所有的dump_SC.txt文件的最后一列是否为0
# 如果不为0，则输出该文件的路径和最后一列的值
check_dump_sc_lastcol() {
    find . -name "dump_SC.txt" -print0 | while IFS= read -r -d '' file; do
        value=$(awk 'NF{last=$2} END{print last}' "$file")

        if [ -z "$value" ]; then
            echo "Empty or invalid file: $file"
            continue
        fi

        awk -v v="$value" 'BEGIN { exit (v == 0 ? 0 : 1) }'
        if [ $? -ne 0 ]; then
            echo "$(dirname "$file")  -> last line second column = $value"
        fi
    done
}


check_dump_hex_lastcol() {
    find . -name "dump_HEX_DIAMOND.txt" -print0 | while IFS= read -r -d '' file; do
        # 提取最后一行有内容的行的第二列值
        value=$(awk 'NF{last=$2} END{print last}' "$file")

        # 处理空文件或无有效内容的情况
        if [ -z "$value" ]; then
            echo "Empty or invalid file (no valid content): $file"
            continue
        fi

        # 检查值是否不等于 1
        awk -v v="$value" 'BEGIN { exit (v == 1 ? 0 : 1) }'
        if [ $? -ne 0 ]; then
            echo "$(dirname "$file")  -> last line second column = $value (expected 1)"
        fi
    done
}


# 该函数将会对thermo.out文件进行画图，分析出平均压强，体积变化率，偏应力等信息
analyze_thermo_out() {
    # 用法:
    #   analyze_thermo_out thermo.out [thermo_every] [thermo_first_step]
    #
    # 例如:
    #   analyze_thermo_out thermo.out 100 100
    #
    # 输出:
    #   thermo_analysis_data.txt
    #   thermo_analysis_3panel.png
    #
    # 说明:
    #   thermo.out 无表头，列定义为：
    #   1:T 2:K 3:U 4:Pxx 5:Pyy 6:Pzz 7:Pyz 8:Pxz 9:Pxy
    #   10:ax 11:ay 12:az 13:bx 14:by 15:bz 16:cx 17:cy 18:cz

    local thermo_file="${1:-thermo.out}"
    local thermo_every="${2:-100}"
    local thermo_first_step="${3:-100}"

    if [[ ! -f "$thermo_file" ]]; then
        echo "Error: File not found: $thermo_file"
        return 1
    fi

    python3 - "$thermo_file" "$thermo_every" "$thermo_first_step" << 'PY'
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

thermo_file = Path(sys.argv[1])
thermo_every = int(sys.argv[2])
thermo_first_step = int(sys.argv[3])

# =========================
# 读取 thermo.out
# =========================
rows = []
with thermo_file.open("r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 18:
            continue
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue

if len(rows) == 0:
    raise ValueError(f"No valid data found in {thermo_file}")

data = np.array(rows, dtype=float)

# =========================
# 取出应力分量
# 注意：
# thermo.out 列从 1 开始计数
# Python 下标从 0 开始
# 4:Pxx 5:Pyy 6:Pzz 7:Pyz 8:Pxz 9:Pxy
# => 3,4,5,6,7,8
# =========================
pxx = data[:, 3]
pyy = data[:, 4]
pzz = data[:, 5]
pyz = data[:, 6]
pxz = data[:, 7]
pxy = data[:, 8]

# =========================
# 计算体积
# 晶胞体积 V = |a · (b × c)|
# a=(10,11,12), b=(13,14,15), c=(16,17,18)
# Python 下标 => a=(9,10,11), b=(12,13,14), c=(15,16,17)
# =========================
a = data[:, 9:12]
b = data[:, 12:15]
c = data[:, 15:18]

volumes = np.abs(np.einsum("ij,ij->i", a, np.cross(b, c)))
v0 = volumes[0]

# 你指定要画的量：(Vt - V0) / V0
rel_volume_change = (volumes - v0) / v0

# =========================
# 平均压强
# 压缩时通常希望画成正值，所以取负号
# P = -(Pxx + Pyy + Pzz)/3
# =========================
p_avg = -(pxx + pyy + pzz) / 3.0

# =========================
# 偏应力
# 先构造偏应力张量:
# s_ij = sigma_ij - mean(sigma)*delta_ij
#
# 这里采用等效偏应力（类似 von Mises）
# q = sqrt( ((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 + 6*(sxy^2+syz^2+sxz^2))/2 )
#
# 对应路径研究中，这是一个很有用的“非静水性”指标
# =========================
mean_normal = (pxx + pyy + pzz) / 3.0
sxx = pxx - mean_normal
syy = pyy - mean_normal
szz = pzz - mean_normal
sxy = pxy
syz = pyz
sxz = pxz

deviatoric_q = np.sqrt(
    (
        (sxx - syy)**2 +
        (syy - szz)**2 +
        (szz - sxx)**2 +
        6.0 * (sxy**2 + syz**2 + sxz**2)
    ) / 2.0
)

# =========================
# 横轴 step
# =========================
n = len(volumes)
steps = thermo_first_step + np.arange(n) * thermo_every

# =========================
# 导出数据
# =========================
out_txt = thermo_file.with_name("thermo_analysis_data.txt")
header = "step volume rel_volume_change pressure_avg deviatoric_q"
np.savetxt(
    out_txt,
    np.column_stack([steps, volumes, rel_volume_change, p_avg, deviatoric_q]),
    fmt=["%d", "%.10f", "%.10f", "%.10f", "%.10f"],
    header=header,
    comments=""
)

# =========================
# 画图：三联图
# =========================
fig, axes = plt.subplots(3, 1, figsize=(7.5, 9.0), sharex=True)

# 图1：相对体积变化
axes[0].plot(steps, rel_volume_change, linewidth=1.2)
axes[0].set_ylabel(r'$(V_t - V_0)/V_0$')
axes[0].set_title(thermo_file.parent.name)

# 图2：平均压强
axes[1].plot(steps, p_avg, linewidth=1.2)
axes[1].set_ylabel('Average pressure')

# 图3：偏应力
axes[2].plot(steps, deviatoric_q, linewidth=1.2)
axes[2].set_ylabel('Deviatoric stress')
axes[2].set_xlabel('Step')

for ax in axes:
    ax.tick_params(direction="in")
    ax.grid(False)

fig.tight_layout()

out_png = thermo_file.with_name("thermo_analysis_3panel.png")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Data saved to: {out_txt}")
print(f"Figure saved to: {out_png}")
PY
}

deal_AlN_diff_axis_deform(){
# 处理AlN不同轴的形变数据
    initpath=$PWD; export -f analyze_mulcrystallinity_fraction; export -f analyze_thermo_out; find ./ -name "dump.xyz" | while read -r i; do
    dir=$(dirname "$i")
    cd "$dir" || continue
    bash -c "analyze_mulcrystallinity_fraction  dump.xyz sc ; analyze_thermo_out" &
    cd "$initpath"
    done

    find ./ -name "dump_mulcrystallinity_fraction.txt" | while read -r file; do
    num=$(awk 'END{print $2}' "$file")
    num=${num:-0}
    # 在文件所在目录生成标记
    touch "$(dirname "$file")/lastsc:$num"
    echo "已生成: $(dirname "$file")/lastsc:$num"
    done
}

# 提交当前目录的 LAMMPS 任务
submit_lmp_matpl() {
    # 1. 获取交互参数：若不传参，则默认使用 lmp.in 和 4 核
    local lmp_input=${1:-lmp.in}
    local corenum=${2:-4}

    # 检查指定的输入文件是否存在，若不存在则提前报错拦截
    if [ ! -f "$lmp_input" ]; then
        echo "❌ 错误: 当前目录下未找到输入文件 [$lmp_input]！"
        return 1
    fi

    # 2. 自动获取当前目录的名称作为作业名
    local JOB_NAME=$(basename "$PWD")
    local SBATCH_FILE="submit_run.sh"

    echo "🚀 正在为当前目录 [$JOB_NAME] 生成并提交 Slurm 任务 (文件: $lmp_input, 核心数: $corenum)..."

    # 3. 动态生成符合当前目录的 Slurm 提交脚本
    cat << EOF > "$SBATCH_FILE"
#!/bin/bash
#SBATCH -p pg_g4J4
#SBATCH -N 1
#SBATCH --ntasks-per-node=$corenum
#SBATCH --cpus-per-task=1
#SBATCH -J ${JOB_NAME}
#SBATCH -D ${PWD}

module purge
module use /public/share/pguan_group/modules
module load python/ase-lammps-kokkos-nep

echo "===== ENV CHECK ====="
module list
echo "Working Directory: \$PWD"
which lmp
echo "====================="

# 运行 LAMMPS 计算（已修改为动态变量）
mpirun -np $corenum lmp -in $lmp_input
EOF

    # 4. 提交任务并删除临时生成的 sbatch 脚本文件（保持目录干净）
    sbatch "$SBATCH_FILE" && rm -f "$SBATCH_FILE"
}


# 提交当前目录的 GPUMD 任务
submit_gpumd() {
    # 1. 获取交互参数：若不传参，则默认输入文件为 run.in
    local input_file=${1:-run.in}
    
    # 检查输入文件是否存在
    if [ ! -f "$input_file" ]; then
        echo "❌ 错误: 当前目录下未找到输入文件 [$input_file]！"
        echo "💡 提示: gpumd 默认需要 run.in，如果是其他名称请指定: submit_gpumd <文件名>"
        return 1
    fi

    # 2. 自动获取当前目录名称作为作业名
    local JOB_NAME=$(basename "$PWD")
    local SBATCH_FILE="submit_gpumd.sh"

    echo "🚀 正在为当前目录 [$JOB_NAME] 准备 GPUMD 任务..."

    # 3. 动态生成 Slurm 提交脚本
    cat << EOF > "$SBATCH_FILE"
#!/bin/bash
#SBATCH -p pg_g4J4
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH -J ${JOB_NAME}
#SBATCH -D ${PWD}

# 加载自定义模块路径
module purge
module use /public/share/pguan_group/modules
module load gpumd/latest

echo "===== ENV CHECK ====="
module list
echo "Working Directory: \$PWD"
which gpumd
echo "====================="

# 运行 GPUMD
# 注意：如果输入文件不是默认的 run.in，gpumd 通常需要重定向或通过特定方式指定
# 但标准 gpumd 默认直接读取当前目录下的 run.in
if [ "$input_file" != "run.in" ]; then
    echo "⚠️ 警告: gpumd 通常默认读取 run.in。当前指定为 $input_file，尝试链接..."
    ln -sf $input_file run.in
fi

gpumd
EOF

    # 4. 提交任务并清理
    sbatch "$SBATCH_FILE" && rm -f "$SBATCH_FILE"
}

alias runlmp=submit_lmp_matpl
alias rungpumd=submit_gpumd
