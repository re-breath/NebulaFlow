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



vaspstart_geo_opt(){
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
 EDIFFG = -0.02
 
 Static Calculation
 ISMEAR =  0            (gaussian smearing method)
 SIGMA  =  0.05         (please check the width of the smearing)
 LORBIT =  11           (PAW radii for projected DOS)
 NEDOS  =  2001         (DOSCAR points)
 #NELM   =  260           (Max electronic SCF steps)
 NSW    =   160
 EDIFF  =  1E-04        (SCF energy convergence, in eV)
 IVDW = 12
 KSPACING = 0.2
EOF
    vasprun_dcu 4

}

vaspstart_single_energy(){
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
NELM   =  260         
EDIFF  =  1E-05        
 
KSPACING = 0.2 
EOF
    echo "task adress : $PWD"
    vasprun
    echo "Single-point energy calculation"
}

