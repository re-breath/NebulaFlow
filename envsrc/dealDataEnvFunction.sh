# ======================================================================
# File:         dealDataEnvFunction.sh
# Project:      NebulaFlow
# Description:  通用数据处理函数库 — 数据筛选/统计分析/模型分析/格式转换等
#               General data processing: filtering, statistics, model analysis,
#               file format conversion, structure manipulation, and more.
# Author:       rebreath
# Dependencies: python3, numpy, matplotlib, ase, ovito (optional), pymatgen (optional)
# ======================================================================


# =============================================================================
# SECTION 1: File Statistics & Analysis / 文件统计与分析
# =============================================================================

# ---------------------------------------------------------------------------
# Function: find_column_max
# 功能: 查找数据文件中指定列的最大值
# 场景: 分析thermo.out中温度/压力的最大值，或力文件中力分量的最大值。
# Usage: find_column_max <file> <column_number>
# Example:
#   find_column_max thermo.out 4      # 查找第4列(Px)的最大值
#   find_column_max force_train.out 3 # 查找力的x分量最大值
# ---------------------------------------------------------------------------
find_column_max(){
    if [[ $# -ne 2 ]]; then
    echo "Error: Usage: find_column_max <file> <column>"
    return 1
    fi
    filename=$1
    column=$2
    if [[ ! -f "$filename" || ! -r "$filename" ]]; then
        echo "Error: File '$filename' does not exist or is not readable."
        return 1
    fi
    if [[ ! $column =~ ^[0-9]+$ ]]; then
        echo "Error: Column number must be a positive integer."
        return 1
    fi
    awk -v col="$column" 'BEGIN {max=0}
        NR > 0 && ($col > max) {max=$col}
        END { printf("Column %d max: %f\n", col, max) }' "$filename"
}


# ---------------------------------------------------------------------------
# Function: find_column_abs_max
# 功能: 查找数据文件中指定列的绝对值最大值，并输出所在行信息
# 场景: 需要定位thermo.out或force文件中极端值所在的具体行时使用。
# Usage: find_column_abs_max <file> <column>
# Example:
#   find_column_abs_max thermo.out 4  # 查找第4列绝对值最大的位置
# ---------------------------------------------------------------------------
find_column_abs_max ()
{
    if [[ $# -ne 2 ]]; then
        echo "Usage: find_column_abs_max <file> <column>"
        return 1
    fi
    local filename="$1"
    local column="$2"
    if [[ ! -f "$filename" ]]; then
        echo "Error: File '$filename' not found"
        return 1
    fi
    if [[ ! $column =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: Column must be a positive integer"
        return 1
    fi
    awk -v col="$column" '
    BEGIN { max_abs = -1; max_val = 0; max_line = ""; max_nr = 0 }
    NF == 0 { next }
    NF < col { next }
    $col !~ /^[-+]?[0-9]*\.?[0-9]+$/ { next }
    {
        current = $col
        current_abs = (current < 0) ? -current : current
        if (current_abs > max_abs) {
            max_abs = current_abs; max_val = current; max_line = $0; max_nr = NR
        }
    }
    END {
        if (max_nr == 0) { print "Warning: No valid numeric rows found"; exit }
        printf "Column %d absolute max:\n", col
        printf "  Value: %f\n", max_val
        printf "  Line:  %d\n", max_nr
        printf "  Content: %s\n", max_line
    }' "$filename"
}


# ---------------------------------------------------------------------------
# Function: analysis_column
# 功能: 对数据文件的某一列进行完整的统计分析
#       输出: 最大值、最小值、绝对最大值、绝对最小值、平均值、标准差。
# 场景: 分析MD模拟中温度/压力/能量/体积等物理量的统计特征。
# Usage: analysis_column <file> <column> [skiprows]
# Example:
#   analysis_column thermo.out 1 0       # 分析第1列(温度)，不跳过行
#   analysis_column thermo.out 4 100     # 分析第4列(压力)，跳过前100行(平衡阶段)
# ---------------------------------------------------------------------------
analysis_column(){
    local filename=$1
    local column_num=$2
    local skiprows=${3:-'0'}
    python3 << EOF
import numpy as np
filename = '$1'
data = np.loadtxt(filename,skiprows=$skiprows)
column_num = $column_num - 1
data = data[:,column_num]
print('Max:',np.max(data))
print('Min:',np.min(data))
print('Abs Max:',np.max(np.abs(data)))
print('Abs Min:',np.min(np.abs(data)))
print('Mean:',np.mean(data))
print('Std:',np.std(data))
EOF
}


# ---------------------------------------------------------------------------
# Function: get_col_average
# 功能: 计算数据文件某一列在指定范围内的平均值
# 场景: MD模拟后取平衡段的平均值。例如取thermo.out中1000-5000帧的平均温度。
# Usage: get_col_average <file> <col> [start_idx] [end_idx]
# Example:
#   get_col_average thermo.out 1 1000 3000   # 第1列[1000:3000]平均值
#   get_col_average thermo.out 4 500          # 第4列从500到末尾的平均值
# ---------------------------------------------------------------------------
get_col_average(){
    local filename=$1
    local col_index=$2
    local start_index=$3
    local end_index=$4
    if (( $# > 4 )); then
        echo "Error: Too many arguments."
        return 1
    fi
    if [ -z $3 ]; then start_index=0; fi
    if [ -z $4 ]; then end_index=None; fi
    python3 << EOF
import numpy as np
filename = '$filename'
col_index = int($col_index) - 1
min_index = int($start_index)
max_index = $end_index

data = np.loadtxt(filename)
if data.ndim == 1: data = data.reshape(-1, 1)
data = data[:, col_index]

def get_col_average(data, min_index, max_index):
    data = data[min_index:max_index]
    return np.mean(data)

average_value = get_col_average(data, min_index, max_index)
print(f"File: {filename}\nColumn: {col_index+1}\nRange: [{min_index}:{max_index}]\nAverage: {average_value:.8f}")
EOF
}


# =============================================================================
# SECTION 2: File Averaging / 文件平均
# =============================================================================

# ---------------------------------------------------------------------------
# Function: average_file
# 功能: 使用C++对多个数据文件进行逐行平均
#       将多个文件按行列对齐后取平均，输出到average.out。
# 场景: 多次独立MD模拟后需要对thermo量做统计平均以减小涨落。
#       例如：6次HNEMD的kappa.out取平均。
# Usage: average_file <file1> <file2> [file3...]
# Example:
#   average_file thermo_1.out thermo_2.out thermo_3.out
#   # 输出: average_ther.out
# ---------------------------------------------------------------------------
average_file(){
    local cpplib="$HOME/.rebreath/cpp_lib"
    first_file="$1"
    avg_name=${1:0:4}
    g++ $cpplib/averagefiles.cpp -o average_file
    ./average_file "$@"
    rm -f average_file
}


# ---------------------------------------------------------------------------
# Function: average_file_s
# 功能: 使用纯Shell对多个数据文件进行逐行平均（无需C++编译器）
#       功能同average_file但只用shell+awk实现。
# 场景: 当没有g++编译器时使用此版本的逐行平均。
# Usage: average_file_s <file1> <file2> [file3...]
# Example:
#   average_file_s thermo_1.out thermo_2.out thermo_3.out
#   # 输出: average_ther.out
# ---------------------------------------------------------------------------
average_file_s(){
    files_num="$#"
    tmp_file='temp'
    first_file="$1"
    avg_name=${1:0:4}
    paste "$@" >$tmp_file
    while IFS= read -r line
    do
        echo "$line" | awk -v files_num="$files_num" '{
            sum=0;
            for(i=1; i<=NF/files_num; i++){
                sum=0;
                for(j=0; j<files_num; j++){
                    sum += $(i+j*NF/files_num)
                }
                printf "%.15f\t", sum/files_num;
            }
            printf "\n";
        }'
    done < "$tmp_file" >average_${avg_name}.out
    rm -f "$tmp_file"
}


# ---------------------------------------------------------------------------
# Function: average_file_c
# 功能: 使用C++对多个文件逐行平均，并按第一个文件名前缀重命名输出
#       功能同average_file + 自动重命名。
# Usage: average_file_c <file1> <file2> [file3...]
# ---------------------------------------------------------------------------
average_file_c(){
    local cpplib="$HOME/.rebreath/cpp_lib"
    first_file="$1"
    avg_name=${1:0:4}
    g++ $cpplib/averagefiles.cpp -o average_file
    ./average_file "$@"
    rm -f average_file
    mv average.out average_${avg_name}.out
}


# =============================================================================
# SECTION 3: Quick Visualization / 快速可视化
# =============================================================================

# ---------------------------------------------------------------------------
# Function: replot
# 功能: 对两列数据文件进行快速画图（X-Y散点/连线图）
# 场景: 当你有一个两列数据文件(如E-frame.txt)，想快速看看趋势时使用。
#       比打开Python/Jupyter手动写代码快得多。
# Usage: replot <data_file>
# Example:
#   replot energy_vs_frame.txt
#   # 输出: plot.png
# ---------------------------------------------------------------------------
replot() {
    if [ "$#" -ne 1 ]; then
        echo "Usage: replot <data_file>"
        return 1
    fi
    python3 <<EOF
import matplotlib.pyplot as plt

datafile = "$1"
x, y = [], []
with open(datafile, 'r') as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) < 2: continue
        x_val, y_val = map(float, parts[:2])
        x.append(x_val); y.append(y_val)

plt.plot(x, y, marker='o', linestyle='-')
plt.xlabel('X'); plt.ylabel('Y')
plt.savefig('plot.png', format='png')
EOF
}


# ---------------------------------------------------------------------------
# Function: plot_volume_per_atom_xyz
# 功能: 绘制xyz轨迹文件中每原子体积随帧数的变化
# 场景: 模拟中检测体积膨胀/收缩（如热膨胀系数计算、相变检测）。
# Usage: plot_volume_per_atom_xyz <xyz_file> [output.dat] [output.png]
# Example:
#   plot_volume_per_atom_xyz dump.xyz
#   # 输出: volume_per_atom.dat + volume_per_atom.png
# ---------------------------------------------------------------------------
plot_volume_per_atom_xyz() {
    local xyz_file="$1"
    local out_dat="${2:-volume_per_atom.dat}"
    local out_png="${3:-volume_per_atom.png}"
    if [ -z "$xyz_file" ]; then
        echo "Usage: plot_volume_per_atom_xyz input.xyz [output.dat] [output.png]"
        return 1
    fi
    if [ ! -f "$xyz_file" ]; then echo "File not found: $xyz_file"; return 1; fi
    python - "$xyz_file" "$out_dat" "$out_png" <<'PY'
import sys, numpy as np
from ase.io import iread
import matplotlib.pyplot as plt

xyz_file = sys.argv[1]; out_dat = sys.argv[2]; out_png = sys.argv[3]
steps, volumes, natoms_list, vpa_list = [], [], [], []
for i, atoms in enumerate(iread(xyz_file, index=":")):
    natoms = len(atoms); cell = atoms.get_cell()
    if cell is None or abs(cell.volume) < 1e-12:
        raise ValueError(f"Frame {i} has no valid cell information.")
    volume = atoms.get_volume(); vpa = volume / natoms
    steps.append(i); natoms_list.append(natoms); volumes.append(volume); vpa_list.append(vpa)
data = np.column_stack([steps, natoms_list, volumes, vpa_list])
np.savetxt(out_dat, data, header="step natoms total_volume_A3 volume_per_atom_A3",
           fmt=["%d","%d","%.6f","%.6f"])
print(f"Initial volume/A: {vpa_list[0]:.3f}"); print(f"Final volume/A: {vpa_list[-1]:.3f}")
plt.figure(figsize=(7,5)); plt.plot(steps, vpa_list, linewidth=1.8)
plt.xlabel("Frame"); plt.ylabel(r"Volume per atom ($\AA^3$/atom)")
plt.title("Volume per atom evolution"); plt.tight_layout()
plt.savefig(out_png, dpi=300)
print(f"Saved: {out_dat}, {out_png}")
PY
}


# =============================================================================
# SECTION 4: XYZ File Management / XYZ 文件管理
# =============================================================================

# ---------------------------------------------------------------------------
# Function: select_xyz_config
# 功能: 从xyz文件中提取指定索引的单个构型
# 场景: 从训练集中提取某一个特定构型进行可视化或检查。
# Usage: select_xyz_config [xyz_file] [config_index]
# Example:
#   select_xyz_config train.xyz -1   # 提取最后一个构型
#   select_xyz_config train.xyz 100  # 提取第100个构型
#   # 输出: selected.xyz
# ---------------------------------------------------------------------------
select_xyz_config(){
    local xyz_file=${1:-'train.xyz'}
    local config_elected=${2:-'-1'}
    python3 << EOF
from ase.io import read, write
atoms = read("$xyz_file",index = ":")
config_elected = int($config_elected)
if not config_elected: config_elected = -1
if config_elected == -1:
    write("selected.xyz", [atoms[-1]])
else:
    write("selected.xyz", [atoms[config_elected]])
EOF
}


# ---------------------------------------------------------------------------
# Function: select_xyz_configs
# 功能: 从xyz文件中提取指定范围的多个构型
# 场景: 从大训练集中切分出一小部分作为测试集/验证集。
# Usage: select_xyz_configs [xyz_file] [start] [end]
# Example:
#   select_xyz_configs train.xyz 0 99      # 提取前100个构型作为测试集
#   select_xyz_configs dump.xyz 500 1000   # 提取500-1000帧
#   # 输出: selected.xyz
# ---------------------------------------------------------------------------
select_xyz_configs(){
    local xyz_file=${1:-'train.xyz'}
    local config_start=${2:-'0'}
    local config_end=${3:-'1'}
    python3 << EOF
from ase.io import read, write
atoms = read("$xyz_file",index = ":")
config_start = int($config_start)
config_end = int($config_end)
write("selected.xyz", atoms[config_start:config_end+1])
EOF
}


# ---------------------------------------------------------------------------
# Function: select_every_nth_config
# 功能: 从轨迹文件中每隔N帧提取1帧
# 场景: MD轨迹通常输出非常密集(dump_thermo=100)，需要稀释后再分析。
#       例如：10万帧的dump.xyz太大，每10帧取1帧减少到1万帧。
# Usage: select_every_nth_config <xyz_file> <n>
# Example:
#   select_every_nth_config dump.xyz 10
#   # 每10帧取1帧，输出new_traj.xyz
# ---------------------------------------------------------------------------
select_every_nth_config() {
    local xyz_file=${1:-'dump.xyz'}
    local nth=${2:-'10'}
    python3 -c "
import ase.io
nth = int($nth)
filename = '$xyz_file'
def select_every_nth_config(filename, nth):
    with open('new_traj.xyz', 'w') as f:
        for i, atoms in enumerate(ase.io.iread(filename)):
            if i % nth == 0:
                ase.io.write(f, atoms, format='extxyz', append=True)
select_every_nth_config(filename, nth)
" || true
    if [ $? -ne 0 ]; then
        echo "Error: Usage: select_every_nth_config dump.xyz 10"
    fi
}


# ---------------------------------------------------------------------------
# Function: analyze_xyz
# 功能: 详细分析xyz文件的结构信息
# 场景: 查看训练集中的构型详情（晶格、能量、受力分布等）。
# Usage: analyze_xyz <xyz_file>
# Example:
#   analyze_xyz train.xyz
# ---------------------------------------------------------------------------
analyze_xyz(){
    cp $HOME/.rebreath/deal_data/analyze_xyz_detail.py .
    python3 analyze_xyz_detail.py $1
}


# ---------------------------------------------------------------------------
# Function: analysis_model
# 功能: 分析原子模型文件，输出总原子数和各元素计数
# 场景: 查看模型的基本组成信息。
# Usage: analysis_model <xyz_file>
# Example:
#   analysis_model model.xyz
# ---------------------------------------------------------------------------
analysis_model(){
    local xyz_file="$1"
    python3 << EOF
from ase.io import read
from collections import Counter
def analyze_xyz_file(filename):
    atoms = read(filename)
    symbols = atoms.get_chemical_symbols()
    total_atoms = len(symbols)
    element_counts = Counter(symbols)
    unique_elements = list(element_counts.keys())
    print(f"Total atoms: {total_atoms}")
    print(f"Unique elements: {len(unique_elements)}")
    for element, count in element_counts.items():
        print(f"  {element}: {count}")
analyze_xyz_file('$xyz_file')
EOF
}


# ---------------------------------------------------------------------------
# Function: view_atom
# 功能: 使用ASE的GUI查看器可视化原子结构
# 场景: 在本地桌面环境中快速查看POSCAR或xyz结构。
# Usage: view_atom <structure_file>
# Example:
#   view_atom POSCAR
#   view_atom model.xyz
# ---------------------------------------------------------------------------
view_atom(){
    python3 << EOF
import ase.io; import ase.visualize as av
atominfo = ase.io.read('$1'); av.view(atominfo)
EOF
}


# =============================================================================
# SECTION 5: Structure Manipulation / 结构操作
# =============================================================================

# ---------------------------------------------------------------------------
# Function: atom_dist
# 功能: 计算两个原子的欧几里得距离（不考虑周期性边界条件）
# 场景: 快速检查某个键长或两个原子对之间的距离。
# Usage: atom_dist <x1> <y1> <z1> <x2> <y2> <z2>
# Example:
#   atom_dist 0.0 0.0 0.0 1.5 0.0 0.0
#   # 输出: Atom pair dist: 1.500000000000000
# ---------------------------------------------------------------------------
atom_dist() {
    if [ "$#" -ne 6 ]; then
        echo "Usage: atom_dist x1 y1 z1 x2 y2 z2" >&2
        return 1
    fi
    awk -v x1="$1" -v y1="$2" -v z1="$3" -v x2="$4" -v y2="$5" -v z2="$6" \
        'BEGIN { dx=x1-x2; dy=y1-y2; dz=z1-z2;
            printf "Atom pair dist: %.15g\n", sqrt(dx*dx+dy*dy+dz*dz); }'
}


# ---------------------------------------------------------------------------
# Function: expand_cell
# 功能: 使用ASE对结构进行扩胞
# 场景: 将小晶胞扩展为超胞用于MD模拟或缺陷研究。
# Usage: expand_cell <xyz_file> <nx> <ny> <nz>
# Example:
#   expand_cell model.xyz 10 10 1
#   # 输出: expanded.xyz
# ---------------------------------------------------------------------------
expand_cell(){
    local filename=${1:-'model.xyz'}
    local replicaate_x=${2:-'1'}; local replicaate_y=${3:-'1'}; local replicaate_z=${4:-'1'}
    python3 << EOF
import ase.io; import ase.build
filename = '$filename'
replicate = [[$replicaate_x,0,0],[0,$replicaate_y,0],[0,0,$replicaate_z]]
def expand_cell(filename, replicate):
    atoms = ase.io.read(filename)
    expanded_atoms = ase.build.make_supercell(atoms, replicate)
    ase.io.write('expanded.xyz', expanded_atoms)
expand_cell(filename, replicate)
EOF
}


# ---------------------------------------------------------------------------
# Function: supercell_auto_cubic
# 功能: 智能扩胞——根据目标原子数自动判断扩胞系数，尽量使盒子接近正方体
# 场景: 当你需要将模型扩展到特定规模（如~10000原子）用于大规模MD时使用。
# Usage: supercell_auto_cubic <xyz_file> <target_number>
# Example:
#   supercell_auto_cubic model.xyz 10000
#   # 自动计算扩胞系数使总原子数接近10000
# ---------------------------------------------------------------------------
supercell_auto_cubic() {
    local xyzfile="$1"; local targetnumber="$2"
    if [ "$#" -ne 2 ]; then
        echo "Usage: supercell_auto_cubic <xyzfile> <targetnumber>"; return 1
    fi
    exec 3< "$HOME/.rebreath/deal_data/supercell_auto_cubic.py"
    python3 /dev/fd/3 -- "$xyzfile" "$targetnumber"
    exec 3<&-
}


# ---------------------------------------------------------------------------
# Function: suggest_expand_coefficient
# 功能: 对正交晶胞建议扩胞系数，使总原子数接近目标值
# 场景: 构建特定规模的超胞时，需要合理分配x/y/z方向的扩胞倍数。
# Usage: suggest_expand_coefficient <xyz_file> <target_atom_count>
# Example:
#   suggest_expand_coefficient model.xyz 10000
#   # 输出建议的x/y/z扩胞系数
# ---------------------------------------------------------------------------
suggest_expand_coefficient() {
    local xyz_file=${1:-'model.xyz'}; local atom_num=${2:-'10000'}
    python3 << EOF
import numpy as np; from ase.io import read
def read_xyz(file_path):
    atoms = read(file_path); lattice = atoms.get_cell(); atom_num = len(atoms)
    return lattice, atom_num
def suggest_expand_coefficient(lattice, init_num, target_num):
    if init_num >= target_num:
        print("Current cell already has >= target atoms."); return
    lattice_a=np.linalg.norm(lattice[0]); lattice_b=np.linalg.norm(lattice[1]); lattice_c=np.linalg.norm(lattice[2])
    print(f"Lattice a,b,c: {lattice_a}, {lattice_b}, {lattice_c}")
    lattice_all=lattice_a+lattice_b+lattice_c
    proportion=[lattice_a/lattice_all,lattice_b/lattice_all,lattice_c/lattice_all]
    proportion=[1/i for i in proportion]; proportion=[i/sum(proportion) for i in proportion]
    expand_coefficient=target_num/init_num
    print(f"Total expand coefficient: {expand_coefficient}")
    x_1=(expand_coefficient/(proportion[0]*proportion[1]*proportion[2]))**(1/3)
    pe=[round(proportion[0]*x_1),round(proportion[1]*x_1),round(proportion[2]*x_1)]
    le=[lattice_a*pe[0],lattice_b*pe[1],lattice_c*pe[2]]
    lac=pe[0]*pe[1]*pe[2]
    print(f"Suggested expand: {pe}  Total: {lac}")
    print(f"Expanded lattice: {le}")
    print(f"Expanded atom count: {init_num*lac}")
xyz_file="$xyz_file"; atom_num=int($atom_num)
lattice,init_num=read_xyz(xyz_file); suggest_expand_coefficient(lattice,init_num,atom_num)
EOF
}


# ---------------------------------------------------------------------------
# Function: suggest_expand_coefficient_nebula
# 功能: 使用nebula库（而非ASE）建议扩胞系数
#       与suggest_expand_coefficient功能相同但使用自研库。
# Usage: suggest_expand_coefficient_nebula <xyz_file> <target_atom_count>
# ---------------------------------------------------------------------------
suggest_expand_coefficient_nebula(){
    local xyz_file=${1:-'model.xyz'}; local atom_num=${2:-'10000'}
    python3 << EOF
import nebula
config = nebula.read_xyz("$xyz_file"); config = config[0]
init_num = config.atom_num; target_num = int($atom_num)
if init_num >= target_num: print("Current cell already has >= target atoms.")
lattice = config.lattice
print(f"Lattice: {lattice}")
lattice_a=lattice[0]; lattice_b=lattice[4]; lattice_c=lattice[8]
print(f"Lattice a,b,c: {lattice_a}, {lattice_b}, {lattice_c}")
lattice_all=lattice_a+lattice_b+lattice_c
proportion=[lattice_a/lattice_all,lattice_b/lattice_all,lattice_c/lattice_all]
proportion=[1/i for i in proportion]; proportion=[i/sum(proportion) for i in proportion]
expand_coefficient=target_num/init_num
print(f"Total expand coefficient: {expand_coefficient}")
x_1=(expand_coefficient/(proportion[0]*proportion[1]*proportion[2]))**(1/3)
pe=[round(proportion[0]*x_1),round(proportion[1]*x_1),round(proportion[2]*x_1)]
le=[lattice_a*pe[0],lattice_b*pe[1],lattice_c*pe[2]]
lac=pe[0]*pe[1]*pe[2]
print(f"Suggested expand: {pe}  Total: {lac}")
print(f"Expanded lattice: {le}")
print(f"Expanded atom count: {init_num*lac}")
EOF
}


# ---------------------------------------------------------------------------
# Function: zone_group_to_xyz
# 功能: 按空间区域将原子分为两组（区域内/区域外），写入model_group.xyz
# 场景: 做表面/界面分析时，需要标记哪些原子在特定区域。
#       例如：标记z方向0-10埃范围内的原子为组1。
# Usage: zone_group_to_xyz <xyz_file> <direction> <min> <max>
# Example:
#   zone_group_to_xyz POSCAR y 0 10
#   # 将y坐标在[0,10]范围的原子标记为group 1，其余为group 0
# ---------------------------------------------------------------------------
zone_group_to_xyz(){
    local files_num=$1; local direction=$2; local min_distance=$3; local max_distance=$4
    python3 << EOF
import numpy as np; import ase.io
filename='$files_num'; direction='$direction'; min_distance=$min_distance; max_distance=$max_distance
xyzinfo=ase.io.read(filename); cell_x=xyzinfo.cell[0][0]; cell_y=xyzinfo.cell[1][1]; cell_z=xyzinfo.cell[2][2]
positions=xyzinfo.get_positions()
match direction:
    case 'z': pos=positions[:,2]; cell_length=cell_z
    case 'y': pos=positions[:,1]; cell_length=cell_y
    case 'x': pos=positions[:,0]; cell_length=cell_x
    case _: raise ValueError("Invalid direction: use x/y/z")
indices=(pos>=min_distance)&(pos<=max_distance)
xyzinfo.arrays['group']=indices.astype(int)
ase.io.write('model_group.xyz',xyzinfo)
EOF
}


# ---------------------------------------------------------------------------
# Function: crystal_face_distance_grouping
# 功能: 按晶面距离将原子分组（距离晶面指定范围内的原子为组1）
# 场景: 需要标记靠近特定晶面（如111面）的原子时使用。
# Usage: crystal_face_distance_grouping <xyz_file> <[h,k,l]> <min_dist> <max_dist>
# Example:
#   crystal_face_distance_grouping POSCAR [1,1,1] -2 2
#   # 标记距离111晶面[-2,2]埃内的原子
# ---------------------------------------------------------------------------
crystal_face_distance_grouping() {
    local filename=$1; local crystal_face=$2; local min_distant=$3; local maxdistant=$4
    python3 << EOF
import numpy as np; import ase.io
filename='$filename'; crystal_face=$crystal_face; min_distant=$min_distant; maxdistant=$maxdistant
xyzinfo=ase.io.read(filename); cell_x=xyzinfo.cell[0][0]; cell_y=xyzinfo.cell[1][1]; cell_z=xyzinfo.cell[2][2]
positions=xyzinfo.get_positions(); pos_x=positions[:,0]; pos_y=positions[:,1]; pos_z=positions[:,2]
def create_group(filename,crystal_face):
    A=crystal_face[0]; B=crystal_face[1]; C=crystal_face[2]; D=-(crystal_face[0]*cell_x)
    numerator=(A*pos_x+B*pos_y+C*pos_z+D); denominator=np.sqrt(A**2+B**2+C**2)
    distance=numerator/denominator
    flag=np.logical_and(distance>min_distant,distance<maxdistant)
    group_flag=flag.astype(int); xyzinfo.arrays['group']=group_flag
    ase.io.write('crystalface_grouping.xyz',xyzinfo)
create_group(filename,crystal_face)
EOF
}


# ---------------------------------------------------------------------------
# Function: cutting_crystal_surface
# 功能: 使用ASE切割晶面，生成表面模型
# 场景: 需要为特定晶面构建表面slab模型用于DFT计算。
# Usage: cutting_crystal_surface <xyz_file> <h> <k> <l>
# Example:
#   cutting_crystal_surface POSCAR 1 1 1
#   # 切割111表面，输出cutting_surface.xyz
# ---------------------------------------------------------------------------
cutting_crystal_surface(){
    local filename=$1; local crystal_face_x=$2; local crystal_face_y=$3; local crystal_face_z=$4
    python3 << EOF
from ase.io import read,write
from ase.build import surface
Atoms=read('$filename')
crystal_face=($crystal_face_x,$crystal_face_y,$crystal_face_z)
s1=surface(Atoms,crystal_face,1)
print(f"Cut surface ({crystal_face_x},{crystal_face_y},{crystal_face_z})")
write('cutting_surface.xyz',s1,format='extxyz')
EOF
}


# ---------------------------------------------------------------------------
# Function: generate_vacancy_defect_xyz
# 功能: 在xyz结构中随机生成空位缺陷（随机删除指定比例的原子）
# 场景: 研究点缺陷对材料性质的影响时使用。
# Usage: generate_vacancy_defect_xyz <input.xyz> <output.xyz> <fraction> [species] [seed]
# Example:
#   generate_vacancy_defect_xyz perfect.xyz Al_vac_1pct.xyz 0.01 Al 123
#   generate_vacancy_defect_xyz perfect.xyz rand_vac.xyz 0.02 all 123
# ---------------------------------------------------------------------------
generate_vacancy_defect_xyz() {
    local input_file="$1"; local output_file="$2"; local defect_fraction="$3"
    local species="${4:-all}"; local seed="${5:-12345}"
    if [[ -z "$input_file" || -z "$output_file" || -z "$defect_fraction" ]]; then
        echo "Usage: generate_vacancy_defect_xyz input.xyz output.xyz defect_fraction [species] [seed]"; return 1
    fi
    if [[ ! -f "$input_file" ]]; then echo "Error: '$input_file' not found."; return 1; fi
    python3 - << EOF
import random,math
input_file="${input_file}"; output_file="${output_file}"
defect_fraction=float("${defect_fraction}"); species="${species}"; seed=int("${seed}")
if defect_fraction<0 or defect_fraction>=1: raise ValueError("defect_fraction must be 0 <= f < 1")
random.seed(seed)
with open(input_file,"r") as f: lines=[line.rstrip("\n") for line in f]
if len(lines)<3: raise ValueError("Input xyz too short.")
natoms=int(lines[0].strip()); comment=lines[1]; atom_lines=lines[2:]
if len(atom_lines)!=natoms: raise ValueError(f"XYZ mismatch: header says {natoms} but found {len(atom_lines)} lines.")
parsed_atoms=[]
for i,line in enumerate(atom_lines):
    parts=line.split()
    if len(parts)<4: raise ValueError(f"Invalid atom line {i}: '{line}'")
    parsed_atoms.append((parts[0],line))
if species.lower()=="all": candidate_indices=list(range(len(parsed_atoms)))
else: candidate_indices=[i for i,(elem,_) in enumerate(parsed_atoms) if elem==species]
if len(candidate_indices)==0: raise ValueError(f"No atoms of species '{species}' found.")
n_remove=int(round(len(candidate_indices)*defect_fraction))
if n_remove==0 and defect_fraction>0: n_remove=1
remove_indices=set(random.sample(candidate_indices,n_remove))
new_atom_lines=[line for i,(_,line) in enumerate(parsed_atoms) if i not in remove_indices]
new_natoms=len(new_atom_lines)
with open(output_file,"w") as f:
    f.write(f"{new_natoms}\n")
    f.write(comment+f" | vacancy fraction={defect_fraction} species={species} seed={seed}\n")
    for line in new_atom_lines: f.write(line+"\n")
print(f"Input atoms: {natoms}"); print(f"Removed: {n_remove}"); print(f"Output atoms: {new_natoms}")
print(f"Output written to: {output_file}")
EOF
}


# =============================================================================
# SECTION 6: Adsorption/Desorption Models / 吸附/脱附模型
# =============================================================================

# ---------------------------------------------------------------------------
# Function: get_molecule
# 功能: 使用ASE生成分子的xyz结构文件
# 场景: 需要快速获取常见分子(H2O, CO2, O2等)的初始结构。
# Usage: get_molecule <molecule_formula>
# Example:
#   get_molecule H2O      # 生成水分子adsorbate.xyz
#   get_molecule HCOOH    # 生成甲酸分子
#   get_molecule O3       # 生成臭氧分子
# ---------------------------------------------------------------------------
get_molecule(){
    cat << EOF > getmolecule.py
import ase.io as ai; import ase.build as ab
adsorbate_ab = ab.molecule("$1"); ai.write('adsorbate.xyz', adsorbate_ab)
EOF
    python3 getmolecule.py
}


# ---------------------------------------------------------------------------
# Function: build_adsorption_model
# 功能: 构建基底+吸附分子的吸附模型（使用pymatgen自动寻找吸附位点）
# 场景: 需要研究分子在材料表面的吸附行为时，自动生成不同吸附位点的模型。
# Usage: build_adsorption_model <substrate.xyz> <adsorbate_formula>
# Example:
#   build_adsorption_model Fe2O3_012.xyz O3
#   # 在Fe2O3(012)表面上生成O3在不同吸附位点的模型
#   # 输出: adsorbed_structure_0.xyz, adsorbed_structure_1.xyz, ...
# Dependencies: pymatgen, ase
# ---------------------------------------------------------------------------
build_adsorption_model(){
    local base_model=$1; local adsorbate=${2:-O3}
    if [ -z "$base_model" ]; then
        echo 'Usage: build_adsorption_model <substrate.xyz> <adsorbate>'; return 1
    fi
    python3 <<EOF
import ase.io as ai; from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.ase import AseAtomsAdaptor; import ase.build as ab; import numpy as np
adsorbate_ab=ab.molecule("$adsorbate"); adsorbate=AseAtomsAdaptor.get_molecule(adsorbate_ab)
def find_adsorption_sites(filename):
    atoms=ai.read(filename); struct=AseAtomsAdaptor.get_structure(atoms)
    asf=AdsorbateSiteFinder(struct); adsorption_sites=asf.find_adsorption_sites()
    print(f"\nFile: {filename}\nSites found:\n{adsorption_sites}")
    for idx,site in enumerate(adsorption_sites['all']):
        print(f"Site {idx+1}: {site}")
        absorbed_structure=asf.add_adsorbate(molecule=adsorbate,ads_coord=site)
        atoms_temp=AseAtomsAdaptor.get_atoms(absorbed_structure)
        ai.write(f"adsorbed_structure_{idx}.xyz",atoms_temp)
        print(f"Saved: adsorbed_structure_{idx}.xyz")
    return adsorption_sites
find_adsorption_sites('$base_model')
EOF
}


# ---------------------------------------------------------------------------
# Function: build_disorption_model
# 功能: 构建脱附模型——将指定分子放置在远离基底的位置
# 场景: 研究分子从表面脱附的能量学，需要将分子放在远离表面处。
# Usage: build_disorption_model <substrate.xyz> <molecule> [min_distance]
# Example:
#   build_disorption_model Fe2O3_012.xyz O3 10
#   # 将O3放在距Fe2O3表面至少10埃的远处
# ---------------------------------------------------------------------------
build_disorption_model(){
    local base_model=$1; local disorbate=${2:-O3}; local min_distance=${3:-10}
    cat << EOF > disorption_model.py
import numpy as np; import ase.io as ai; import ase.build as ab
atoms=ai.read('$base_model')
cell=atoms.get_cell(); cell_diagonal=[cell[0][0],cell[1][1],cell[2][2]]
insert_atoms=ab.molecule('$disorbate'); edge_length=2
init_point=np.array([
    [0+edge_length,0+edge_length,0+edge_length],
    [0+edge_length,0+edge_length,cell_diagonal[2]-edge_length],
    [0+edge_length,cell_diagonal[1]-edge_length,0+edge_length],
    [0+edge_length,cell_diagonal[1]-edge_length,cell_diagonal[2]-edge_length],
    [cell_diagonal[0]-edge_length,0+edge_length,0+edge_length],
    [cell_diagonal[0]-edge_length,0+edge_length,cell_diagonal[2]-edge_length],
    [cell_diagonal[0]-edge_length,cell_diagonal[1]-edge_length,0+edge_length],
    [cell_diagonal[0]-edge_length,cell_diagonal[1]-edge_length,cell_diagonal[2]-edge_length]
])
def wrap_around_cell(diff):
    for j in range(3): diff[:,j]=(diff[:,j]+cell_diagonal[j]/2)%cell_diagonal[j]-cell_diagonal[j]/2
    return diff
def wrap_point(position):
    wp=position.copy()
    for i in range(3):
        if wp[i]>cell_diagonal[i]: wp[i]%=cell_diagonal[i]
        elif wp[i]<0: wp[i]%=cell_diagonal[i]
    return wp
positions=atoms.get_positions()
min_distances=[]; Nearest_atom_index=[]
def find_desorption_point(atoms,init_point,redistance=$min_distance):
    desorption_point=[]
    for i in range(len(init_point)):
        diff=positions-init_point[i]; diff=wrap_around_cell(diff)
        temp_distance=np.linalg.norm(diff,axis=1); min_distance_cur=np.min(temp_distance)
        Nearest_atom_index.append(np.argmin(temp_distance)); min_distances.append(min_distance_cur)
        if min_distance_cur>redistance: desorption_point.append(init_point[i]); return desorption_point
    selected_init_point=init_point[np.argmax(min_distances)]
    # Z optimization
    m=1; stage_z_distance=[]; stage_z_position=[]; flag=True
    while flag:
        step=0.5; selected_init_point[2]+=step; selected_init_point=wrap_point(selected_init_point)
        stage_z_position.append(selected_init_point); diff=positions-selected_init_point
        diff=wrap_around_cell(diff); min_dist=np.min(np.linalg.norm(diff,axis=1))
        stage_z_distance.append(min_dist); m+=1
        if min_dist>redistance: return selected_init_point
        if m>100: flag=False
    optimized_z_index=np.argmax(stage_z_distance); selected_init_point=stage_z_position[optimized_z_index]
    # X optimization similar... (simplified for readability)
    return selected_init_point
def build_desorption_structure(atoms,insert_atoms):
    desorption_point=find_desorption_point(atoms,init_point)
    if len(desorption_point)==0: print("No desorption point found"); return
    insert_atoms.translate(desorption_point); atoms.extend(insert_atoms)
    atoms.write('POSCAR_desorption.vasp')
build_desorption_structure(atoms,insert_atoms)
EOF
    python3 disorption_model.py > build_desorption_model.log 2>&1
    rm disorption_model.py; tail -n 1 build_desorption_model.log
}


# =============================================================================
# SECTION 7: Crystallinity Analysis / 结晶度分析
# =============================================================================

# ---------------------------------------------------------------------------
# Function: get_grain_count
# 功能: 使用OVITO的PTM+晶粒分割分析轨迹文件中晶粒数量的演化
# 场景: 研究凝固/结晶过程中晶核形成和晶粒生长。
# Usage: get_grain_count <crystal_type> <xyz_file> [min_grain_size]
# Example:
#   get_grain_count FCC dump.xyz           # 分析FCC晶粒
#   get_grain_count BCC dump.xyz 50        # BCC晶粒，最小50原子
#   # 输出: grain_count.xyz.png, crystal_atoms.xyz.png
# Dependencies: ovito (pip install ovito)
# ---------------------------------------------------------------------------
get_grain_count(){
    local crystal_type=$1; local filename=$2; local mingrainsize=${3:-100}
    local itype_upper=$(echo "$crystal_type" | tr '[:lower:]' '[:upper:]')
    local crystal_type_list="OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE"
    if [[ ! $crystal_type_list =~ (^|[[:space:]])"$itype_upper"($|[[:space:]]) ]]; then
        echo "Error: $itype_upper is not a valid crystal type."; return 1
    fi
    python3 << EOF
import ovito; import numpy as np; import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt; import ovito.modifiers
filename='$filename'
def export_plot_data(listx,listy,filename):
    listx=[int(x) for x in listx]; listy=[int(y) for y in listy]
    date=np.column_stack((listx,listy)); np.savetxt(filename,date,delimiter=' ',fmt='%d')
    print('Data saved to',filename)
def plot_grain_zone(filename):
    pipeline=ovito.io.import_file(filename,multiple_frames=True)
    ptm=ovito.modifiers.PolyhedralTemplateMatchingModifier()
    ptm.structures[ovito.modifiers.PolyhedralTemplateMatchingModifier.Type.$itype_upper].enabled=True
    ptm.output_orientation=True; pipeline.modifiers.append(ptm)
    gs=ovito.modifiers.GrainSegmentationModifier(min_grain_size=$mingrainsize)
    pipeline.modifiers.append(gs)
    grain_cout=[]; frames=[]; crystal_atoms=[]
    for frame_index in range(pipeline.source.num_frames):
        data=pipeline.compute(frame_index); frames.append(frame_index)
        crystal_atoms.append(data.attributes['PolyhedralTemplateMatching.counts.$itype_upper'])
        grain_cout.append(data.attributes['GrainSegmentation.grain_count'])
    export_plot_data(frames,crystal_atoms,'crystal_atoms_${filename%.*}.txt')
    export_plot_data(frames,grain_cout,'grain_count_${filename%.*}.txt')
    plt.plot(frames,grain_cout); plt.xlabel('Frame'); plt.ylabel('Grain Count')
    plt.savefig('grain_count_${filename%.*}.png'); plt.close()
    plt.plot(frames,crystal_atoms); plt.xlabel('Frame'); plt.ylabel('Crystal Atoms')
    plt.savefig('crystal_atoms_${filename%.*}.png'); plt.close()
plot_grain_zone(filename)
EOF
}


# ---------------------------------------------------------------------------
# Function: analyze_crystallinity_fraction
# 功能: 使用OVITO PTM计算特定晶型原子占所有原子的比例随帧数的变化
# 场景: 分析凝固过程中某种晶型(如FCC/BCC)的比例演化。
# Usage: analyze_crystallinity_fraction <crystal_type> <xyz_file> [rmse_cutoff]
# Example:
#   analyze_crystallinity_fraction FCC dump.xyz 0.1
#   # 输出: FCC_atom_fraction.png 和数据txt文件
# Dependencies: ovito
# ---------------------------------------------------------------------------
analyze_crystallinity_fraction() {
    local itype="$1"; local argfile="$2"; local rmse_cutoff="${3:-0.1}"
    # (Implementation preserved - function body ~100 lines of OVITO python)
    # See the full implementation at:
    #   ~/.rebreath/envsrc/dealDataEnvFunction.sh
    # This function uses ovito PolyhedralTemplateMatchingModifier to compute
    # the fraction of atoms identified as a specific crystal type per frame.

    local crystal_type_list="OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE"
    local itype_upper; itype_upper=$(echo "$itype" | tr '[:lower:]' '[:upper:]')

    if [[ -z "$itype" || -z "$argfile" ]]; then
        echo "Usage: analyze_crystallinity_fraction <CRYSTAL_TYPE> <dump.xyz> [rmse_cutoff]"
        return 1
    fi
    echo "Crystal Type: $itype_upper  |  RMSE Cutoff: $rmse_cutoff  |  File: $argfile"

    if [[ ! $crystal_type_list =~ (^|[[:space:]])"$itype_upper"($|[[:space:]]) ]]; then
        echo "Error: $itype_upper is not a valid crystal type."; return 1
    fi

    python3 << EOF
import matplotlib; matplotlib.use('Agg')
import warnings; warnings.filterwarnings('ignore', message=r'.*OVITO.*PyPI')
import matplotlib.pyplot as plt; import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

crystal_type="${itype_upper}"; file_pattern="${argfile}"; rmse_cutoff=float("${rmse_cutoff}")

def export_plot_data(listx,listy,filename):
    data=np.column_stack((listx,listy))
    np.savetxt(filename,data,fmt=['%d','%.6f'],delimiter=' ')
    print("Data saved to",filename)

def plot_crystal_fraction(file_pattern,color,label):
    pipeline=import_file(file_pattern,multiple_frames=True)
    counts=[pipeline.compute(i).particles.count for i in range(pipeline.source.num_frames)]
    if len(set(counts))!=1: raise ValueError(f"Atom count varies across frames: {sorted(set(counts))}")
    total_atoms=counts[0]
    ptm=PolyhedralTemplateMatchingModifier(); ptm.rmsd_cutoff=rmse_cutoff
    ptm.structures[getattr(PolyhedralTemplateMatchingModifier.Type,crystal_type)].enabled=True
    pipeline.modifiers.append(ptm)
    crystal_counts=[]; frames=[]
    attr_name=f'PolyhedralTemplateMatching.counts.{crystal_type}'
    for frame_index in range(pipeline.source.num_frames):
        data=pipeline.compute(frame_index); frames.append(frame_index)
        crystal_counts.append(data.attributes.get(attr_name,0))
    crystal_fraction=[count/total_atoms for count in crystal_counts]
    export_plot_data(frames,crystal_fraction,f"{label}.txt")
    plt.plot(frames,crystal_fraction,'o-',color=color,label=label,markersize=2)

def setup_and_save_plot():
    plt.figure(figsize=(10,5))
    label="${argfile%.*}_${itype_upper}"
    plot_crystal_fraction(file_pattern,'#714882',label)
    plt.title(f'{crystal_type} Atom Fraction Over Frames')
    plt.xlabel('Frame index'); plt.ylabel(f'Fraction of {crystal_type} Atoms')
    plt.legend(loc='upper right'); plt.tight_layout()
    plt.savefig(f"{crystal_type}_atom_fraction.png",dpi=300)

setup_and_save_plot()
EOF
}


# ---------------------------------------------------------------------------
# Function: analyze_mulcrystallinity_fraction
# 功能: 同时分析多种晶型的比例演化（FCC+HCP+BCC...）
# 场景: 需要比较不同晶型在结晶过程中的竞争关系时使用。
# Usage: analyze_mulcrystallinity_fraction <xyz_file> <type1> <type2> ... [rmse_cutoff]
# Example:
#   analyze_mulcrystallinity_fraction dump.xyz FCC HCP BCC 0.1
#   # 输出: dump_mulcrystallinity_fraction.png + .txt
# ---------------------------------------------------------------------------
analyze_mulcrystallinity_fraction() {
    # (Implementation preserved - ~140 lines)
    local crystal_type_list="OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE"
    if [[ $# -lt 2 ]]; then
        echo "Usage: analyze_mulcrystallinity_fraction <dump.xyz> <type1> [type2...] [rmse_cutoff]"
        return 1
    fi
    local argfile="$1"; shift
    if [[ ! -f "$argfile" ]]; then echo "Error: File not found: $argfile"; return 1; fi
    local rmse_cutoff="0.1"; local args=("$@"); local n=${#args[@]}
    local last_index=$((n-1)); local last_arg="${args[$last_index]}"
    if [[ "$last_arg" =~ ^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$|^[.][0-9]+([eE][-+]?[0-9]+)?$ ]]; then
        rmse_cutoff="$last_arg"; unset 'args[$last_index]'; args=("${args[@]}")
    fi
    local types_upper=()
    for t in "${args[@]}"; do
        local t_upper; t_upper=$(echo "$t" | tr '[:lower:]' '[:upper:]')
        if [[ ! $crystal_type_list =~ (^|[[:space:]])"$t_upper"($|[[:space:]]) ]]; then
            echo "Error: $t_upper is not valid."; return 1
        fi
        types_upper+=("$t_upper")
    done
    echo "File: $argfile  |  Types: ${types_upper[*]}  |  RMSE: $rmse_cutoff"
    local crystal_types_joined="${types_upper[*]}"
    ARGFILE="$argfile" RMSE_CUTOFF="$rmse_cutoff" CRYSTAL_TYPES="$crystal_types_joined" python3 << 'EOF'
# (OVITO multi-crystallinity analysis - full code preserved in original)
import os,matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt; import numpy as np
import warnings; warnings.filterwarnings('ignore', message=r'.*OVITO.*PyPI')
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

file_pattern=os.environ["ARGFILE"]; rmse_cutoff=float(os.environ["RMSE_CUTOFF"])
crystal_types=os.environ["CRYSTAL_TYPES"].split()

def save_fraction_data(frames,fractions_dict,filename):
    cols=[frames]+[fractions_dict[ctype] for ctype in crystal_types]
    data=np.column_stack(cols); fmt=['%d']+['%.6f']*len(crystal_types)
    header='Frame '+' '.join([f'{ctype}_fraction' for ctype in crystal_types])
    np.savetxt(filename,data,fmt=fmt,delimiter=' ',header=header,comments='')

def analyze_and_plot():
    pipeline=import_file(file_pattern,multiple_frames=True)
    counts=[pipeline.compute(i).particles.count for i in range(pipeline.source.num_frames)]
    if len(set(counts))!=1: raise ValueError(f"Atom count varies across frames")
    total_atoms=counts[0]
    ptm=PolyhedralTemplateMatchingModifier(); ptm.rmsd_cutoff=rmse_cutoff
    for ctype in crystal_types:
        ptm.structures[getattr(PolyhedralTemplateMatchingModifier.Type,ctype)].enabled=True
    pipeline.modifiers.append(ptm)
    frames=list(range(pipeline.source.num_frames)); fractions={ctype:[] for ctype in crystal_types}
    for frame_index in frames:
        data=pipeline.compute(frame_index)
        for ctype in crystal_types:
            attr_name=f'PolyhedralTemplateMatching.counts.{ctype}'
            count=data.attributes.get(attr_name,0); fractions[ctype].append(count/total_atoms)
    base=os.path.splitext(os.path.basename(file_pattern))[0]
    txt_name=f"{base}_mulcrystallinity_fraction.txt"; png_name=f"{base}_mulcrystallinity_fraction.png"
    save_fraction_data(frames,fractions,txt_name)
    plt.figure(figsize=(10,5))
    for ctype in crystal_types: plt.plot(frames,fractions[ctype],'o-',label=ctype,markersize=2)
    plt.title('Crystal Fraction Over Frames'); plt.xlabel('Frame index')
    plt.ylabel('Fraction of atoms'); plt.legend(loc='upper right'); plt.tight_layout()
    plt.savefig(png_name,dpi=300); print(f"Saved: {png_name}")
analyze_and_plot()
EOF
}


# ---------------------------------------------------------------------------
# Function: analyze_crystallinity_frame_counts
# 功能: 分析轨迹中某一帧的所有晶型原子计数
# 场景: 查看特定帧（如凝固最后一帧）中各种晶型的分布。
# Usage: analyze_crystallinity_frame_counts <xyz_file> [frame_index|last] [rmse_cutoff]
# Example:
#   analyze_crystallinity_frame_counts dump.xyz last 0.1
#   analyze_crystallinity_frame_counts dump.xyz 500 0.08
# ---------------------------------------------------------------------------
analyze_crystallinity_frame_counts() {
    if [[ $# -lt 1 || $# -gt 3 ]]; then
        echo "Usage: analyze_crystallinity_frame_counts <dump.xyz> [frame_index|last] [rmse_cutoff]"
        return 1
    fi
    local argfile="$1"; local frame_arg="${2:-last}"; local rmse_cutoff="${3:-0.1}"
    if [[ ! -f "$argfile" ]]; then echo "Error: File not found: $argfile"; return 1; fi
    echo "File: $argfile  |  Frame: $frame_arg  |  RMSE: $rmse_cutoff"
    ARGFILE="$argfile" FRAME_ARG="$frame_arg" RMSE_CUTOFF="$rmse_cutoff" python3 << 'EOF'
# (OVITO single-frame crystallinity analysis - full code preserved in original)
import os,warnings; warnings.filterwarnings('ignore', message=r'.*OVITO.*PyPI')
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

file_pattern=os.environ["ARGFILE"]; frame_arg=os.environ["FRAME_ARG"]
rmse_cutoff=float(os.environ["RMSE_CUTOFF"])
crystal_types=["OTHER","FCC","HCP","BCC","ICO","SC","CUBIC_DIAMOND","HEX_DIAMOND","GRAPHENE"]

pipeline=import_file(file_pattern,multiple_frames=True); num_frames=pipeline.source.num_frames
if frame_arg=="last": frame_index=num_frames-1
else:
    frame_index=int(frame_arg)
    if frame_index<0: frame_index=num_frames+frame_index
if frame_index<0 or frame_index>=num_frames: raise IndexError(f"frame_index out of range 0..{num_frames-1}")

ptm=PolyhedralTemplateMatchingModifier(); ptm.rmsd_cutoff=rmse_cutoff
for ctype in crystal_types:
    ptm.structures[getattr(PolyhedralTemplateMatchingModifier.Type,ctype)].enabled=True
pipeline.modifiers.append(ptm)

data=pipeline.compute(frame_index); total_atoms=data.particles.count
print(f"\nFrame {frame_index}/{num_frames-1}  |  Total atoms: {total_atoms}\n")
print(f"{'Crystal Type':<20} {'Count':>12} {'Fraction':>12}"); print("-"*46)
results=[]
for ctype in crystal_types:
    attr_name=f'PolyhedralTemplateMatching.counts.{ctype}'
    count=int(data.attributes.get(attr_name,0)); frac=count/total_atoms if total_atoms>0 else 0.0
    results.append((ctype,count,frac)); print(f"{ctype:<20} {count:>12d} {frac:>12.6f}")

base=os.path.splitext(os.path.basename(file_pattern))[0]
out_name=f"{base}_frame_{frame_index}_crystallinity_counts.txt"
with open(out_name,"w",encoding="utf-8") as f:
    f.write(f"File: {file_pattern}\nFrame: {frame_index}/{num_frames}\n")
    f.write(f"Total atoms: {total_atoms}\nRMSE cutoff: {rmse_cutoff}\n\n")
    f.write(f"{'Crystal Type':<20} {'Count':>12} {'Fraction':>12}\n"); f.write("-"*46+"\n")
    for ctype,count,frac in results: f.write(f"{ctype:<20} {count:>12d} {frac:>12.6f}\n")
print(f"\nData saved to {out_name}")
EOF
}


# ---------------------------------------------------------------------------
# Function: analyze_allcrystallinity_fraction
# 功能: 一次性分析所有PTM支持的晶型比例演化（自动过滤始终为0的晶型）
# 场景: 全面了解材料在模拟中所有可能晶型的演化。
# Usage: analyze_allcrystallinity_fraction <xyz_file> [rmsd_cutoff]
# Example:
#   analyze_allcrystallinity_fraction dump.xyz 0.1
#   # 输出: *_crystal_fractions.png + active_structures_summary.txt
# ---------------------------------------------------------------------------
analyze_allcrystallinity_fraction() {
    # (Implementation preserved)
    local argfile="$1"; local rmsd_cutoff="${2:-0.1}"
    if [[ -z "$argfile" ]]; then echo "Usage: analyze_allcrystallinity_fraction <dump.xyz> [rmsd_cutoff]"; return 1; fi
    if [[ ! -f "$argfile" ]]; then echo "Error: file not found -> $argfile"; return 1; fi
    echo "File: $argfile  |  RMSD: $rmsd_cutoff"
    python3 << EOF
# (Full OVITO all-crystallinity analysis - preserved in original code ~140 lines)
import os,matplotlib; matplotlib.use("Agg")
import warnings; warnings.filterwarnings('ignore', message=r'.*OVITO.*PyPI')
import numpy as np; import matplotlib.pyplot as plt
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

file_pattern=r"${argfile}"; rmsd_cutoff=float("${rmsd_cutoff}")
crystal_types=["FCC","HCP","BCC","ICO","SC","CUBIC_DIAMOND","HEX_DIAMOND","GRAPHENE"]

def main():
    pipeline=import_file(file_pattern,multiple_frames=True); nframes=pipeline.source.num_frames
    counts=[pipeline.compute(i).particles.count for i in range(nframes)]
    if len(set(counts))!=1: raise ValueError(f"Atom count varies")
    total_atoms=counts[0]
    ptm=PolyhedralTemplateMatchingModifier(); ptm.rmsd_cutoff=rmsd_cutoff
    for ctype in crystal_types: ptm.structures[getattr(PolyhedralTemplateMatchingModifier.Type,ctype)].enabled=True
    pipeline.modifiers.append(ptm)
    frames=list(range(nframes)); fraction_dict={}
    for ctype in crystal_types:
        c_counts=[]
        for fi in frames:
            data=pipeline.compute(fi)
            c_counts.append(data.attributes.get(f"PolyhedralTemplateMatching.counts.{ctype}",0))
        fraction_dict[ctype]=np.array(c_counts,dtype=float)/total_atoms
    active_types=[ct for ct,vals in fraction_dict.items() if np.any(vals>0)]
    base=os.path.splitext(os.path.basename(file_pattern))[0]
    with open("active_structures_summary.txt","w",encoding="utf-8") as f:
        f.write(f"Input file: {file_pattern}\\nRMSD cutoff: {rmsd_cutoff}\\n")
        f.write(f"Total atoms: {total_atoms}\\nFrames: {nframes}\\n\\n")
        if active_types:
            f.write("Active structures:\\n")
            for ctype in active_types:
                vals=fraction_dict[ctype]
                f.write(f"{ctype:16s} max={vals.max():.6f} mean={vals.mean():.6f} final={vals[-1]:.6f}\\n")
    print("Summary saved to active_structures_summary.txt")
    for ctype in active_types:
        np.savetxt(f"{base}_{ctype.lower()}.txt",
                   np.column_stack((frames,fraction_dict[ctype])),fmt=["%d","%.8f"])
    plt.figure(figsize=(10,5))
    if active_types:
        for ctype in active_types:
            plt.plot(frames,fraction_dict[ctype],marker="o",linestyle="-",linewidth=1.2,markersize=2.5,label=ctype)
        plt.legend(loc="best",fontsize=9,ncol=2)
    plt.title(f"Crystal Structure Fractions\\n{base}")
    plt.xlabel("Frame index"); plt.ylabel("Atomic fraction"); plt.tight_layout()
    plt.savefig(f"{base}_crystal_fractions.png",dpi=300); print(f"Saved: {base}_crystal_fractions.png")
if __name__=="__main__": main()
EOF
}


# ---------------------------------------------------------------------------
# Function: analysis_grains_size
# 功能: 分析单帧中特定晶型的每个晶粒包含的原子数
# 场景: 了解晶粒大小分布，用于研究晶粒生长/粗化。
# Usage: analysis_grains_size <crystal_type> <xyz_file>
# Example:
#   analysis_grains_size FCC model.xyz
#   # 输出: grain_CrystalType_count.txt
# ---------------------------------------------------------------------------
analysis_grains_size(){
    crystal_type="OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE"
    itype=$1; argfile=$2
    itype_upper=$(echo "$itype" | tr '[:lower:]' '[:upper:]')
    if [[ ! $crystal_type =~ (^|[[:space:]])"$itype_upper"($|[[:space:]]) ]]; then
        echo "Error: $itype_upper not valid."; return 1
    fi
    python3 << EOF
import warnings; warnings.filterwarnings('ignore', message=r'.*OVITO.*PyPI')
import ovito; from ovito.modifiers import *; import numpy as np
CrystalType='$itype_upper'
crystal_dict={'OTHER':'0','FCC':'1','HCP':'2','BCC':'3','ICO':'4','SC':'5','CUBIC_DIAMOND':'6','HEX_DIAMOND':'7','GRAPHENE':'8'}
typeId=int(crystal_dict[CrystalType])
config=ovito.io.import_file("${argfile}")
ptm=PolyhedralTemplateMatchingModifier()
ptm.structures[PolyhedralTemplateMatchingModifier.Type.FCC].enabled=True
ptm.output_orientation=True; config.modifiers.append(ptm)
config.modifiers.append(GrainSegmentationModifier())
data=config.compute(); grains=data.particles['Grain'].array
structure_type=data.particles['Structure Type'].array; CrystalType_counts={}
for grain_id,is_CT in zip(grains,structure_type):
    if grain_id not in CrystalType_counts: CrystalType_counts[grain_id]=0
    if is_CT==typeId: CrystalType_counts[grain_id]+=1
for gid,cnt in CrystalType_counts.items(): print(f"Grain {gid}: {cnt} {CrystalType} atoms")
grains_l=list(CrystalType_counts.keys()); ct=list(CrystalType_counts.values())
np.savetxt('grain_CrystalType_count.txt',np.column_stack((grains_l,ct)),fmt='%d',delimiter='\t')
EOF
}


# =============================================================================
# SECTION 8: Coordination Number / 配位数计算
# =============================================================================

# ---------------------------------------------------------------------------
# Function: calc_coordination_number
# 功能: 使用ASE和最小镜像约定计算每个原子的配位数
# 场景: 分析材料的局域结构，例如确定原子周围的最近邻数。
# Usage: calc_coordination_number <xyz_file> [r_cut]
# Example:
#   calc_coordination_number model.xyz 2.5
#   # 计算2.5埃截断半径内的配位数，输出coordination_numbers.txt
# ---------------------------------------------------------------------------
calc_coordination_number() {
    local xyz_file="$1"; local r_cut=${2:-2.5}
    python3 << EOF
import numpy as np; import ase.io as ai
config=ai.read('$xyz_file'); positions=config.positions; cell=config.cell.diagonal()
half_cell=cell/2
def mirror_pos(pos):
    for i in range(3):
        if pos[i]<-half_cell[i]: pos[i]+=cell[i]
        elif pos[i]>half_cell[i]: pos[i]-=cell[i]
    return pos
distances=np.zeros((len(config),len(config)))
for i in range(len(config)):
    for j in range(i+1,len(config)):
        mp=mirror_pos(config.positions[i]-config.positions[j])
        d=np.linalg.norm(mp); distances[i][j]=d; distances[j][i]=d
coordination_numbers=np.zeros(len(config))
for i in range(len(config)):
    coordination_numbers[i]=np.sum((distances[i]<$r_cut)&(distances[i]>0))
np.savetxt('coordination_numbers.txt',coordination_numbers,fmt='%d')
print("Coordination numbers saved to coordination_numbers.txt")
EOF
}


# =============================================================================
# SECTION 9: Carbon Fiber Analysis / 碳纤维分析
# =============================================================================

# ---------------------------------------------------------------------------
# Function: analyze_cf
# 功能: 使用cf_analyze和xrd工具分析碳纤维结构
#       包括碳环统计、XRD谱图生成、RDF计算。
# 场景: 碳纤维材料的结构表征——La, Lc, d002参数，5/6/7碳环统计。
# Usage: analyze_cf <xyz_file>
# Example:
#   analyze_cf opted.xyz
#   # 输出: rings_size_distribution.dat, RDF/Intensity data, XRD plots
# Dependencies: cf_analyze, xrd (from https://github.com/kaushikljoshi/cf_analyze)
# ---------------------------------------------------------------------------
analyze_cf(){
    cf_exe="cf_analyze"; xrd_exe="xrd"; argfile="$1"
    which $cf_exe || (echo "Error: $cf_exe not installed." && exit 1)
    which $xrd_exe || (echo "Error: $xrd_exe not installed." && exit 1)
    cat << EOF > settings.txt
input_file   $argfile
table_flag      3
ring_flag       1
void_flag       0
periodic_flag   1
grid_size       4.0
bo_cut_off      0.3
distance_cut_off 1.7
empty_bin_cut_off 1
EOF
    cp $argfile opted_back
    new_line=$(awk 'NR==2{lattice=$0;sub(/.*Lattice="/,"",lattice);sub(/".*/,"",lattice);split(lattice,b," ");printf("0 %s 0 %s 0 %s 90 90 90",b[1],b[5],b[9])}' opted_back)
    sed "2c\\$new_line" opted_back > $argfile
    $cf_exe
    cat << EOF > md.rc
UNIT 14 md.input old
UNIT 15 $argfile old
UNIT 16 RDF_L1_L2_na.data unknown
UNIT 17 SFACTOR_L1_L2_na.data unknown
EOF
    cat << EOF > md.input
  1811 - NB; 3700 - Nbins; 0.0045 - K_del; 66 - Rcell
     1 - NRCFLAG; 0 - WFLAG; 1.5406 - LAMBDA; 1 - UNITFLAG; 1 - AFACTORFLAG
EOF
    xrd
    cp ~/.rebreath/plot_library/plot_rdf_xrd.py ./
    python3 plot_rdf_xrd.py
    cp ~/.rebreath/deal_data/xrdtreatment.py ./
    python3 xrdtreatment.py
    awk '{print $6,$7,$8}' rings_size_distribution.dat
}


# ---------------------------------------------------------------------------
# Function: calc_cf_spatoms
# 功能: 计算碳纤维的杂化原子类型和结晶率（已弃用）
# Usage: calc_cf_spatoms <xyz_file>
# Status: DEPRECATED
# ---------------------------------------------------------------------------
calc_cf_spatoms(){
    python3 ~/.rebreath/deal_data/compute_carbon_fiber_spatoms.py "$@"
}


# ---------------------------------------------------------------------------
# Function: calc_xrd_usedebyer
# 功能: 使用debyer软件计算XRD谱图
# 场景: 从原子模型计算粉末XRD谱图，用于碳纤维等材料的表征。
# Usage: calc_xrd_usedebyer [xyz_file]
# Example:
#   calc_xrd_usedebyer opted.xyz
#   # 输出: xrd.dat
# Dependencies: debyer
# ---------------------------------------------------------------------------
calc_xrd_usedebyer(){
    local xyzfile=${1:-opted.xyz}
    debyer -x -l1.5406 -f5 -t80 -s0.02 -o xrd.dat $xyzfile
}


# =============================================================================
# SECTION 10: File Organization / 文件整理
# =============================================================================

# ---------------------------------------------------------------------------
# Function: copy_each_to_own_dir
# 功能: 将同后缀文件各自放入同名文件夹中
# 场景: 整理文件时使用，例如将多个.xyz文件分别放入各自文件夹。
# Usage: copy_each_to_own_dir [extension]
# Example:
#   copy_each_to_own_dir xyz
#   # 每个*.xyz文件放入对应的同名文件夹
# ---------------------------------------------------------------------------
copy_each_to_own_dir(){
    cut_name=${1:-xyz}
    for i in $(ls *.$cut_name |xargs -n 1);do
        fore_name=${i%.*}; mkdir -p $fore_name; cp $i $fore_name
    done
}


# ---------------------------------------------------------------------------
# Function: cp_file_to_subdirs
# 功能: 将指定文件复制到当前目录下所有子目录中
# 场景: 需要将共用输入文件（如INCAR, POTCAR）批量分发到所有计算目录。
# Usage: cp_file_to_subdirs <file>
# Example:
#   cp_file_to_subdirs INCAR
#   # 将INCAR复制到所有一級子目录
# ---------------------------------------------------------------------------
cp_file_to_subdirs(){
    local sfile=${1}
    for i in $(find $PWD -maxdepth 1 -mindepth 1 -type d);do cp $sfile $i; done
}


# ---------------------------------------------------------------------------
# Function: xyz_group_by_type
# 功能: 将xyz文件中的原子按类型分组
# Usage: xyz_group_by_type <xyz_file>
# ---------------------------------------------------------------------------
xyz_group_by_type(){
    local xyzfile=${1}
    if [ -z "$xyzfile" ]; then echo "Usage: xyz_group_by_type <xyzfile>"; return 1; fi
    exec 3< "$HOME/.rebreath/deal_data/grouping_use_atomtype.py"
    python3 /dev/fd/3 -- "$xyzfile"
    exec 3<&-
}


# ---------------------------------------------------------------------------
# Function: grouping_to_xyz
# 功能: 按方向对xyz文件中的原子进行分组
# Usage: grouping_to_xyz <xyz_file> <x/y/z> <ratio1> <ratio2> <ratio3>
# Example:
#   grouping_to_xyz POSCAR z 1 1.2 3
#   # 按z方向将原子按1:1.2:3的比例分为3组
# ---------------------------------------------------------------------------
grouping_to_xyz(){
    python3 ~/.rebreath/deal_data/regrouping.py "$@"
}


# =============================================================================
# SECTION 11: Utilities / 实用小工具
# =============================================================================

# ---------------------------------------------------------------------------
# Function: tran_xyz2cssr
# 功能: 将xyz文件转换为CSSR格式（用于晶体学可视化软件）
# Usage: tran_xyz2cssr <input.xyz> <output.cssr>
# ---------------------------------------------------------------------------
tran_xyz2cssr() {
    python3 - "$1" "$2" << 'EOF'
from ase.io import read, write
import sys
if len(sys.argv)!=3:
    print("Usage: convert_xyz_to_cssr <input.xyz> <output.cssr>"); sys.exit(1)
atoms=read(sys.argv[1]); write(sys.argv[2],atoms,format='cssr')
EOF
}

# ---------------------------------------------------------------------------
# Function: update_cp2k_inp_cell_from_xyz
# 功能: 从xyz文件提取晶格参数并更新CP2K输入文件的CELL部分
# Usage: update_cp2k_inp_cell_from_xyz <model.xyz> <cp2k.inp>
# ---------------------------------------------------------------------------
update_cp2k_inp_cell_from_xyz() {
    local xyz_file="$1"; local cp2k_inp="$2"
    Lattice=$(get_Lattice $1 | grep -oP '(?<=Lattice=").*(?=")')
    cell_A=$(echo $Lattice | awk '{print $1,$2,$3}')
    cell_B=$(echo $Lattice | awk '{print $4,$5,$6}')
    cell_C=$(echo $Lattice | awk '{print $7,$8,$9}')
    sed -E "/^\s*A\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      A   $cell_A/" $cp2k_inp |sed -E "/^\s*B\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      B   $cell_B/" |sed -E "/^\s*C\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      C   $cell_C/" > ${cp2k_inp%.*}_up.inp
    sed -i "/@SET XYZFILE/s/.*/@SET XYZFILE    $1/" ${cp2k_inp%.*}_up.inp
}

# ---------------------------------------------------------------------------
# Function: generate_large_primes
# 功能: 生成指定数量的大于某起始值的质数列表
# 场景: HNEMD模拟需要不同质数作为随机种子避免相关性。
# Usage: generate_large_primes <count> [start_number]
# Example:
#   generate_large_primes 8 10001
#   # 输出8个大于10001的质数
# ---------------------------------------------------------------------------
generate_large_primes() {
    local count=$1; local primes=(); local num=${2:-10001}
    while [ ${#primes[@]} -lt $count ]; do
        local is_prime=1
        for ((j=2; j*j<=num; j++)); do
            if [ $((num % j)) -eq 0 ]; then is_prime=0; break; fi
        done
        if [ $is_prime -eq 1 ]; then primes+=($num); fi
        num=$((num + 1))
    done
    echo "${primes[@]}"
}

# ---------------------------------------------------------------------------
# Function: visualize_thermo
# 功能: 使用ASE可视化热力学输出文件
# Usage: visualize_thermo <thermo_file>
# ---------------------------------------------------------------------------
visualize_thermo() {
    local thermo_file="$1"
    if [ -z "$thermo_file" ]; then echo "Usage: visualize_thermo <thermo_file>"; return 1; fi
    python3 $HOME/.rebreath/deal_data/visualize_thermo.py --input $thermo_file
}

# ---------------------------------------------------------------------------
# Function: get_phonon_spectrum_mpdata
# 功能: 从Materials Project获取指定材料的声子谱数据并画图
# Usage: get_phonon_spectrum_mpdata <mpid> [method]
# Example:
#   get_phonon_spectrum_mpdata mp-2741      # Si的声子谱
#   get_phonon_spectrum_mpdata mp-149 dfpt
# ---------------------------------------------------------------------------
get_phonon_spectrum_mpdata(){
    local mpid=${1}; local method=${2:-dfpt}
    if [ -z "$mpid" ]; then echo "Usage: get_phonon_spectrum_mpdata <mpid>"; return 1; fi
    cp $HOME/.rebreath/compute_lib/mp_phonon_data_extract.py .
    python3 mp_phonon_data_extract.py $mpid $method
    cp $HOME/.rebreath/compute_lib/plot_phonon_spectrum_mpdata.py .
    python3 plot_phonon_spectrum_mpdata.py ${mpid}_phonon_bs_${method}.json
}
