# 该分库使用来存放处理数据的函数，包括对于xyz类文件的处理等，较为通用的函数

atom_dist() {
# 该函数使用来快速计算两个原子对之间的距离
# 使用方法：atom_dist x1 y1 z1 x2 y2 z2 输出为两个原子对的距离，不考虑周期性边界条件
    if [ "$#" -ne 6 ]; then
        echo "Usage: atom_dist x1 y1 z1 x2 y2 z2" >&2
        return 1
    fi

    # 用 awk 做浮点计算，避免 bash 整数截断
    awk -v x1="$1" -v y1="$2" -v z1="$3" \
        -v x2="$4" -v y2="$5" -v z2="$6" \
        'BEGIN {
            dx = x1 - x2;
            dy = y1 - y2;
            dz = z1 - z2;
            printf "Atom pair dist : %.15g\n", sqrt(dx*dx + dy*dy + dz*dz);
        }'
}

view_atom(){
# 该函数使用来快速查看原子位置
# 该函数可以快速查看原子结构，使用方法为 view_atom POSCAR（或者类似的文件）
    python3 << EOF
import ase.io
import ase.visualize as av
atominfo = ase.io.read('$1')
av.view(atominfo)
EOF
}


zone_group_to_xyz(){
# idea from bessel
# 本脚本的分组方式为按照指定的区域来进行分组，可以将模型文件按照指定的区域分为区域内与区域外两组
# 使用方法为 zone_group_to_xyz POSCAR y 0 10 将会按照y方向的距离来进行分组，距离在0-10之间为组1,其他的为组0
    local files_num=$1
    local direction=$2
    local min_distance=$3
    local max_distance=$4
    python3 << EOF
import numpy as np
import ase.io
filename = '$files_num'
direction = '$direction'  
min_distance = $min_distance  
max_distance = $max_distance
xyzinfo = ase.io.read(filename)
cell_x = xyzinfo.cell[0][0]
cell_y = xyzinfo.cell[1][1]
cell_z = xyzinfo.cell[2][2]
positions = xyzinfo.get_positions()
match direction:
    case 'z':
        pos = positions[:, 2]
        cell_length = cell_z
    case 'y':
        pos = positions[:, 1]
        cell_length = cell_y
    case 'x':
        pos = positions[:, 0]
        cell_length = cell_x
    case _:
        print("Invalid direction")
        raise ValueError("Invalid direction specified.")
indices = (pos >= min_distance)&(pos <= max_distance)
xyzinfo.arrays['group'] = indices.astype(int)
ase.io.write('model_group.xyz', xyzinfo)
EOF
}

crystal_face_distance_grouping() {
# 本脚本的分组方式为按照晶面与晶面的距离来进行分组，可以将模型文件按照晶面与晶面的距离分为距离内与距离外两组
# 使用方法为 crystal_face_distance_grouping POSCAR [1,1,1] -2 2 将会按照原子距离晶面111面之间的距离进行分组，将POSCAR距离晶面111在-2，2之间为组1,其他的为组0
local filename=$1
local crystal_face=$2
local min_distant=$3
local maxdistant=$4
python3 << EOF
import numpy as np
import ase.io
filename = '$filename'
crystal_face = $crystal_face
min_distant = $min_distant
maxdistant = $maxdistant
xyzinfo = ase.io.read(filename)
cell_x = xyzinfo.cell[0][0]
cell_y = xyzinfo.cell[1][1]
cell_z = xyzinfo.cell[2][2]
positions = xyzinfo.get_positions()
pos_x = positions[:, 0]
pos_y = positions[:, 1]
pos_z = positions[:, 2]
def create_group(filename, crystal_face):
    A = crystal_face[0]
    B = crystal_face[1]
    C = crystal_face[2]
    D = -(crystal_face[0] * cell_x)
    numerator = (A*pos_x + B*pos_y + C*pos_z + D)
    denominator = np.sqrt(A**2 + B**2 + C**2)
    distance = numerator / denominator
    flag = np.logical_and(distance > min_distant, distance < maxdistant)
    group_flag = flag.astype(int)
    xyzinfo.arrays['group'] = group_flag
    ase.io.write('crystalface_grouping.xyz', xyzinfo)
create_group(filename, crystal_face)
EOF
}

cutting_crystal_surface(){
# 该脚本用于将晶面上的原子进行切割，生成晶面切割后的模型文件
# 使用方法为 cutting_crystal_surface POSCAR [1,1,1] 将会按照晶面111面的距离进行切割，并且将111面的原子放在z轴的方向
    local filename=$1
    local crystal_face_x=$2
    local crystal_face_y=$3
    local crystal_face_z=$4
    python3 << EOF
from ase.io import read,write
from ase.visualize import view
from ase.build import surface

filename = '$filename'
Atoms = read(filename)
crystal_face = ($crystal_face_x, $crystal_face_y, $crystal_face_z)
s1 = surface(Atoms, crystal_face , 1)
print(f"已完成对于晶面{crystal_face}的切割")
write('cutting_surface.xyz',s1, format='extxyz')
EOF
}



get_grain_count(){
#该函数能计算出指定的轨迹文件（例如dump.xyz）中的晶粒个数，并输出到文件中，并画出晶粒个数随时间变化的图
#使用案例：get_dump_grain_count FCC dump_900.xyz
    local crystal_type=$1
    local filename=$2
    local mingrainsize=${3:-100}
    local itype_upper=$(echo "$crystal_type" | tr '[:lower:]' '[:upper:]')
    local crystal_type_list="OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE"
    if [[ ! $crystal_type_list =~ (^|[[:space:]])"$itype_upper"($|[[:space:]]) ]]; then
        echo "Error: $itype_upper is not a valid crystal type."
        echo "Valid crystal types are: $crystal_type_list"
        return 1
    fi

    python3 << EOF
import ovito
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import ovito.modifiers
matplotlib.use('Agg')

filename = '$filename'

def export_plot_data(listx,listy,filename:str):
    #该函数接受列表x与列表y,将列表x与列表y写入文件第一列与第二列
    listx = [int(x) for x in listx]
    listy = [int(y) for y in listy]
    date = np.column_stack((listx,listy))
    #print(date)
    np.savetxt(filename,date,delimiter=' ',fmt='%d')
    print('数据已保存至',filename)

def plot_grain_zone(filename):
    """
    读取一个轨迹文件，并进行晶区个数的识别。
    """
    pipeline = ovito.io.import_file(filename, multiple_frames=True)

    ptm_modifier = ovito.modifiers.PolyhedralTemplateMatchingModifier()
    ptm_modifier.structures[ovito.modifiers.PolyhedralTemplateMatchingModifier.Type.$itype_upper].enabled = True
    ptm_modifier.output_orientation = True
    pipeline.modifiers.append(ptm_modifier)
    
    grain_segmentation_modifier = ovito.modifiers.GrainSegmentationModifier(min_grain_size=$mingrainsize)
    pipeline.modifiers.append(grain_segmentation_modifier)
    
    grain_cout = []
    frames = []
    crystal_atoms = []

    for frame_index in range(pipeline.source.num_frames):
        data = pipeline.compute(frame_index)
        frames.append(frame_index)
        crystal_atoms.append(data.attributes['PolyhedralTemplateMatching.counts.$itype_upper'])
        grain_cout.append(data.attributes['GrainSegmentation.grain_count'])
    
    export_plot_data(frames,crystal_atoms,'crystal_atoms_${filename%.*}.txt')
    export_plot_data(frames,grain_cout,'grain_count_${filename%.*}.txt')
    # 绘制晶粒数量变化的图
    plt.plot(frames, grain_cout)
    plt.xlabel('Frame')
    plt.ylabel('Grain Count')
    plt.savefig('grain_count_${filename%.*}.png')
    plt.close()

    plt.plot(frames, crystal_atoms)
    plt.xlabel('Frame')
    plt.ylabel('Crystal Atoms')
    plt.savefig('crystal_atoms_${filename%.*}.png')
    plt.close()

plot_grain_zone(filename)
EOF

}

analysis_column(){
    # 该函数用来分析一个文件中的某一列数据，并输出其最大值，最小值，绝对最大值，绝对最小值，平均值，标准差
    local filename=$1
    local column_num=$2
    local skiprows=${3:-'0'}
    python3 << EOF
import numpy as np
filename = '$1'
data = np.loadtxt(filename,skiprows=$skiprows)
column_num = $column_num - 1 
data = data[:,column_num] 
print('最大值:',np.max(data))
print('最小值:',np.min(data))
print('绝对最大值:',np.max(np.abs(data)))
print('绝对最小值:',np.min(np.abs(data)))
print('平均值:',np.mean(data))
print('标准差:',np.std(data))
EOF
}

tran_xyz2cssr() {
# 内嵌的Python脚本
    python3 - "$1" "$2" << 'EOF'
from ase.io import read, write

# 转换函数
def convert(file_in, file_out):
    # 读取XYZ文件
    atoms = read(file_in)
    
    # 写出CSSR文件
    write(file_out, atoms, format='cssr')

# 主函数
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print("用法: convert_xyz_to_cssr <input.xyz> <output.cssr>")
        sys.return(1)

    # 执行转换
    convert(sys.argv[1], sys.argv[2])

EOF
}



select_xyz_config(){
#选择xyz文件中的某一构型(默认选择最后一个构型)，并输出到一个新的xyz文件中
    local xyz_file=${1:-'train.xyz'}
    local config_elected=${2:-'-1'}
    python3 << EOF
from  ase.io import read, write
atoms = read("$xyz_file",index = ":")
config_elected = int($config_elected)
if not config_elected:
    config_elected = -1
if config_elected == -1:
    write("selected.xyz", [atoms[-1]])
else:
    write("selected.xyz", [atoms[config_elected]])
EOF
}

# 该函数用来专门解析碳纤维材料，解析碳纤维的La,Lc,d002,5，6，7碳圆环的数量
analyze_cf(){

    #需要使用到的软件
    cf_exe="cf_analyze"
    xrd_exe="xrd"
    argfile="$1"

    which $cf_exe || (echo "Error: $cf_exe 未安装，请先安装 from https://github.com/kaushikljoshi/cf_analyze" && exit 1) 
    which $xrd_exe || (echo "Error: $xrd_exe 未安装，请先安装 from https://github.com/kaushikljoshi/cf_analyze" && exit 1)

    # 生成cf_analyze需要的参数文件
    cat << EOF > settings.txt
input_file   $argfile  #Name of the input xyz structure
table_flag      3    #1 = lammps connection table, 2 = reax connection table, 3 = build_distance_based table
ring_flag       1    #1 means identify carbon rings
void_flag       0    #1 means identify voids.
periodic_flag   1    #1 means system is periodic. Dims should be in input_file
grid_size       4.0  #Grid size for binning. Will be used only if table_flag is 1
bo_cut_off      0.3  # bo cut-off for identifying molecules for lammps/reax connection table
distance_cut_off 1.7 #This cut_off will be used to identify bonds based on interatomic distace. Default is 1.7
empty_bin_cut_off 1   #If any bin contains atoms less than or equal to this number, then it will be tagged as empty
EOF

    cp $argfile opted_back

    # 提取第二行的 Lattice 信息并生成替换行
    new_line=$(awk 'NR==2 {
    lattice=$0
    sub(/.*Lattice="/,"",lattice)
    sub(/".*/,"",lattice)
    split(lattice,b," ")
    printf("0 %s 0 %s 0 %s 90 90 90", b[1], b[5], b[9])
    }' opted_back)

    # 用 sed 替换第二行
    sed "2c\\$new_line" opted_back > $argfile


    # 运行cf_analyze
    $cf_exe 

    # 生成xrd需要的参数文件
    cat << EOF > md.rc
UNIT 14 md.input old
UNIT 15 $argfile old
UNIT 16 RDF_L1_L2_na.data unknown
UNIT 17 SFACTOR_L1_L2_na.data unknown
EOF

    cat << EOF > md.input
  1811      - NB (Number of number of bins to calculate Q)
  3700      - Nbins (Number of number of bins to calculate g(r)
0.0045      - K_del (step for Q integration)
    66      - Rcell (Size of computational cell for g(r))	
     1      - NRCFLAG (!Nearest image convention flag (0 = Not used, 1 = used))
     0      - WFLAG (If window function is used (0 = Not used, 1 = sinc function is used, 2 = Hann function))	
1.5406      - LAMBDA (Wavelength to calculate in 2Theta, ANGSTROMS)
     1      - UNITFLAG (Units to calculate scattering angle (0 = Inverse angstroms, 1 = 2Theta))
     1      - AFACTORFLAG (Account for atomic scattering factor (0 = Calculate structure factor, 1 = Calculate intensity))
EOF
    

   
    xrd
    cp ~/.rebreath/plot_library/plot_rdf_xrd.py ./
    python3 plot_rdf_xrd.py 
    cp ~/.rebreath/deal_data/xrdtreatment.py ./
    python3 xrdtreatment.py 

    # 打印5，6，7圆环原子数量 rings_size_distribution.dat
    awk '{print $6,$7,$8}' rings_size_distribution.dat
}

# 处理XRD.exe程序输出的文件
# deal_xrd(){
#     cp ~/.rebreath/plot_library/plot_rdf_xrd.py ./
#     python3 plot_rdf_xrd.py 
# }

select_xyz_configs(){
#选择xyz文件中的某一组构型，并输出到一个新的xyz文件中
    local xyz_file=${1:-'train.xyz'}
    local config_start=${2:-'0'}
    local config_end=${3:-'1'}

    python3 << EOF
from  ase.io import read, write
atoms = read("$xyz_file",index = ":")
config_start = int($config_start)
config_end = int($config_end)
write("selected.xyz", atoms[config_start:config_end+1])
EOF
}

select_every_nth_config() {
    # 处理轨迹类型的文件，按照n取一来保存到一个新的文件中
    # 使用的方法 ： select_every_nth_config dump.xyz 10
    # 注意：该函数需要python3.8以上版本，并且需要安装ase包
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
        echo -e "Error: 使用方法应为 select_every_nth_config dump.xyz 10\n表示将dump.xyz文件中的每10帧保存到一个新的文件中"
    fi
}

grouping_to_xyz(){
#该函数用来对xyz类型的文件进行分组
#使用方式为grouping_to_xyz POSCAR x/y/z 1 1.2 3 表示将POSCAR文件中的原子按照x/y/z方向的距离进行分组，分成三份，按照1:1.2:3的比例进行分组
    python3  ~/.rebreath/deal_data/regrouping.py "$@"

}

calc_cf_spatoms(){
# 该函数用来计算碳纤维的杂化原子类型,并计算出结晶率
# 使用方式为 calc_cf_spatoms train.xyz
# 该函数弃用
    python3 ~/.rebreath/deal_data/compute_carbon_fiber_spatoms.py "$@"
}

analyze_crystallinity_fraction() {
    # 用法:
    # plot_crystallinity_fraction FCC dump.xyz [rmse_cutoff]

    local crystal_type_list="OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE"
    local itype="$1"
    local argfile="$2"
    local rmse_cutoff="${3:-0.1}"

    if [[ -z "$itype" || -z "$argfile" ]]; then
        echo "Usage: plot_crystallinity_fraction <CRYSTAL_TYPE> <dump.xyz> [rmse_cutoff]"
        return 1
    fi

    local itype_upper
    itype_upper=$(echo "$itype" | tr '[:lower:]' '[:upper:]')

    echo "Analyzing Crystal Type: $itype_upper"
    echo "Rmse Cutoff = $rmse_cutoff"
    echo "Reading File: $argfile"

    if [[ ! $crystal_type_list =~ (^|[[:space:]])"$itype_upper"($|[[:space:]]) ]]; then
        echo "Error: $itype_upper is not a valid crystal type."
        echo "Valid crystal types are: $crystal_type_list"
        return 1
    fi

    python3 << EOF
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

crystal_type = "${itype_upper}"
file_pattern = "${argfile}"
rmse_cutoff = float("${rmse_cutoff}")

def export_plot_data(listx, listy, filename: str):
    data = np.column_stack((listx, listy))
    np.savetxt(filename, data, fmt=['%d', '%.6f'], delimiter=' ')
    print("数据已保存至", filename)

def plot_crystal_fraction(file_pattern, color, label):
    pipeline = import_file(file_pattern, multiple_frames=True)

    counts = [pipeline.compute(i).particles.count for i in range(pipeline.source.num_frames)]
    if len(set(counts)) != 1:
        raise ValueError(f"不同帧的原子总数不一致，唯一值为: {sorted(set(counts))}")
    total_atoms = counts[0]

    ptm_modifier = PolyhedralTemplateMatchingModifier()
    ptm_modifier.rmsd_cutoff = rmse_cutoff
    ptm_modifier.structures[getattr(PolyhedralTemplateMatchingModifier.Type, crystal_type)].enabled = True
    pipeline.modifiers.append(ptm_modifier)

    crystal_counts = []
    frames = []

    attr_name = f'PolyhedralTemplateMatching.counts.{crystal_type}'

    for frame_index in range(pipeline.source.num_frames):
        data = pipeline.compute(frame_index)
        frames.append(frame_index)
        crystal_counts.append(data.attributes.get(attr_name, 0))

    crystal_fraction = [count / total_atoms for count in crystal_counts]

    export_plot_data(frames, crystal_fraction, f"{label}.txt")
    plt.plot(frames, crystal_fraction, 'o-', color=color, label=label, markersize=2)

def setup_and_save_plot():
    plt.figure(figsize=(10, 5))
    label = "${argfile%.*}_${itype_upper}"
    plot_crystal_fraction(file_pattern, '#714882', label)

    plt.title(f'{crystal_type} Atom Fraction Over Frames')
    plt.xlabel('Frame index')
    plt.ylabel(f'Fraction of {crystal_type} Atoms')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f"{crystal_type}_atom_fraction.png", dpi=300)

setup_and_save_plot()
EOF
}


analysis_grains_size(){
# 该函数使用来分析模型文件中的晶粒大小，其将会分析某种特定的晶粒类型（如FCC）的每个晶粒中该类型的原子数量，并将其保存到文件中
# 该函数的使用方法为 analysis_grains_size FCC model.xyz 会出书每个FCC晶粒中的FCC的原子的数量。
    crystal_type="OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE"

    itype=$1
    argfile=$2
    # 将 itype 转换为大写
    itype_upper=$(echo "$itype" | tr '[:lower:]' '[:upper:]')

    # 检查 itype_upper 是否在 crystal_type 列表中
    if [[ ! $crystal_type =~ (^|[[:space:]])"$itype_upper"($|[[:space:]]) ]]; then
        echo "Error: $itype_upper is not a valid crystal type."
        echo "Valid crystal types are: $crystal_type"
        return 1
    fi

    python3 << EOF
import ovito
from ovito.modifiers import *
import numpy as np
CrystalType = '$itype_upper'
print(f"Analysing {CrystalType} grains in ${argfile}")
crystal_dict = {'OTHER': '0','FCC' : '1', 'HCP' : '2', 'BCC' : '3', 'ICO' : '4', 'SC' : '5', 'CUBIC_DIAMOND' : '6', 'HEX_DIAMOND' : '7', 'GRAPHENE' : '8'}
typeId = int(crystal_dict[CrystalType])
print(f"Type ID: {typeId}")
config = ovito.io.import_file("${argfile}")

ptm_modifier = PolyhedralTemplateMatchingModifier()
ptm_modifier.structures[PolyhedralTemplateMatchingModifier.Type.FCC].enabled = True
ptm_modifier.output_orientation = True
config.modifiers.append(ptm_modifier)

cluster_modifier = GrainSegmentationModifier()
config.modifiers.append(cluster_modifier)
data = config.compute()

grains = data.particles['Grain'].array
structure_type = data.particles['Structure Type'].array
CrystalType_counts = {}

for grain_id, is_CrystalType in zip(grains, structure_type):
    if grain_id not in CrystalType_counts:
        CrystalType_counts[grain_id] = 0
    if is_CrystalType == typeId:
        CrystalType_counts[grain_id] += 1
for grain_id, count in CrystalType_counts.items():
    print(f"Grain ID: {grain_id}, {CrystalType} atoms: {count}")
total_CrystalType_atoms = sum(CrystalType_counts.values())
print(f"Total FCC atoms: {total_CrystalType_atoms}")
grains = list(CrystalType_counts.keys())
CrystalType_counts = list(CrystalType_counts.values())
context = np.column_stack((grains,CrystalType_counts))
np.savetxt('grain_CrystalType_count.txt',context,fmt='%d',delimiter='\t')
EOF
}

expand_cell(){
# 该函数可以使用ase进行扩胞
# 输入参数为xyz文件名，以及扩胞系数，使用的方法为expand_cell model.xyz 10 10 10
    local filename=${1:-'model.xyz'}
    local replicaate_x=${2:-'1'}
    local replicaate_y=${3:-'1'}
    local replicaate_z=${4:-'1'}
    python3 << EOF
import ase.io
import ase.build

filename = '$filename'
replicate = [[$replicaate_x, 0, 0], 
             [0, $replicaate_y, 0],
             [0, 0, $replicaate_z]]
def expand_cell(filename, replicate):
    atoms = ase.io.read(filename)
    expanded_atoms = ase.build.make_supercell(atoms, replicate)
    ase.io.write('expanded.xyz', expanded_atoms)
expand_cell(filename, replicate)
EOF
}

suggest_expand_coefficient() {
    # 对给定的正交的晶胞进行建议，计算出特定的原子个数左右需要进行扩胞的扩胞系数
    # 使用方法：suggest_expand_coefficient model.xyz 10000,会输出建议的扩胞系数
    local xyz_file=${1:-'model.xyz'}
    local atom_num=${2:-'10000'}

    python3 << EOF
import numpy as np
from ase import Atoms
from ase.io import read

def read_xyz(file_path):
    atoms = read(file_path)
    lattice = atoms.get_cell()
    atom_num = len(atoms)
    return lattice, atom_num

def suggest_expand_coefficient(lattice, init_num, target_num):
    if init_num >= target_num:
        print("当前晶胞中原子数目已经大于或等于目标数目，不需要扩胞。")
        return
    
    lattice_a = np.linalg.norm(lattice[0])
    lattice_b = np.linalg.norm(lattice[1])
    lattice_c = np.linalg.norm(lattice[2])
    
    print("当前晶胞的晶格参数为：", lattice)
    print("当前晶胞的a,b,c分别为：", lattice_a, lattice_b, lattice_c)
    
    lattice_all = lattice_a + lattice_b + lattice_c
    proportion = [lattice_a / lattice_all, lattice_b / lattice_all, lattice_c / lattice_all]
    print("当前晶胞的各个方向的占比为：", proportion)
    
    proportion = [1 / i for i in proportion]
    print("各个晶胞方向的权重为：", proportion)
    proportion = [i / sum(proportion) for i in proportion]
    
    expand_coefficient = target_num / init_num
    print("当前晶胞需要扩胞的总系数为：", expand_coefficient)
    x_1 = (expand_coefficient / (proportion[0] * proportion[1] * proportion[2])) ** (1 / 3)
    
    proportion_expand = [round(proportion[0] * x_1), round(proportion[1] * x_1), round(proportion[2] * x_1)]
    print("当前晶胞需要扩胞的各个方向的系数为：", proportion_expand)
    
    lattice_xyz = [lattice_a * proportion_expand[0], lattice_b * proportion_expand[1], lattice_c * proportion_expand[2]]
    last_all_expand_coefficient = proportion_expand[0] * proportion_expand[1] * proportion_expand[2]
    
    print("\n——————>建议结果如下：")
    print(f"扩胞系数为：{proportion_expand}   总的扩胞系数为：{last_all_expand_coefficient}")
    print("扩胞后的晶格参数为：", lattice_xyz)
    print("扩胞后的晶胞中原子数目为：", init_num * last_all_expand_coefficient)

# 主程序
xyz_file = "$xyz_file"
atom_num = int($atom_num)
lattice, init_num = read_xyz(xyz_file)
suggest_expand_coefficient(lattice, init_num, atom_num)
EOF
}



suggest_expand_coefficient_nebula(){
#对给定的正交的晶胞进行建议，计算出特定的原子个数左右需要进行扩胞的扩胞系数
#使用方法：suggest_expand_coefficient model.xyz 10000,会输出建议的扩胞系数
    local xyz_file=${1:-'model.xyz'}
    local atom_num=${2:-'10000'}
    python3 << EOF
import nebula
config = nebula.read_xyz("$xyz_file")
config = config[0]
init_num = config.atom_num
target_num = int($atom_num)
if init_num >= target_num:
    print("当前晶胞中原子数目已经大于或等于目标数目，不需要扩胞。")

lattice = config.lattice
print("当前晶胞的晶格参数为：",lattice)
lattice_a = lattice[0]
lattice_b = lattice[4]
lattice_c = lattice[8]
print("当前晶胞的a,b,c分别为：",lattice_a,lattice_b,lattice_c)
lattice_all = lattice_a + lattice_b + lattice_c
proportion = [lattice_a/lattice_all, lattice_b/lattice_all, lattice_c/lattice_all]
print("当前晶胞的各个方向的占比为：",proportion)
proportion = [1/i for i in proportion]
print("各个晶胞方向的权重为：",proportion)
proportion = [i/sum(proportion) for i in proportion]

expand_coefficient = target_num/init_num
print("当前晶胞需要扩胞的总系数为：",expand_coefficient)
x_1 = (expand_coefficient/(proportion[0]*proportion[1]*proportion[2]))**(1/3)

proportion_expand = [round(proportion[0]*x_1), round(proportion[1]*x_1), round(proportion[2]*x_1)]
print("当前晶胞需要扩胞的各个方向的系数为：",proportion_expand)
lattice_xyz = [lattice_a*proportion_expand[0], lattice_b*proportion_expand[1], lattice_c*proportion_expand[2]]
last_all_expand_coefficient = proportion_expand[0]*proportion_expand[1]*proportion_expand[2]
print("\n——————>建议结果如下：")
print(f"扩胞系数为：{proportion_expand}   总的扩胞系数为：{last_all_expand_coefficient}")
print("扩胞后的晶格参数为：",lattice_xyz)
print("扩胞后的晶胞中原子数目为：",init_num*last_all_expand_coefficient)
EOF
}


find_column_max(){
#检查文件指定列的最大值 使用方式为 find_column_max 4 thermo.out 4 (检查文件第四列的最大值) 
    if [[ $# -ne 2 ]]; then
    echo "错误：需要提供文件名和列号作为参数。"
    return 1
    fi

    filename=$1
    column=$2

    # 检查文件是否存在且可读
    if [[ ! -f "$filename" || ! -r "$filename" ]]; then
        echo "错误：文件 '$filename' 不存在或不可读。"
        return 1
    fi
    # 验证列号是否为正整数
    if [[ ! $column =~ ^[0-9]+$ ]]; then
        echo "错误：列号 '$column' 必须是正整数。"
        return 1
    fi
    # 使用双引号包含变量，并在 awk 脚本中安全地使用变量
    awk -v col="$column" 'BEGIN {max=0} 
        NR > 0 && ($col > max) {max=$col} 
        END { printf("第 %d 列的最大值为 ：%f\n", col, max) }' "$filename"
}

find_column_abs_max(){
#检查文件指定列的绝对值最大值 使用方式为 find_column_max thermo.out 4 (检查文件第四列的最大值) 
    if [[ $# -ne 2 ]]; then
    echo "错误：需要提供文件名和列号作为参数。"
    return 1
    fi

    filename=$1
    column=$2

    # 检查文件是否存在且可读
    if [[ ! -f "$filename" || ! -r "$filename" ]]; then
        echo "错误：文件 '$filename' 不存在或不可读。"
        return 1
    fi
    # 验证列号是否为正整数
    if [[ ! $column =~ ^[0-9]+$ ]]; then
        echo "错误：列号 '$column' 必须是正整数。"
        return 1
    fi
    # 使用双引号包含变量，并在 awk 脚本中安全地使用变量
    awk -v col="$column" 'BEGIN {max=0} 
        NR > 0 && (sqrt(($col)^2) > sqrt((max)^2)) {max=$col} 
        END { printf("第 %d 列的绝对值最大值为：%f\n", col, max) }' "$filename"
}


replot() {
#该函数将指定文件的前两列进行画图
    if [ "$#" -ne 1 ]; then
        echo "用法: replot <数据文件>"
        return 1
    fi

    python3 <<EOF
import matplotlib.pyplot as plt
import sys

datafile = "$1"

x, y = [], []
with open(datafile, 'r') as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) < 2:
            continue  # 跳过不正确的行
        x_val, y_val = map(float, parts[:2])
        x.append(x_val)
        y.append(y_val)

plt.plot(x, y, marker='o', linestyle='-')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig('plot.png', format='png')
EOF
}


get_col_average(){
#对文件的某一列进行平均，并输出为1列
#使用的方式为 get_col_average thermo.out 10 1 20 表示从thermo.out文件中取第10列，索引1到索引20的平均值,注意最大值最小值接受的为数列切片的最大最小值
#使用案例2：get_col_average thermo.out 10 1
    local filename=$1
    local col_index=$2
    local start_index=$3
    local end_index=$4
    local flag=0
    if (( $# > 4 )); then
        echo "错误：参数过多，请检查命令。"
        return 1
    fi
    if [ -z $3 ]; then
        start_index=0
    fi
    if [ -z $4 ]; then
        end_index=None
    fi
    python3 << EOF
import numpy as np
filename = '$filename'
col_index = int($col_index) - 1
min_index = int($start_index)
max_index = $end_index

data = np.loadtxt(filename)
if data.ndim == 1:
    data = data.reshape(-1, 1)

data = data[:, col_index]
def get_col_average(data, min_index, max_index):
    data = data[min_index:max_index]
    return np.mean(data)

average_value = get_col_average(data, min_index, max_index)

print(f"average file : {filename}\nThe column : {col_index+1}\nRange : [{min_index}:{max_index}]\nAverage value : {average_value:.8f}")
EOF
}

average_file(){
#将输入文件进行平均,使用C++实现，平均后输入到average.out文件中
    local cpplib="$HOME/.rebreath/cpp_lib"
    first_file="$1"
    avg_name=${1:0:4}
    g++ $cpplib/averagefiles.cpp -o average_file
    ./average_file "$@" 
    rm -f average_file   
}

average_file_s(){
#将输入文件进行平均,使用shell实现
    files_num="$#"
    tmp_file='temp'
    first_file="$1"
    avg_name=${1:0:4}
    paste "$@" >$tmp_file
    # 逐行读取处理，无论多少行都能处理
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

average_file_c(){
#将输入文件进行平均,使用C++实现，平均后输入到average.out文件中,后将文件改名
    local cpplib="$HOME/.rebreath/cpp_lib"
    first_file="$1"
    avg_name=${1:0:4}
    g++ $cpplib/averagefiles.cpp -o average_file
    ./average_file "$@" 
    rm -f average_file
    mv average.out average_${avg_name}.out
}



generate_large_primes() {
# 辅助函数：生成大于10000的质数列表
    local count=$1
    local primes=()
    local num=${2:-10001}
    while [ ${#primes[@]} -lt $count ]; do
        local is_prime=1
        for ((j=2; j*j<=num; j++)); do
            if [ $((num % j)) -eq 0 ]; then
                is_prime=0
                break
            fi
        done
        if [ $is_prime -eq 1 ]; then
            primes+=($num)
        fi
        num=$((num + 1))
    done
    echo "${primes[@]}"
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
    sed -E "/^\s*A\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      A   $cell_A/"  $cp2k_inp |sed -E "/^\s*B\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      B   $cell_B/" |sed -E "/^\s*C\s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+ \s*[0-9]*\.[0-9]+/s/.*/      C   $cell_C/" > ${cp2k_inp%.*}_up.inp

     sed -i "/@SET XYZFILE/s/.*/@SET XYZFILE    $1/"  ${cp2k_inp%.*}_up.inp
}

calc_coordination_number() {
# 该函数将会计算一个给定模型的配位数
# 使用方式：calc_coordination_number model.xyz 2.5 表示将会计算model.xyz模型的配位数，配位半径为2.5中每个原子有多少个邻居,该脚本使用了ase库来读取模型文件，并且使用了最小镜像约定来计算周期性情况的配位数
    local xyz_file="$1"
    local r_cut=${2:-2.5}
    python3 << EOF
# 该脚本为使用python来计算配位数
import numpy as np
import ase.io as ai
from ase.geometry import get_distances
config = ai.read('$xyz_file')
positions = config.positions
cell = config.cell.diagonal()  # 获取晶胞对角线上的元素
half_cell = cell / 2
def mirror_pos(pos):
    """
    实现最小化镜像约定，原子间的距离必须要小于半晶胞的距离
    """
    for i in range(3):
        if pos[i] < -half_cell[i]:
            pos[i] = pos[i] + cell[i]
        elif pos[i] > half_cell[i]:
            pos[i] = pos[i] - cell[i]
    return pos

positions = config.positions
distances = np.zeros((len(config), len(config)))
for i in range(len(config)):
    for j in range(i+1, len(config)):
        mirrpos = mirror_pos(config.positions[i]-config.positions[j])
        distances[i][j] = np.linalg.norm(mirrpos)   
        distances[j][i] = distances[i][j]
coordination_numbers = np.zeros(len(config))
for i in range(len(config)):
    coordination_numbers[i] = np.sum((distances[i] < $r_cut)&(distances[i] > 0))
np.savetxt('coordination_numbers.txt', coordination_numbers, fmt='%d')
print("\n配位数计算完成，结果保存在coordination_numbers.txt文件中 OVO")
EOF
}

calc_coordination_number_cpp() {
    # 该函数将会计算一个给定模型的配位数
    # 使用方式：calc_coordination_number model.xyz 2.5 表示将会计算model.xyz模型的配位数，配位半径为2.5中每个原子有多少个邻居,该脚本使用了ase库来读取模型文件，并且使用了最小镜像约定来计算周期性情况的配位数
    # 该函数可以使用了C++来实现，并且使用了最小镜像约定来计算周期性情况的配位数，并且可以计算斜胞的配位数
    local xyz_file="$1"
    local r_cut=${2:-2.5}
    cat > coordination_numbers.cpp << EOF
// 该程序使用来计算模型中每个原子的配位数
// 该程序将会试图处理斜胞的情况
#include <iostream>
#include <vector>
#include <fstream>
#include <string>s
#include <sstream>
#include <regex>
#include <cmath>
using namespace std;
struct Atoms{
    int atom_num;
    double box[9];
    double invBox[9];
    vector<string> atomType;
    vector<double> coordX, coordY, coordZ;
};
enum class LineType{
    FirstLine,
    SecondLine,
    OtherLine
};
double getDet(double (&box)[9]){
    //获得矩阵的模
    return box[0] * (box[4] * box[8] - box[5] * box[7])
         - box[1] * (box[3] * box[8] - box[5] * box[6])
         + box[2] * (box[3] * box[7] - box[4] * box[6]);
}
void InverseBox(Atoms& atoms) {
    //获得矩阵的逆矩阵
    double (&box)[9] = atoms.box;
    double det = getDet(atoms.box);
    atoms.invBox[0] = (box[4] * box[8] - box[5] * box[7]) / det;
    atoms.invBox[1] = (box[2] * box[7] - box[1] * box[8]) / det;
    atoms.invBox[2] = (box[1] * box[5] - box[2] * box[4]) / det;
    atoms.invBox[3] = (box[5] * box[6] - box[3] * box[8]) / det;
    atoms.invBox[4] = (box[0] * box[8] - box[2] * box[6]) / det;
    atoms.invBox[5] = (box[2] * box[3] - box[0] * box[5]) / det;
    atoms.invBox[6] = (box[3] * box[7] - box[4] * box[6]) / det;
    atoms.invBox[7] = (box[1] * box[6] - box[0] * box[7]) / det;
    atoms.invBox[8] = (box[0] * box[4] - box[1] * box[3]) / det;
}
string tolow(string& str) {
    for (auto& c : str) {
        c = tolower(c);
    }
    return str;
}
Atoms read_xyz(string filename){
    //该函数将会试图读取一个xyz文件，并将其中的原子坐标和盒子信息读入到Atoms结构体中
    //使用正则表达式来处理第二行可能会更好
    Atoms atoms;
    ifstream infile(filename);
    string line;
    LineType state {LineType::FirstLine};
    while (getline(infile, line)) {
        if (state == LineType::FirstLine) {
            atoms.atom_num = stoi(line);
            state = LineType::SecondLine;
        }
        else if (state == LineType::SecondLine) {
                line = tolow(line);
                regex re(R"(\blattice=.*\"\s*([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+).*\")");
                smatch match;
                
                if (regex_search(line, match, re) && match.size() > 1) {

                    for (size_t i = 1; i < match.size(); ++i) {
                        atoms.box[i - 1] = std::stod(match[i].str());
                    }
            }
            state = LineType::OtherLine;
        }
        else if (state == LineType::OtherLine) {
            istringstream iss(line);
            double x, y, z;
            string atom_type;
            iss >> atom_type >> x >> y >> z;
            atoms.atomType.push_back(atom_type);
            atoms.coordX.push_back(x);
            atoms.coordY.push_back(y);
            atoms.coordZ.push_back(z);
            //cout << atom_type << " " << x << " " << y << " " << z << endl;
        }
        }
    infile.close();
    return atoms;
}

double Mircone(double& xyz){
    if (xyz < -0.5) {
        return xyz + 1.0;
    }
    else if (xyz > 0.5) {
        return xyz - 1.0;
    }
    return xyz;
}

void applyMirc(double (&box)[9],double (&InverseBox)[9], double& x, double& y, double& z){
    // 该函数将会应用镜像约定法，让坐标满足周期性
    // 输入的为两个原子间的距离的坐标
    double fracX , fracY, fracZ;
    fracX = InverseBox[0] * x + InverseBox[1] * y + InverseBox[2] * z;
    fracY = InverseBox[3] * x + InverseBox[4] * y + InverseBox[5] * z;
    fracZ = InverseBox[6] * x + InverseBox[7] * y + InverseBox[8] * z;
    x = Mircone(fracX);
    y = Mircone(fracY);
    z = Mircone(fracZ);
    x = box[0] * x + box[1] * y + box[2] * z;
    y = box[3] * x + box[4] * y + box[5] * z;
    z = box[6] * x + box[7] * y + box[8] * z;
}

int main()
{
    Atoms atoms = read_xyz("$xyz_file");

    //转换坐标到分数坐标
    InverseBox(atoms);
    for (auto i:atoms.invBox){
        cout<<i<<" ";
    }
    vector<int> coordination_number(atoms.atom_num, 0);
    vector<vector<double>> distances(atoms.atom_num, vector<double>(atoms.atom_num));
    for (int i = 0; i < atoms.atom_num; i++) {
        for (int j = i + 1; j < atoms.atom_num; j++) {
            double dx = atoms.coordX[i] - atoms.coordX[j];
            double dy = atoms.coordY[i] - atoms.coordY[j];
            double dz = atoms.coordZ[i] - atoms.coordZ[j];
            applyMirc(atoms.box, atoms.invBox, dx, dy, dz);
            distances[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
            distances[j][i] = distances[i][j];
        }
    }
    for (int i = 0; i < atoms.atom_num; i++) {
        for (int j = 0; j < atoms.atom_num; j++) {
            //cout << i << " " << j << " " << distances[i][j] << endl;
            if (distances[i][j] < $r_cut && distances[i][j] > 0) {
                coordination_number[i]++;
            }
        }
    }
    fstream outfile("coordination_number.txt", ios::out);
    for (int i = 0; i < atoms.atom_num; i++) {
        outfile << atoms.atomType[i] << " " << coordination_number[i] << endl;
    }
    outfile.close();
    return 0;
}
EOF
    g++ -O3 -o coordination_number coordination_numbers.cpp
    ./coordination_number
    rm -f coordination_numbers.cpp coordination_number
}


analysis_model(){
# 该函数使用来分析一个原子model文件输出其原子的个数，与每种原子的种类
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
    print(f"Total number of atoms: {total_atoms}")
    print(f"Number of unique elements: {len(unique_elements)}")
    print("Element counts:")
    for element, count in element_counts.items():
        print(f"{element}: {count}")
    return total_atoms, unique_elements, element_counts
filename = '$xyz_file'  
analyze_xyz_file(filename)
EOF
}

build_adsorption_model(){
# 该函数可以指定基底模型与需要吸附的模型，而后生成各种吸附位点的吸附模型
# 使用方法： build_adsorption_model alphaFe2O3_012.xyz O3 表示将会构建一个alphaFe2O3的基底模型，并吸附O3
    local base_model=$1
    local adsorbate=${2:-O3}
    if [ -z "$base_model" ]; then
        echo 'Use method: build_adsorption_model $base_model $adsorbate'
        return 1
    fi

    python3 <<EOF
import ase.io as ai
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab
from pymatgen.io.ase import AseAtomsAdaptor
import ase.build as ab
import matplotlib.pyplot as plt
import numpy as np

# 构建吸附分子
adsorbate_ab = ab.molecule("$adsorbate")
adsorbate = AseAtomsAdaptor.get_molecule(adsorbate_ab)

def calc_key_lengthandangle(atoms) -> tuple:
    """计算分子的键长和键角"""
    key_length = atoms.get_distance(0, 1)
    key_angle = atoms.get_angle(0, 1, 2)
    return key_length, key_angle

def find_adsorption_sites(filename) -> list:
    """查找吸附点并生成吸附结构"""
    atoms = ai.read(filename)
    struct = AseAtomsAdaptor.get_structure(atoms)  # ✅ 直接转换，无需写入 .cif

    # 创建吸附查找器
    asf = AdsorbateSiteFinder(struct)

    # 查找吸附点
    adsorption_sites = asf.find_adsorption_sites()
    print(f"\nFile: {filename}\nAdsorption sites found:\n {adsorption_sites}")

    # 遍历吸附点并生成结构
    for idx, site in enumerate(adsorption_sites['all']):
        print(f"Adsorption site {idx+1}: {site}")
        absorbed_structure = asf.add_adsorbate(molecule=adsorbate, ads_coord=site)
        atoms_temp = AseAtomsAdaptor.get_atoms(absorbed_structure)
        output_filename = f"adsorbed_structure_{idx}.xyz"
        ai.write(output_filename, atoms_temp)
        print(f"✅ Saved adsorbed structure to {output_filename}")

    return adsorption_sites

# 执行吸附点查找
adsorption_sites = find_adsorption_sites('$base_model')
EOF
}

get_molecule(){
    # 该函数用来获取一个分子的xyz文件
    # 使用方法：get_molecule HCOOH 表示将会获取一个HCOOH的分子的xyz文件
    cat << EOF > getmolecule.py
import ase.io as ai
import ase.build as ab

# 构建吸附分子
adsorbate_ab = ab.molecule("$1")
ai.write('adsorbate.xyz', adsorbate_ab)
EOF
    python3 getmolecule.py
}

build_disorption_model(){
    # 该函数用来构建脱附模型，可以指定一个模型，并指定需要脱附的原子，生成脱附模型，脱附模型会尽可能的离其他的原子更远
    # 使用方法：build_disorption_model Fe2O3_012.xyz O3 表示将会构建一个Fe2O3的脱附模型，并脱附O3
    local base_model=$1
    local disorbate=${2:-O3}
    local min_distance=${3:-10}

    cat << EOF > disorption_model.py
import numpy as np
import ase.io as ai
#import ase.visualize as av
import ase.build as ab

atoms = ai.read('$base_model')
#av.view(atoms)
cell = atoms.get_cell()
cell_diagonal = [cell[0][0], cell[1][1], cell[2][2]]
print(cell_diagonal)

insert_atoms = ab.molecule('$disorbate')
edge_length = 2

# 初始的点为八个顶点
init_point = np.array([
    [0 + edge_length, 0 + edge_length, 0 + edge_length],
    [0 + edge_length, 0 + edge_length, cell_diagonal[2] - edge_length],
    [0 + edge_length, cell_diagonal[1] - edge_length, 0 + edge_length],
    [0 + edge_length, cell_diagonal[1] - edge_length, cell_diagonal[2] - edge_length],
    [cell_diagonal[0] - edge_length, 0 + edge_length, 0 + edge_length],
    [cell_diagonal[0] - edge_length, 0 + edge_length, cell_diagonal[2] - edge_length],
    [cell_diagonal[0] - edge_length, cell_diagonal[1] - edge_length, 0 + edge_length],
    [cell_diagonal[0] - edge_length, cell_diagonal[1] - edge_length, cell_diagonal[2] - edge_length]
])

def wrap_around_cell(diff):
    """
    该函数使用来应用周期性边界条件，将超出一定范围的原子进行抓取回来
    """
    for j in range(3):
        diff[:, j] = (diff[:, j] + cell_diagonal[j] / 2) % cell_diagonal[j] - cell_diagonal[j] / 2
    return diff



def wrap_point(position):
    """
    该函数用于应用周期性边界条件，将超出一定范围的单个点进行折回
    """
    wrapped_position = position.copy()  # 确保不修改原始位置
    for i in range(3):
        if wrapped_position[i] > cell_diagonal[i]:
            wrapped_position[i] = wrapped_position[i] % cell_diagonal[i]
        elif wrapped_position[i] < 0:
            wrapped_position[i] = wrapped_position[i] % cell_diagonal[i]
    return wrapped_position

min_distances = []
Nearest_atom_index = []

def find_desorption_point(atoms, init_point, redistance=$min_distance):
    """
    该函数用于初步寻找脱附位置，并返回脱附位置的坐标
    """
    positions = atoms.get_positions()
    desorption_point = []
    for i in range(len(init_point)):
        diff = positions - init_point[i]
        diff = wrap_around_cell(diff)
        temp_distance = np.linalg.norm(diff, axis=1)
        min_distance = np.min(temp_distance)
        Nearest_atom_index.append(np.argmin(temp_distance))
        min_distances.append(min_distance)

        if min_distance > redistance:
            desorption_point.append(init_point[i])
            return desorption_point

    selected_init_point = init_point[np.argmax(min_distances)]

    # Z轴的方向优化    
    m = 1
    stage_z_distance = []
    stage_z_position = []
    flag = True
    while flag:
        step = 0.5
        selected_init_point[2] = selected_init_point[2] + step
        selected_init_point = wrap_point(selected_init_point)
        print(f"Z轴第{m}次优化，当前位置为{selected_init_point}")
        stage_z_position.append(selected_init_point)
        diff = positions - selected_init_point
        diff = wrap_around_cell(diff)
        temp_distance = np.linalg.norm(diff, axis=1)
        min_distance = np.min(temp_distance)
        stage_z_distance.append(min_distance)
        m += 1
        print(f"最近的原子的距离为{min_distance}")
        if min_distance > redistance:
            desorption_point = selected_init_point
            return desorption_point
        if m > 100:
            flag = False

    optimized_z_index = np.argmax(stage_z_distance)
    selected_init_point = stage_z_position[optimized_z_index]
    print(f"------->经过Z轴优化后的点为{selected_init_point},最近的距离为{stage_z_distance[optimized_z_index]}")
    print("start进行X轴优化")

    # X轴进行优化
    m = 1
    stage_x_distance = []
    stage_x_position = []
    stage_x_distance.append(stage_z_distance[optimized_z_index])
    stage_x_position.append(selected_init_point)
    flag = True
    while flag:
        step = 0.5
        selected_init_point[0] = selected_init_point[0] + step
        selected_init_point = wrap_point(selected_init_point)
        stage_x_position.append(selected_init_point)
        print(f"X轴第{m}次优化，当前位置为{selected_init_point}")
        diff = positions - selected_init_point
        diff = wrap_around_cell(diff)
        temp_distance = np.linalg.norm(diff, axis=1)
        min_distance = np.min(temp_distance)
        stage_x_distance.append(min_distance)
        m += 1
        print(f"最近的原子的距离为{min_distance}")
        if min_distance > redistance:
            desorption_point = selected_init_point
            return desorption_point
        if m > 50:
            flag = False

    optimal_index = np.argmax(stage_x_distance)
    optimal_point = stage_x_position[optimal_index]
    print(f"------->经过X轴的优化后的点为{optimal_point},最近的距离为{stage_x_distance[optimal_index]}")

    # Y轴进行优化
    m = 1
    stage_y_distance = []
    stage_y_position = []
    stage_y_distance.append(stage_x_distance[optimal_index])
    stage_y_position.append(optimal_point)
    flag = True
    while flag:
        step = 0.5
        selected_init_point[1] = selected_init_point[1] + step
        selected_init_point = wrap_point(selected_init_point)
        stage_y_position.append(selected_init_point)
        print(f"Y轴第{m}次优化，当前位置为{selected_init_point}")
        diff = positions - selected_init_point
        diff = wrap_around_cell(diff)
        temp_distance = np.linalg.norm(diff, axis=1)
        min_distance = np.min(temp_distance)
        stage_y_distance.append(min_distance)
        m += 1
        print(f"最近的原子的距离为{min_distance}")
        if min_distance > redistance:
            desorption_point = selected_init_point
            return desorption_point
        if m > 50:
            flag = False

    optimal_index = np.argmax(stage_y_distance)
    optimal_distance = stage_y_distance[optimal_index]
    optimal_point = stage_y_position[optimal_index]
    print(f"------->经过Y轴的优化后的点为{optimal_point},最近的距离为{optimal_distance}")

    return optimal_point
    
        
def build_desorption_structure(atoms, insert_atoms):
    """
    该函数使用来构建脱附结构
    """
    desorption_point = find_desorption_point(atoms, init_point)
    if len(desorption_point) == 0:
        print("没有找到脱附位置，请重新运行程序")
        return
    insert_atoms.translate(desorption_point)
    insert_positions = insert_atoms.get_positions()
    #wrapped_positions = np.array([wrap_point(pos) for pos in insert_positions])
    insert_atoms.set_positions(insert_positions)
    
    #print(f"{insert_atoms}, type: {type(insert_atoms)}")
    atoms.extend(insert_atoms)
    atoms.write('POSCAR_desorption.vasp')
    #av.view(atoms)

build_desorption_structure(atoms, insert_atoms)
EOF

    python3 disorption_model.py > build_desorption_model.log 2>&1
    rm disorption_model.py
    tail -n 1 build_desorption_model.log
}

analyze_xyz(){
# 该函数用来分析xyz文件的详细信息
    cp $HOME/.rebreath/deal_data/analyze_xyz_detail.py .
    python3 analyze_xyz_detail.py $1
}


copy_each_to_own_dir(){
# 该函数可以读取一类指定后缀的文件，将会创建一系列文件夹而后将文件分类到文件夹中
# 例如指定xyz类型的文件，则会创建一系列去除xyz后缀的文件夹，将xyz文件放到指定的xyz文件夹中

    cut_name=${1:-xyz}
    for i in $(ls *.$cut_name |xargs -n 1);do 
        fore_name=${i%.*}
        mkdir -p $fore_name
        cp $i $fore_name
    done
}

cp_file_to_subdirs(){
#该函数将会将一个指定的文件送到当前目录下所有的文件夹中
    local sfile=${1}
    for i in $(find $PWD -maxdepth 1 -mindepth 1  -type d);do
        cp $sfile $i
    done
}

calc_xrd_usedebyer(){
# 该函数使用debyer计算xrd,目前主要计算碳纤维的xrd
    local xyzfile=${1:-opted.xyz}
    debyer -x -l1.5406 -f5 -t80 -s0.02 -o xrd.dat $xyzfile
}

xyz_group_by_type(){
# 该函数用来将原子按照原子类型分组
# 使用方法：xyz_group_by_type opted.xyz
    local xyzfile=${1}
    if [ -z "$xyzfile" ]; then
        echo "Error: xyzfile is not set."
        echo "Usage: xyz_group_by_type <xyzfile>"
        return 1
    fi
    exec 3< "$HOME/.rebreath/deal_data/grouping_use_atomtype.py"
    python3 /dev/fd/3 -- "$xyzfile"
    exec 3<&-
}


supercell_auto_cubic() {
# 该函数用来智能扩胞，根据原子的分布情况，自动判断如何扩胞才能达到目标原子数，并且尽可能使box为正方体形状
# 使用方法：supercell_auto_cubic opted.xyz 10000
    local xyzfile="$1"
    local targetnumber="$2"

    if [ "$#" -ne 2 ]; then
    echo "Usage: supercell_auto_cubic <xyzfile> <targetnumber>"
    return 1
    fi

    exec 3< "$HOME/.rebreath/deal_data/supercell_auto_cubic.py"
    python3 /dev/fd/3 -- "$xyzfile" "$targetnumber"
    exec 3<&-
}

get_phonon_spectrum_mpdata(){
# 该函数用来获取material project 中指定mpid的 phonon 谱图
# 使用方法：get_phonon_spectrum_mpdata mp-2741 
# 可选参数：method, 默认为dfpt
# 原理：
# 1. 从material project 中下载指定mpid的 phonon 谱图
# 2. 调用plot_phonon_spectrum_mpdata.py 绘制 phonon 谱图

    local mpid=${1}
    local method=${2:-dfpt}

    if [ -z "$mpid" ]; then
        echo "Error: mpid is not set."
        echo "Usage: get_phonon_spectrum_mpdata <mpid>"
        return 1
    fi
    cp $HOME/.rebreath/compute_lib/mp_phonon_data_extract.py  .
    python3 mp_phonon_data_extract.py $mpid $method
    cp $HOME/.rebreath/compute_lib/plot_phonon_spectrum_mpdata.py  .
    python3 plot_phonon_spectrum_mpdata.py ${mpid}_phonon_bs_${method}.json
}


plot_volume_per_atom_xyz() {
    local xyz_file="$1"
    local out_dat="${2:-volume_per_atom.dat}"
    local out_png="${3:-volume_per_atom.png}"

    if [ -z "$xyz_file" ]; then
        echo "Usage: plot_volume_per_atom_xyz input.xyz [output.dat] [output.png]"
        return 1
    fi

    if [ ! -f "$xyz_file" ]; then
        echo "File not found: $xyz_file"
        return 1
    fi

    python - "$xyz_file" "$out_dat" "$out_png" <<'PY'
import sys
import numpy as np
from ase.io import iread
import matplotlib.pyplot as plt

xyz_file = sys.argv[1]
out_dat = sys.argv[2]
out_png = sys.argv[3]

steps = []
volumes = []
natoms_list = []
vpa_list = []

for i, atoms in enumerate(iread(xyz_file, index=":")):
    natoms = len(atoms)
    cell = atoms.get_cell()

    if cell is None or abs(cell.volume) < 1e-12:
        raise ValueError(
            f"Frame {i} has no valid cell information. "
            "Your xyz trajectory must contain lattice/cell info for volume calculation."
        )

    volume = atoms.get_volume()
    vpa = volume / natoms

    steps.append(i)
    natoms_list.append(natoms)
    volumes.append(volume)
    vpa_list.append(vpa)

data = np.column_stack([steps, natoms_list, volumes, vpa_list])
header = "step natoms total_volume_A3 volume_per_atom_A3"
np.savetxt(out_dat, data, header=header, fmt=["%d", "%d", "%.6f", "%.6f"])

print(f"prim volume_per_atom_A3 = {vpa_list[0]:.3f}")
print(f"last volume_per_atom_A3 = {vpa_list[-1]:.3f}")

plt.figure(figsize=(7, 5))
plt.plot(steps, vpa_list, linewidth=1.8)
plt.xlabel("Frame")
plt.ylabel(r"Volume per atom ($\AA^3$/atom)")
plt.title("Volume per atom evolution")
plt.tight_layout()
plt.savefig(out_png, dpi=300)
print(f"Saved data to: {out_dat}")
print(f"Saved figure to: {out_png}")
PY
}


generate_vacancy_defect_xyz() {
    # 用法:
    # generate_vacancy_defect_xyz input.xyz output.xyz defect_fraction [species] [seed]
    #
    # 例子:
    # generate_vacancy_defect_xyz perfect.xyz Al_vac_1pct.xyz 0.01 Al 123
    # generate_vacancy_defect_xyz perfect.xyz N_vac_2pct.xyz 0.02 N 123
    # generate_vacancy_defect_xyz perfect.xyz rand_vac_1pct.xyz 0.01 all 123
    #
    # 参数:
    #   input.xyz         输入xyz文件
    #   output.xyz        输出xyz文件
    #   defect_fraction   缺陷比例，例如 0.01 表示删除1%
    #   species           删除的元素类型: Al / N / all，默认 all
    #   seed              随机种子，默认 12345

    local input_file="$1"
    local output_file="$2"
    local defect_fraction="$3"
    local species="${4:-all}"
    local seed="${5:-12345}"

    if [[ -z "$input_file" || -z "$output_file" || -z "$defect_fraction" ]]; then
        echo "Usage: generate_vacancy_defect_xyz input.xyz output.xyz defect_fraction [species] [seed]"
        return 1
    fi

    if [[ ! -f "$input_file" ]]; then
        echo "Error: input file '$input_file' not found."
        return 1
    fi

    python3 - << EOF
import random
import math

input_file = "${input_file}"
output_file = "${output_file}"
defect_fraction = float("${defect_fraction}")
species = "${species}"
seed = int("${seed}")

if defect_fraction < 0 or defect_fraction >= 1:
    raise ValueError("defect_fraction must satisfy 0 <= defect_fraction < 1")

random.seed(seed)

with open(input_file, "r") as f:
    lines = [line.rstrip("\n") for line in f]

if len(lines) < 3:
    raise ValueError("Input xyz file is too short.")

try:
    natoms = int(lines[0].strip())
except Exception:
    raise ValueError("First line of xyz must be the atom count.")

comment = lines[1]
atom_lines = lines[2:]

if len(atom_lines) != natoms:
    raise ValueError(f"XYZ format mismatch: header says {natoms} atoms, but found {len(atom_lines)} atom lines.")

# 解析原子行
parsed_atoms = []
for i, line in enumerate(atom_lines):
    parts = line.split()
    if len(parts) < 4:
        raise ValueError(f"Invalid atom line at index {i}: '{line}'")
    elem = parts[0]
    parsed_atoms.append((elem, line))

# 选择可删除的原子索引
if species.lower() == "all":
    candidate_indices = list(range(len(parsed_atoms)))
else:
    candidate_indices = [i for i, (elem, _) in enumerate(parsed_atoms) if elem == species]

if len(candidate_indices) == 0:
    raise ValueError(f"No atoms of species '{species}' found in file.")

n_remove = int(round(len(candidate_indices) * defect_fraction))

if n_remove == 0 and defect_fraction > 0:
    n_remove = 1

if n_remove > len(candidate_indices):
    raise ValueError("Requested removal count exceeds available candidate atoms.")

remove_indices = set(random.sample(candidate_indices, n_remove))

new_atom_lines = [line for i, (_, line) in enumerate(parsed_atoms) if i not in remove_indices]
new_natoms = len(new_atom_lines)

with open(output_file, "w") as f:
    f.write(f"{new_natoms}\n")
    f.write(comment + f" | vacancy_defect fraction={defect_fraction} species={species} seed={seed}\n")
    for line in new_atom_lines:
        f.write(line + "\n")

print(f"Input atoms      : {natoms}")
print(f"Candidate atoms  : {len(candidate_indices)}")
print(f"Removed atoms    : {n_remove}")
print(f"Output atoms     : {new_natoms}")
print(f"Output written to: {output_file}")
EOF
}