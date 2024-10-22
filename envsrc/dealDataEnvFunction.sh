# 该分库使用来存放处理数据的函数，包括对于xyz类文件的处理等，较为通用的函数

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

#调用ovito转换文件的格式
xyz_to_poscar() {
    python3 -c "from ovito.io import import_file, export_file; pipeline = import_file('$1'); export_file(pipeline, 'POSCAR_convered', 'vasp')"
}

poscar_to_xyz() {     
       	python3 -c "from ovito.io import import_file, export_file; pipeline = import_file('$1'); export_file(pipeline, 'model_conversed.xyz', 'xyz',columns=['Particle Type', 'Position.X', 'Position.Y', 'Position.Z'])"
}


get_grain_count(){
#该函数能计算出指定的轨迹文件（例如dump.xyz）中的晶粒个数，并输出到文件中，并画出晶粒个数随时间变化的图
#使用案例：get_dump_grain_count FCC dump_900.xyz
    local crystal_type=$1
    local filename=$2
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
    np.savetxt(filename,date,delimiter=' ')
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
    
    grain_segmentation_modifier = ovito.modifiers.GrainSegmentationModifier()
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

xyz_to_cssr() {
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

xyz_to_cif (){ 
    #其实不只是xyz到cif，很多其他的格式也可以到达cif的格式，ase，让我看看你的实力了
    #不过注意xyz文件如果想要使用ovito来进行转化的时候需要注意，最好只输出xyz方向的坐标，多了其他的东西可能会无法识别。
    cifname={2:-'file.cif'}
    python3  <<'EOF'
from ase.io import read, write

def convert(file_in, file_out):
    # 读取XYZ文件
    atoms = read(file_in)
    # 输出为CIF文件
    write(file_out, atoms)
infile = '$1'
outfile = '$cifname
convert(infile, outfile)
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

plot_crysralinity_fraction(){
# 本脚本的功能为对dump.xyz文件的每一帧进行晶型分析，分析其蜂窝状的原子的个数并画出图形进行可视化
# 输入参数为dump.xyz文件名，以及晶型类型，如plot_crysralinity_fraction FCC dump.xyz 
# 能处理的晶型有：OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE
    crystal_type="OTHER FCC HCP BCC ICO SC CUBIC_DIAMOND HEX_DIAMOND GRAPHENE"
    # 输入的晶型
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
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier
from ovito.data import *

crystal_type = '$itype_upper'
# 设置matplotlib不显示图形界面
matplotlib.use('Agg')
def export_plot_data(listx,listy,filename:str):
    #该函数接受列表x与列表y，将列表x与列表y写入文件第一列与第二列
    listx = [int(x) for x in listx]
    listy = [int(y) for y in listy]
    date = np.column_stack((listx,listy))
    #print(date)
    np.savetxt(filename,date,delimiter=' ')
    print('数据已保存至',filename)

def plot_graphene_counts(file_pattern, color, label):
    #计算并绘制每个文件的石墨烯原子计数
    pipeline = import_file(file_pattern, multiple_frames=True)

    ptm_modifier = PolyhedralTemplateMatchingModifier()
    ptm_modifier.structures[PolyhedralTemplateMatchingModifier.Type.$itype_upper].enabled = True
    pipeline.modifiers.append(ptm_modifier)

    graphene_counts = []
    frames = []

    for frame_index in range(pipeline.source.num_frames):
        data = pipeline.compute(frame_index)
        frames.append(frame_index)
        graphene_counts.append(data.attributes['PolyhedralTemplateMatching.counts.$itype_upper'])

    export_plot_data(frames, graphene_counts, f"{label}.txt")
    plt.plot(frames, graphene_counts, 'o-', color=color, label=label,markersize=2)
    


def setup_and_save_plot(file_list):
    #设置图表并保存为PNG文件
    plt.figure(figsize=(10, 5))
    for file_data in file_list:
        plot_graphene_counts(file_data['pattern'], file_data['color'], file_data['label'])

    plt.title('$itype_upper Atom Count Over Time')
    plt.xlabel('Time/10$^{-2}ns$')
    plt.ylabel('Count of $itype_upper Atoms')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.savefig("${itype_upper}_atom_count.png")


file_list = [
    {'pattern': '${argfile}', 'color': '#714882', 'label': '${argfile%.*}_$itype_upper'}
]  # 更多文件可以添加到列表

# 调用函数
setup_and_save_plot(file_list)
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