# 该分库使用来存放处理cp2k数据的函数

cp2krestart2cif() {
# 使用了Multiwfn软件，将cp2k的restart文件转换为cif文件
    local cifilename=${1:-"opt.cif"}
    Multiwfn *.restart << EOF > /dev/null
    100
    2
    33
    $cifilename
    0
    q
EOF
}

cp2krestart2xyz() {
# 使用了Multiwfn软件，将cp2k的restart文件转换为xyz文件
    local xyzfilename=${1:-"opt.xyz"}
    Multiwfn *.restart << EOF  > /dev/null
    100
    2
    2
    $xyzfilename
    0
    q
EOF
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
    sed -i "/@SET PROJECT_NAME/s/.*/@SET PROJECT_NAME    ${xyz_file%.*}/"  ${cp2k_inp%.*}_up.inp
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

get_cp2k_energy() {
# 该函数从CP2K的log文件中提取能量
# 使用方式：get_cp2k_energy cp2k.log
# 达成目的：提取CP2K的能量
    local logfile=${1:-"cp2k.log"}
    grep "ENERGY|" $logfile |tail -n 5
}