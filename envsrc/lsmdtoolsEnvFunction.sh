# 该文件库包含了测试自研软件LSMD的各种环境函数

# ---------------- Test lmp model ----------------
# 测试lmp模型
# 使用方法：testlmpmod lmpmodelfile
# 将测试脚本放到~/.comparewithlmp/lmptest目录下，而后运行测试脚本
testlmpmod(){
    lmpmodelfile=${1}
    lmpmodelfile_address=$PWD
    mkdir -p ~/.comparewithlmp/lmptest
    rm -f ~/.comparewithlmp/* || true
    cp ${lmpmodelfile} ~/.comparewithlmp/lmptest
    cp ~/.rebreath/inp_lib/Lammps/testall.inp ~/.comparewithlmp/lmptest
    sed -E -i '/^[[:space:]]*read_data/s/.*/read_data   '"${lsmdmodelfile}"'/' ~/.comparewithlmp/lmptest/testall.inp
    conda activate deepmd || true
    cd ~/.comparewithlmp/lmptest/
    mpirun -np 4 lmp -i testall.inp
    cd ${lmpmodelfile_address}
    echo "${lmpmodelfile} has been tested."
}

# 在当前目录下测试lmp模型
testlmpdir(){
    lmpmodelfile=${1}
    lmptemplate=${2:-300}  # 默认为300K
    lmpmodelfile_address=$PWD
    cp ~/.rebreath/inp_lib/Lammps/testall.inp .
    sed -E -i '/^[[:space:]]*read_data/s/.*/read_data   '"${lsmdmodelfile}"'/' testall.inp
    #velocity all create 300.0 12345 mom yes rot yes dist gaussian
    sed -E -i '/^([[:space:]]*)velocity all create/s~.*~\1velocity all create '"${lmptemplate}"' 12345 mom yes rot yes dist gaussian~' testall.inp
    conda activate deepmd || true
    mpirun -np 1 lmp -i testall.inp
    echo "${lmpmodelfile} has been tested."
}


# ---------------- Test lsmd model ----------------
# 测试lsmd模型
# 使用方法：testlsmdmod lsmdmodelfile
testlsmdmod(){

    lsmdmodelfile=${1}
    lsmdmodelfile_address=$PWD
    mkdir -p ~/.comparewithlmp/lsmdtest
    rm -f ~/.comparewithlmp/* || true
    cp ${lsmdmodelfile} ~/.comparewithlmp/lsmdtest

    # 开启全部的opls力的计算
    cat << EOF > ~/.comparewithlmp/lsmdtest/md.inp
ZOOM         1
OPLSlmp  ${lsmdmodelfile}

Potential         opls
nobond_constant  0 0 0.5  10  0.25   
opls_start    1 1 1 0 1    
UnitStyle         real
T_init            300
outvirial  1  virialall.txt   

ensemble       NVE
Timestep        0
writedata   1   xyz   allforce.xyz      
MDRUN           1

Build_NEIGHBOR  2            
build_neighbor_interval 10   
USEOMP   0                   
NEIGHBORSKIN  1              
Maxbonding   6               
compute_virial  1            
TIMELOG     0                
OUTMOLECULEVIRIAL   1        
SKIPCOULOMB  1               
OUTRAJ      2                
EOF
    cd ~/.comparewithlmp/lsmdtest/
    lsmd

    #只开启bond的计算
    sed -E -i '/^opls_start/s/.*/opls_start    1 0 0 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialbond.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   bondforce.xyz      /' md.inp
    lsmd

    #只开启angle的计算
    sed -E -i '/^opls_start/s/.*/opls_start    0 1 0 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialangle.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   angleforce.xyz      /' md.inp
    lsmd

    #只开启dihedral的计算
    sed -E -i '/^opls_start/s/.*/opls_start    0 0 1 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialdihedral.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   dihedralforce.xyz      /' md.inp
    lsmd

    #只开启pair的计算
    sed -E -i '/^opls_start/s/.*/opls_start    0 0 0 0 1/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialpair.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   pairforce.xyz      /' md.inp

    # 开启全部的opls力的计算
    sed -E -i '/^opls_start/s/.*/opls_start    1 1 1 0 1/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialall.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   allforce.xyz      /' md.inp
    lsmd

    cd ${lsmdmodelfile_address}
    echo "${lsmdmodelfile} has been tested."
}

clearcompare(){
    rm -rf ~/.comparewithlmp/lmptest || true
    rm -rf ~/.comparewithlmp/lsmdtest || true
}

cdcomparedir(){
    cd ~/.comparewithlmp/
}

cdlmpdir(){
    if [ -z "$lmpmodelfile_address" ]; then
        echo "Error: lmpmodelfile_address is not set."
        return 1
    fi
    cd ${lmpmodelfile_address}
}

# 在当前的文件夹下测试lsmd
testlsmdir(){

    lsmdmodelfile=${1}
    lsmdtemplate=${2:-300}  # 默认为300K
    lsmdmodelfile_address=$PWD

    # 开启全部的opls力的计算
    cat << EOF > md.inp
ZOOM         1
OPLSlmp  ${lsmdmodelfile}

Potential         opls
nobond_constant  0 0 0.5  10  0.25   
opls_start    1 1 1 0 1    
UnitStyle         real
T_init            ${lsmdtemplate}
outvirial  1  virialall.txt   

ensemble       NVE
Timestep        0
writedata   1   xyz   allforce.xyz      
MDRUN           1

Build_NEIGHBOR  2            
build_neighbor_interval 10   
USEOMP   0                   
NEIGHBORSKIN  1              
Maxbonding   6               
compute_virial  1            
TIMELOG     0                
OUTMOLECULEVIRIAL   1        
SKIPCOULOMB  1               
OUTRAJ      2                
EOF
    lsmd

    #只开启bond的计算
    sed -E -i '/^opls_start/s/.*/opls_start    1 0 0 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialbond.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   bondforce.xyz      /' md.inp
    lsmd

    #只开启angle的计算
    sed -E -i '/^opls_start/s/.*/opls_start    0 1 0 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialangle.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   angleforce.xyz      /' md.inp
    lsmd

    #只开启dihedral的计算
    sed -E -i '/^opls_start/s/.*/opls_start    0 0 1 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialdihedral.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   dihedralforce.xyz      /' md.inp
    lsmd

    #只开启pair的计算
    sed -E -i '/^opls_start/s/.*/opls_start    0 0 0 0 1/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialpair.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   pairforce.xyz      /' md.inp
    lsmd

    # 开启全部的opls力的计算
    sed -E -i '/^opls_start/s/.*/opls_start    1 1 1 0 1/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialall.txt/' md.inp
    sed -E -i '/^writedata/s/.*/writedata   1   xyz   allforce.xyz      /' md.inp
    lsmd

    cd ${lsmdmodelfile_address}
    echo "${lsmdmodelfile} has been tested."
}

clearcompare(){
    rm -rf ~/.comparewithlmp/lmptest || true
    rm -rf ~/.comparewithlmp/lsmdtest || true
}

cdcomparedir(){
    cd ~/.comparewithlmp/
}

cdlmpdir(){
    if [ -z "$lmpmodelfile_address" ]; then
        echo "Error: lmpmodelfile_address is not set."
        return 1
    fi
    cd ${lmpmodelfile_address}
}

# 将lammps的metal单位转化为real单位，单位转换因子为23.060547830619
# 用法示例
# metal2real 125.3725 154.4014 25.6099 0.0000 -132.9088 -1429.0067        
# 2891.0898 3560.4963 590.5643 0.0000 -3064.8769 -32952.8945 
metal2real() {
    echo
    echo "metal2real:"
    for val in "$@"; do
        awk -v v="$val" 'BEGIN{printf "%.4f ", v*23.060547830619}'
    done
}
