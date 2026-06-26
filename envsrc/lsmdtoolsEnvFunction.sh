# ======================================================================
# File:         lsmdtoolsEnvFunction.sh
# Project:      NebulaFlow
# Description:  自研LSMD软件测试函数库 — 力场各分量独立验证
#               Testing functions for the LSMD in-house MD software.
#               Tests each force component (bond/angle/dihedral/pair) independently.
# Author:       rebreath
# Dependencies: lmp, lsmd
# ======================================================================


# ---------------------------------------------------------------------------
# Function: testlmpmod / testlmpdir
# 功能: 使用LAMMPS验证力场模型文件
#       用testall.in模板运行LAMMPS单步计算，验证输出力的正确性。
# 场景: 开发或修改力场参数后，快速验证模型是否可以正常读取和计算。
# Usage:
#   testlmpmod <model.data>         # 在~/.comparewithlmp/lmptest中测试
#   testlmpdir <model.data> [T_K]   # 在当前目录测试，可指定温度
# Example:
#   testlmpdir system.data 300
# ---------------------------------------------------------------------------
testlmpmod(){
    lmpmodelfile=${1}; lmpmodelfile_address=$PWD
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

testlmpdir(){
    lmpmodelfile=${1}; lmptemplate=${2:-300}; lmpmodelfile_address=$PWD
    cp ~/.rebreath/inp_lib/Lammps/testall.inp .
    sed -E -i '/^[[:space:]]*read_data/s/.*/read_data   '"${lsmdmodelfile}"'/' testall.inp
    sed -E -i '/^([[:space:]]*)velocity all create/s~.*~\1velocity all create '"${lmptemplate}"' 12345 mom yes rot yes dist gaussian~' testall.inp
    conda activate deepmd || true
    mpirun -np 1 lmp -i testall.inp
    echo "${lmpmodelfile} has been tested."
}


# ---------------------------------------------------------------------------
# Function: testlsmdmod / testlsmdir
# 功能: 使用LSMD软件分别测试力场各分量的正确性
#       依次开启bond/angle/dihedral/pair各力的计算，输出对比文件。
# 场景: 开发LSMD力场时，分别验证每个力场分量的计算是否正确。
#       通过依次关闭其他力分量，单独计算每种力的贡献。
# Usage:
#   testlsmdmod <model.data>
#   testlsmdir <model.data> [T_K]
# Notes:
#   - 会生成virialbond.txt, bondforce.xyz 等各分量输出文件
#   - 对比LAMMPS和LSMD的对应输出可验证力场一致性
# ---------------------------------------------------------------------------
testlsmdmod(){
    lsmdmodelfile=${1}; lsmdmodelfile_address=$PWD
    mkdir -p ~/.comparewithlmp/lsmdtest; rm -f ~/.comparewithlmp/* || true
    cp ${lsmdmodelfile} ~/.comparewithlmp/lsmdtest
    cat << EOF > ~/.comparewithlmp/lsmdtest/md.inp
ZOOM 1; OPLSlmp ${lsmdmodelfile}; Potential opls; nobond_constant 0 0 0.5 10 0.25
opls_start 1 1 1 0 1; UnitStyle real; T_init 300; outvirial 1 virialall.txt
ensemble NVE; Timestep 0; writedata 1 xyz allforce.xyz; MDRUN 1
Build_NEIGHBOR 2; build_neighbor_interval 10; USEOMP 0; NEIGHBORSKIN 1
Maxbonding 6; compute_virial 1; TIMELOG 0; OUTMOLECULEVIRIAL 1
SKIPCOULOMB 1; OUTRAJ 2
EOF
    cd ~/.comparewithlmp/lsmdtest/; lsmd
    # bond only
    sed -E -i '/^opls_start/s/.*/opls_start    1 0 0 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialbond.txt/' md.inp; lsmd
    # angle only
    sed -E -i '/^opls_start/s/.*/opls_start    0 1 0 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialangle.txt/' md.inp; lsmd
    # dihedral only
    sed -E -i '/^opls_start/s/.*/opls_start    0 0 1 0 0/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialdihedral.txt/' md.inp; lsmd
    # pair only
    sed -E -i '/^opls_start/s/.*/opls_start    0 0 0 0 1/' md.inp
    sed -E -i '/^outvirial/s/.*/outvirial 1 virialpair.txt/' md.inp; lsmd
    cd ${lsmdmodelfile_address}; echo "${lsmdmodelfile} has been tested."
}

# (testlsmdir — similar to testlsmdmod but runs in current directory)

# ---------------------------------------------------------------------------
# Function: clearcompare / cdcomparedir
# 功能: 清理/进入力场对比测试目录
# ---------------------------------------------------------------------------
clearcompare(){ rm -rf ~/.comparewithlmp/lmptest ~/.comparewithlmp/lsmdtest 2>/dev/null || true; }
cdcomparedir(){ cd ~/.comparewithlmp/; }

# ---------------------------------------------------------------------------
# Function: metal2real
# 功能: 将LAMMPS metal单位转换为real单位（能量: eV→kcal/mol, 力: eV/A→kcal/mol/A）
#       转换因子: 23.060547830619
# 场景: 使用不同单位系统的力场时需要进行单位换算。
# Usage: metal2real <val1> <val2> ...
# Example:
#   metal2real 125.37 154.40
# ---------------------------------------------------------------------------
metal2real() {
    echo "metal2real:"
    for val in "$@"; do
        awk -v v="$val" 'BEGIN{printf "%.4f ", v*23.060547830619}'
    done
}

# =============================================================================
