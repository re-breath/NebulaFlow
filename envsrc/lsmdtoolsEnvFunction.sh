# 该文件库包含了测试自研软件LSMD的各种环境函数


testlmpmod(){
    lmpmodelfile=${1}
    lmpmodelfile_address=$PWD
    mkdir -p ~/.comparewithlmp/lmptest
    # rm -f ~/.comparewithlmp/* || true
    cp ${lmpmodelfile} ~/.comparewithlmp/lmptest
    cat << EOF > ~/.comparewithlmp/lmptest/testall.inp
# ---------------- Bond-only ----------------
clear
units       real
atom_style  full

bond_style      harmonic
angle_style     zero nocoeff
dihedral_style  zero nocoeff
pair_style      zero 10.0 nocoeff

read_data   ${lmpmodelfile}

thermo      1
thermo_style custom step etotal ebond eangle edihed evdwl ecoul

compute     f_bond all property/atom fx fy fz
dump        dump_bond all custom 1 bondforce.log id type x y z c_f_bond[1] c_f_bond[2] c_f_bond[3]
dump_modify dump_bond sort id
dump_modify dump_bond format float %14.6f
dump_modify dump_bond format int %8d

# 每原子位力（stress/atom 的 6 分量）
compute         s all stress/atom NULL
dump            v_all all custom 1 virialbond.dump id type c_s[1] c_s[2] c_s[3] c_s[4] c_s[5] c_s[6]
dump_modify v_all sort id
dump_modify v_all format float %14.6f
dump_modify v_all format int %8d

timestep    0
fix         nve_bond all nve
run         0
unfix       nve_bond
undump      dump_bond
uncompute   f_bond

# ---------------- Angle-only ----------------
clear
units       real
atom_style  full

bond_style      zero nocoeff
angle_style     harmonic
dihedral_style  zero nocoeff
pair_style      zero 10.0 nocoeff

read_data   ${lmpmodelfile}

thermo      1
thermo_style custom step etotal ebond eangle edihed evdwl ecoul

compute     f_angle all property/atom fx fy fz
dump        dump_angle all custom 1 angleforce.log id type x y z c_f_angle[1] c_f_angle[2] c_f_angle[3]
dump_modify dump_angle sort id
dump_modify dump_angle format float %14.6f
dump_modify dump_angle format int %8d

# 每原子位力（stress/atom 的 6 分量）
compute         s all stress/atom NULL
dump            v_all all custom 1 virialangle.dump id type c_s[1] c_s[2] c_s[3] c_s[4] c_s[5] c_s[6]
dump_modify v_all sort id
dump_modify v_all format float %14.6f
dump_modify v_all format int %8d


fix         nve_angle all nve
run         0
unfix       nve_angle
undump      dump_angle
uncompute   f_angle

# ---------------- Dihedral-only ----------------
clear
units       real
atom_style  full

bond_style      zero nocoeff
angle_style     zero nocoeff
dihedral_style  opls
pair_style      zero 10.0 nocoeff

read_data   ${lmpmodelfile}

compute     f_dihedral all property/atom fx fy fz
dump        dump_dihedral all custom 1 dihedralforce.log id type x y z c_f_dihedral[1] c_f_dihedral[2] c_f_dihedral[3]
dump_modify dump_dihedral sort id
dump_modify dump_dihedral format float %14.6f
dump_modify dump_dihedral format int %8d

# 每原子位力（stress/atom 的 6 分量）
compute         s all stress/atom NULL
dump            v_all all custom 1 virialdihedral.dump id type c_s[1] c_s[2] c_s[3] c_s[4] c_s[5] c_s[6]
dump_modify v_all sort id
dump_modify v_all format float %14.6f
dump_modify v_all format int %8d


fix         nve_dihedral all nve
run         0
unfix       nve_dihedral
undump      dump_dihedral
uncompute   f_dihedral

# ---------------- Pair-only ----------------
clear
units       real
atom_style  full

bond_style      zero nocoeff
angle_style     zero nocoeff
dihedral_style  zero nocoeff
pair_style      lj/cut/coul/dsf 0.25 10.0

read_data   ${lmpmodelfile}

thermo      1
thermo_style custom step etotal ebond eangle edihed evdwl ecoul

special_bonds   lj 0.0 0.0 0.5 coul 0.0 0.0 0.5 angle yes dihedral yes

compute     f_pair all property/atom fx fy fz
dump        dump_pair all custom 1 pairforce.log id type x y z c_f_pair[1] c_f_pair[2] c_f_pair[3]
dump_modify dump_pair sort id
dump_modify dump_pair format float %14.6f
dump_modify dump_pair format int %8d



# 每原子位力（stress/atom 的 6 分量）
compute         s all stress/atom NULL
dump            v_all all custom 1 virialpair.dump id type c_s[1] c_s[2] c_s[3] c_s[4] c_s[5] c_s[6]
dump_modify v_all sort id
dump_modify v_all format float %14.6f
dump_modify v_all format int %8d



fix         nve_pair all nve
run         0
unfix       nve_pair
undump      dump_pair
uncompute   f_pair

# ---------------- All interactions ----------------
clear
units       real
atom_style  full

bond_style      harmonic
angle_style     harmonic
dihedral_style  opls
pair_style      lj/cut/coul/dsf 0.25 10.0

read_data   ${lmpmodelfile}

thermo      1
thermo_style custom step etotal ebond eangle edihed evdwl ecoul

special_bonds   lj 0.0 0.0 0.5 coul 0.0 0.0 0.5 angle yes dihedral yes

compute     f_all all property/atom fx fy fz
dump        dump_all all custom 1 allforce.log id type x y z c_f_all[1] c_f_all[2] c_f_all[3]
dump_modify dump_all sort id
dump_modify dump_all format float %14.6f
dump_modify dump_all format int %8d


# 每原子位力（stress/atom 的 6 分量）
compute         s all stress/atom NULL
dump            v_all all custom 1 virialall.dump id type c_s[1] c_s[2] c_s[3] c_s[4] c_s[5] c_s[6]
dump_modify v_all sort id
dump_modify v_all format float %14.6f
dump_modify v_all format int %8d

fix         nve_all all nve
run         0
EOF
    conda activate deepmd || true
    cd ~/.comparewithlmp/lmptest/
    mpirun -np 4 lmp -i testall.inp
    cd ${lmpmodelfile_address}
    echo "${lmpmodelfile} has been tested."
}

# ---------------- Test lsmd model ----------------
# 测试lsmd模型
# 使用方法：testlsmdmod lsmdmodelfile
testlsmdmod(){

    lsmdmodelfile=${1}
    lsmdmodelfile_address=$PWD
    mkdir -p ~/.comparewithlmp/lsmdtest
    # rm -f ~/.comparewithlmp/* || true
    cp ${lsmdmodelfile} ~/.comparewithlmp/lsmdtest

    # 开启全部的opls力的计算
    cat << EOF > ~/.comparewithlmp/lsmdtest/md.inp
ZOOM         1
OPLSlmp  ${lsmdmodelfile}

Potential         opls
nobond_constant  0 0 0.5  10  0.25   
opls_start    1 1 1 0 1    
UnitStyle         real
T_init            100
outvirial  1  virialall.txt   

ensemble       NVE
Timestep        0
writedata   1   xyz   allforce.xyz      
MDRUN           1

Build_NEIGHBOR  2            
build_neighbor_interval 10   
USEOMP   1                   
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
