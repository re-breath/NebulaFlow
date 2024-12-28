#!/bin/bash
# 该脚本使用来完成vasp的三个阶段的结构优化，并且进行单点能的计算

source ~/.rebreath/rebreath-env-function


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
 SIGMA  =  0.1        (please check the width of the smearing)
 LORBIT =  11           (PAW radii for projected DOS)
 NEDOS  =  2001         (DOSCAR points)
 #NELM   =  260           (Max electronic SCF steps)
 NSW    =   30
 EDIFF  =  1E-02      (SCF energy convergence, in eV)
 IVDW = 12
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
 SIGMA  =  0.1        (please check the width of the smearing)
 LORBIT =  11           (PAW radii for projected DOS)
 NEDOS  =  2001         (DOSCAR points)
 #NELM   =  260           (Max electronic SCF steps)
 NSW    =   100
 EDIFF  =  1E-04     (SCF energy convergence, in eV)
 IVDW = 12
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
 SIGMA  =  0.1        (please check the width of the smearing)
 LORBIT =  11           (PAW radii for projected DOS)
 NEDOS  =  2001         (DOSCAR points)
 #NELM   =  260           (Max electronic SCF steps)
 NSW    =   100
 EDIFF  =  1E-04     (SCF energy convergence, in eV)
 IVDW = 12
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
IVDW = 12
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
SIGMA  =  0.1        
LORBIT =  11           
NEDOS  =  2001        
NELM   =  120         
EDIFF  =  1E-05        
IVDW = 12 
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

vaspstart_geo_optstage1
wait_complete_vasp_run  'pwd'
mkdir -p ../stage2

cp CONTCAR  ../stage2/POSCAR

cd ../stage2
vaspstart_geo_optstage2
wait_complete_vasp_run  'pwd'

mkdir -p ../stage3
cp CONTCAR  ../stage3/POSCAR

cd ../stage3
vaspstart_geo_optstage3
wait_complete_vasp_run 'pwd'

vaspstart_single_energy_after_relaxation

