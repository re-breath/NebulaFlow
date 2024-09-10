#!/bin/bash
#本脚本致力于二维材料与三维材料进行高精度优化

#基础设定控制台（只需要对此进行更改，其他模块无需更改）
init_address=$(pwd)   #保存初始位置，无需更该
codenum=${codenum:=1}     #核数（cpu跑就填cpu核数，gpu跑就填gpu核数）
dim=${dim:=3}        #维数：块状材料（bulk）就填3，二维材料（2D）填2


set -e

######————————下面为工作区——————#################

#1.自动生成KPOINTS与POTCAR文件
#echo -e "102 \n1 \n0.04 " |vaspkit

#2.第一次优化
sed -i '/EDIFFG/ s/.*/#EDIFFG=1e-2/' INCAR
case $dim in
  2)
    sed -i '/ISIF/ s/.*/ISIF=4/' INCAR
    ;;
  3)
    sed -i '/ISIF/ s/.*/ISIF=3/' INCAR
    ;;
  *)
    echo "Error: dim must be 2 or 3"
    exit 1
    ;;
esac
sed -i '/IBRION/ s/.*/IBRION=1/' INCAR

mpirun -np $codenum vasp
mv POSCAR POSCAR-init
cp CONTCAR POSCAR

 #3.第二次优化
 sed -i '/EDIFFG/ s/.*/EDIFFG=1e-3/' INCAR
 case $dim in
  2)
    sed -i '/ISIF/ s/.*/ISIF=2/' INCAR
    ;;
  3)
    sed -i '/ISIF/ s/.*/ISIF=3/' INCAR
    ;;
  *)
    echo "Error: dim must be 2 or 3"
    exit 1
    ;;
esac
 sed -i '/IBRION/ s/.*/IBRION=1/' INCAR

 mpirun -np $codenum vasp
cp CONTCAR POSCAR

 #4.第三次优化
 sed -i '/EDIFFG/ s/.*/EDIFFG=1e-4/' INCAR
 sed -i '/ISIF/ s/.*/ISIF=2/' INCAR
 sed -i '/IBRION/ s/.*/IBRION=1/' INCAR


mpirun -np $codenum vasp
cp CONTCAR POSCAR


  #5.第四次优化
 sed -i '/EDIFFG/ s/.*/EDIFFG=1e-5/' INCAR
 sed -i '/ISIF/ s/.*/ISIF=2/' INCAR
 sed -i '/IBRION/ s/.*/IBRION=1/' INCAR

 mpirun -np $codenum vasp
 cp CONTCAR POSCAR

  #6.第五次优化
 sed -i '/EDIFFG/ s/.*/EDIFFG=1e-6/' INCAR
 sed -i '/ISIF/ s/.*/ISIF=2/' INCAR
 sed -i '/IBRION/ s/.*/IBRION=1/' INCAR


 mpirun -np $codenum vasp
 cp CONTCAR POSCAR
    #7.第六次优化
  sed -i '/EDIFFG/ s/.*/EDIFFG=1e-7/' INCAR
  sed -i '/ISIF/ s/.*/ISIF=2/' INCAR
  sed -i '/IBRION/ s/.*/IBRION=1/' INCAR


  mpirun -np $codenum vasp
cp CONTCAR POSCAR
   #8.第七次优化
  sed -i '/EDIFFG/ s/.*/EDIFFG=1e-8/' INCAR
  sed -i '/ISIF/ s/.*/ISIF=2/' INCAR
  sed -i '/IBRION/ s/.*/IBRION=1/' INCAR


  mpirun -np $codenum vasp
cp CONTCAR POSCAR

   #9.第八次优化
  sed -i '/EDIFFG/ s/.*/EDIFFG=1e-14/' INCAR
  sed -i '/ISIF/ s/.*/ISIF=2/' INCAR
  sed -i '/IBRION/ s/.*/IBRION=1/' INCAR


  mpirun -np $codenum vasp
cp CONTCAR POSCAR

echo " >>hmidfy.sh已完成十层楼高的精度优化O.<"

