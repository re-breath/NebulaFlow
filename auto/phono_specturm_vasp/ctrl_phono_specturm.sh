#!/bin/bash
#from rebreath
#本脚本为计算声子谱控制台，需要提前准备好INCAR（结构优化的），INCAR-new（计算声子谱需要的），POSCAR
#本脚本的使用需要对整个计算过程有一定了解，能将所有脚本的文件参数进行修改
#如果你做到上面的，只需运行本脚本，即可一键使用vasp的有限位移法计算声子谱
#注意初始的INCAR文件为结构优化的，INCAR-new为计算声子谱的输入文件，band.conf需要进行一些修改

#**注意：扩胞默认为3 3 3，需要对本文件与band.conf进行修改

#----------------------controls------------------------
dim=3
codenum=1
deform_cell_num=0.998   #该参数为phonopy扩胞后再对其扩胞后的POSCAR进行拉伸，防止出现虚频
#-------------------------------------------------------
set e
init_address=$PWD
library_address='/$HOME/.rebreath/auto/phono_specturm_vasp'
incarfiles="INCAR INCAR-new band.conf"
for i in $incarfiles
do
  if [ ! -f "$PWD/$i" ];then
    echo "当前文件夹没有找到$i，正在复制脚本文件..."
    cp $library_address/$i  $init_address
  fi
done

all_sh_files="hmidfy.sh  mk_doc.sh  name_mul.sh  name.sh  sbatch.sh"
for i in $all_sh_files
do
  if [ ! -f "$init_address/$i" ];then
    echo "当前文件夹没有找到$i，正在复制脚本文件..."
    cp $library_address/$i  $init_address
  fi
done

#function zone
prevention(){
    dos2unix * || true
    conda init
    conda activate 4090 || true
}

contend_cell_and_rename(){
#检查phonopy是否存在并进行扩胞
bela=$(which phonopy)
if [ -z "$bela" ];then
  echo "当前没有找到phonopy，正在安装..."
  pip install phonopy
  if [ $? -ne 0 ];then
    echo "安装失败，请手动安装phonopy"
    exit 1
  fi
  echo "安装完成"
fi
phonopy --amplitude=0.005  -d --dim='3 3 3'
echo "已完成扩胞"
mv INCAR INCAR-old
mv INCAR-new INCAR
mv POSCAR POSCAR-unitcell
echo "已完成job需要的文件"
}

reverse_incar(){
  mv INCAR INCAR-new
  mv INCAR-old INCAR
}

deform_cell(){
  for i in $(find $init_address -type d -regex ".*job-[0-9]+");do  sed -i "2s/.*/$deform_cell_num/" $i/POSCAR;done
}

prevention
#source ./hmidfy.sh
echo "已完成高精度优化"
contend_cell_and_rename
source ./mk_doc.sh
echo "已完成分发job需要的输入文件"
source ./name.sh
deform_cell
source ./sbatch.sh
echo "已成功完成声子谱的计算任务"
reverse_incar
mv vasprun.xml vasprun..xml
source ./name_mul.sh
mv vasprun..xml vasprun.xml
phonopy -f vasprun.xml-0*
phonopy -c CONTCAR  $init_address/band.conf -s -p

phonopy-bandplot --gnuplot > phonon.out

echo "完美完成任务>.<"

