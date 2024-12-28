#!/bin/bash
# 该脚本使用来完成vasp的三个阶段的结构优化，并且进行单点能的计算

source ~/.rebreath/rebreath-env-function

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

