#!/bin/bash
#from rebreath
#本脚本作用：自动进行三个方向的应力应变曲线的整合，自动识别拉伸方向与该文件夹的工作
#使用方法：source auto_plot_xyz_strain_stress_curve.sh为了使用环境变量中的函数

init_address=$PWD
init_dir_name=$(basename $PWD)

mkdir -p strain_stress_curve

#Identify the stretching direction
idfy_stretch_direction(){
    #识别拉伸方向并将拉伸的文件进行改名
    dirname=0
    dir_name=$(basename $PWD)
    if [[ $dir_name =~ _([xyz]) ]];then
        stretch_direction=${BASH_REMATCH[1]}
    fi
    cp  stress_strain_curve.txt  strain_stress_curve_${stretch_direction}.txt
}

bulid_plot_stretch_warehouse(){
    for i in $(find $init_address -type d -regex ".*deform.*_[xyz].*")
    do
        cd $i
        echo "当前工作目录为$i"
        #使用rebreath函数库的函数
        plot_stress_strain_curve
        idfy_stretch_direction
        cp strain_stress_curve_*.txt  $init_address/strain_stress_curve/
        cd $init_address
    done
}

deform_string(){
    #进行字符串的规整化
    init_str=$1
    step1=${init_str/run_/ }
    step2=${step1/change/change_}
    step3=${step2/dumpto/_to_}
    step4=${step3/4800/from_4800k}
    finish_str=${step4//_/ }
    echo $finish_str
}



plot_xyz_strain_stress_curve(){
    #进行三图合一
    python3 <<EOF
import matplotlib.pyplot as plt
import numpy as np

aw = 2
fs = 16
font = {'size'   : fs}
plt.rc('font', **font)
plt.rc('axes' , linewidth=aw)


def nonuniform_sampling(file,start,end,step):
    #进行非均匀的采样
    data=np.loadtxt(file)
    date1=data[:start]
    date2=data[start:end:step]
    date3=data[end:]
    return np.concatenate((date1,date2,date3))

def plot_strain_stress(data,pngname='strain_stress_curve.png'):
    #对二维列表进行画图
    plt.figure(figsize=(10,6))
    plt.title('strain-stress curve')
    plt.xlabel('strain')
    plt.ylabel('stress')
    plt.plot(data[:,0],data[:,1],'o-')
    plt.savefig(pngname)

def export_data(data,file):
    np.savetxt(file,data)


data_x=np.loadtxt('strain_stress_curve_x.txt')
data_y=np.loadtxt('strain_stress_curve_y.txt')
data_z=np.loadtxt('strain_stress_curve_z.txt')


plt.figure(figsize=(10,6))
plt.title('strain-stress curve ${deformed_name}')
plt.xlabel('strain')
plt.ylabel('stress(Gpa)')

plt.plot(data_x[:,0],data_x[:,1],'o-',label='X-direction stretching')
plt.plot(data_y[:,0],data_y[:,1],'o-',label='Y-direction stretching')
plt.plot(data_z[:,0],data_z[:,1],'o-',label='Z-direction stretching')
plt.legend()
plt.savefig('strain-stress curve ${deformed_name}.png')
EOF
}


deformed_name="$(deform_string $init_dir_name)"
bulid_plot_stretch_warehouse
cd strain_stress_curve/
plot_xyz_strain_stress_curve
cd $init_address
