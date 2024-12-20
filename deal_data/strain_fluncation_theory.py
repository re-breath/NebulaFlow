# 该脚本使用应变波动法计算弹性模量
# rebreath练习使用版
# 计算晶系为 cublic晶系

import numpy as np
import pandas as pd

T = 1800
def read_thermo(filename):
    """
    读取gpumd的输出文件thermo.out
    """
    data = np.loadtxt(filename)
    thermo = dict()
    if data.shape[1] == 12:
        labels = ['T','K','U','Px','Py','Pz','Pyz','Pxz','Pxy','Lx','Ly','Lz']
        for i in range(12):
            thermo[labels[i]] = data[:,i]
    if data.shape[1] == 18:
        labels = ['T','K','U','Px','Py','Pz','Pyz','Pxz','Pxy','ax','ay','az','bx','by','bz','cx','cy','cz']
        for i in range(18):
            thermo[labels[i]] = data[:,i]
    return thermo
def calc_angle(v1,v2):
    """
    该函数使用来计算两个矢量之间的夹角
    """
    angle = np.arccos(np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2))
    return angle
def calc_elastic(strain_ij,strain_kl):
    """
    该函数通过应变波动法来计算弹性模量
    """
    return np.average(strain_ij*strain_kl)-np.average(strain_ij)*np.average(strain_kl)  
def get_Lattice_from_file(filename):
    """
    该函数可以通过读取其他的文件来获得晶格常数
    该文件的格式需要为9列，分别为ax,ay,az,bx,by,bz,cx,cy,cz
    """
    data = np.loadtxt(filename)
    label = ['ax','ay','az','bx','by','bz','cx','cy','cz']
    lattice = dict()
    for i in range(9):
        lattice[label[i]] = data[i]
    return lattice
def get_Lattice_from_thermo():
    """
    该函数使用来获取thermo字典中的晶格常数
    """
    thermo = read_thermo('thermo.out')
    label = ['ax','ay','az','bx','by','bz','cx','cy','cz']
    lattice = dict()
    for i in label:
        lattice[i] = thermo[i]
    return lattice

lattice = get_Lattice_from_thermo()
#print(lattice)
def compute_cubic_elastic(lattice):
    """
    通过应变波动法计算弹性模量
    """
    slice_num = 10
    length = len(lattice['ax'])
    #sliced = int(length/3)
    for i in lattice.keys():
        lattice[i] = lattice[i][1000:]
    lattice_all = lattice.copy()
    Cij = np.zeros((slice_num,3))
    for time in range(slice_num):
        for key in lattice.keys():
            lattice[key] = lattice_all[key][time*1000:(time+1)*1000]
        # 计算二面角
        va = np.column_stack((lattice['ax'],lattice['ay'],lattice['az']))
        vb = np.column_stack((lattice['bx'],lattice['by'],lattice['bz']))
        vc = np.column_stack((lattice['cx'],lattice['cy'],lattice['cz']))

        alpha = np.zeros(len(va))
        beta = np.zeros_like(alpha)
        gamma = np.zeros_like(alpha)
        for i in range(len(va)):
            alpha[i] = calc_angle(vb[i],vc[i]) #yz面
            beta[i] = calc_angle(vc[i],va[i])  #zx面
            gamma[i] = calc_angle(va[i],vb[i]) #xy面
        
        #print(f"alpha:{alpha[0]}")
        # 计算应变
        strain = dict()
        strain['11'] = lattice['ax']/np.average(lattice['ax'])-1
        strain['22'] = lattice['by']/np.average(lattice['by'])-1
        strain['33'] = lattice['cz']/np.average(lattice['cz'])-1
        strain['23'] = strain['32'] = (alpha - np.pi/2)/2 
        strain['13'] = (beta - np.pi/2)/2
        strain['12'] = (gamma - np.pi/2)/2
        # df = pd.DataFrame(strain)
        # df.to_csv(f"strain.csv",mode = 'a',sep='\t',header=True)
        # strain_array = df.values
        # with open(f"strain.txt",mode = 'a') as f:
        #     np.savetxt(f,strain_array)

        # 计算应变波动
        volum = np.average(lattice_all['ax']*lattice_all['by']*lattice_all['cz']) # 这里的体积使用的是总的体积的平均值
        #print(f"V:{volum}")
        KB = 1.38064852
        scale = 100/(T*KB)*volum
        # 计算弹性模量的逆矩阵
        S1111 = scale*calc_elastic(strain["11"], strain["11"])
        S2222 = scale*calc_elastic(strain["22"], strain["22"])
        S3333 = scale*calc_elastic(strain["33"], strain["33"])
        S1122 = scale*calc_elastic(strain["11"], strain["22"])
        S1133 = scale*calc_elastic(strain["11"], strain["33"])
        S2233 = scale*calc_elastic(strain["22"], strain["33"])
        S2323 = scale*calc_elastic(strain["23"], strain["23"])
        S1313 = scale*calc_elastic(strain["13"], strain["13"])
        S1212  = scale*calc_elastic(strain["12"], strain["12"])

        S1123 = scale*calc_elastic(strain["11"], strain["23"])
        S1113 = scale*calc_elastic(strain["11"], strain["13"])
        S1112 = scale*calc_elastic(strain["11"], strain["12"])
        S1112 = scale*calc_elastic(strain["11"], strain["12"])

        S11 = (S1111+S2222+S3333)/3
        S12 = (S1122+S1133+S2233)/3
        S44 = 4*(S2323+S1313+S1212)/3
        Spq = np.array([[S11, S12, S12,   0,   0,   0], 
                        [S12, S11, S12,   0,   0,   0],
                        [S12, S12, S11,   0,   0,   0],
                        [0,     0,   0, S44,   0,   0],
                        [0,     0,   0,   0, S44,   0],
                        [0,     0,   0,   0,   0, S44]])
            
        # Convert Spq to Cpq and 
        Cpq = np.linalg.inv(Spq)
        #print(Cpq)
        C11 = Cpq[1,1]
        C12 = Cpq[1,2]
        C44 = Cpq[4,4]
        Cij[time] = np.array([C11, C12, C44])
    return Cij


def get_cubic_modulus(cij):
    """
    计算bulk modulus，shear modulus，Young's modulus
    """
    B=(cij[0]+2*cij[1])/3
    print("The bulk modulus of Si estimated by gpumd with npt_scr for 1ns:")
    print(B)
    G= cij[2]
    print("The shear modulus of Si estimated by gpumd with npt_scr for 1ns:")
    print(G)
    E=9*B*G/(3*B+G)
    print("The Young's modulus of Si estimated by gpumd with npt_scr for 1ns:")
    print(E)



# 注意该处计算的为立方晶系
Cij = compute_cubic_elastic(lattice)
average_cij = np.average(Cij,axis=0)
print(f"Cij : {average_cij}")

