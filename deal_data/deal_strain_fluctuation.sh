#!/bin/bash
#本脚本用于辅助处理应变波动法计算弹性模量
python3 << EOF
from thermo.gpumd.data import load_thermo
from pylab import *
import math
import numpy as np
def thermo_slic(filepath1,filepath2):
    thermo=np.loadtxt(filepath1)
    thermo=np.delete(thermo, [6,7,8], axis=1)
    np.savetxt(filepath2, thermo)
    
def calc_angle(v1, v2):
    alpha = math.acos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    return alpha

def calc_elstic(strain_ij, strain_kl):
    s_ijkl = average(strain_ij*strain_kl)-average(strain_ij)*average(strain_kl)
    return s_ijkl

def output_elastic_gpumd(slice_num):
    #Obtain the strain tensor     directory=path,
    thermo = load_thermo(filename = 'thermo_slic.out')
    for keys in thermo:
        thermo[keys]=thermo[keys][1000:] #Discard the first 100 ps data
    Cij = np.zeros((slice_num,3))
    for i in range(slice_num): #Split to 10 slices
        # print(i)
        thermo_slice = dict()
        for keys in thermo:
            thermo_slice[keys]=thermo[keys][1000*i:1000*(i+1)]
        alpha = np.zeros(len(thermo_slice["ax"]))
        beta  = np.zeros_like(alpha)
        gamma = np.zeros_like(alpha)
        for j in range(len(thermo_slice["ax"])):
            va = [thermo_slice["ax"][j], thermo_slice["ay"][j], thermo_slice["az"][j]]
            vb = [thermo_slice["bx"][j], thermo_slice["by"][j], thermo_slice["bz"][j]]
            vc = [thermo_slice["cx"][j], thermo_slice["cy"][j], thermo_slice["cz"][j]]
            alpha[j] = calc_angle(vb, vc)
            beta[j]  = calc_angle(va, vc)
            gamma[j] = calc_angle(va, vb)
           
        strain = dict()
        strain["11"] = thermo_slice["ax"]/average(thermo_slice["ax"])-1 #xx 
        strain["22"] = thermo_slice["by"]/average(thermo_slice["by"])-1 #yy
        strain["33"] = thermo_slice["cz"]/average(thermo_slice["cz"])-1 #zz
        strain["23"] = (alpha-pi/2)/2 #yz
        strain["13"] = (beta-pi/2)/2 #xz
        strain["12"] = (gamma-pi/2)/2 #xy
        
        #Calculate the cubic compliance tensor
        #S11=S1111=S2222=S3333,S12=S1122=S1133=S2233; S44=S2323=S1313=S1212
        V = average(thermo["ax"]*thermo["by"]*thermo["cz"]) #volue, unit in angstrom^{3}
        T = $1 #temperature, unit in K
        kB = 1.38064852 #unit in e-23 J/K
        scale = 100/(T*kB)*V #unit in GPa^{-1}, scale=V/(KB*T)
        S1111 = scale*calc_elstic(strain["11"], strain["11"])
        S2222 = scale*calc_elstic(strain["22"], strain["22"])
        S3333 = scale*calc_elstic(strain["33"], strain["33"])
        S1122 = scale*calc_elstic(strain["11"], strain["22"])
        S1133 = scale*calc_elstic(strain["11"], strain["33"])
        S2233 = scale*calc_elstic(strain["22"], strain["33"])
        S2323 = scale*calc_elstic(strain["23"], strain["23"])
        S1313 = scale*calc_elstic(strain["13"], strain["13"])
        S1212  = scale*calc_elstic(strain["12"], strain["12"])
        
        # All the below value should be very close to zero
        S1123 = scale*calc_elstic(strain["11"], strain["23"])
        S1113 = scale*calc_elstic(strain["11"], strain["13"])
        S1112 = scale*calc_elstic(strain["11"], strain["12"])
        
        # Convert Sijkl to  Spq 
        # NOTED that the Spq = Sijkl when p q = 1, 2, 3 BUT Spq = 4*Sijkl when p q = 4, 5, 6.
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
        # print(np.array([C11, C12, C44]))
        Cij[i] = np.array([C11, C12, C44])
    return Cij

def calc_ste(array): # calculate the standard error
    ste = np.zeros(3)
    for i in range(array.shape[1]):
        ste[i] = sqrt(sum(abs(array[:,i] - array[:,i].mean())**2))/len(array[:,i])
    return ste
        
thermo_slic("./thermo.out", "./thermo_slic.out")
Cij_gpumd_1ns = output_elastic_gpumd(10)

print("[C11, C12, C44] (unit in GPa) of Si estimated by gpumd with npt_scr for 1ns:")
cij=np.average(Cij_gpumd_1ns, axis=0)
print(cij)
#print(np.average(Cij_gpumd_1ns, axis=0))
print("With standard error (unit in GPa) of: ")
print(calc_ste(Cij_gpumd_1ns))
B=(cij[0]+2*cij[1])/3
print("The bulk modulus of Si estimated by gpumd with npt_scr for 1ns:")
print(B)
G= cij[2]
print("The shear modulus of Si estimated by gpumd with npt_scr for 1ns:")
print(G)
E=9*B*G/(3*B+G)
print("The Young's modulus of Si estimated by gpumd with npt_scr for 1ns:")
print(E)
EOF
